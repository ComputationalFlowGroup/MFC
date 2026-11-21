!>
!! @file
!! @brief Contains module m_hyperelastic

#:include 'macros.fpp'

!> @brief Computes the left Cauchy-Green deformation tensor and hyperelastic stress source terms, from the reference map (Reference
!! Map Technique, rmt) advected by m_riemann_solver_hll.fpp

module m_hyperelastic

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_finite_differences
    use m_graded, only: s_grade_Ca_inv, Ca_inv_graded_flag

    implicit none

    private; public :: s_hyperelastic_rmt_stress_update, s_initialize_hyperelastic_module, s_finalize_hyperelastic_module, &
        & s_hyperelastic_update_fd_coeffs

    ! The left Cauchy-Green tensor at cell centers. Computed once per RHS call and reused for both the Cauchy stress and
    ! the elastic strain energy, so it is kept at module scope rather than recomputed.
    type(vector_field) :: btensor
    $:GPU_DECLARE(create='[btensor]')

    real(wp), allocatable, dimension(:,:) :: fd_coeff_x_hyper
    real(wp), allocatable, dimension(:,:) :: fd_coeff_y_hyper
    real(wp), allocatable, dimension(:,:) :: fd_coeff_z_hyper
    $:GPU_DECLARE(create='[fd_coeff_x_hyper, fd_coeff_y_hyper, fd_coeff_z_hyper]')
    real(wp), allocatable, dimension(:) :: Ca_invs_hyper
    $:GPU_DECLARE(create='[Ca_invs_hyper]')
    real(wp), allocatable, dimension(:) :: mr_c2s_hyper  !< Mooney-Rivlin's second material constant (C2), per fluid
    $:GPU_DECLARE(create='[mr_c2s_hyper]')

contains

    !> Initialize the hyperelastic module
    impure subroutine s_initialize_hyperelastic_module

        integer :: i  !< generic iterator

        ! Sized to m_alloc/n_alloc/p_alloc (the largest AMR fine-block batch extent, not the coarse grid's own m/n/p): the
        ! AMR fine advance swaps m/n/p to a block's (or batch's) own values and reuses this same module-level workspace.

        @:ALLOCATE(btensor%vf(1:b_size))
        do i = 1, b_size
            @:ALLOCATE(btensor%vf(i)%sf(0:m_alloc, 0:n_alloc, 0:p_alloc))
        end do
        @:ACC_SETUP_VFs(btensor)

        @:ALLOCATE(Ca_invs_hyper(1:num_fluids), mr_c2s_hyper(1:num_fluids))
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            Ca_invs_hyper(i) = fluid_pp(i)%Ca_inv
            mr_c2s_hyper(i) = fluid_pp(i)%mr_c2
        end do
        $:GPU_UPDATE(device='[Ca_invs_hyper, mr_c2s_hyper]')

        ! s_compute_finite_difference_coefficients always extends fd_number beyond the interior on each side
        @:ALLOCATE(fd_coeff_x_hyper(-fd_number:fd_number, -fd_number:m_alloc + fd_number))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y_hyper(-fd_number:fd_number, -fd_number:n_alloc + fd_number))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z_hyper(-fd_number:fd_number, -fd_number:p_alloc + fd_number))
        end if

        ! Computing centered finite difference coefficients
        call s_hyperelastic_update_fd_coeffs()

    end subroutine s_initialize_hyperelastic_module

    !> (Re)compute the centered finite-difference coefficients from the current grid globals (m/n/p, x/y/z_cc) and push them to the
    !! device. Called at init and by the AMR fine-block swap/restore (m_amr_frame.fpp), where the grid globals flip between the
    !! coarse grid and a fine block: the coefficients are spacing-dependent, so reusing coarse ones on the fine grid would silently
    !! corrupt every deformation-gradient component the stress and wave-speed correction are built from.
    impure subroutine s_hyperelastic_update_fd_coeffs()

        ! the AMR fine-IB setup swaps grids before this module initializes (s_initialize_modules order); no RHS runs until
        ! after init, so skipping is correct -- init then computes the coarse coefficients and every subsequent swap/restore
        ! recomputes for the active grid
        if (.not. allocated(fd_coeff_x_hyper)) return

        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x_hyper, buff_size, fd_number, fd_order)
        $:GPU_UPDATE(device='[fd_coeff_x_hyper]')
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y_hyper, buff_size, fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_y_hyper]')
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z_hyper, buff_size, fd_number, fd_order)
            $:GPU_UPDATE(device='[fd_coeff_z_hyper]')
        end if

    end subroutine s_hyperelastic_update_fd_coeffs

    !> Compute the left Cauchy-Green deformation tensor and update the hyperelastic stress
    subroutine s_hyperelastic_rmt_stress_update(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf

        #:if USING_AMD
            real(wp), dimension(10) :: tensora, tensorb
        #:else
            real(wp), dimension(tensor_size) :: tensora, tensorb
        #:endif

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_k, alpha_rho_k
            real(wp), dimension(3) :: Ca_inv_graded_loc
        #:else
            real(wp), dimension(num_fluids) :: alpha_k, alpha_rho_k
            real(wp), dimension(num_fluids) :: Ca_inv_graded_loc
        #:endif
        real(wp) :: rho, gamma, pi_inf, qv
        real(wp) :: Ca_inv_local, mr_c2_local
        real(wp) :: xi_x, xi_y, xi_z
        integer  :: j, k, l, i, r

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l, r, alpha_K, alpha_rho_K, Ca_inv_graded_loc, rho, gamma, pi_inf, qv, &
                            & Ca_inv_local, mr_c2_local, tensora, tensorb, xi_x, xi_y, xi_z]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_species_fraction(q_cons_vf, j, k, l, alpha_rho_k, alpha_k)

                    xi_x = q_prim_vf(eqn_idx%xi%beg)%sf(j, k, l)
                    xi_y = 0._wp; xi_z = 0._wp
                    if (num_dims >= 2) xi_y = q_prim_vf(eqn_idx%xi%beg + 1)%sf(j, k, l)
                    if (num_dims == 3) xi_z = q_prim_vf(eqn_idx%xi%end)%sf(j, k, l)

                    ! Graded Ca_inv: default every entry to the static (non-graded) value, then override only the
                    ! fluids flagged graded_Ca_inv. A no-op when graded = F, since Ca_inv_graded_flag is then all-false.
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        Ca_inv_graded_loc(i) = Ca_invs_hyper(i)
                        if (graded) then
                            if (Ca_inv_graded_flag(i)) call s_grade_Ca_inv(xi_x, xi_y, xi_z, i, Ca_inv_graded_loc(i))
                        end if
                    end do

                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha_k, alpha_rho_k, &
                        & Ca_inv_K=Ca_inv_local, Ca_inv=Ca_inv_graded_loc)
                    rho = max(rho, sgm_eps)
                    Ca_inv_local = max(Ca_inv_local, sgm_eps)

                    ! Mooney-Rivlin's second material constant: a plain alpha-weighted mixture average, same as
                    ! Ca_inv's own static (non-graded) mixing -- grading C2 independently is not supported.
                    mr_c2_local = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        mr_c2_local = mr_c2_local + alpha_k(i)*mr_c2s_hyper(i)
                    end do

                    if (Ca_inv_local > verysmall) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, tensor_size
                            tensora(i) = 0._wp
                        end do
                        ! STEP 1-2b: grad_xi (tensora) and its determinant (tensorb(tensor_size)), branched on num_dims at
                        ! runtime (num_dims is a module variable here, not an fypp constant -- it is only available at Fypp
                        ! preprocessing time under case_optimization builds, see toolchain/mfc/case.py's __get_sim_fpp).
                        if (num_dims == 1) then
                            ! grad_xi is the scalar d(xi_x)/dx
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                tensora(1) = tensora(1) + q_prim_vf(eqn_idx%xi%beg)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                            end do
                            ! the "determinant" of a 1x1 matrix is itself
                            tensorb(tensor_size) = tensora(1)
                        else if (num_dims == 2) then
                            ! grad_xi is 2x2: tensora organization 1-2: dxix_dx, dxiy_dx; 3-4: dxix_dy, dxiy_dy
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                tensora(1) = tensora(1) + q_prim_vf(eqn_idx%xi%beg)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                                tensora(2) = tensora(2) + q_prim_vf(eqn_idx%xi%beg + 1)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                                tensora(3) = tensora(3) + q_prim_vf(eqn_idx%xi%beg)%sf(j, k + r, l)*fd_coeff_y_hyper(r, k)
                                tensora(4) = tensora(4) + q_prim_vf(eqn_idx%xi%beg + 1)%sf(j, k + r, l)*fd_coeff_y_hyper(r, k)
                            end do
                            ! adjoint of a 2x2 [[a,b],[c,d]] is [[d,-b],[-c,a]]
                            tensorb(1) = tensora(4)
                            tensorb(2) = -tensora(2)
                            tensorb(3) = -tensora(3)
                            tensorb(4) = tensora(1)
                            tensorb(tensor_size) = tensora(1)*tensora(4) - tensora(3)*tensora(2)
                        else
                            ! grad_xi is 3x3: tensora organization 1-3: dxix_dx, dxiy_dx, dxiz_dx; 4-6: dxix_dy, dxiy_dy,
                            ! dxiz_dy; 7-9: dxix_dz, dxiy_dz, dxiz_dz
                            $:GPU_LOOP(parallelism='[seq]')
                            do r = -fd_number, fd_number
                                ! derivatives in the x-direction
                                tensora(1) = tensora(1) + q_prim_vf(eqn_idx%xi%beg)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                                tensora(2) = tensora(2) + q_prim_vf(eqn_idx%xi%beg + 1)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                                tensora(3) = tensora(3) + q_prim_vf(eqn_idx%xi%end)%sf(j + r, k, l)*fd_coeff_x_hyper(r, j)
                                ! derivatives in the y-direction
                                tensora(4) = tensora(4) + q_prim_vf(eqn_idx%xi%beg)%sf(j, k + r, l)*fd_coeff_y_hyper(r, k)
                                tensora(5) = tensora(5) + q_prim_vf(eqn_idx%xi%beg + 1)%sf(j, k + r, l)*fd_coeff_y_hyper(r, k)
                                tensora(6) = tensora(6) + q_prim_vf(eqn_idx%xi%end)%sf(j, k + r, l)*fd_coeff_y_hyper(r, k)
                                ! derivatives in the z-direction
                                tensora(7) = tensora(7) + q_prim_vf(eqn_idx%xi%beg)%sf(j, k, l + r)*fd_coeff_z_hyper(r, l)
                                tensora(8) = tensora(8) + q_prim_vf(eqn_idx%xi%beg + 1)%sf(j, k, l + r)*fd_coeff_z_hyper(r, l)
                                tensora(9) = tensora(9) + q_prim_vf(eqn_idx%xi%end)%sf(j, k, l + r)*fd_coeff_z_hyper(r, l)
                            end do
                            ! adjoint of the 3x3 grad_xi tensor, for the inverse
                            tensorb(1) = tensora(5)*tensora(9) - tensora(6)*tensora(8)
                            tensorb(2) = -(tensora(2)*tensora(9) - tensora(3)*tensora(8))
                            tensorb(3) = tensora(2)*tensora(6) - tensora(3)*tensora(5)
                            tensorb(4) = -(tensora(4)*tensora(9) - tensora(6)*tensora(7))
                            tensorb(5) = tensora(1)*tensora(9) - tensora(3)*tensora(7)
                            tensorb(6) = -(tensora(1)*tensora(6) - tensora(4)*tensora(3))
                            tensorb(7) = tensora(4)*tensora(8) - tensora(5)*tensora(7)
                            tensorb(8) = -(tensora(1)*tensora(8) - tensora(2)*tensora(7))
                            tensorb(9) = tensora(1)*tensora(5) - tensora(2)*tensora(4)
                            tensorb(tensor_size) = tensora(1)*(tensora(5)*tensora(9) - tensora(6)*tensora(8)) - tensora(2) &
                                    & *(tensora(4)*tensora(9) - tensora(6)*tensora(7)) + tensora(3)*(tensora(4)*tensora(8) &
                                    & - tensora(5)*tensora(7))
                        end if

                        if (tensorb(tensor_size) > verysmall) then
                            ! STEP 2c-3: F = (grad_xi)^-1 (tensora becomes F), J = det(F), B = F F^T -- same per-dimension
                            ! branch as STEP 1-2b above
                            if (num_dims == 1) then
                                ! F = 1/grad_xi; J = det(F) = F (1x1 again); B_xx = F**2
                                tensora(1) = 1._wp/tensorb(tensor_size)
                                tensorb(tensor_size) = tensora(1)
                                tensorb(1) = tensora(1)**2
                            else if (num_dims == 2) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, tensor_size - 1
                                    tensora(i) = tensorb(i)/tensorb(tensor_size)
                                end do
                                tensorb(tensor_size) = 1._wp/tensorb(tensor_size)
                                ! B = F F^T: B_xx, B_yy, B_xy (tensora organization 1-2: F_xx, F_yx; 3-4: F_xy, F_yy)
                                tensorb(1) = tensora(1)**2 + tensora(3)**2
                                tensorb(3) = tensora(2)**2 + tensora(4)**2
                                tensorb(2) = tensora(1)*tensora(2) + tensora(3)*tensora(4)
                            else
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, tensor_size - 1
                                    tensora(i) = tensorb(i)/tensorb(tensor_size)
                                end do
                                tensorb(tensor_size) = 1._wp/tensorb(tensor_size)
                                ! B = F F^T, the left Cauchy-Green tensor
                                tensorb(1) = tensora(1)**2 + tensora(2)**2 + tensora(3)**2
                                tensorb(5) = tensora(4)**2 + tensora(5)**2 + tensora(6)**2
                                tensorb(9) = tensora(7)**2 + tensora(8)**2 + tensora(9)**2
                                tensorb(2) = tensora(1)*tensora(4) + tensora(2)*tensora(5) + tensora(3)*tensora(6)
                                tensorb(3) = tensora(1)*tensora(7) + tensora(2)*tensora(8) + tensora(3)*tensora(9)
                                tensorb(6) = tensora(4)*tensora(7) + tensora(5)*tensora(8) + tensora(6)*tensora(9)
                            end if

                            ! STEP 4: update the btensor, in the same symmetric order (xx, xy, yy, xz, yz, zz) eqn_idx%stress
                            ! uses, truncated to b_size - 1 entries for num_dims < 3
                            if (num_dims == 1) then
                                btensor%vf(1)%sf(j, k, l) = tensorb(1)
                            else if (num_dims == 2) then
                                btensor%vf(1)%sf(j, k, l) = tensorb(1)
                                btensor%vf(2)%sf(j, k, l) = tensorb(2)
                                btensor%vf(3)%sf(j, k, l) = tensorb(3)
                            else
                                btensor%vf(1)%sf(j, k, l) = tensorb(1)
                                btensor%vf(2)%sf(j, k, l) = tensorb(2)
                                btensor%vf(3)%sf(j, k, l) = tensorb(5)
                                btensor%vf(4)%sf(j, k, l) = tensorb(3)
                                btensor%vf(5)%sf(j, k, l) = tensorb(6)
                                btensor%vf(6)%sf(j, k, l) = tensorb(9)
                            end if
                            ! store the determinant at the last entry of the btensor
                            btensor%vf(b_size)%sf(j, k, l) = tensorb(tensor_size)
                            ! STEP 5a: updating the Cauchy stress primitive scalar field
                            if (hyper_model == 1) then
                                call s_neoHookean_cauchy_solver(btensor%vf, q_prim_vf, Ca_inv_local, j, k, l)
                            else if (hyper_model == 2) then
                                call s_Mooney_Rivlin_cauchy_solver(btensor%vf, q_prim_vf, Ca_inv_local, mr_c2_local, j, k, l)
                            end if
                            ! STEP 5b: updating the pressure field (subtracting the elastic strain energy). Both Cauchy
                            ! solvers now stash the FULL strain energy (already scaled by their own material
                            ! constant(s)) at xi%end+1, not energy-per-unit-modulus, since Mooney-Rivlin has two
                            ! independent constants that do not scale together.
                            q_prim_vf(eqn_idx%E)%sf(j, k, l) = q_prim_vf(eqn_idx%E)%sf(j, k, &
                                      & l) - q_prim_vf(eqn_idx%xi%end + 1)%sf(j, k, l)/gamma
                            ! STEP 5c: updating the Cauchy stress conservative scalar field
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, b_size - 1
                                q_cons_vf(eqn_idx%stress%beg + i - 1)%sf(j, k, &
                                          & l) = rho*q_prim_vf(eqn_idx%stress%beg + i - 1)%sf(j, k, l)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_hyperelastic_rmt_stress_update

    !> Trace of B (the full 3D invariant I1, with (3 - num_dims) added for the directions frozen at identity stretch -- see
    !! misc/hyperelastic_acoustic_tensor.py's plane_strain_reduction() for the derivation) and the in-place deviatoric projection
    !! dev(B) = B - tr(B)/3 * I. Shared by both Cauchy-stress solvers below. The diagonal positions in the symmetric (xx, xy, yy,
    !! xz, yz, zz) storage order are the triangular numbers k(k+1)/2 for k = 1..num_dims, so a single runtime loop over num_dims (a
    !! module variable, not an Fypp constant -- see the note in s_hyperelastic_rmt_stress_update) covers every dimensionality with
    !! no branching at all.
    subroutine s_hyperelastic_trace_and_deviator(btensor_in, j, k, l, trace)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(b_size), intent(inout) :: btensor_in
        integer, intent(in)                                  :: j, k, l
        real(wp), intent(out)                                :: trace
        real(wp), parameter                                  :: f13 = 1._wp/3._wp
        integer                                              :: dk, diag_idx

        trace = real(3 - num_dims, wp)
        $:GPU_LOOP(parallelism='[seq]')
        do dk = 1, num_dims
            trace = trace + btensor_in(dk*(dk + 1)/2)%sf(j, k, l)
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do dk = 1, num_dims
            diag_idx = dk*(dk + 1)/2
            btensor_in(diag_idx)%sf(j, k, l) = btensor_in(diag_idx)%sf(j, k, l) - f13*trace
        end do

    end subroutine s_hyperelastic_trace_and_deviator

    !> Compute the neo-Hookean Cauchy stress from the left Cauchy-Green tensor
    subroutine s_neoHookean_cauchy_solver(btensor_in, q_prim_vf, Ca_inv_param, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(inout)   :: btensor_in
        real(wp), intent(in)                                   :: Ca_inv_param
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: trace
        integer                                                :: i  !< Generic loop iterators

        call s_hyperelastic_trace_and_deviator(btensor_in, j, k, l, trace)

        ! dividing by the jacobian for neo-Hookean model setting the tensor to the stresses for riemann solver
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, b_size - 1
            q_prim_vf(eqn_idx%stress%beg + i - 1)%sf(j, k, l) = Ca_inv_param*btensor_in(i)%sf(j, k, l)/btensor_in(b_size)%sf(j, &
                      & k, l)
        end do
        ! First invariant strain energy: W = Ca_inv/2 * (I1 - 3), neo-Hookean model. Stashed as the FULL energy
        ! (already scaled by Ca_inv_param), not energy-per-unit-modulus -- see the STEP 5b comment in
        ! s_hyperelastic_rmt_stress_update for why.
        q_prim_vf(eqn_idx%xi%end + 1)%sf(j, k, l) = 0.5_wp*Ca_inv_param*(trace - 3.0_wp)/btensor_in(b_size)%sf(j, k, l)
        ! tr(B) and J = det(F), stashed for the Riemann solver's direct wave speed (see m_riemann_solver_hll.fpp)
        q_prim_vf(eqn_idx%xi%end + 2)%sf(j, k, l) = trace
        q_prim_vf(eqn_idx%xi%end + 3)%sf(j, k, l) = btensor_in(b_size)%sf(j, k, l)

    end subroutine s_neoHookean_cauchy_solver

    !> Compute the Mooney-Rivlin Cauchy stress from the left Cauchy-Green tensor: tau = (1/J)*(2*C1*dev(B) - 2*C2*dev(B^-1)), C1 =
    !! Ca_inv_param/2, C2 = mr_c2_param (reduces exactly to neo-Hookean at C2 = 0). Derived and verified in
    !! misc/hyperelastic_acoustic_tensor.py's mooney_rivlin() (the continuum form) and mooney_rivlin_implementation_check() (the
    !! discrete identities used below: B^-1 computed via its own adjoint/determinant rather than a second grad_xi inversion, the
    !! second invariant I2 = 0.5*(tr(B)^2 - tr(B^2)), and the 1D/2D frozen-direction convention for tr(B^2)).
    subroutine s_Mooney_Rivlin_cauchy_solver(btensor_in, q_prim_vf, Ca_inv_param, mr_c2_param, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), dimension(b_size), intent(inout)   :: btensor_in
        real(wp), intent(in)                                   :: Ca_inv_param, mr_c2_param
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: trB, trB2, trBinv, I2, C1, C2, J_det, detB
        real(wp)                                               :: sum_all, sum_diag
        real(wp), parameter                                    :: f13 = 1._wp/3._wp
        real(wp), dimension(6)                                 :: Binv
        integer                                                :: i, dk, diag_idx

        C1 = 0.5_wp*Ca_inv_param
        C2 = mr_c2_param
        J_det = btensor_in(b_size)%sf(j, k, l)
        detB = J_det*J_det

        ! tr(B^2), via tr(B^2) = 2*sum(all stored components^2) - sum(diagonal components^2) + (3 - num_dims) for the
        ! frozen directions' unit diagonal entries (each off-diagonal component appears twice in the full symmetric
        ! matrix, each diagonal component once) -- avoids a per-component diagonal/off-diagonal branch.
        sum_all = 0._wp; sum_diag = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, b_size - 1
            sum_all = sum_all + btensor_in(i)%sf(j, k, l)**2
        end do
        $:GPU_LOOP(parallelism='[seq]')
        do dk = 1, num_dims
            sum_diag = sum_diag + btensor_in(dk*(dk + 1)/2)%sf(j, k, l)**2
        end do
        trB2 = 2._wp*sum_all - sum_diag + real(3 - num_dims, wp)

        ! B^-1 via its own adjoint/determinant (determinant of B is J**2), the same cofactor pattern already used to
        ! invert grad_xi in STEP 2c above -- computed from B directly (not grad_xi) so this routine does not need
        ! grad_xi's pre-inversion value, which is no longer available by this point in the caller.
        if (num_dims == 1) then
            Binv(1) = 1._wp/btensor_in(1)%sf(j, k, l)
        else if (num_dims == 2) then
            Binv(1) = btensor_in(3)%sf(j, k, l)/detB
            Binv(3) = btensor_in(1)%sf(j, k, l)/detB
            Binv(2) = -btensor_in(2)%sf(j, k, l)/detB
        else
            Binv(1) = (btensor_in(3)%sf(j, k, l)*btensor_in(6)%sf(j, k, l) - btensor_in(5)%sf(j, k, l)**2)/detB
            Binv(3) = (btensor_in(1)%sf(j, k, l)*btensor_in(6)%sf(j, k, l) - btensor_in(4)%sf(j, k, l)**2)/detB
            Binv(6) = (btensor_in(1)%sf(j, k, l)*btensor_in(3)%sf(j, k, l) - btensor_in(2)%sf(j, k, l)**2)/detB
            Binv(2) = (btensor_in(4)%sf(j, k, l)*btensor_in(5)%sf(j, k, l) - btensor_in(2)%sf(j, k, l)*btensor_in(6)%sf(j, k, &
                 & l))/detB
            Binv(4) = (btensor_in(2)%sf(j, k, l)*btensor_in(5)%sf(j, k, l) - btensor_in(3)%sf(j, k, l)*btensor_in(4)%sf(j, k, &
                 & l))/detB
            Binv(5) = (btensor_in(2)%sf(j, k, l)*btensor_in(4)%sf(j, k, l) - btensor_in(1)%sf(j, k, l)*btensor_in(5)%sf(j, k, &
                 & l))/detB
        end if

        trBinv = real(3 - num_dims, wp)
        $:GPU_LOOP(parallelism='[seq]')
        do dk = 1, num_dims
            trBinv = trBinv + Binv(dk*(dk + 1)/2)
        end do

        ! dev(B^-1) in place in the local Binv array (dev(B) itself is taken care of by the trace-and-deviator call
        ! below, which mutates btensor_in into dev(B) in place -- after that call, both arrays hold the deviatoric
        ! tensor, so the stress loop further down needs no diagonal/off-diagonal branch either)
        $:GPU_LOOP(parallelism='[seq]')
        do dk = 1, num_dims
            diag_idx = dk*(dk + 1)/2
            Binv(diag_idx) = Binv(diag_idx) - f13*trBinv
        end do

        call s_hyperelastic_trace_and_deviator(btensor_in, j, k, l, trB)
        I2 = 0.5_wp*(trB*trB - trB2)

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, b_size - 1
            q_prim_vf(eqn_idx%stress%beg + i - 1)%sf(j, k, l) = (2._wp*C1*btensor_in(i)%sf(j, k, l) - 2._wp*C2*Binv(i))/J_det
        end do

        ! Strain energy W = C1*(I1 - 3) + C2*(I2 - 3), stashed as the full energy (see the STEP 5b comment above)
        q_prim_vf(eqn_idx%xi%end + 1)%sf(j, k, l) = (C1*(trB - 3._wp) + C2*(I2 - 3._wp))/J_det
        ! tr(B) and J = det(F), stashed for the Riemann solver's direct wave speed (see m_riemann_solver_hll.fpp)
        q_prim_vf(eqn_idx%xi%end + 2)%sf(j, k, l) = trB
        q_prim_vf(eqn_idx%xi%end + 3)%sf(j, k, l) = J_det

    end subroutine s_Mooney_Rivlin_cauchy_solver

    !> Finalize the hyperelastic module
    impure subroutine s_finalize_hyperelastic_module()

        integer :: i  !< iterator

        do i = 1, b_size
            @:DEALLOCATE(btensor%vf(i)%sf)
        end do
        @:DEALLOCATE(Ca_invs_hyper, mr_c2s_hyper)
        @:DEALLOCATE(fd_coeff_x_hyper)
        if (n > 0) then
            @:DEALLOCATE(fd_coeff_y_hyper)
            if (p > 0) then
                @:DEALLOCATE(fd_coeff_z_hyper)
            end if
        end if

    end subroutine s_finalize_hyperelastic_module

end module m_hyperelastic
