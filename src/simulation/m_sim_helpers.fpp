!>
!! @file
!! @brief Contains module m_sim_helpers

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Simulation helper routines for cell state, CFL calculation, and stability checks
module m_sim_helpers

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_riemann_state, only: f_elastic_signal_speed, f_hyperelastic_signal_speed, Ca_inv_2s_rs, Ca_invs_rs, Res_gs
    use m_graded, only: s_grade_Ca_inv, Ca_inv_graded_flag, s_grade_Re, Re_graded_flag

    implicit none

    private; public :: s_compute_cell_state, s_compute_cell_diffusivity, s_compute_stability_from_dt, s_compute_dt_from_cfl, &
        & s_amr_scale_stability, s_compute_elastic_signal_speeds, dt_limiter, dt_limiter_names

    !> Criterion currently limiting the adaptive time step (ICFL, VCFL, CCFL, TCFL, the collision cap, or the ramp limiter)
    character(len=4)                          :: dt_limiter = 'none'
    character(len=4), dimension(6), parameter :: dt_limiter_names = (/'ICFL', 'VCFL', 'CCFL', 'TCFL', 'COLL', 'WCOL'/)

contains

    !> Brackbill's capillary density (rho_1 + rho_2)/2 from the phase densities of an interface cell. Zero where a phase is absent,
    !! so single-phase cells, which carry no capillary force, impose no capillary limit
    subroutine s_compute_capillary_rho(alpha, alpha_rho, rho_c)

        $:GPU_ROUTINE(function_name='s_compute_capillary_rho', parallelism='[seq]', cray_inline=True)
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha, alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha, alpha_rho
        #:endif
        real(wp), intent(out) :: rho_c
        integer               :: i

        rho_c = 0._wp
        if (minval(alpha(1:num_fluids)) < capillary_alpha_min) return
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_c = rho_c + alpha_rho(i)/alpha(i)
        end do
        rho_c = rho_c/real(num_fluids, wp)

    end subroutine s_compute_capillary_rho

    !> Computes the modified dtheta for Fourier filtering in azimuthal direction
    function f_compute_filtered_dtheta(k, l) result(fltr_dtheta)

        $:GPU_ROUTINE(parallelism='[seq]')
        integer, intent(in) :: k, l
        real(wp)            :: fltr_dtheta
        integer             :: Nfq

        if (grid_geometry == 3) then
            if (k == 0) then
                fltr_dtheta = 2._wp*pi*y_cb(0)/3._wp
            else if (k <= fourier_rings) then
                Nfq = min(floor(2._wp*real(k, wp)*pi), (p + 1)/2 + 1)
                fltr_dtheta = 2._wp*pi*y_cb(k - 1)/real(Nfq, wp)
            else
                fltr_dtheta = y_cb(k - 1)*dz(l)
            end if
        else
            fltr_dtheta = 0._wp
        end if

    end function f_compute_filtered_dtheta

    !> Sum over directions of 1/spacing^2 at cell (j, k, l). An explicit diffusion with diffusivity D is RK-stable for D*dt*sum <=
    !! 2.51/4 (RK3; 2/4 for RK1 and RK2), so the diffusive CFL numbers are D*dt*sum
    subroutine s_compute_inv_dx2_sum(j, k, l, s)

        $:GPU_ROUTINE(function_name='s_compute_inv_dx2_sum', parallelism='[seq]', cray_inline=True)
        integer, intent(in)   :: j, k, l
        real(wp), intent(out) :: s

        s = 1._wp/dx(j)**2
        if (n > 0) s = s + 1._wp/dy(k)**2
        if (p > 0) then
            if (grid_geometry == 3) then
                s = s + 1._wp/f_compute_filtered_dtheta(k, l)**2
            else
                s = s + 1._wp/dz(l)**2
            end if
        end if

    end subroutine s_compute_inv_dx2_sum

    !> Momentum diffusivity of the viscous stress: normal stresses diffuse with (4/3 mu + mu_b)/rho, the largest coefficient
    subroutine s_compute_visc_diffusivity(Re_l, rho, nu)

        $:GPU_ROUTINE(function_name='s_compute_visc_diffusivity', parallelism='[seq]', cray_inline=True)
        real(wp), dimension(2), intent(in) :: Re_l
        real(wp), intent(in)               :: rho
        real(wp), intent(out)              :: nu

        nu = (4._wp/(3._wp*Re_l(1)) + 1._wp/Re_l(2))/rho

    end subroutine s_compute_visc_diffusivity

    !> Computes the mixture coefficients, velocity and pressure of one cell
    subroutine s_compute_cell_state(q_prim_vf, pres, rho, gamma, pi_inf, Re, alpha, alpha_rho, vel, vel_sum, qv, j, k, l, Ca_inv_K)

        $:GPU_ROUTINE(function_name='s_compute_cell_state',parallelism='[seq]', cray_inline=True)

        type(scalar_field), intent(in), dimension(sys_size) :: q_prim_vf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), intent(inout), dimension(3) :: alpha, alpha_rho
            real(wp), intent(inout), dimension(3) :: vel
        #:else
            real(wp), intent(inout), dimension(num_fluids) :: alpha, alpha_rho
            real(wp), intent(inout), dimension(num_vels)   :: vel
        #:endif
        real(wp), intent(inout)               :: rho, gamma, pi_inf, vel_sum, pres
        real(wp), intent(out)                 :: qv
        integer, intent(in)                   :: j, k, l
        real(wp), dimension(2), intent(inout) :: Re
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: Ca_inv_graded_loc
        #:else
            real(wp), dimension(num_fluids) :: Ca_inv_graded_loc
        #:endif
        real(wp), dimension(2, Re_size_max) :: Res_graded_loc
        !> Mixture inverse Cauchy number, for callers that need the elastic signal speed (see s_compute_elastic_signal_speeds).
        !! Graded where the case grades it, so the automatic-CFL bound built from it is correct for a spatially varying modulus too.
        real(wp), optional, intent(out) :: Ca_inv_K
        real(wp)                        :: Ca_inv_local
        real(wp)                        :: xi_x, xi_y, xi_z
        integer                         :: i, j_re

        call s_compute_species_fraction(q_prim_vf, j, k, l, alpha_rho, alpha)

        ! Graded material properties are evaluated at the reference map, so build this cell's overrides first -- the same
        ! default-then-override loops m_hyperelastic.fpp and m_riemann_solver_hll.fpp use, so the dt bound this feeds sees
        ! the same moduli and viscosities the fluxes do. Absent, the conversion kernel falls back to its own static tables.
        if (graded) then
            xi_x = q_prim_vf(eqn_idx%xi%beg)%sf(j, k, l)
            xi_y = 0._wp; xi_z = 0._wp
            if (num_dims >= 2) xi_y = q_prim_vf(eqn_idx%xi%beg + 1)%sf(j, k, l)
            if (num_dims == 3) xi_z = q_prim_vf(eqn_idx%xi%end)%sf(j, k, l)

            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                Ca_inv_graded_loc(i) = Ca_invs_rs(i)
                if (Ca_inv_graded_flag(i)) call s_grade_Ca_inv(xi_x, xi_y, xi_z, i, Ca_inv_graded_loc(i))
            end do

            if (viscous) then
                $:GPU_LOOP(parallelism='[seq]')
                do i = 1, 2
                    $:GPU_LOOP(parallelism='[seq]')
                    do j_re = 1, Re_size_max
                        Res_graded_loc(i, j_re) = Res_gs(i, j_re)
                    end do
                    $:GPU_LOOP(parallelism='[seq]')
                    do j_re = 1, Re_size(i)
                        if (Re_graded_flag(Re_idx(i, j_re))) then
                            call s_grade_Re(xi_x, xi_y, xi_z, Re_idx(i, j_re), i, Res_graded_loc(i, j_re))
                        end if
                    end do
                end do
            end if
        end if

        Ca_inv_local = 0._wp
        if (hypoelasticity .or. hyperelasticity) then
            if (graded .and. viscous) then
                call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, Ca_inv_local, &
                    & Ca_inv_graded_loc, Res_graded_loc)
            else if (graded) then
                call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, Ca_inv_local, &
                    & Ca_inv_graded_loc)
            else
                call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, Ca_inv_local)
            end if
        else if (graded .and. viscous) then
            call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re, Res=Res_graded_loc)
        else
            call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv, alpha, alpha_rho, Re)
        end if
        if (present(Ca_inv_K)) Ca_inv_K = Ca_inv_local

        if (igr) then
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel(i) = q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)/rho
            end do
        else
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_vels
                vel(i) = q_prim_vf(eqn_idx%cont%end + i)%sf(j, k, l)
            end do
        end if

        vel_sum = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_vels
            vel_sum = vel_sum + vel(i)**2._wp
        end do

        if (igr) then
            pres = (q_prim_vf(eqn_idx%E)%sf(j, k, l) - pi_inf - qv - 5.e-1_wp*rho*vel_sum)/gamma
        else
            pres = q_prim_vf(eqn_idx%E)%sf(j, k, l)
        end if

    end subroutine s_compute_cell_state

    !> Per-direction signal speed for the advective CFL criterion. For an elastic material the fastest characteristic is not the
    !! acoustic speed c but c stiffened by the elastic tangent modulus, and that stiffening is DIRECTION-DEPENDENT, so each
    !! direction gets its own speed: an isotropic worst case would be needlessly restrictive in some directions and wrong in others.
    !! Entries beyond num_dims are filled with c and never read. Reduces to c exactly when neither elastic model is active, so
    !! callers can pass the result unconditionally.
    subroutine s_compute_elastic_signal_speeds(q_prim_vf, c, Ca_inv_K, alpha, rho, j, k, l, c_dir)

        $:GPU_ROUTINE(function_name='s_compute_elastic_signal_speeds', parallelism='[seq]', cray_inline=True)

        type(scalar_field), intent(in), dimension(sys_size) :: q_prim_vf
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha
        #:endif
        real(wp), intent(in)                :: c, Ca_inv_K, rho
        integer, intent(in)                 :: j, k, l
        real(wp), dimension(3), intent(out) :: c_dir
        real(wp)                            :: trB, trBinv, J_det, Ca_inv_2_K
        integer                             :: d, diag_idx

        c_dir = c
        if (.not. (hypoelasticity .or. hyperelasticity)) return

        if (hypoelasticity) then
            ! The transported stress is already the right linearization point for a rate law: tau_nn for the
            ! direction in question is the diagonal entry at the triangular position d(d+1)/2.
            do d = 1, num_dims
                diag_idx = d*(d + 1)/2
                c_dir(d) = f_elastic_signal_speed(c, Ca_inv_K, q_prim_vf(eqn_idx%stress%beg - 1 + diag_idx)%sf(j, k, l), rho)
            end do
        else
            ! Finite strain: read the cell-centered stash m_hyperelastic.fpp wrote this step -- tr(B), J, tr(B^-1)
            ! and, per direction, the aligned diagonal entries of B and B^-1 (see m_global_parameters_common.fpp).
            trB = q_prim_vf(eqn_idx%trB)%sf(j, k, l)
            J_det = q_prim_vf(eqn_idx%Jdet)%sf(j, k, l)
            trBinv = q_prim_vf(eqn_idx%trBinv)%sf(j, k, l)
            Ca_inv_2_K = 0._wp
            do d = 1, num_fluids
                Ca_inv_2_K = Ca_inv_2_K + alpha(d)*Ca_inv_2s_rs(d)
            end do
            do d = 1, num_dims
                c_dir(d) = f_hyperelastic_signal_speed(c, Ca_inv_K, Ca_inv_2_K, trB, q_prim_vf(eqn_idx%Bdiag%beg - 1 + d)%sf(j, &
                      & k, l), trBinv, q_prim_vf(eqn_idx%Binvdiag%beg - 1 + d)%sf(j, k, l), J_det, rho)
            end do
        end if

    end subroutine s_compute_elastic_signal_speeds

    !> Thermal diffusivity k/(rho cv) of cell state (alpha, alpha_rho) for the explicit conduction step limit, zero without it
    subroutine s_compute_cell_diffusivity(alpha, alpha_rho, Dth)

        $:GPU_ROUTINE(function_name='s_compute_cell_diffusivity', parallelism='[seq]', cray_inline=True)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha, alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha, alpha_rho
        #:endif
        real(wp), intent(out) :: Dth
        real(wp)              :: k_mix, rho_cv
        integer               :: i

        Dth = 0._wp
        if (heat_conduction) then
            k_mix = 0._wp
            rho_cv = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                k_mix = k_mix + alpha(i)*fluid_k_therm(i)
                rho_cv = rho_cv + alpha_rho(i)*cvs(i)
            end do
            Dth = k_mix/max(rho_cv, sgm_eps)
        end if

    end subroutine s_compute_cell_diffusivity

    !> Computes stability criterion for a specified dt
    subroutine s_compute_stability_from_dt(vel, c, rho, Re_l, alpha, alpha_rho, Dth, j, k, l, icfl, vcfl, Rc, ccfl, tcfl, c_dir)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in), dimension(num_vels) :: vel
        real(wp), intent(in)                      :: c, rho
        real(wp), intent(inout)                   :: icfl
        real(wp), intent(inout)                   :: vcfl, Rc, ccfl, tcfl
        real(wp), dimension(2), intent(in)        :: Re_l
        real(wp), intent(in)                      :: Dth
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha, alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha, alpha_rho
        #:endif
        integer, intent(in) :: j, k, l
        !> Per-direction advective signal speed (s_compute_elastic_signal_speeds). Absent, the acoustic c is used in every
        !! direction, which is the correct ICFL only when no elastic model is active. Rc below deliberately keeps c: it reports
        !! viscous-against-advective scaling, not a wave-speed bound.
        real(wp), dimension(3), optional, intent(in) :: c_dir
        real(wp)                                     :: fltr_dtheta
        real(wp)                                     :: rho_c, nu, inv_dx2
        real(wp), dimension(3)                       :: a

        call s_compute_inv_dx2_sum(j, k, l, inv_dx2)

        a = c
        if (present(c_dir)) a = c_dir

        ! Inviscid CFL calculation
        ! The multi-dimensional CFL terms are written out here rather than
        ! obtained from a shared helper procedure: NVHPC 25.5's fort2 segfaults
        ! when a routine containing a call to that helper is cross-file inlined
        ! by -Minline (the IPO setup in cmake/MFCTargets.cmake).
        if (p > 0) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                if (grid_geometry == 3) then
                    fltr_dtheta = f_compute_filtered_dtheta(k, l)
                    icfl = dt/min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)), fltr_dtheta/(abs(vel(3)) + a(3)))
                else
                    icfl = dt/min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)), dz(l)/(abs(vel(3)) + a(3)))
                end if
            #:endif
        else if (n > 0) then
            icfl = dt/min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)))
        else
            icfl = (dt/dx(j))*(abs(vel(1)) + a(1))
        end if

        ! Viscous calculations
        if (viscous) then
            call s_compute_visc_diffusivity(Re_l, rho, nu)
            vcfl = dt*nu*inv_dx2
            if (p > 0) then
                #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                    if (grid_geometry == 3) then
                        fltr_dtheta = f_compute_filtered_dtheta(k, l)
                        Rc = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), fltr_dtheta*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                    else
                        Rc = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c), dz(l)*(abs(vel(3)) + c))/maxval(1._wp/Re_l)
                    end if
                #:endif
            else if (n > 0) then
                Rc = min(dx(j)*(abs(vel(1)) + c), dy(k)*(abs(vel(2)) + c))/maxval(1._wp/Re_l)
            else
                Rc = dx(j)*(abs(vel(1)) + c)/maxval(1._wp/Re_l)
            end if
        end if

        ! Capillary CFL calculation
        if (surface_tension) then
            ccfl = 0._wp
            call s_compute_capillary_rho(alpha, alpha_rho, rho_c)
            if (rho_c > 0._wp) then
                if (p > 0) then
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (grid_geometry == 3) then
                            fltr_dtheta = f_compute_filtered_dtheta(k, l)
                            ccfl = dt*sqrt(2._wp*pi*sigma/(rho_c*min(dx(j), dy(k), fltr_dtheta)**3._wp))
                        else
                            ccfl = dt*sqrt(2._wp*pi*sigma/(rho_c*min(dx(j), dy(k), dz(l))**3._wp))
                        end if
                    #:endif
                else if (n > 0) then
                    ccfl = dt*sqrt(2._wp*pi*sigma/(rho_c*min(dx(j), dy(k))**3._wp))
                else
                    ccfl = dt*sqrt(2._wp*pi*sigma/(rho_c*dx(j)**3._wp))
                end if
            end if
        end if

        ! Thermal diffusion CFL
        tcfl = dt*Dth*inv_dx2

    end subroutine s_compute_stability_from_dt

    !> The stability numbers of the finest cell over a coarse cell, from the coarse cell's: rf is its spacing ratio (amr_cfl_rf) and
    !! pw1, pw2, pw3 the powers of it for the inviscid, diffusive and capillary criteria (amr_cfl_pow).
    subroutine s_amr_scale_stability(rf, pw1, pw2, pw3, icfl, vcfl, Rc, ccfl, tcfl)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)    :: rf, pw1, pw2, pw3
        real(wp), intent(inout) :: icfl, vcfl, Rc, ccfl, tcfl

        icfl = icfl*rf**pw1
        vcfl = vcfl*rf**pw2
        ccfl = ccfl*rf**pw3
        tcfl = tcfl*rf**pw2
        Rc = Rc/rf

    end subroutine s_amr_scale_stability

    !> Computes the candidate dts for a specified CFL number: max_dt(1) from the inviscid, max_dt(2) the viscous, max_dt(3) the
    !! capillary, and max_dt(4) the thermal diffusion criterion (huge where the criterion is inactive)
    subroutine s_compute_dt_from_cfl(vel, c, max_dt, rho, Re_l, alpha, alpha_rho, Dth, j, k, l, c_dir)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), dimension(num_vels), intent(in) :: vel
        real(wp), intent(in)                      :: c, rho
        real(wp), dimension(4), intent(out)       :: max_dt
        real(wp), dimension(2), intent(in)        :: Re_l
        real(wp), intent(in)                      :: Dth
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: alpha, alpha_rho
        #:else
            real(wp), dimension(num_fluids), intent(in) :: alpha, alpha_rho
        #:endif
        integer, intent(in) :: j, k, l
        !> Per-direction advective signal speed (s_compute_elastic_signal_speeds). Absent, the acoustic c is used in every
        !! direction, which is the correct bound only when no elastic model is active.
        real(wp), dimension(3), optional, intent(in) :: c_dir
        real(wp)                                     :: ccfl_dt, rho_c, hmin, nu, inv_dx2
        real(wp)                                     :: fltr_dtheta
        real(wp), dimension(3)                       :: a

        max_dt(2) = huge(1._wp)
        max_dt(3) = huge(1._wp)
        max_dt(4) = huge(1._wp)

        call s_compute_inv_dx2_sum(j, k, l, inv_dx2)

        a = c
        if (present(c_dir)) a = c_dir

        ! Inviscid CFL calculation
        ! The multi-dimensional CFL terms are written out here rather than
        ! obtained from a shared helper procedure: NVHPC 25.5's fort2 segfaults
        ! when a routine containing a call to that helper is cross-file inlined
        ! by -Minline (the IPO setup in cmake/MFCTargets.cmake).
        if (p > 0) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                if (grid_geometry == 3) then
                    fltr_dtheta = f_compute_filtered_dtheta(k, l)
                    max_dt(1) = cfl_target*min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)), &
                           & fltr_dtheta/(abs(vel(3)) + a(3)))
                else
                    max_dt(1) = cfl_target*min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)), dz(l)/(abs(vel(3)) + a(3)))
                end if
            #:endif
        else if (n > 0) then
            max_dt(1) = cfl_target*min(dx(j)/(abs(vel(1)) + a(1)), dy(k)/(abs(vel(2)) + a(2)))
        else
            max_dt(1) = cfl_target*(dx(j)/(abs(vel(1)) + a(1)))
        end if

        ! Viscous calculations
        if (chem_params%transport_dt) then
            ! Re from s_compute_transport_dt_re already sums the num_dims directions, so it keeps the min-width form
            hmin = dx(j)
            if (n > 0) hmin = min(hmin, dy(k))
            if (p > 0) then
                if (grid_geometry == 3) then
                    hmin = min(hmin, f_compute_filtered_dtheta(k, l))
                else
                    hmin = min(hmin, dz(l))
                end if
            end if
            max_dt(2) = cfl_target*hmin**2/maxval(1/(rho*Re_l))
        else if (viscous) then
            call s_compute_visc_diffusivity(Re_l, rho, nu)
            max_dt(2) = min(cfl_target, visc_dt_cap)/(nu*inv_dx2)
        end if

        ! Capillary CFL calculations
        if (surface_tension) then
            ccfl_dt = huge(1._wp)
            call s_compute_capillary_rho(alpha, alpha_rho, rho_c)
            if (rho_c > 0._wp) then
                if (p > 0) then
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (grid_geometry == 3) then
                            fltr_dtheta = f_compute_filtered_dtheta(k, l)
                            ccfl_dt = cfl_target*sqrt(rho_c*min(dx(j), dy(k), fltr_dtheta)**3._wp/(2._wp*pi*sigma))
                        else
                            ccfl_dt = cfl_target*sqrt(rho_c*min(dx(j), dy(k), dz(l))**3._wp/(2._wp*pi*sigma))
                        end if
                    #:endif
                else if (n > 0) then
                    ccfl_dt = cfl_target*sqrt(rho_c*min(dx(j), dy(k))**3._wp/(2._wp*pi*sigma))
                else
                    ccfl_dt = cfl_target*sqrt(rho_c*dx(j)**3._wp/(2._wp*pi*sigma))
                end if
            end if
            max_dt(3) = ccfl_dt
        end if

        ! Thermal diffusion CFL: dt <= cfl/(D sum(1/dx^2))
        if (Dth > 0._wp) max_dt(4) = cfl_target/(Dth*inv_dx2)

    end subroutine s_compute_dt_from_cfl

end module m_sim_helpers
