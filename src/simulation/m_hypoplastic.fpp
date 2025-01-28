!>
!! @file m_hypoplastic.f90
!! @brief Contains module m_hypoplastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoplastic model
module m_hypoplastic

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_finite_differences

    implicit none

    private; public :: s_initialize_hypoplastic_module, &
 s_finalize_hypoplastic_module, &
 s_compute_hypoplastic_rhs

    real(wp), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(wp), allocatable, dimension(:, :, :) :: du_dx, du_dy
    real(wp), allocatable, dimension(:, :, :) :: dv_dx, dv_dy
    !$acc declare create(du_dx,du_dy,dv_dx,dv_dy)

    real(wp), allocatable, dimension(:, :) :: fd_coeff_x, fd_coeff_y, fd_coeff_z
    !$acc declare create(fd_coeff_x,fd_coeff_y,fd_coeff_z)

contains
    !>   The following subroutine handles the hypoelastic evolution
        !! equation for Johnson-Cook plasticity model, which requires
        !! the Jaumann-Zaremba rate of the Kirchhoff stress deviator
    subroutine s_initialize_hypoplastic_module

        integer :: i, k, r

        @:ALLOCATE(Gs(1:num_fluids))
        @:ALLOCATE(du_dx(0:m,0:n,0:p))
        @:ALLOCATE(du_dy(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p))

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        @:ALLOCATE(fd_coeff_x(-fd_number:fd_number, 0:m))
        if (n > 0) then
            @:ALLOCATE(fd_coeff_y(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then
            @:ALLOCATE(fd_coeff_z(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                      fd_number, fd_order)
        !$acc update device(fd_coeff_x)
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_y)
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_z)
        end if

    end subroutine s_initialize_hypoplastic_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoplastic_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(wp) :: rho_K, G_K, wtensor
        real(wp), dimension(2*num_dims) :: stensor, tensora, devdtensor, Dp

        integer :: i, k, l, p, r, q !< Loop variables

        real(wp) :: energy, alf, dyn_p, pi_inf
        real(wp) :: gamma, rho, pres, stress, mom, temp, G
        real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(wp) :: theta_m, tempref, theta_hat, sigma_bar, &
                    dp_JC, d_p, equiv_tens_stress

        real(wp), dimension(num_species) :: rhoYks

        du_dx(:, :, :) = 0._wp
        du_dy(:, :, :) = 0._wp
        dv_dx(:, :, :) = 0._wp
        dv_dy(:, :, :) = 0._wp

        if (num_dims == 1) then
            ! For quasi-1D case
            l = 0
            q = 0
            !$acc parallel loop collapse(1) gang vector default (present)
            do k = 0, m
                do r = -fd_number, fd_number
                    du_dx(k, l, q) = du_dx(k, l, q) + &
                                     q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
                end do
            end do
            !$acc end parallel loop
            do k = 0, m
                energy = q_cons_vf(E_idx)%sf(k, l, q)
                dyn_p = 0._wp
                do i = momxb, momxe
                    dyn_p = dyn_p + 5d-1*q_cons_vf(i)%sf(k, l, q)*q_prim_vf(i)%sf(k, l, q)
                end do

                rho_K = 0._wp; G_K = 0._wp; 
                ! STEP 3.2 : Compute mixtures in preparation for pressure and temperature
                do i = 1, num_fluids
                    rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)
                    G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)
                    alpha_rho_K(i) = q_prim_vf(i)%sf(k, l, q)
                    alpha_K(i) = q_prim_vf(advxb + i - 1)%sf(k, l, q)
                end do

                ! STEP 3.3: TODO MIRELYS
                if (G_K > sgm_eps) then
                    !STEP 3.4 : Compute mixture pressure and temperature
                    call s_compute_pressure(energy, 0._wp, dyn_p, pi_inf, 0._wp, rho, 0._wp, &
                                            rhoYks, pres, temp, 0._wp, 0._wp, 0._wp, alpha_K, alpha_rho_K)

                    call s_compute_temperature(energy, dyn_p, temp, alpha_K, alpha_rho_K)
                    ! STEP 3.5 : Compute theta_m, theta_hat, and sigma_bar
                    ! compute theta_m from equation 4.10
                    ! jcook(6) = theta_m0, jcook(8) = pres_init, jcook(9) = d, assuming presref = 0
                    theta_m = jcook6(1)*(1_wp + (pres/jcook8(1)))**(1_wp/jcook9(1))
                    ! compute theta_hat from equation 4.9
                    tempref = jcook11(1)
                    !This line is here because temperature subroutine gives
                    !me the increase is temperature from reference
                    !temperature and not the absolute temperature
                    temp = temp + tempref

                    if (temp < tempref) then
                        theta_hat = 0._wp
                    elseif (temp <= theta_m) then
                        theta_hat = (temp - tempref)/(theta_m - tempref)
                    else
                        theta_hat = 1._wp + verysmall
                    end if

                    !could alternatively compute subtract tempref in both temp subroutine and theta_m
                    ! compute sigma_bar = sqrt(3/2) * | S |
                    sigma_bar = sqrt(1.5_wp)*abs(q_prim_vf(strxb)%sf(k, l, q))

                    ! STEP 3.6 : Compute d^p and update rhs
                    ! compute d^p_JC from equation 4.7
                    ! _wp = 1 s^-1, jcook(4) = C, jcook(1) = A, jcook(2) = B,
                    ! jcook(10) = _wp = R_tilde nondimensionally
                    equiv_tens_stress = (jcook1(1) + &
                                         jcook2(1)*q_prim_vf(plasidx)%sf(k, l, q)**jcook3(1))*(1_wp - theta_hat**jcook5(1))
                    dp_JC = jcook10(1)*exp((1_wp/jcook4(1))*(sigma_bar/equiv_tens_stress - 1_wp))
                    !jcook2(1)*q_prim_vf(plasidx)%sf(k, l, q)**jcook3(1),&
                    !(1_wp - theta_hat**jcook5(1))
                    !if (dp_JC .gt. sgm_eps) then
                    ! print *, dp_JC
                    !end if
                    ! compute d^p from equation 4.6
                    ! jcook(7) = d^p_lim
                    !if (sigma_bar .gt. sgm_eps) then
                    d_p = ((1_wp/dp_JC) + (1_wp/jcook7(1)))**(-1_wp)
                    ! compute D^p using equation 4.5
                    do i = strxb, strxe
                        Dp(i - strxb + 1) = 1.5_wp*(d_p/sigma_bar)*q_prim_vf(i)%sf(k, l, q)
                    end do
                    !print *, 'I got here F'
                    !else
                    !    d_p   = 0._wp
                    !    Dp(:) = 0._wp
!                 print *, 'I got here G'
                    !end if

                    ! STEP 4: Compute rhs source terms
                    devdtensor(1) = 0.5_wp*du_dx(k, l, q)

                    do i = strxb, strxe
                        rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) + &
                                                2_wp*rho_K*G_K*(devdtensor(i - strxb + 1) - Dp(i - strxb + 1))
                    end do
                    if (d_p > sgm_eps) then
                        print *, d_p
                    end if
                    !print *, rho_K
                    ! STEP 5: Compute hardening rhs term
                    rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q) + rho_K*d_p
                end if
            end do

        else if (num_dims == 2) then
            ! compute velocity gradients and rho_K and G_K
            q = 0
            !$acc parallel loop collapse(2) gang vector default(present)
            do l = 0, n
                do k = 0, m
                    do r = -fd_number, fd_number
                        du_dx(k, l, q) = du_dx(k, l, q) + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
                        du_dy(k, l, q) = du_dy(k, l, q) + q_prim_vf(momxb)%sf(k, l + r, q)*fd_coeff_y(r, l)
                        dv_dx(k, l, q) = dv_dx(k, l, q) + q_prim_vf(momxb + 1)%sf(k + r, l, q)*fd_coeff_x(r, k)
                        dv_dy(k, l, q) = dv_dy(k, l, q) + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y(r, l)
                    end do
                end do
            end do
            !$acc end parallel loop

            tensora(:) = 0._wp
            stensor(:) = 0._wp
            !$acc parallel loop collapse(2) gang vector default(present) &
            !$acc private(rho_K,G_K,alpha_rho_K,alpha_K,Dp)
            do l = 0, n
                do k = 0, m
                    ! STEP 1 : Compute the first additional term in rhs: -SW + WS
                    ! Let wtensor = W12, tensora = -SW + WS
                    wtensor = 0.5_wp*(du_dy(k, l, q) - dv_dx(k, l, q))
                    stensor(1) = 2._wp*q_prim_vf(strxb + 1)%sf(k, l, q) !2*S12
                    stensor(2) = q_prim_vf(strxb + 2)%sf(k, l, q) - q_prim_vf(strxb)%sf(k, l, q) !S22 - S11
                    stensor(3) = -stensor(1) !-2*S12
                    tensora(1) = wtensor*stensor(1)
                    tensora(2) = wtensor*stensor(2)
                    tensora(3) = wtensor*stensor(3)
                    tensora(4) = 0._wp

                    ! STEP 2: Compute the deviatoric part of D, symmetric part of velocity gradient
                    ! dtrace = du_dx(k, l, q) + dv_dy(k, l, q)
                    devdtensor(1) = du_dx(k, l, q) - (1._wp/3._wp)*(du_dx(k, l, q) + dv_dy(k, l, q))
                    devdtensor(2) = 0.5_wp*(du_dy(k, l, q) + dv_dx(k, l, q))
                    devdtensor(3) = dv_dy(k, l, q) - (1._wp/3._wp)*(du_dx(k, l, q) + dv_dy(k, l, q))
                    devdtensor(4) = -(1._wp/3._wp)*(du_dx(k, l, q) + dv_dy(k, l, q))

                    ! STEP 3: Compute the equivalent plastic strain rate, d^p
                    ! STEP 3.1 : Compute mixtures variables for computing
                    ! pressure and temperature
                    energy = q_cons_vf(E_idx)%sf(k, l, q)
                    dyn_p = 0._wp
                    do i = momxb, momxe
                        dyn_p = dyn_p + 0.5_wp*q_cons_vf(i)%sf(k, l, q)*q_prim_vf(i)%sf(k, l, q)
                    end do

                    rho_K = 0._wp
                    G_K = 0._wp
                    ! STEP 3.2 : Compute mixtures in preparation for pressure and temperature
                    do i = 1, num_fluids
                        rho_K = rho_K + q_prim_vf(i)%sf(k, l, q)
                        G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)
                        alpha_rho_K(i) = q_prim_vf(i)%sf(k, l, q)
                        alpha_K(i) = q_prim_vf(advxb + i - 1)%sf(k, l, q)
                    end do

                    call s_compute_pressure(energy, 0._wp, dyn_p, pi_inf, 0._wp, rho, 0._wp, &
                                            rhoYks, pres, temp, 0._wp, 0._wp, 0._wp, alpha_K, alpha_rho_K)

                    call s_compute_temperature(energy, dyn_p, temp, alpha_K, alpha_rho_K)

                    ! STEP 3.5 : Compute theta_m, theta_hat, and sigma_bar
                    ! compute theta_m from equation 4.10
                    ! jcook(6) = theta_m0, jcook(8) = pres_init, jcook(9) = d, assuming presref = 0
                    theta_m = jcook6(1)*(1._wp + (pres/jcook8(1)))**(1._wp/jcook9(1))
                    ! compute theta_hat from equation 4.9
                    tempref = jcook11(1)
                    if (temp < tempref) then
                        theta_hat = 0
                    elseif (temp <= theta_m) then
                        theta_hat = (temp - tempref)/(theta_m - tempref)
                    else
                        theta_hat = 1
                    end if
                    !could alternatively compute subtract tempref in both temp subroutine and theta_m
                    ! compute sigma_bar = sqrt(3/2) * | S |
                    sigma_bar = sqrt(1.5_wp)*(q_prim_vf(strxb)%sf(k, l, q)**2._wp + &
                                              2._wp*q_prim_vf(strxb + 1)%sf(k, l, q)**2._wp + q_prim_vf(strxe - 1)%sf(k, l, q)**2._wp + &
                                              q_prim_vf(strxe)%sf(k, l, q)**2._wp)**(0.5_wp)

                    ! STEP 3.6 : Compute d^p and update rhs
                    ! compute d^p_JC from equation 4.7
                    ! _wp = 1 s^-1, jcook(4) = C, jcook(1) = A, jcook(2) = B,
                    ! jcook(10) = _wp = R_tilde nondimensionally
                    dp_JC = jcook10(1)*exp((1._wp/jcook4(1))*(sigma_bar/ &
                                                              ((jcook1(1) + jcook2(1)*q_prim_vf(plasidx)%sf(k, l, q)**jcook3(1)) &
                                                               *(1._wp - theta_hat**jcook5(1))) - 1._wp))
                    ! compute d^p from equation 4.6
                    ! jcook(7) = d^p_lim
                    if (sigma_bar > 1.e-16_wp) then
                        d_p = ((1._wp/dp_JC) + (jcook10(1)/jcook7(1)))**(-1._wp)
                        ! compute D^p using equation 4.5
                        do i = strxb, strxe
                            Dp(i - strxb + 1) = 1.5_wp*(d_p/sigma_bar)*q_prim_vf(i)%sf(k, l, q)
                        end do
                    else
                        d_p = 0._wp
                        Dp(:) = 0._wp
                    end if

                    ! STEP 4: Compute rhs source terms
                    do i = strxb, strxe
                        rhs_vf(i)%sf(k, l, q) = rhs_vf(i)%sf(k, l, q) + rho_K*tensora(i - strxb + 1) + &
                                                2._wp*rho_K*G_K*(devdtensor(i - strxb + 1) - Dp(i - strxb + 1))
                    end do

                    ! STEP 5: Compute hardening rhs term
                    rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q) + rho_K*d_p
                end do
            end do
            !$acc end parallel loop
        end if
    end subroutine s_compute_hypoplastic_rhs

    subroutine s_finalize_hypoplastic_module()

        @:DEALLOCATE(Gs)
        @:DEALLOCATE(du_dx)
        @:DEALLOCATE(fd_coeff_x)
        @:DEALLOCATE(du_dy,dv_dx,dv_dy)
        @:DEALLOCATE(fd_coeff_y)

    end subroutine s_finalize_hypoplastic_module

end module m_hypoplastic
