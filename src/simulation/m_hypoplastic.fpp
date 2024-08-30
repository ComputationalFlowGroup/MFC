!>
!! @file m_hypoplastic.f90
!! @brief Contains module m_hypoplastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoplastic model
module m_hypoplastic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_hypoplastic_module, &
 s_finalize_hypoplastic_module, &
 s_compute_hypoplastic_rhs

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), Gs)
    !$acc declare link(Gs)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), du_dx, du_dy)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dv_dx, dv_dy)
    !$acc declare link(du_dx,du_dy,dv_dx,dv_dy)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), rho_K_field, G_K_field)
    !$acc declare link(rho_K_field, G_K_field)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), allocatable, dimension(:, :), fd_coeff_x, fd_coeff_y)
    !$acc declare link(fd_coeff_x,fd_coeff_y)

#else
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: du_dx, du_dy
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dv_dx, dv_dy
    !$acc declare create(du_dx,du_dy,dv_dx,dv_dy)

    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_x, fd_coeff_y
    !$acc declare create(fd_coeff_x,fd_coeff_y)
#endif

contains
   !>   The following subroutine handles the hypoelastic evolution
        !! equation for Johnson-Cook plasticity model, which requires
        !! the Jaumann-Zaremba rate of the Kirchhoff stress deviator
    subroutine s_initialize_hypoplastic_module

        integer :: i, k, r

        @:ALLOCATE_GLOBAL(Gs(1:num_fluids))
        @:ALLOCATE_GLOBAL(du_dx(0:m,0:n,0:p))
        @:ALLOCATE_GLOBAL(du_dy(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p))

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        @:ALLOCATE_GLOBAL(fd_coeff_x(-fd_number:fd_number, 0:m))
        @:ALLOCATE_GLOBAL(fd_coeff_y(-fd_number:fd_number, 0:n))
   
        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                      fd_number, fd_order)
        !$acc update device(fd_coeff_x)
        call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                          fd_number, fd_order)
        !$acc update device(fd_coeff_y)
    
    end subroutine s_initialize_hypoplastic_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoplastic_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rho_K, G_K, wtensor
        real(kind(0d0)), dimension(num_dims*(num_dims + 1)/2) :: stensor, tensora, devdtensor, Dp

        integer :: i, k, l, p, r, q !< Loop variables

        real(kind(0d0)) :: energy, alf, dyn_p, pi_inf
        real(kind(0d0)) ::  gamma, rho, qv, pres, stress, mom, temp, G
        real(kind(0d0)), dimension(num_fluids) ::  alpha_K, alpha_rho_K
        real(kind(0d0)) :: theta_m, tempref, theta_hat, sigma_bar, dp_JC, d_p

        ! compute velocity gradients and rho_K and G_K        
        du_dx(:, :, :) = 0d0; du_dy(:, :, :) = 0d0
        dv_dx(:, :, :) = 0d0; dv_dy(:, :, :) = 0d0

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

        tensora(:) = 0d0
        !$acc parallel loop collapse(2) gang vector default(present) &
        !$acc private(rho_K,G_K,alpha_rho_K,alpha_K)
        do l = 0, n
          do k = 0, m
             ! STEP 1 : Compute the first additional term in rhs: -SW + WS
             ! Let wtensor = W12, tensora = -SW + WS
             wtensor = 5d-1*(du_dy(k, l, q) - dv_dx(k, l, q))
             stensor(1) = 2d0*q_prim_vf(strxe - 1)%sf(k, l, q) !2*S12
             stensor(2) = q_prim_vf(strxe)%sf(k, l, q) - q_prim_vf(strxb)%sf(k, l, q) !S22 - S11
             stensor(3) = -stensor(1) !-2*S12
             tensora(1) = wtensor*stensor(1)
             tensora(2) = wtensor*stensor(2)
             tensora(3) = wtensor*stensor(3)
           
             ! STEP 2: Compute the deviatoric part of D, symmetric part of velocity gradient
             ! dtrace = du_dx(k, l, q) + dv_dy(k, l, q)
             devdtensor(1) = du_dx(k, l, q) - (1d0/3d0)*(du_dx(k, l, q) + dv_dy(k, l, q))
             devdtensor(2) = 5d-1*(du_dy(k, l, q) + dv_dx(k, l, q))
             devdtensor(3) = dv_dy(k, l, q) - (1d0/3d0)*(du_dx(k, l, q) + dv_dy(k, l, q))
 
!             print *, 'I got here A' 

             ! STEP 3: Compute the equivalent plastic strain rate, d^p 
             ! STEP 3.1 : Compute mixtures variables for computing
             ! pressure and temperature
             energy = q_cons_vf(E_idx)%sf(k, l, q) 
             dyn_p = 0d0          
             do i = momxb, momxe
                dyn_p = dyn_p + 5d-1*q_cons_vf(i)%sf(k, l, q)*q_prim_vf(i)%sf(k, l, q)
             end do
!              print *, 'I got here B' 
             
             rho_K = 0d0; G_K = 0d0;
             ! STEP 3.2 : Compute mixtures in preparation for pressure and temperature
             do i = 1, num_fluids
                rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) 
                G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i) 
                alpha_rho_K(i) = q_prim_vf(i)%sf(k, l, q)
                alpha_K(i) = q_prim_vf(advxb + i - 1)%sf(k, l, q)
             end do
!              print *, 'I got here C' 

             ! STEP 3.3: TODO MIRELYS
             if (G_K .gt. verysmall) then 
!               print *, 'I got here D' 
              ! STEP 3.4 : Compute mixture pressure and temperature
!                print *, 'energy ::', energy, 'alf ::', alf, 'dyn_p ::',&
!dyn_p, 'pi_inf ::', pi_inf, 'gamma ::', gamma, 'rho ::', rho, 'qv ::', &
!qv, 'stress ::', stress, 'mom ::', mom, 'G ::', G, 'alpha_K ::', &
!alpha_K, 'alpha_rho_K ::', alpha_rho_K
                call s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, & 
                                        pres, stress, mom, G, alpha_K, alpha_rho_K)
                call s_compute_temperature(energy, dyn_p, pi_inf, gamma, rho, qv, & 
                                          temp, alpha_K, alpha_rho_K)
!                print *, 'pressure :: ', pres, 'temperature ::', temp
                ! STEP 3.5 : Compute theta_m, theta_hat, and sigma_bar
                ! compute theta_m from equation 4.10
                ! jcook(6) = theta_m0, jcook(8) = pres_init, jcook(9) = d, assuming presref = 0
                theta_m = jcook6(1)*(1d0 + (pres/jcook8(1)))**(1d0/jcook9(1))
                ! compute theta_hat from equation 4.9
                tempref = jcook11(1)
                if (temp .lt. tempref) then
                   theta_hat = 0
                elseif (temp .le. theta_m) then
                   theta_hat = (temp - tempref)/(theta_m - tempref)
                else
                   theta_hat = 1
                end if
!                print *, 'I got here E' 
             !could alternatively compute subtract tempref in both temp subroutine and theta_m
                ! compute sigma_bar = sqrt(3/2) * | S | 
                sigma_bar = dsqrt(1.5d0) * (q_prim_vf(strxb)%sf(k, l, q)**2d0 + & 
                            2d0*q_prim_vf(strxb + 1)%sf(k, l, q)**2d0 + q_prim_vf(strxe - 1)%sf(k, l, q)**2d0 + &
                            q_prim_vf(strxe)%sf(k, l, q)**2d0)**(5d-1)
!                print *, 'sigma_bar ::', sigma_bar

                ! STEP 3.6 : Compute d^p and update rhs
                ! compute d^p_JC from equation 4.7
                ! d0 = 1 s^-1, jcook(4) = C, jcook(1) = A, jcook(2) = B,
                ! jcook(10) = d0 = R_tilde nondimensionally
                 dp_JC = jcook10(1) * dexp( (1d0/jcook4(1)) * (sigma_bar / &
                     ((jcook1(1) + jcook2(1)*q_prim_vf(plasidx)%sf(k, l, q))**(jcook3(1)) &
                     *(1d0 - theta_hat**jcook5(1)))) - 1d0)
                ! compute d^p from equation 4.6
                ! jcook(7) = d^p_lim
                if (sigma_bar .gt. 1d0*verysmall) then
                    d_p = ((1d0/dp_JC) + (1d0/jcook7(1)))**(-1d0)
                    ! compute D^p using equation 4.5
                    Dp(1) = 1.5d0*(d_p / sigma_bar) * q_prim_vf(strxb)%sf(k, l, q)
                    Dp(2) = 1.5d0*(d_p / sigma_bar) * q_prim_vf(strxb + 1)%sf(k, l, q)
                    Dp(3) = 1.5d0*(d_p / sigma_bar) * q_prim_vf(strxe)%sf(k, l, q)
!                 print *, 'I got here F' 
               else
                    d_p = 0d0
                    Dp(:) = 0d0
!                  print *, 'I got here G' 
               end if

                ! STEP 4: Compute rhs source terms
                rhs_vf(strxb + 0)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K*tensora(1) + & 
                                                2d0*rho_K*G_K*(devdtensor(1) - Dp(1))
                      
                rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K*tensora(2) + &
                                                2d0*rho_K*G_K*(devdtensor(2) - Dp(2))
                                                     
                rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K*tensora(3) + &
                                                2d0*rho_K*G_K*(devdtensor(3) - Dp(3))
               
               ! STEP 5: Compute hardening rhs term
               rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q) + rho_K*d_p
             end if
           end do
         end do
         !$acc end parallel loop

    end subroutine s_compute_hypoplastic_rhs

    subroutine s_finalize_hypoplastic_module() ! --------------------

        @:DEALLOCATE_GLOBAL(Gs)
        @:DEALLOCATE_GLOBAL(du_dx)
        @:DEALLOCATE_GLOBAL(fd_coeff_x)
        @:DEALLOCATE_GLOBAL(du_dy,dv_dx,dv_dy)
        @:DEALLOCATE_GLOBAL(fd_coeff_y)

    end subroutine s_finalize_hypoplastic_module

end module m_hypoplastic
