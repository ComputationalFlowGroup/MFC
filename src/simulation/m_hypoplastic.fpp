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

    real(kind(0d0)), allocatable, dimension(:, :, :) :: rho_K_field, G_K_field
    !$acc declare create(rho_K_field, G_K_field)

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
        @:ALLOCATE_GLOBAL(rho_K_field(0:m,0:n,0:p), G_K_field(0:m,0:n,0:p))
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
        !!  @param idir Dimension splitting index
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoplastic_rhs(idir, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rho_K, G_K
        real(kind(0d0)), dimension(num_dims**2) :: atensor, tensora, devdtensor, Dp

        integer :: i, k, l, p, r, q !< Loop variables

        real(kind(0d0)), dimension(10) :: jcook
        real(kind(0d0)) :: energy, alf, dyn_p, pi_inf
        real(kind(0d0)) :: gamma, rho, qv, pres, stress, mom, temp, G
        real(kind(0d0)), dimension(num_fluids) ::  alpha_K, alpha_rho_K
        real(kind(0d0)) :: theta_m, tempref, theta_hat, sigma_bar, dp_JC, d_p


        ! compute velocity gradients and rho_K and G_K        
        !$acc parallel loop collapse(2) gang vector default(present)
        do l = 0, n
          do k = 0, m
            du_dx(k, l, q) = 0d0
            du_dy(k, l, q) = 0d0
            dv_dx(k, l, q) = 0d0
            dv_dy(k, l, q) = 0d0
            do r = -fd_number, fd_number
               du_dx(k, l, q) = du_dx(k, l, q) &
               + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
               du_dy(k, l, q) = du_dy(k, l, q) &
               + q_prim_vf(momxb)%sf(k, l + r, q)*fd_coeff_y(r, l)
               dv_dx(k, l, q) = dv_dx(k, l, q) &
               + q_prim_vf(momxb + 1)%sf(k + r, l, q)*fd_coeff_x(r, k)
               dv_dy(k, l, q) = dv_dy(k, l, q) &
               + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y(r, l)
            end do
          end do
        end do
        !$acc end parallel loop

        !$acc parallel loop collapse(2) gang vector default(present)
        do l = 0, n
           do k = 0, m
             rho_K = 0d0; G_K = 0d0
             do i = 1, num_fluids
                rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) !alpha_rho_K(1)
                G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)  !alpha_K(1) * Gs(1)
             end do
                rho_K_field(k, l, q) = rho_K
                G_K_field(k, l, q) = G_K
                !TODO: take this out if not needed
                if (G_K < verysmall) then
                    G_K_field(k, l, q) = 0
                end if
           end do
        end do
        !$acc end parallel loop


        tensora(1) = 0d0;  tensora(4) = 0d0
        !$acc parallel loop collapse(2) gang vector default(present)
        do l = 0, n
          do k = 0, m
             ! STEP 1 : Compute the first additional term in rhs: -rho((SW)^T - SW)
             ! Let atensor = SW, attensor = (SW)^T, tensora = attensor - atensor
             atensor(1) = (1d0/4d0)*(dv_dx(k, l, q)**2 - du_dy(k, l, q)**2)
             atensor(2) = (1d0/2d0)*(du_dy(k, l, q)*du_dx(k, l, q) - &
                  dv_dx(k, l, q)*du_dx(k, l, q))
             atensor(3) = (1d0/2d0)*(dv_dx(k, l, q)*dv_dy(k, l, q) - &
                  du_dy(k, l, q)*dv_dy(k, l, q))
             atensor(4) = (1d0/4d0)*(du_dy(k, l, q)**2 - dv_dx(k, l, q)**2)
             tensora(2) = atensor(3) - atensor(2)
             tensora(3) = atensor(2) - atensor(3)
            
             ! STEP 2: Compute the deviatoric part of D, symmetric part of velocity gradient
             ! dtrace = du_dx(k, l, q) + dv_dy(k, l, q)
             devdtensor(1) = du_dx(k, l, q) - (1d0/3d0)*(du_dx(k, l, q) + dv_dy(k, l, q))
             devdtensor(2) = (1d0/2d0)*(du_dy(k, l, q) + dv_dx(k, l, q))
             devdtensor(3) = devdtensor(2)
             devdtensor(4) = dv_dy(k, l, q) - (1d0/3d0)*(du_dx(k, l, q) + dv_dy(k, l, q))
            
             ! STEP 3: Compute the equivalent plastic strain rate, d^p 
             ! STEP 3.1 : Compute mixture pressure and temperature
             call s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, pres, stress, mom, G, alpha_K, alpha_rho_K)
             call s_compute_temperature(energy, dyn_p, pi_inf, gamma, rho, qv, temp, alpha_K, alpha_rho_K)

             ! STEP 3.2 : Compute theta_m, theta_hat, and sigma_bar
             ! compute theta_m from equation 4.10
	     ! jcook(6) = theta_m0, jcook(8) = pres_init, jcook(9) = d, assuming presref = 0
	     theta_m = jcook(6)*(1d0 + (pres/jcook(8)))**(1d0/jcook(9))
             ! compute theta_hat from equation 4.9
             tempref = 298 ! DO NOT DO: HARDCODED REFERENCE TEMPERATURE
             theta_hat = (temp - tempref)/(theta_m - tempref) 
             !could alternatively compute subtract tempref in both temp subroutine and theta_m
             ! compute sigma_bar = sqrt(3/2) * | S | 
             sigma_bar = sqrt(3d0/2d0) * (du_dx(k, l, q)*dv_dy(k, l, q) - &
                         (1d0/2d0)*du_dy(k, l, q)*dv_dx(k, l, q) - &
                         (1d0/40)*(du_dy(k, l, q)**2 * dv_dx(k, l, q)**2))
             ! STEP 3.3 : Compute d^p and update rhs
             ! compute d^p_JC from equation 4.7
             ! d0 = 1 s^-1, jcook(4) = C, jcook(1) = A, jcook(2) = B, 
             dp_JC = exp( (1d0/jcook(4)) * (sigma_bar / &
                    ((jcook(1) + jcook(2)*q_prim_vf(plasidx)%sf(k, l,q)) * &
                    (1d0 - theta_hat))) - 1d0)
             ! compute d^p from equation 4.6
             ! jcook(7) = d^p_lim
             d_p = ((1d0/dp_JC) + (1d0/jcook(7)))**(-1d0)
             ! compute D^p using equation 4.5
             Dp(1) = ((3d0*d_p) / (2d0*sigma_bar)) * du_dx(k, l, q)
             Dp(2) = ((3d0*d_p) / (2d0*sigma_bar)) * (1d0/2d0)*(du_dy(k, l, q) + dv_dx(k, l, q))
             Dp(3) = Dp(2)
             Dp(4) = ((3d0*d_p) / (2d0*sigma_bar)) * dv_dy(k, l, q)

             ! STEP 4: Compute rhs source terms
             rhs_vf(strxb + 0)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l ,q)*tensora(1) + & 
               2d0*rho_K_field(k, l, q)*G_K_field(k, l, q)*(devdtensor(1) - Dp(1))
                      
             rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* tensora(2) + &
               2d0*rho_K_field(k, l, q)*G_K_field(k, l, q)*(devdtensor(2) - Dp(2))
                                                     
             rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)*tensora(3) + &
               2d0*rho_K_field(k, l, q)*G_K_field(k, l, q)*(devdtensor(3) - Dp(3))
               
             rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + rho_K_field(k, l, q)*tensora(4) + &
               2d0*rho_K_field(k, l, q)*G_K_field(k, l, q)*(devdtensor(4) - Dp(4))             
             ! TODO: IS THIS RIGHT?
             rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q)*d_p
            end do
         end do
         !$acc end parallel loop

    end subroutine s_compute_hypoplastic_rhs

    subroutine s_finalize_hypoplastic_module() ! --------------------

        @:DEALLOCATE_GLOBAL(Gs)
        @:DEALLOCATE_GLOBAL(rho_K_field, G_K_field)
        @:DEALLOCATE_GLOBAL(du_dx)
        @:DEALLOCATE_GLOBAL(fd_coeff_x)
        @:DEALLOCATE_GLOBAL(du_dy,dv_dx,dv_dy)
        @:DEALLOCATE_GLOBAL(fd_coeff_y)

    end subroutine s_finalize_hypoplastic_module

end module m_hypoplastic
