!>
!! @file m_miegruneisen.f90
!! @brief Contains module m_miegruneisen

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for miegruneisen model
module m_miegruneisen

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_helper

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_miegruneisen_module, &
 s_finalize_miegruneisen_module, &
 s_compute_miegruneisen_rhs

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), Gs)
    !$acc declare link(Gs)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), du_dx)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dv_dy)
    !$acc declare link(du_dx,dv_dy)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), rho_K_field, G_K_field)
    !$acc declare link(rho_K_field, G_K_field)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), allocatable, dimension(:, :), fd_coeff_x, fd_coeff_y)
    !$acc declare link(fd_coeff_x,fd_coeff_y)

#else
    real(kind(0d0)), allocatable, dimension(:, :, :) :: du_dx
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dv_dy
    !$acc declare create(du_dx,dv_dy)

    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_x, fd_coeff_y
    !$acc declare create(fd_coeff_x,fd_coeff_y)
#endif

contains
   !>   The following subroutine handles the extra evolution
        !! equations required for handling Miegruneisen EoS
    subroutine s_initialize_miegruneisen_module

        integer :: i, k, r

        @:ALLOCATE_GLOBAL(du_dx(0:m,0:n,0:p))
        @:ALLOCATE_GLOBAL(fd_coeff_x(-fd_number:fd_number, 0:m))
        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                      fd_number, fd_order)

        !$acc update device(fd_coeff_x)
        if (num_dims /= 1) then
        @:ALLOCATE_GLOBAL(dv_dy(0:m,0:n,0:p))
        @:ALLOCATE_GLOBAL(fd_coeff_y(-fd_number:fd_number, 0:n))
        call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                          fd_number, fd_order)

        !$acc update device(fd_coeff_y)
        end if
    
    end subroutine s_initialize_miegruneisen_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_miegruneisen_rhs(q_prim_vf, q_cons_vf, rhs_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rho_K, G_K, mg_exp, rhs_mgidx2_mix, rhs_mgidx3_mix
        real(kind(0d0)) :: A_cv, rho0_mix, theta_E, gamma_inf, gamma0, phi_mix, dummy 

        integer :: i, k, l, p, r, q !< Loop variables

        real(kind(0d0)) :: energy, alf, dyn_p, pi_inf
        real(kind(0d0)) :: gamma, rho, qv, pres, mom, temp, G 
        real(kind(0d0)), dimension(num_fluids) ::  alpha_K, alpha_rho_K
        if (num_dims == 1) then
        ! Writing code for quasi-1D not true 1D
        du_dx(:, :, :) = 0d0; 
        !$acc parallel loop collapse(3) gang vector default(present)
        do q = 0, p
            do l = 0, n
                do k = 0, m
                    do r = -fd_number, fd_number
                        du_dx(k, l, q) = du_dx(k, l, q) + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

        !$acc parallel loop collapse(2) gang vector default(present) &
        !$acc private(rho_K,G_K,alpha_rho_K,alpha_K,&
        !$acc mg_exp,rhs_mgidx2_mix, rhs_mgidx3_mix,&
        !$acc A_cv, rho0_mix, theta_E, gamma_inf, gamma0)
        l=0
        q=0
        do k = 0, m
             ! STEP 3.1 : Compute mixtures variables for computing
             ! pressure and temperature
             energy = q_cons_vf(E_idx)%sf(k, l, q) 
             dyn_p  = 0d0          
             do i = momxb, momxe
                dyn_p = dyn_p + 5d-1*q_cons_vf(i)%sf(k, l, q)*q_prim_vf(i)%sf(k, l, q)
             end do
             !print *, 'I got here B' 
             
             rho_K          = 0d0
             mg_exp         = 0d0
             rhs_mgidx2_mix = 0d0
             rhs_mgidx3_mix = 0d0
             A_cv           = 0d0
             rho0_mix       = 0d0
             theta_E        = 0d0
             gamma_inf      = 0d0
             gamma0         = 0d0
             ! STEP 3.2 : Compute mixtures in preparation for pressure and temperature
             do i = 1, num_fluids
                rho_K          = rho_K    + q_prim_vf(i)%sf(k, l, q) 
                alpha_rho_K(i) = q_prim_vf(i)%sf(k, l, q)
                alpha_K(i)     = q_prim_vf(advxb + i - 1)%sf(k, l, q)
                mg_exp         = mg_exp   + alpha_K(i)*mg_b(i)
                rho0_mix       = rho0_mix + alpha_K(i)*rho0(i)
                A_cv           = A_cv     + alpha_K(i)*ein_cv1(i)
                theta_E        = theta_E  + alpha_K(i)*ein_cv2(i)
                gamma_inf      = gamma_inf+ alpha_K(i)*mg_a(i)
                gamma0         = gamma0   + alpha_K(i)*gammas(i)
             end do
!              print *, 'I got here C'
             !phi_mix = ((rho0_mix/rho_K)**(-gamma_inf))*&
             !               dexp((gamma0 - gamma_inf)*&
             !               (1d0- (rho0_mix/rho_K)**mg_exp))
             
             do i = 1, num_fluids
                !rhs_mgidx2_mix = rhs_mgidx2_mix +&
                !                 pi_infs(i)*alpha_K(i)*rho_K/rho0(i) +&
                !                 dlog(rho_K/rho0_mix)*alpha_K(i)&
                !                 *pi_infs(i)*(qvs(i)-2d0)*rho_K/rho0(i)
                dummy = alpha_K(i)/alpha_rho_K(i)
               phi_mix = &
               (alpha_K(i)*rho0(i)/alpha_rho_K(i)**(-mg_a(i)))*&
                   dexp((gammas(i)-mg_a(i))*(1d0- &
               (alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i)))
                !rhs_mgidx2_mix = rhs_mgidx2_mix +&
                !                 pi_infs(i)*alpha_rho_K(i)/rho0(i) +&
                !                 dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i)))*alpha_rho_K(i)&
                !                 *pi_infs(i)*(qvs(i)-2d0)/rho0(i)
                rhs_mgidx2_mix = rhs_mgidx2_mix + &
                    ((qvs(i)*dummy)**2d0/(1d0/rho0(i)-1.51d0*(1d0/rho0(i)-dummy))**2d0+&
                    2d0*((qvs(i)*dummy)**2d0)*1.51*(1d0/rho0(i)-dummy)/&
                    (1d0/rho0(i)-1.51d0*(1d0/rho0(i)-dummy))**3d0)
                                  
                !rhs_mgidx3_mix = rhs_mgidx3_mix +&
                !                 (1d0/q_prim_vf(mgidxb)%sf(k, l, q))*alpha_rho_K(i)*A_cv*((theta_E*phi_mix)**2d0)&
                !                 *dexp(theta_E*phi_mix)/((dexp(theta_E*phi_mix)-1d0)**2d0)                 
                rhs_mgidx3_mix = rhs_mgidx3_mix +&
                                 (1d0/q_prim_vf(mgidxb)%sf(k, l, q))*alpha_rho_K(i)*ein_cv1(i)*((ein_cv2(i)*phi_mix)**2d0)&
                                 *dexp(ein_cv2(i)*phi_mix)/((dexp(ein_cv2(i)*phi_mix)-1d0)**2d0)                 
             end do
             
             rhs_mgidx2_mix = rhs_mgidx2_mix*q_prim_vf(mgidxb)%sf(k, l, q)
             ! STEP 3.3: TODO MIRELYS
!              if (G_K .gt. verysmall) then 
!               print *, 'I got here D' 
              ! STEP 3.4 : Compute mixture pressure and temperature
!                print *, 'energy ::', energy, 'alf ::', alf, 'dyn_p ::',&
!                dyn_p, 'pi_inf ::', pi_inf, 'gamma ::', gamma, 'rho ::', rho, 'qv ::', &
!                qv, 'stress ::', stress, 'mom ::', mom, 'G ::', G, 'alpha_K ::', &
!                alpha_K, 'alpha_rho_K ::', alpha_rho_K
               ! print *,q_cons_vf(mgidxe)%sf(k, l, q)
               ! if (q_cons_vf(mgidxe)%sf(k, l, q) /= q_cons_vf(mgidxe)%sf(k,l,q)) then
               !     print *,'rho_eref is NaN in miegruneisen'
               ! end if
                call s_compute_pressure(energy, 0d0, dyn_p, pi_inf, q_cons_vf(mgidxb)%sf(k, l, q), rho_K, 0d0, & 
                                        pres, 0d0, 0d0, G, &
                                        q_cons_vf(mgidxb+1)%sf(k, l, q), &
                                        q_cons_vf(mgidxe)%sf(k, l, q))
!               print *, 'pressure :: ', pres, 'temperature ::', temp
!               rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q) + rho_K*d_p
!             Compute the three rhs for the mie-gruneisen eos as written
!             in the overleaf
                rhs_vf(mgidxb)%sf(k, l, q)   = rhs_vf(mgidxb)%sf(k, l, q)&
                                            + (q_prim_vf(mgidxb)%sf(k, l, q)*&
                                            (alpha_K(1)*(1d0-mg_b(1))+alpha_K(2)*(1d0-mg_b(2))))*du_dx(k, l, q)
                rhs_vf(mgidxb+1)%sf(k, l, q) = rhs_vf(mgidxb+1)%sf(k, l, q)+ &
                                            (q_cons_vf(mgidxb+1)%sf(k, l, q)*&
                                            (alpha_K(2)*(1d0-mg_b(2))+(alpha_K(1)*(1d0-mg_b(1)))) -&
                                            q_prim_vf(mgidxb)%sf(k,l,q)*alpha_rho_K(1)*rhs_mgidx2_mix)&
                                            *du_dx(k, l, q)
                !if (rhs_vf(mgidxe)%sf(k,l,q) /= rhs_vf(mgidxe)%sf(k,l,q)) then
                !  print *,'k',k, rhs_vf(mgidxe)%sf(k, l, q)      
                !end if
                rhs_vf(mgidxe)%sf(k, l, q) = rhs_vf(mgidxe)%sf(k, l, q)&
                                            -0.5d0*(q_prim_vf(mgidxb+1)%sf(k,l,q)+&
                                            alpha_K(1)*rhs_mgidx2_mix*(((alpha_rho_K(1)/alpha_K(1))**2d0)/rho0(1)&
                                        - alpha_rho_K(1)/alpha_K(1)))*du_dx(k, l, q)
           end do
         !$acc end parallel loop
        else if (num_dims == 2) then 
        ! compute velocity gradients and rho_K and G_K        
        du_dx(:, :, :) = 0d0; dv_dy(:, :, :) = 0d0

        !$acc parallel loop collapse(2) gang vector default(present)
        do l = 0, n
          do k = 0, m
            do r = -fd_number, fd_number
               du_dx(k, l, q) = du_dx(k, l, q) + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
               dv_dy(k, l, q) = dv_dy(k, l, q) + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y(r, l)
            end do
          end do
        end do
        !$acc end parallel loop
        q = 0
        !$acc parallel loop collapse(2) gang vector default(present) &
        !$acc private(rho_K,G_K,alpha_rho_K,alpha_K,&
        !$acc mg_exp,rhs_mgidx2_mix, rhs_mgidx3_mix,&
        !$acc A_cv, rho0_mix, theta_E, gamma_inf, gamma0)
        do l = 0, n
          do k = 0, m
             ! STEP 3.1 : Compute mixtures variables for computing
             ! pressure and temperature
             energy = q_cons_vf(E_idx)%sf(k, l, q) 
             dyn_p  = 0d0          
             do i = momxb, momxe
                dyn_p = dyn_p + 5d-1*q_cons_vf(i)%sf(k, l, q)*q_prim_vf(i)%sf(k, l, q)
             end do
!             print *, 'I got here B' 
             
             rho_K          = 0d0
             mg_exp         = 0d0
             rhs_mgidx2_mix = 0d0
             rhs_mgidx3_mix = 0d0
             A_cv           = 0d0
             rho0_mix       = 0d0
             theta_E        = 0d0
             gamma_inf      = 0d0
             gamma0         = 0d0
             ! STEP 3.2 : Compute mixtures in preparation for pressure and temperature
             do i = 1, num_fluids
                rho_K          = rho_K    + q_prim_vf(i)%sf(k, l, q) 
                alpha_rho_K(i) = q_prim_vf(i)%sf(k, l, q)
                alpha_K(i)     = q_prim_vf(advxb + i - 1)%sf(k, l, q)
                mg_exp         = mg_exp   + alpha_K(i)*mg_b(i)
                rho0_mix       = rho0_mix + alpha_K(i)*rho0(i)
                A_cv           = A_cv     + alpha_K(i)*ein_cv1(i)
                theta_E        = theta_E  + alpha_K(i)*ein_cv2(i)
                gamma_inf      = gamma_inf+ alpha_K(i)*mg_a(i)
                gamma0         = gamma0   + alpha_K(i)*gammas(i)
             end do
!              print *, 'I got here C'

                phi_mix = ((rho0_mix/rho_K)**(-gamma_inf))*&
                            dexp((gamma0 - gamma_inf)*&
                            (1d0- (rho0_mix/rho_K)**mg_exp))
 
             do i = 1, num_fluids
                rhs_mgidx2_mix = rhs_mgidx2_mix +&
                                 pi_infs(i)*alpha_K(i)*rho_K/rho0(i) +&
                                 dlog(rho_K/rho0_mix)*alpha_K(i)&
                                 *pi_infs(i)*(qvs(i)-2d0)*rho_K/rho0(i)
               
                rhs_mgidx3_mix = rhs_mgidx3_mix +&
                                 (1d0/q_prim_vf(mgidxb)%sf(k, l, q))*alpha_rho_K(i)*A_cv*((theta_E*phi_mix)**2d0)&
                                 *dexp(theta_E*phi_mix)/((dexp(theta_E*phi_mix)-1d0)**2d0)
                                 
             end do
             rhs_mgidx2_mix = rhs_mgidx2_mix*q_prim_vf(mgidxb)%sf(k, l, q)
             ! STEP 3.3: TODO MIRELYS
!              if (G_K .gt. verysmall) then 
!               print *, 'I got here D' 
              ! STEP 3.4 : Compute mixture pressure and temperature
!                print *, 'energy ::', energy, 'alf ::', alf, 'dyn_p ::',&
!                dyn_p, 'pi_inf ::', pi_inf, 'gamma ::', gamma, 'rho ::', rho, 'qv ::', &
!                qv, 'stress ::', stress, 'mom ::', mom, 'G ::', G, 'alpha_K ::', &
!                alpha_K, 'alpha_rho_K ::', alpha_rho_K
                call s_compute_pressure(energy, alf, dyn_p, pi_inf, q_cons_vf(mgidxb)%sf(k, l, q), rho_K, 0d0, & 
                                        pres, 0d0, 0d0, G, &
                                        q_cons_vf(mgidxb+1)%sf(k, l, q), &
                                        q_cons_vf(mgidxe)%sf(k, l, q))
!               print *, 'pressure :: ', pres, 'temperature ::', temp
!               rhs_vf(plasidx)%sf(k, l, q) = rhs_vf(plasidx)%sf(k, l, q) + rho_K*d_p
!             Compute the three rhs for the mie-gruneisen eos as written
!             in the overleaf
                rhs_vf(mgidxb)%sf(k, l, q)   = rhs_vf(mgidxb)%sf(k, l, q)&
                                            + q_prim_vf(mgidxb)%sf(k, l, q)*&
                                            (1d0- mg_exp)*(du_dx(k, l, q)+ dv_dy(k, l, q))
                rhs_vf(mgidxb+1)%sf(k, l, q) = rhs_vf(mgidxb+1)%sf(k, l, q) &
                                            -(q_cons_vf(mgidxb+1)%sf(k, l, q)*mg_exp + rhs_mgidx2_mix)&
                                            *(du_dx(k, l, q)+dv_dy(k, l, q))
                rhs_vf(mgidxe)%sf(k, l, q)   = rhs_vf(mgidxe)%sf(k, l, q)&
                                            +(-q_prim_vf(mgidxb+1)%sf(k, l, q)+&
                                            rhs_mgidx3_mix)*(du_dx(k, l, q)+ dv_dy(k, l, q))
                                            
!             end if
           end do
         end do
         !$acc end parallel loop
     end if

    end subroutine s_compute_miegruneisen_rhs

    subroutine s_finalize_miegruneisen_module() ! --------------------

        @:DEALLOCATE_GLOBAL(du_dx)
        @:DEALLOCATE_GLOBAL(fd_coeff_x)
        if (num_dims /=1) then
            @:DEALLOCATE_GLOBAL(dv_dy)
            @:DEALLOCATE_GLOBAL(fd_coeff_y)
        end if

    end subroutine s_finalize_miegruneisen_module

end module m_miegruneisen
