#:def s_compute_speed_of_sound()
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c, alpha_rho_K, pref)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_compute_speed_of_sound
#else
        !$acc routine seq
#endif
        real(kind(0d0)), intent(IN) :: pres
        real(kind(0d0)), intent(IN) :: rho, gamma, pi_inf
        real(kind(0d0)), intent(IN) :: H
        real(kind(0d0)), dimension(num_fluids), intent(IN) :: adv
        real(kind(0d0)), intent(IN) :: vel_sum
        real(kind(0d0)), optional, dimension(num_fluids), intent(IN) :: alpha_rho_K
        real(kind(0d0)), optional, intent(IN) :: pref
        real(kind(0d0)), intent(OUT) :: c
        real(kind(0d0)) :: blkmod1, blkmod2
       

        !Local variables used for computation only
        real(kind(0d0)) :: rho0_mix, gamma_inf, gamma0, mg_exp, A_cv,&
                            theta_E, logrho, phi_mix, gamma_inv, pref1, dummy
        integer :: q, r

        if (alt_soundspeed) then
            blkmod1 = ((gammas(1) + 1d0)*pres + &
                       pi_infs(1))/gammas(1)
            blkmod2 = ((gammas(2) + 1d0)*pres + &
                       pi_infs(2))/gammas(2)
            c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))

        elseif (model_eqns == 3) then
            c = 0d0
            !$acc loop seq
            do q = 1, num_fluids
                c = c + adv(q)*(1d0/gammas(q) + 1d0)* &
                    (pres + pi_infs(q)/(gammas(q) + 1d0))
            end do
            c = c/rho

        elseif (((model_eqns == 4) .or. (model_eqns == 2 .and. bubbles))) then
            ! Sound speed for bubble mmixture to order O(\alpha)

            if (mpp_lim .and. (num_fluids > 1)) then
                c = (1d0/gamma + 1d0)* &
                    (pres + pi_inf/(gamma + 1d0))/rho
            else
                c = &
                    (1d0/gamma + 1d0)* &
                    (pres + pi_inf/(gamma + 1d0))/ &
                    (rho*(1d0 - adv(num_fluids)))
            end if

        elseif (model_eqns == 5) then
            !Note that pref and gamma are primitive state
            !variables for Mie-Gruneisen EoS and gamma = (1/Gamma(rho)) 
            rho0_mix      = 0d0
            gamma_inf     = 0d0
            gamma0        = 0d0
            mg_exp        = 0d0
            A_cv          = 0d0
            theta_E       = 0d0
            c             = 0d0
            do q = 1, num_fluids
               rho0_mix = rho0_mix   + adv(q)*rho0(q)
               gamma_inf= gamma_inf  + adv(q)*mg_a(q)
               gamma0   = gamma0     + adv(q)*gammas(q)
               mg_exp   = mg_exp     + adv(q)*mg_b(q)
               A_cv     = A_cv       + adv(q)*ein_cv1(q)
               theta_E  = theta_E    + adv(q)*ein_cv2(q)
            end do
            phi_mix = ((rho0_mix/rho)**(-gamma_inf))*&
                        dexp((gamma0 - gamma_inf)*&
                        (1d0- (rho0_mix/rho)**mg_exp))
            
            do q = 1, num_fluids
               phi_mix = &
               ((adv(q)*rho0(q)/alpha_rho_K(q))**(-mg_a(q)))*&
               dexp((gammas(q)-mg_a(q))*(1.0d0-(adv(q)*rho0(q)/alpha_rho_K(q))**mg_b(q)))
               
               dummy = adv(q)/alpha_rho_K(q) 
               pref1 = pi_infs(q) + &
               qvs(q)*(1d0/rho0(q)-dummy)/(1d0/rho0(q)-1.51d0*(1d0/rho0(q)-dummy))**2d0
               !pi_infs(q)*dlog(alpha_rho_K(q)/(adv(q)*rho0(q)))*(1d0+0.5d0*(qvs(q)-2d0)*dlog(alpha_rho_K(q)/(adv(q)*rho0(q))))

               gamma_inv = &
               1d0/(mg_a(q)+(gammas(q)-mg_a(q))*(dummy*rho0(q))**mg_b(q)) 
               c = c + &
               (alpha_rho_K(q)/(alpha_rho_K(1)+alpha_rho_K(2)))*(dummy*pres*(gamma_inv + &
                    1d0-mg_b(q)*gamma_inv) + & !(mg_b(q)-1d0)*pref1*dummy*gamma_inv + &
                    gamma_inv*((qvs(q)/rho0(q))**2d0/(1d0/rho0(q)-1.51d0*(1d0/rho0(q)-dummy))**2d0+&
                    2d0*((qvs(q)/rho0(q))**2d0)*1.51*(1d0/rho0(q)-dummy)/&
                    (1d0/rho0(q)-1.51d0*(1d0/rho0(q)-dummy))**3d0)-&
                    (pref1*dummy)) !+&
                   ! (1d0/gamma_inv)*ein_cv1(q)*((phi_mix*ein_cv2(q))**2d0)*dexp(phi_mix*ein_cv2(q))/&
                   ! (dexp(phi_mix*ein_cv2(q))-1.0d0)**2d0)
!             c = c + pres*gamma*adv(q)*(1d0-mg_b(q))&
!                    + pref*gamma*adv(q)*mg_b(q)&
!                    + gamma*alpha_rho_K(q)*pi_infs(q)/rho0(q)&
!                    + gamma*dlog(rho/rho0_mix)*alpha_rho_K(q)*pi_infs(q)*(qvs(q)-2d0)/rho0(q)&
!                    + (1d0/gamma)*alpha_rho_K(q)*A_cv*((theta_E*phi_mix)**2d0)&
!                    *dexp(theta_E*phi_mix)/((dexp(theta_E*phi_mix)-1d0)**2d0)
               
            end do 
!            c = c + pres - pref
!            c = c/(rho*gamma)
            c = c/gamma 
           ! if (c<0d0) then
           !  print *,&
           !  'c',c,'rho1',alpha_rho_K(1)/adv(1),'rho2',alpha_rho_K(2)/adv(2),'sum',adv(1)+adv(2) 
           ! end if
!            print *,c
           ! if (c /=c) then
           !     print *,'c is NaN',c,alpha_rho_K(1),alpha_rho_K(2),adv(1),adv(2),pres 
!          !      call s_MPI_abort()
           ! end if
        else           
            c = ((H - 5d-1*vel_sum)/gamma)
        end if

        if (mixture_err .and. c < 0d0) then
            c = 100.d0*sgm_eps
        else
            c = sqrt(c)
            !print *, 'c ::', c
        end if
    end subroutine s_compute_speed_of_sound
#:enddef
