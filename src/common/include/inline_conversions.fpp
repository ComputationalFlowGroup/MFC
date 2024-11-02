#:def s_compute_speed_of_sound()
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c, alpha_rho_K)
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
        real(kind(0d0)), intent(OUT) :: c
        real(kind(0d0)) :: blkmod1, blkmod2
       

        !Local variables used for computation only
        real(kind(0d0)) :: pref, gamma_inv, xi, rho_K, pref_prime,&
        rho_eref_prime, gamma_avg, gam
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

        elseif ((model_eqns == 5) .and. (MGEoS_model == 1)) then
            !Note that pref and gamma are primitive state
            !variables for Mie-Gruneisen EoS and gamma = (1/Gamma(rho))
            gamma_avg = 0d0
            c = 0d0
            do q = 1, num_fluids
               !rho_K = alpha_rho_K(q)/adv(q)
               if (adv(q) .gt. 1d-8) then
               xi    = 1d0 - rho0(q)*adv(q)/alpha_rho_K(q)
            !   if (xi .lt. 0d0) then 
            !       print *,'xi :', xi
            !   end if
               if (xi .le. 1d-16) then 
                   xi = 1d-16
               end if
               pref  = pi_infs(q) + &
               rho0(q)*(mg_a(q)**2d0)*xi/(1d0-mg_b(q)*xi)**2d0
                
               gamma_avg = gamma_avg + &
               adv(q)*((alpha_rho_K(q))/(adv(q)*rho0(q)))**qvps(q)/gammas(q)
                
               gamma_inv = &
               (alpha_rho_K(q)/(adv(q)*rho0(q)))**qvps(q)/gammas(q)
               
               gam = &
               gammas(q)*(adv(q)*rho0(q)/(alpha_rho_K(q)))**qvps(q)

               pref_prime  = &
               (mg_a(q)**2d0)*((rho0(q)*adv(q)/alpha_rho_K(q))**2d0)*(1d0+mg_b(q)*xi)/(1d0-mg_b(q)*xi)**3d0
               
               if ((1d0 - mg_b(q)*xi) .lt. 1d-16) then 
                   pref_prime = 0d0
                   pref = pi_infs(q)
               end if

               rho_eref_prime = &
               0.5d0*(pref*adv(q)/alpha_rho_K(q) + &
               pref_prime*(alpha_rho_K(q)/(adv(q)*rho0(q))-1d0))
               
               !Mie-gruneisen sound-speed mixture
               !c = c + &
               !(alpha_rho_K(q)/rho)*((1d0+(1d0-qvps(q))*gamma_inv)*(pres-pref)*(adv(q)/(alpha_rho_K(q)+sgm_eps)) &
               !+pref*adv(q)/(alpha_rho_K(q)+sgm_eps) + pref_prime*gamma_inv - rho_eref_prime)
              
              !Bubbles sound-speed mixture
              !c = c + &
              ! (adv(q)**2d0/(alpha_rho_K(q)+sgm_eps))/(((1d0+(1d0-qvps(q))*gamma_inv)*(pres-pref)*(adv(q)/(alpha_rho_K(q)+sgm_eps)) &
              ! +pref*adv(q)/(alpha_rho_K(q)+sgm_eps) + &
              !  pref_prime*gamma_inv - rho_eref_prime+sgm_eps)/gamma_inv)
               
               !Frozen speed of sound
              c = c + &
               alpha_rho_K(q)*((gam+(1d0-qvps(q)))*(pres-pref)*(adv(q)/(alpha_rho_K(q))) &
               +gam*pref*adv(q)/alpha_rho_K(q) + &
                pref_prime - rho_eref_prime*gam)
            if (c /= c .or. (c .lt. -1d-16)) then
                print &
                *,'c:',c,'alpha_rho_K(i)',alpha_rho_K(q),'pref_prime',pref_prime,'rho_eref_prime',rho_eref_prime,'pres &
                (MPa)',pres/1d6,&
                    'gamma',gam,'adv(i)',adv(q),'q',q,'xi',xi
            end if
            else
                c = c + 0d0
            end if

            end do 
            !c = c/gamma_avg
            !c = 1d0/(rho*c)
            c = c/rho
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
