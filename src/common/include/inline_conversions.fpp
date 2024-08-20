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
        real(kind(0d0)) :: log_rho_mix_ratio, deno_rho_sq, phi_mix, theta_E, rho0_mix, term1, term2, term3 
       

        !Local variables used for computation only
       

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
            deno_rho_sq = 0.d0
            phi_mix     = 0.d0
            theta_E     = 0.d0
            rho0_mix    = 0.d0
            term1       = 0.d0
            rho0_mix    = 0.d0
            term2       = 0.d0
            term3       = 0.d0
            do r = 1, num_fluids
                !print*, 'rho0_mix' , rho0_mix
                rho0_mix    = rho0_mix + adv(r)*rho0(r)

                deno_rho_sq = deno_rho_sq+mg_a(r)*adv(r)*&
                                rho0(r)**2.d0+mg_b(r)*&
                               alpha_rho_K(r)*rho0(r)
                phi_mix     = phi_mix+adv(r)*gammas(r)-&
                                adv(r)*gammas(r)*rho0(r)/rho
                theta_E     = theta_E+adv(r)*ein_cv2(r)
                term1       = term1 + adv(r)*mg_a(r)*rho0(r)+mg_b(r)*alpha_rho_K(r)
                term2       = term2 + adv(r)*gammas(r)*mg_a(r)*&
                                (rho0(r)**2.d0)*ein_cv1(r)*ein_cv2(r)**2.d0
                term3       = term3+adv(r)*gammas(r)*mg_b(r)*rho0(r)*ein_cv1(r)*ein_cv2(r)**2.d0
            end do

                log_rho_mix_ratio = log(rho/rho0_mix)
           !print*,'adv1',adv(1),'adv2',adv(2),'rho0_mix',rho0_mix,'log',log_rho_mix_ratio,'rho0(1)',rho0(1),'rho0(2)',rho0(2)                 
                phi_mix = exp(phi_mix)
                c = 0.d0
            do r = 1,num_fluids     
                c = c   +alpha_rho_K(r)*pi_infs(r)/rho0(r)&
                        +log_rho_mix_ratio*(alpha_rho_K(r)*pi_infs(r)*(qvs(r)-1.d0)/rho0(r))&
                        +(log_rho_mix_ratio**2.d0)*(alpha_rho_K(r)*pi_infs(r)*0.5d0*(qvs(r)-2.d0)/rho0(r)) &
                        +pres*(alpha_rho_K(r)*mg_b(r))/term1 &
                        -log_rho_mix_ratio*alpha_rho_K(r)**2.d0*mg_b(r)*pi_infs(r)/deno_rho_sq &
                        -(log_rho_mix_ratio**2.d0)*(alpha_rho_K(r)**2.d0)*mg_b(r)*pi_infs(r)*0.5d0*(qvs(r)-2.d0)/deno_rho_sq &
                        +pres*(gammas(r)*mg_a(r)*adv(r)*rho0(r)/rho0_mix+mg_b(r)*gammas(r)*adv(r)) &
                        -log_rho_mix_ratio*(gammas(r)*adv(r)*mg_a(r)*pi_infs(r) &
                        +gammas(r)*alpha_rho_K(r)*mg_b(r)*pi_infs(r)/rho0_mix ) &
                        -(log_rho_mix_ratio**2.d0)*(gammas(r)*adv(r)*mg_a(r)*pi_infs(r)*0.5d0*(qvs(r)-2.d0) &
                        +gammas(r)*alpha_rho_K(r)*mg_b(r)*pi_infs(r)*0.5d0*(qvs(r)-2.d0)/rho0_mix) 
            end do
            c = c + (phi_mix**2.d0)*(exp(phi_mix*theta_E)/(exp(phi_mix*theta_E)-1.d0)**2.d0)*(term2/rho+term3)
            c = c/rho
            !c = 1d0
            !print *, 'print c before the grinder :: c :: ',sqrt(c)
        else           
            c = ((H - 5d-1*vel_sum)/gamma)
        end if

        if (mixture_err .and. c < 0d0) then
            c = 100.d0*sgm_eps
        else
            c = sqrt(c)
        end if
    end subroutine s_compute_speed_of_sound
#:enddef
