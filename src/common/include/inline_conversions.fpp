#:def s_compute_speed_of_sound()
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c, G, alpha_rho_K)
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
        real(kind(0d0)), optional, dimension(num_fluids), intent(IN) :: G, alpha_rho_K
        real(kind(0d0)), intent(OUT) :: c
        real(kind(0d0)) :: blkmod1, blkmod2
        real(kind(0d0)) :: log_rho_mix_ratio, deno_rho_sq, phi_mix, theta_E, rho0_mix 
       

        !Local variables used for computation only
       

        integer :: q

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
            log_rho_mix_ratio = log(rho/sum(adv(:)*fluid_pp(:)%rho0))
            deno_rho_sq = sum(fluid_pp(:)%mg_a*adv(:)*&
                                fluid_pp(:)%rho0**2+fluid_pp(:)%mg_b*&
                                alpha_rho_K(:)*fluid_pp(:)%rho0)
            phi_mix = exp(sum(adv(:)*fluid_pp(:)%gamma-&
                       adv(:)*fluid_pp(:)%gamma*fluid_pp(:)%rho0/rho))
            theta_E = sum(adv(:)*fluid_pp(:)%ein_cv(2))
            rho0_mix = sum(adv(:)*fluid_pp(:)%rho0)
            c = sum(alpha_rho_K(:)*fluid_pp(:)%pi_inf/fluid_pp(:)%rho0)&
                   +log_rho_mix_ratio*sum(alpha_rho_K(:)*fluid_pp(:)%pi_inf*(fluid_pp(:)%qv-1.d0)/fluid_pp(:)%rho0)&
                   +(log_rho_mix_ratio**2)*sum(alpha_rho_K(:)*fluid_pp(:)%pi_inf*0.5d0*(fluid_pp(:)%qv-2.d0)/fluid_pp(:)%rho0)&
                   +pres*sum(alpha_rho_K(:)*fluid_pp(:)%mg_b)/sum(adv(:)*fluid_pp(:)%mg_a*fluid_pp(:)%rho0+fluid_pp(:)%mg_b*alpha_rho_K(:))&
                   -log_rho_mix_ratio*sum((alpha_rho_K(:)**2)*fluid_pp(:)%mg_b*fluid_pp(:)%pi_inf)/deno_rho_sq&
                   -(log_rho_mix_ratio**2)*sum((alpha_rho_K(:)**2)*fluid_mg(:)%mg_b*fluid_pp(:)%pi_inf*0.5d0*(fluid_pp(:)%qv-2.d0))/deno_rho_sq&
                   +pres*(sum(fluid_pp(:)%gamma*fluid_pp(:)%mg_a*adv(:)*fluid_pp(:)%rho0)/rho0_mix+sum(fluid_pp(:)%mg_b*fluid_pp(:)%gamma*adv(:)))&
                   -log_rho_mix_ratio*(fluid_pp(:)%gamma*adv(:)*fluid_pp(:)%mg_a*fluid_pp(:)%pi_inf&
                   +sum(fluid_pp(:)%gamma*alpha_rho_K(:)*fluid_pp(:)%mg_b*fluid_pp(:)%pi_inf)/rho0_mix)&
                   -(log_rho_mix_ratio**2)*(sum(fluid_pp(:)%gamma*adv(:)*fluid_pp(:)%mg_a*fluid_pp(:)%pi_inf*0.5d0*(fluid_pp(:)%qv-2.d0))&
                   +sum(fluid_pp(:)%gamma*alpha_rho_K(:)*fluid_pp(:)%mg_b*fluid_pp(:)%pi_inf*0.5d0*(fluid_pp(:)%qv-2.d0))/rho0_mix)&
                   +(phi_mix**2)*(exp(phi_mix*theta_E)/(exp(phi_mix*theta_E)-1)**2)*(sum(adv(:)*fluid_pp(:)%gamma*fluid_pp(:)%mg_a*&
                   (fluid_pp(:)%rho0**2)*fluid_pp(:)%ein_cv(1)*fluid_pp(:)%ein_cv(2))&
                   +sum(adv(:)*fluid_pp(:)%gamma*fluid_pp(:)%mg_b*fluid_pp(:)%rho0*fluid_pp(:)%ein_cv(1)*fluid_pp(:)%ein_cv(2)**2))
            c = c/rho
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

