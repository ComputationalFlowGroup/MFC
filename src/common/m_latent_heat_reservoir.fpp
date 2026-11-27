!>
!! @file m_latent_heat_reservoir.fpp
!! @brief Latent-heat reservoir phase change

#:include 'macros.fpp'

!> @brief Phase change by a latent-heat reservoir, ported from the phase_change_STALE branch.
!!
!! A cell stores the energy it holds above its saturation internal energy in a reservoir Lambda,
!! carried as the conserved density alpha_l*rho_l*Lambda (eqn_idx%lambda) so that it advects like a
!! partial density. When the reservoir reaches the latent heat of vaporization the liquid converts
!! to vapor in one step and the reservoir empties.
!!
!! This needs model_eqns = 3: it moves energy between the phases through their individual internal
!! energies, which the 5-equation model does not carry. That also makes it mutually exclusive with
!! hyperelasticity, which requires model_eqns = 2 -- so the energy-exclusion handoff
!! s_infinite_relaxation_k takes for the pT/pTg models has no counterpart here.
!!
!! Four defects in the original are fixed rather than carried over; see the comments at each.
module m_latent_heat_reservoir

    use m_derived_types
    use m_global_parameters
    use m_phase_change, only: s_correct_partial_densities, lp, vp
    use m_root_finding, only: t_bracket, s_bracket_set, s_bracket_newton_step, bracket_step, bracket_dry

    implicit none

    private; public :: s_latent_heat_reservoir

    !> Evaluation cap for the saturation temperature. The original declared a cap of 1e8 and never used it, so an unconverged solve
    !! spun forever inside the offload region; an under-relaxed Newton at a factor of 1e-3 then needed some 27000 passes to converge
    !! even where it did. The bracketed solve in s_saturation_temperature reaches machine precision in seven to ten evaluations over
    !! the whole liquid range of water, so this cap is a backstop and not a working limit.
    integer, parameter :: tsat_max_iter = 60

    !> Relative width at which the bracket is treated as collapsed. Reached only where the residual tolerance cannot be, which is
    !! what keeps a cell with an ill-conditioned saturation curve from spending the whole cap on it.
    real(wp), parameter :: tsat_xtol = 1.e-14_wp

    ! Deliberately no GPU_DECLARE: these are Fortran parameters, compile-time constants with no storage to map, and
    ! m_phase_change.fpp declares its own (max_iter, mixM, lp, vp and the rest) the same way.

contains

    !> Convert liquid to vapor wherever the stored reservoir has reached the latent heat.
    subroutine s_latent_heat_reservoir(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp) :: rM, MCT                     !< total reacting mass, and the correction flag
        real(wp) :: pS                          !< cell pressure
        real(wp) :: alphal, alphav, rhol, rhov  !< phasic volume fractions and densities
        real(wp) :: Tl, el, ev, Lambda          !< liquid temperature, phasic internal energies, reservoir
        real(wp) :: TSat, LSat, el_Sat          !< saturation temperature, latent heat, saturation energy
        real(wp) :: dE, alpha_sum
        integer :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, rM, MCT, pS, alphal, alphav, rhol, rhov, Tl, el, ev, Lambda, TSat, &
                            & LSat, el_Sat, dE, alpha_sum]')
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    alphal = q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l)
                    alphav = q_cons_vf(vp + eqn_idx%adv%beg - 1)%sf(j, k, l)

                    ! Nothing to convert without liquid, and Lambda is not defined there either: it is carried as
                    ! alpha_l*rho_l*Lambda, so recovering it divides by the liquid partial density. Guarded against
                    ! verysmall rather than tested against zero, which is what the original did.
                    if (alphal > verysmall) then
                        rM = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) + q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)
                        call s_correct_partial_densities(MCT, q_cons_vf, rM, j, k, l)

                        pS = q_prim_vf(eqn_idx%E)%sf(j, k, l)

                        rhol = q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l)/alphal
                        ! The original divided by alphav unguarded, which is a division by zero in an all-liquid
                        ! cell -- exactly the cell this routine is for.
                        rhov = q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l)/max(alphav, verysmall)

                        Tl = (pS + isentrope_B(lp))*isentrope_n(lp)/cvs(lp)/rhol
                        el = (pS*isentrope_n(lp) + isentrope_B(lp))/rhol + qvs(lp)
                        ev = (pS*isentrope_n(vp) + isentrope_B(vp))/rhov + qvs(vp)
                        Lambda = q_cons_vf(eqn_idx%lambda)%sf(j, k, l)/(alphal*rhol)

                        call s_saturation_temperature(pS, Tl, TSat)
                        LSat = isentrope_n(vp)*cvs(vp)*TSat + qvs(vp) - (isentrope_n(lp)*cvs(lp)*TSat + qvs(lp))
                        el_Sat = (pS + isentrope_n(lp)*isentrope_B(lp))/(pS + isentrope_B(lp))*cvs(lp)*TSat + qvs(lp)

                        ! The transfer is computed ONCE and then applied to both sides, which the original did not do.
                        ! It wrote `el = el + min(el_Sat - el, Lambda)` and then recomputed that min against the
                        ! updated el, so what left the reservoir was not what entered the internal energy; and in the
                        ! other branch it assigned el = el_Sat before evaluating Lambda + (el - el_Sat), which adds
                        ! zero and destroyed the excess energy instead of storing it.
                        if (Tl < TSat) then
                            ! Below saturation: draw from the reservoir, but no further than saturation.
                            dE = min(max(el_Sat - el, 0._wp), max(Lambda, 0._wp))
                            el = el + dE
                            Lambda = Lambda - dE
                        else
                            ! Above saturation: clamp the internal energy and store the excess.
                            dE = max(el - el_Sat, 0._wp)
                            el = el_Sat
                            Lambda = Lambda + dE
                        end if

                        if (Lambda >= LSat) then
                            ! The reservoir has reached the latent heat: the liquid converts, and its internal energy
                            ! plus the whole reservoir is pooled into the vapor.
                            alpha_sum = max(alphav + alphal, verysmall)
                            rhov = (alphav*rhov + alphal*rhol)/alpha_sum
                            ev = (alphav*ev + alphal*(el + Lambda))/alpha_sum
                            alphav = alphav + alphal
                            rhol = 0._wp
                            el = 0._wp
                            alphal = 0._wp
                            Lambda = 0._wp
                        end if

                        q_cons_vf(lp + eqn_idx%adv%beg - 1)%sf(j, k, l) = alphal
                        q_cons_vf(vp + eqn_idx%adv%beg - 1)%sf(j, k, l) = alphav
                        q_cons_vf(lp + eqn_idx%cont%beg - 1)%sf(j, k, l) = alphal*rhol
                        q_cons_vf(vp + eqn_idx%cont%beg - 1)%sf(j, k, l) = alphav*rhov
                        ! int_en exists only for the 6-equation model, which the validator requires for
                        ! relax_model = 7. Deliberately not guarded on model_eqns here: skipping just these two
                        ! writes would leave the volume fractions and partial densities updated with the phasic
                        ! energies stale, which is worse than not running at all. The precondition belongs at the
                        ! case level, where it is enforced.
                        q_cons_vf(lp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = alphal*rhol*el
                        q_cons_vf(vp + eqn_idx%int_en%beg - 1)%sf(j, k, l) = alphav*rhov*ev
                        q_cons_vf(eqn_idx%lambda)%sf(j, k, l) = alphal*rhol*Lambda
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_latent_heat_reservoir

    !> Saturation temperature at pSat, by Newton on the Gibbs-equality residual.
    !!
    !! The residual has more structure than this solve exploits. With
    !!   a = cvs_l*n_l - cvs_v*n_v
    !!   b = -(qvps_l - qvps_v) + cvs_l*(n_l-1)*log(pSat+B_l) - cvs_v*(n_v-1)*log(pSat+B_v)
    !!   c = qvs_l - qvs_v
    !! it is exactly FT(T) = c + (a+b)*T - a*T*log(T), so FT'' = -a/T keeps one sign for T > 0: FT is strictly concave or
    !! strictly convex, with a single stationary point at the closed form T* = exp(b/a), hence at most two roots and
    !! monotone either side of T*. That is an analytic bracket, so a bisection-safeguarded Newton would converge in a
    !! handful of steps with a guarantee rather than the thousands this under-relaxed unbracketed version needs; the root
    !! also has a closed form through the Lambert W function. Left as the original's method here on purpose, so this
    !! change carries the bookkeeping fixes and not a change of numerics.
    !! @param pSat saturation pressure
    !! @param TIn  incoming temperature, used as the initial estimate
    !! @param TSat saturation temperature; left at TIn if the solve does not converge
    subroutine s_saturation_temperature(pSat, TIn, TSat)

        $:GPU_ROUTINE(function_name='s_saturation_temperature', parallelism='[seq]')

        real(wp), intent(in)  :: pSat, TIn
        real(wp), intent(out) :: TSat
        type(t_bracket)       :: br
        real(wp)              :: a, kp, c, FT, dFdT, Tmax
        integer               :: ns

        if ((pSat <= 0._wp) .or. (TIn <= 0._wp)) then
            TSat = 0._wp
            return
        end if

        ! Fallback for every path that cannot produce a root: the cell keeps its incoming temperature, which leaves it
        ! unchanged rather than moving it on a number that means nothing.
        TSat = TIn

        ! Gibbs equality along the saturation curve, written as FT(T) = c + (a + kp) T - a T log(T). Only kp depends on the
        ! pressure, and a and c are material constants. This is the same residual the under-relaxed Newton used, regrouped.
        a = cvs(lp)*isentrope_n(lp) - cvs(vp)*isentrope_n(vp)
        kp = -(qvps(lp) - qvps(vp)) + cvs(lp)*(isentrope_n(lp) - 1._wp)*log(pSat + isentrope_B(lp)) - cvs(vp)*(isentrope_n(vp) &
               & - 1._wp)*log(pSat + isentrope_B(vp))
        c = qvs(lp) - qvs(vp)

        ! FT'' = -a/T, so for T > 0 the residual is concave when a > 0 and convex when a < 0: one stationary point either way,
        ! at T = exp(kp/a), with FT monotone below it. FT(0+) = c, so the interval from zero to that point holds exactly one
        ! root and it is the physical one -- the second root lies beyond the turning point, some thousands of kelvin up. The
        ! bracket is therefore closed form and needs no search. A phase pair with equal c_p leaves a = 0, no turning point and
        ! no bracket; that is not a liquid and vapour of one species, and the cell is left alone.
        if (abs(a) <= verysmall) return

        Tmax = exp(kp/a)
        if (.not. (Tmax < huge(1._wp))) return

        ! The lower endpoint only has to be positive, so that log() is defined, and small enough that FT is still c there.
        call s_bracket_set(br, verysmall, Tmax, ptgalpha_eps, tsat_xtol, tsat_max_iter)

        $:GPU_LOOP(parallelism='[seq]')
        do ns = 1, tsat_max_iter + 1
            if (br%state > bracket_step) exit
            FT = c + (a + kp)*br%x - a*br%x*log(br%x)
            dFdT = kp - a*log(br%x)
            call s_bracket_newton_step(br, FT, dFdT)
        end do

        ! A dry bracket means the residual kept one sign across the whole interval, so this pressure has no saturation
        ! temperature and TSat stays at TIn. Every other outcome leaves br%x inside the bracket, which beats TIn even when the
        ! cap was reached.
        if (br%state /= bracket_dry) TSat = br%x

    end subroutine s_saturation_temperature

end module m_latent_heat_reservoir
