!>
!! @file
!! @brief Contains module m_pressure_relaxation

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Pressure relaxation for the six-equation multi-component model via Newton--Raphson equilibration and volume-fraction
!! correction
module m_pressure_relaxation

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort, s_mpi_reduce_int_sum, s_mpi_allreduce_max
    use m_variables_conversion, only: s_convert_species_to_mixture_variables_kernel
    use m_eos, only: f_pressure, s_phase_internal_energy, s_phase_coefficients, s_phase_density_on_isentrope, f_is_state_dependent
    use m_root_finding, only: t_bracket, s_bracket_set, s_bracket_newton_step, bracket_step, bracket_done

    implicit none

    private; public :: s_pressure_relaxation_procedure, s_report_pressure_relaxation

    !> Cell updates in which the solve below stopped on its iteration cap, or on a bracket with no sign change, rather than on the
    !! tolerance. Counted because the answer then depends on that cap, and nothing else in the code says so.
    integer  :: n_hit_cap = 0
    real(wp) :: worst_residual = 0._wp  !< Largest |f| left behind when it did

contains

    !> The main pressure relaxation procedure
    subroutine s_pressure_relaxation_procedure(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer                                                :: i, j, k, l

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: alpha_rho, alpha
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho, alpha
        #:endif
        real(wp) :: rho, gamma, pi_inf, qv_mix
        integer  :: hit_cap, hit_cap_sum, unusable, unusable_sum
        real(wp) :: resid, resid_max

        ! Formed here, not one call deeper: CCE OpenACC accepts a num_fluids-sized array passed to a device routine from a
        ! parallel-loop body, and rejects the same call from inside another acc routine seq.
        hit_cap_sum = 0
        unusable_sum = 0
        resid_max = 0._wp
        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, alpha_rho, alpha, rho, gamma, pi_inf, qv_mix, hit_cap, resid, unusable]', &
                            & collapse=3, reduction='[[hit_cap_sum, unusable_sum], [resid_max]]', reductionOp='[+, max]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (mpp_lim) call s_correct_volume_fractions(q_cons_vf, j, k, l)

                    hit_cap = 0
                    resid = 0._wp
                    unusable = 0
                    if (s_needs_pressure_relaxation(q_cons_vf, j, k, l)) then
                        call s_equilibrate_pressure(q_cons_vf, j, k, l, hit_cap, resid, unusable)
                    end if
                    hit_cap_sum = hit_cap_sum + hit_cap
                    unusable_sum = unusable_sum + unusable
                    resid_max = max(resid_max, resid)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        alpha_rho(i) = q_cons_vf(i)%sf(j, k, l)
                        alpha(i) = q_cons_vf(eqn_idx%E + i)%sf(j, k, l)
                    end do

                    call s_convert_species_to_mixture_variables_kernel(rho, gamma, pi_inf, qv_mix, alpha, alpha_rho)

                    call s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        n_hit_cap = n_hit_cap + hit_cap_sum
        worst_residual = max(worst_residual, resid_max)

        ! Equilibration produced a density that is not a usable number. Continuing means advancing a cell that has
        ! no physical state, which is the kind of quiet wrongness that only shows up as a bad answer much later.
        if (unusable_sum > 0) then
            call s_mpi_abort('Pressure relaxation produced a non-physical phasic density. Exiting.')
        end if

    end subroutine s_pressure_relaxation_procedure

    !> One line at the end of a run if the equilibration ever stopped without reaching its tolerance, whether on the iteration cap
    !! or on a bracket whose ends did not straddle a root. Silence means every cell converged.
    impure subroutine s_report_pressure_relaxation

        integer  :: n_hit_glb
        real(wp) :: resid_glb

        call s_mpi_reduce_int_sum(n_hit_cap, n_hit_glb)
        call s_mpi_allreduce_max(worst_residual, resid_glb)
        if (n_hit_glb > 0 .and. proc_rank == 0) then
            print '(A,I0,A,ES10.3,A)', ' Pressure relaxation did not converge in ', n_hit_glb, &
                & ' cell updates; worst residual left behind ', resid_glb, ' against a 1e-10 tolerance.'
        end if

    end subroutine s_report_pressure_relaxation

    !> Check if pressure relaxation is needed for this cell
    logical function s_needs_pressure_relaxation(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: j, k, l
        integer                                             :: i

        s_needs_pressure_relaxation = .true.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > (1._wp - sgm_eps)) then
                s_needs_pressure_relaxation = .false.
            end if
        end do

    end function s_needs_pressure_relaxation

    !> A phase takes part in the equilibration if it is present with a positive mass. A vanishing phase transported to round-off can
    !! carry a negative partial density, which has no density to equilibrate; it keeps its state instead
    logical function f_phase_active(q_cons_vf, i, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: i, j, k, l

        f_phase_active = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > sgm_eps .and. q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, &
                                   & k, l) > 0._wp

    end function f_phase_active

    !> Correct volume fractions to physical bounds
    subroutine s_correct_volume_fractions(q_cons_vf, j, k, l)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp)                                               :: sum_alpha
        integer                                                :: i

        sum_alpha = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if ((q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) < 0._wp) .or. (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                & l) < 0._wp)) then
                q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 0._wp
                q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = 0._wp
            end if
            if (q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) > 1._wp) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = 1._wp
            sum_alpha = sum_alpha + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
        end do

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)/sum_alpha
        end do

    end subroutine s_correct_volume_fractions

    !> Main pressure equilibration using Newton-Raphson
    subroutine s_equilibrate_pressure(q_cons_vf, j, k, l, hit_cap, resid, unusable)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        integer, intent(out)                                   :: hit_cap, unusable
        real(wp), intent(out)                                  :: resid
        real(wp)                                               :: pres_relax, f_pres, df_pres
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: pres_K_init, rho_K_init, rho_K_s
        #:else
            real(wp), dimension(num_fluids) :: pres_K_init, rho_K_init, rho_K_s
        #:endif
        real(wp)           :: gamma_K, pi_inf_K, dpi_K, dgamma_K, c2_K, alpha_i, alpha_rho_i, rho_i, p_i, rho_s_i
        real(wp)           :: p_lo, p_hi, p_floor
        type(t_bracket)    :: br
        integer, parameter :: MAX_ITER = 50
        ! Pressure relaxation convergence tolerance
        real(wp), parameter :: TOLERANCE = 1.e-10_wp
        ! Relative bracket width at which the interval counts as collapsed, for a cell whose residual cannot reach TOLERANCE
        real(wp), parameter :: XTOLERANCE = 1.e-14_wp
        integer             :: iter, i
        logical             :: usable

        ! Initialize pressures
        pres_relax = 0._wp
        p_floor = -huge(1._wp)
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            rho_K_s(i) = 0._wp
            if (f_phase_active(q_cons_vf, i, j, k, l)) then
                ! Phasic internal energy carries the formation energy: alpha_rho_k*qv_k must be
                ! removed before inverting the stiffened-gas EOS, or a nonzero qv inflates the
                ! phasic pressure by rho_k*qv_k/gamma_k (this is what breaks the reactive burn).
                alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
                alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
                call s_phase_coefficients(alpha_rho_i, alpha_i, i, rho_i, gamma_K, pi_inf_K, dpi_K, dgamma_K)
                rho_K_init(i) = rho_i
                ! Seeded with the incoming density so that a cell whose phases already share a pressure -- where the bracket
                ! below is empty and the loop never runs -- commits the volume fractions it came in with, rather than whatever
                ! happened to be left in this array.
                rho_K_s(i) = rho_i
                pres_K_init(i) = ((q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, &
                            & k, l)*qvs(i))/q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l) - pi_inf_K)/gamma_K
                if (.not. f_is_state_dependent(i)) then
                    ! The lowest pressure at which this phase's isentrope is still defined. Each phasic pressure is held above
                    ! its own, exactly as before; the largest over the phases is kept separately and bounds the bracket from
                    ! below, because at it p + isentrope_B goes to zero for the limiting phase, its density with it, and the
                    ! residual to plus infinity.
                    if (pres_K_init(i) <= -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp) pres_K_init(i) = -(1._wp - 1.e-8_wp) &
                        & *isentrope_B(i) + 1.e-8_wp
                    p_floor = max(p_floor, -(1._wp - 1.e-8_wp)*isentrope_B(i) + 1.e-8_wp)
                end if
            else
                pres_K_init(i) = 0._wp
            end if
            pres_relax = pres_relax + q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)*pres_K_init(i)
        end do

        ! Bracket, in closed form. Every phasic density on its isentrope rises with pressure, so f(p) = -1 + sum m_k/rho_k^s(p)
        ! falls with it. At the lowest of the phasic pressures each rho_k^s is at or below rho_k^0, so the sum is at least
        ! sum alpha_k = 1 and f is non-negative; at the highest each is at or above it, so f is non-positive. The equilibrium
        ! pressure therefore lies between the lowest and the highest phasic pressure, which is where it has to be physically,
        ! and nothing has to be searched for to find the ends. The lower end is additionally held at the floor above, since
        ! below it a phase with a larger isentrope_B would be asked for the root of a negative number.
        p_lo = huge(1._wp)
        p_hi = -huge(1._wp)
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (f_phase_active(q_cons_vf, i, j, k, l)) then
                p_lo = min(p_lo, pres_K_init(i))
                p_hi = max(p_hi, pres_K_init(i))
            end if
        end do
        p_lo = max(p_lo, p_floor)

        ! The old iteration was an unsafeguarded Newton run for a fixed 50 passes from the alpha-weighted mean pressure, with
        ! the step clamped to the floor whenever it fell through it. It reached its cap often enough that the routine counts
        ! the cells where it did. The residual is evaluated here exactly as it was; only the choice of the next pressure
        ! changes, to a Newton step that is replaced by the bisection of the bracket whenever it would leave it.
        ! Only so that the report below reads a defined value if the bracket is empty and the loop never runs.
        f_pres = 0._wp
        call s_bracket_set(br, p_lo, p_hi, TOLERANCE, XTOLERANCE, MAX_ITER)

        $:GPU_LOOP(parallelism='[seq]')
        do iter = 1, MAX_ITER + 1
            if (br%state > bracket_step) exit
            pres_relax = br%x

            f_pres = -1._wp
            df_pres = 0._wp
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                if (f_phase_active(q_cons_vf, i, j, k, l) .and. any_state_dependent_eos) then
                    rho_i = rho_K_init(i)
                    p_i = pres_K_init(i)
                    call s_phase_density_on_isentrope(i, rho_i, p_i, pres_relax, rho_s_i, c2_K)
                    rho_K_s(i) = rho_s_i
                    f_pres = f_pres + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
                    df_pres = df_pres - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/(rho_K_s(i)**2*c2_K)
                else if (f_phase_active(q_cons_vf, i, j, k, l)) then
                    ! Isentropic relation: rho = rho0 * (p/p0)^(1/gamma), Saurel et al. JFM (2009)
                    rho_K_s(i) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/max(q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                            & l), sgm_eps)*((pres_relax + isentrope_B(i))/(pres_K_init(i) + isentrope_B(i)))**(1._wp/isentrope_n(i))
                    f_pres = f_pres + q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
                    df_pres = df_pres - q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, &
                                                  & l)/(isentrope_n(i)*rho_K_s(i)*(pres_relax + isentrope_B(i)))
                end if
            end do

            call s_bracket_newton_step(br, f_pres, df_pres)
        end do

        ! Nothing reads pres_relax past this point: what the cell commits below is rho_K_s, which belongs to the last pressure
        ! evaluated -- and for a fixed bracket that is br%x in every outcome, so the densities go with the pressure the solve
        ! settled on.
        hit_cap = 0
        resid = 0._wp
        unusable = 0
        ! p_lo >= p_hi is a cell whose phases already share a pressure: the bracket is empty, the solve does not run, and the
        ! cell is left as it came in. That is the right answer, not a miss, so it is not counted.
        if (br%state /= bracket_done .and. p_lo < p_hi) then
            hit_cap = 1
            ! A NaN residual has to be mapped, not passed on: max() with a NaN keeps the other operand, so the
            ! reduction would report the miss as zero.
            resid = abs(f_pres)
            if (f_pres /= f_pres) resid = huge(1._wp)
        end if

        ! An unconverged-but-physical density is still used: the Newton often stops on the iteration cap and the
        ! answer then depends on that cap, so refusing it would move every six-equation result. A density that is
        ! not a usable number is different -- the equilibration has failed, and the cell has no physical state to
        ! continue from. Flag it and let the caller stop the run rather than quietly leaving the cell unrelaxed.
        ! Written as .not. (> 0) so a NaN takes the same path as a non-positive value.
        ! No early return: Cray OpenACC rejects a RETURN inside an accelerator loop.
        usable = .true.
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            if (f_phase_active(q_cons_vf, i, j, k, l) .and. .not. (rho_K_s(i) > 0._wp)) usable = .false.
        end do

        if (usable) then
            $:GPU_LOOP(parallelism='[seq]')
            do i = 1, num_fluids
                if (f_phase_active(q_cons_vf, i, j, k, l)) q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, &
                    & l) = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)/rho_K_s(i)
            end do
        else
            unusable = 1
        end if

    end subroutine s_equilibrate_pressure

    !> Correct internal energies using equilibrated pressure
    subroutine s_correct_internal_energies(q_cons_vf, j, k, l, rho, gamma, pi_inf, qv_mix)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        integer, intent(in)                                    :: j, k, l
        real(wp), intent(in)                                   :: rho, gamma, pi_inf, qv_mix
        real(wp)                                               :: dyn_pres, pres_relax, alpha_i, alpha_rho_i, e_i
        integer                                                :: i

        dyn_pres = 0._wp
        $:GPU_LOOP(parallelism='[seq]')
        do i = eqn_idx%mom%beg, eqn_idx%mom%end
            dyn_pres = dyn_pres + 5.e-1_wp*q_cons_vf(i)%sf(j, k, l)*q_cons_vf(i)%sf(j, k, l)/max(rho, sgm_eps)
        end do

        pres_relax = f_pressure(q_cons_vf(eqn_idx%E)%sf(j, k, l) - dyn_pres, gamma, pi_inf, qv_mix)

        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            alpha_i = q_cons_vf(i + eqn_idx%adv%beg - 1)%sf(j, k, l)
            alpha_rho_i = q_cons_vf(i + eqn_idx%cont%beg - 1)%sf(j, k, l)
            call s_phase_internal_energy(pres_relax, alpha_i, alpha_rho_i, i, e_i)
            q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(j, k, l) = e_i
        end do

    end subroutine s_correct_internal_energies

end module m_pressure_relaxation
