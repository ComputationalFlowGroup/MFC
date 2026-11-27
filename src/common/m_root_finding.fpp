!>
!! @file m_root_finding.fpp
!! @brief Contains module m_root_finding

#:include 'macros.fpp'

!> @brief Scalar root finding by Newton with a bisection safeguard, shared by the relaxation solvers.
!!
!! The solvers in m_latent_heat_reservoir, m_pressure_relaxation and m_phase_change all look for the single root of a scalar
!! residual that is monotone over the interval of interest. Each one carried its own loop: a bare Newton with no bracket, no
!! fallback, and nothing to do when the step left the physical range except take it anyway. The iteration logic is the same in
!! every case, so it lives here once and the residuals stay with the physics.
!!
!! The interface is reverse communication: this module never calls the residual. The caller evaluates f and f' at t_bracket%x in
!! its own loop body and passes them back. That keeps procedure pointers and callbacks out of device code -- a GPU_ROUTINE cannot
!! call through one -- and lets each caller keep whatever state its residual needs in its own scope. It also means the residual
!! is written out once per caller: the endpoints are probed by the same loop body as the iterations, so no call site has to
!! evaluate its residual a second time to set the bracket up.
!!
!! Call pattern, identical at every site:
!!
!!     call s_bracket_set(br, lo, hi, ftol, xtol, max_iter)      ! or s_bracket_seek for an open far end
!!     do it = 1, br%max_iter + 1
!!         if (br%state > bracket_step) exit
!!         <evaluate f and df at br%x>
!!         call s_bracket_newton_step(br, f, df)
!!     end do
!!     if (br%state /= bracket_dry) <accept br%x>
!!
!! The first two passes evaluate the two endpoints; from the third the point is chosen by Newton, or by bisection where Newton
!! would leave the bracket. The outer `do` is bounded rather than a `do while` so the loop has a trip count the compiler can
!! see, and it exits rather than returns: Cray OpenACC rejects a RETURN inside an accelerator loop.
!!
!! Nothing here is allocated, and no array crosses the interface, so there is nothing to map to the device and none of the Cray
!! restrictions on passing device-global-bound arrays into a seq routine apply.
module m_root_finding

    use m_precision_select, only: wp

    implicit none

    private
    public :: t_bracket, s_bracket_set, s_bracket_seek, s_bracket_newton_step
    public :: bracket_probe, bracket_seek, bracket_step, bracket_done, bracket_capped, bracket_dry

    !> Evaluating the near endpoint, to learn which sign the residual has there.
    integer, parameter :: bracket_probe = 1
    !> Evaluating the far endpoint, and pushing it outward while the sign has not changed.
    integer, parameter :: bracket_seek = 2
    !> Bracketed; taking safeguarded Newton steps. Every state above this one is terminal, which is what `state > bracket_step`
    !! tests for at the call sites.
    integer, parameter :: bracket_step = 3
    !> Converged: the residual met ftol, or the interval collapsed to xtol.
    integer, parameter :: bracket_done = 4
    !> Stopped on the iteration cap with the residual still above ftol. The answer is still inside the bracket.
    integer, parameter :: bracket_capped = 5
    !> No sign change was found, so there is no root to bracket and br%x means nothing. The caller keeps its own estimate.
    integer, parameter :: bracket_dry = 6

    !> State of one scalar solve. Scalars only: the whole object is a local variable of the calling device routine, so it is private
    !! per thread with nothing to map.
    type :: t_bracket
        real(wp) :: lo        !< Near endpoint, whose residual has the sign recorded in slo
        real(wp) :: hi        !< Far endpoint; still moving while the state is bracket_seek
        real(wp) :: x         !< Point at which the caller must evaluate the residual next, and the answer once terminal
        real(wp) :: slo       !< Sign of the residual at lo, +1 or -1; set when the near endpoint is probed
        real(wp) :: grow      !< Factor by which the probe span is multiplied while seeking; at or below one the bracket is fixed
        real(wp) :: ftol      !< Converged when abs(f) falls to this
        real(wp) :: xtol      !< Converged when the interval width falls to xtol*(abs(lo) + abs(hi)); zero disables the test
        integer  :: state     !< One of the bracket_* parameters above
        integer  :: iter      !< Residual evaluations consumed so far
        integer  :: max_iter  !< Cap on that count
    end type t_bracket

contains

    !> Solve inside a bracket that is known in closed form. lo and hi must straddle the root; if they do not, the solve ends in
    !! bracket_dry after evaluating both, and the caller keeps its own estimate.
    subroutine s_bracket_set(br, lo, hi, ftol, xtol, max_iter)

        $:GPU_ROUTINE(function_name='s_bracket_set', parallelism='[seq]')

        type(t_bracket), intent(out) :: br
        real(wp), intent(in)         :: lo, hi, ftol, xtol
        integer, intent(in)          :: max_iter

        ! A growth factor of one is what tells the seek below that the far endpoint is fixed.
        call s_bracket_seek(br, lo, hi, 1._wp, ftol, xtol, max_iter)

    end subroutine s_bracket_set

    !> Solve where only the near endpoint is known. lo is that endpoint, x0 the first probe beyond it, and the span x0 - lo is
    !! multiplied by grow until the residual changes sign. Used where the far end of the physical range is a pressure with no finite
    !! bound, only the knowledge that the residual eventually turns over. grow must exceed one; the iteration cap bounds the search
    !! either way, so a factor too close to one ends in bracket_dry rather than running on.
    subroutine s_bracket_seek(br, lo, x0, grow, ftol, xtol, max_iter)

        $:GPU_ROUTINE(function_name='s_bracket_seek', parallelism='[seq]')

        type(t_bracket), intent(out) :: br
        real(wp), intent(in)         :: lo, x0, grow, ftol, xtol
        integer, intent(in)          :: max_iter

        br%lo = lo
        br%hi = x0
        br%x = lo
        br%slo = 0._wp
        br%grow = grow
        br%ftol = ftol
        br%xtol = xtol
        br%iter = 0
        br%max_iter = max_iter
        br%state = bracket_probe

        ! An empty or inverted interval has no root to find. Written as .not. (>) so a NaN endpoint takes this path too.
        if (.not. (x0 > lo)) br%state = bracket_dry

    end subroutine s_bracket_seek

    !> Consume the residual f and its derivative df at br%x, and leave the next trial point in br%x.
    !!
    !! In every terminal state but one, br%x is the point the caller evaluated last, so anything else the residual produced on
    !! the way -- densities, temperatures, whatever the caller kept -- belongs to br%x and can be committed with it. The
    !! exception is a bracket_dry reached while seeking, where br%x is a probe that was never evaluated; a bracket_dry from a
    !! fixed bracket still holds the last point. The cap is not tested while the near endpoint is being probed, so that this
    !! holds for any cap of three or more.
    subroutine s_bracket_newton_step(br, f, df)

        $:GPU_ROUTINE(function_name='s_bracket_newton_step', parallelism='[seq]')

        type(t_bracket), intent(inout) :: br
        real(wp), intent(in)           :: f, df
        real(wp)                       :: span, xn

        if (br%state > bracket_step) return

        br%iter = br%iter + 1

        ! Tested before anything else, so a point that lands on the root while the endpoints are still being probed is kept
        ! rather than stepped over.
        if (abs(f) <= br%ftol) then
            br%state = bracket_done
            return
        end if

        if (br%state == bracket_probe) then
            ! The near endpoint only has to give up its sign. Note that a NaN there yields slo = +-1 from sign(), which is
            ! harmless: the seek below then compares against it and either brackets or runs out, and no NaN is propagated.
            br%slo = sign(1._wp, f)
            br%x = br%hi
            br%state = bracket_seek
            return
        end if

        if (br%state == bracket_seek) then
            ! Same sign as the near endpoint: the root, if there is one, lies further out. A NaN residual fails this test and
            ! is taken as the sign change, which ends the seek and hands the point to the bisection below; that is the safe
            ! direction, since pushing the probe outward on a NaN would run the span up to infinity.
            if (f*br%slo > 0._wp) then
                ! A fixed bracket has nowhere to grow to, so the two endpoints not straddling the root is the whole answer.
                if (br%grow <= 1._wp) then
                    br%state = bracket_dry
                    return
                end if

                ! The interval already searched is not worth revisiting, so the near endpoint moves up to the probe.
                span = br%x - br%lo
                br%lo = br%x
                br%x = br%lo + br%grow*span

                ! The span is geometric and reaches the overflow threshold in a few tens of steps. Written as .not. (<) to
                ! catch an infinite or NaN probe as well.
                if (.not. (br%x < huge(1._wp))) br%state = bracket_dry
                if (br%iter >= br%max_iter) br%state = bracket_dry
                return
            end if

            ! Sign change: the probe is the far endpoint, and the solve proceeds from here as a bracketed one.
            br%hi = br%x
            br%state = bracket_step
        else
            ! Ordinary step: the evaluated point replaces whichever endpoint shares its sign.
            if (f*br%slo > 0._wp) then
                br%lo = br%x
            else
                br%hi = br%x
            end if
        end if

        ! Convergence on a collapsed interval. br%x is left where it is: it is the last point evaluated, and it lies inside.
        if (br%hi - br%lo <= br%xtol*(abs(br%lo) + abs(br%hi))) then
            br%state = bracket_done
            return
        end if

        if (br%iter >= br%max_iter) then
            br%state = bracket_capped
            return
        end if

        ! Newton where it stays inside the bracket, bisection where it does not. The derivative is tested rather than divided
        ! by blindly, both because a build with floating-point trapping would stop on the division and because an infinite or
        ! NaN quotient would otherwise have to be caught afterwards. A trial point that is outside, infinite or NaN fails the
        ! test below -- NaN compares false either way -- and falls back to the midpoint.
        if (df == 0._wp) then
            xn = 0.5_wp*(br%lo + br%hi)
        else
            xn = br%x - f/df
            if (.not. (xn > br%lo .and. xn < br%hi)) xn = 0.5_wp*(br%lo + br%hi)
        end if
        br%x = xn

    end subroutine s_bracket_newton_step

end module m_root_finding
