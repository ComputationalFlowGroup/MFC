!>
!! @file m_graded.fpp
!! @brief Contains module m_graded

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Setup for graded capability

module m_graded

    use m_derived_types
    use m_global_parameters_common

    implicit none

    private; public :: s_initialize_graded, s_grade_Ca_inv, s_grade_Re, s_finalize_graded

    ! flat arrays
    logical, allocatable, dimension(:)    :: Ca_inv_graded_flag, Re_graded_flag
    integer, allocatable, dimension(:)    :: graded_type_arr, graded_profile_arr
    real(wp), allocatable, dimension(:)   :: Ca_inv_init_arr, Ca_inv_end_arr
    real(wp), allocatable, dimension(:)   :: Re_init_arr, Re_end_arr
    real(wp), allocatable, dimension(:)   :: graded_k, graded_n, graded_a
    real(wp), allocatable, dimension(:,:) :: graded_beg_arr, graded_end_arr, graded_center_arr
    real(wp), allocatable, dimension(:)   :: graded_r_beg_arr, graded_r_end_arr

#ifdef MFC_SIMULATION
    $:GPU_DECLARE(create='[Ca_inv_graded_flag, Re_graded_flag]')
    $:GPU_DECLARE(create='[graded_type_arr, graded_profile_arr]')
    $:GPU_DECLARE(create='[graded_k, graded_n, graded_a]')
    $:GPU_DECLARE(create='[Ca_inv_init_arr, Ca_inv_end_arr]')
    $:GPU_DECLARE(create='[Re_init_arr, Re_end_arr]')
    $:GPU_DECLARE(create='[graded_beg_arr, graded_end_arr, graded_center_arr]')
    $:GPU_DECLARE(create='[graded_r_beg_arr, graded_r_end_arr]')
#endif

contains

    !> Initialize the graded module
    impure subroutine s_initialize_graded

        integer :: i, j

        @:ALLOCATE(Ca_inv_graded_flag(1:num_fluids), Re_graded_flag(1:num_fluids), graded_type_arr(1:num_fluids), &
                   & graded_profile_arr(1:num_fluids), Ca_inv_init_arr(1:num_fluids), Ca_inv_end_arr(1:num_fluids), &
                   & Re_init_arr(1:num_fluids), Re_end_arr(1:num_fluids), graded_k(1:num_fluids), graded_n(1:num_fluids), &
                   & graded_a(1:num_fluids), graded_beg_arr(1:3, 1:num_fluids), graded_end_arr(1:3, 1:num_fluids), &
                   & graded_center_arr(1:3, 1:num_fluids), graded_r_beg_arr(1:num_fluids), graded_r_end_arr(1:num_fluids))

        do i = 1, num_fluids
            Ca_inv_graded_flag(i) = fluid_pp(i)%graded_Ca_inv
            Re_graded_flag(i) = fluid_pp(i)%graded_Re
            graded_type_arr(i) = fluid_pp(i)%graded_type
            graded_profile_arr(i) = fluid_pp(i)%graded_profile
            Ca_inv_init_arr(i) = fluid_pp(i)%graded_Ca_inv_init
            Ca_inv_end_arr(i) = fluid_pp(i)%graded_Ca_inv_end
            Re_init_arr(i) = fluid_pp(i)%graded_Re_init
            Re_end_arr(i) = fluid_pp(i)%graded_Re_end
            graded_r_beg_arr(i) = fluid_pp(i)%graded_r_beg
            graded_r_end_arr(i) = fluid_pp(i)%graded_r_end
            graded_k(i) = fluid_pp(i)%graded_pf_coeff
            graded_n(i) = fluid_pp(i)%graded_exp
            graded_a(i) = fluid_pp(i)%graded_scaling
            do j = 1, 3
                graded_beg_arr(j, i) = fluid_pp(i)%graded_beg_loc(j)
                graded_end_arr(j, i) = fluid_pp(i)%graded_end_loc(j)
                graded_center_arr(j, i) = fluid_pp(i)%graded_center_loc(j)
            end do
        end do

        $:GPU_UPDATE(device='[Ca_inv_graded_flag, Re_graded_flag, graded_type_arr, graded_profile_arr, Ca_inv_init_arr, &
                     & Ca_inv_end_arr, Re_init_arr, Re_end_arr, graded_beg_arr, graded_end_arr, graded_center_arr, &
                     & graded_r_beg_arr, graded_r_end_arr, graded_n, graded_k, graded_a]')

    end subroutine s_initialize_graded

    !> Helper function for graded projection calculationn
    function f_calc_projection(xi_x, xi_y, xi_z, fluid_idx) result(proj)

        real(wp), intent(in) :: xi_x, xi_y, xi_z
        integer, intent(in)  :: fluid_idx
        real(wp)             :: proj
        real(wp)             :: dx, dy, dz
        real(wp)             :: rx, ry, rz
        real(wp)             :: denom

        $:GPU_ROUTINE(parallelism='[seq]')

        if (graded_type_arr(fluid_idx) == 1) then
            ! cartesian, project onto segment from graded_beg to graded_end
            dx = graded_end_arr(1, fluid_idx) - graded_beg_arr(1, fluid_idx)
            dy = graded_end_arr(2, fluid_idx) - graded_beg_arr(2, fluid_idx)
            dz = graded_end_arr(3, fluid_idx) - graded_beg_arr(3, fluid_idx)

            rx = xi_x - graded_beg_arr(1, fluid_idx)
            ry = xi_y - graded_beg_arr(2, fluid_idx)
            rz = xi_z - graded_beg_arr(3, fluid_idx)

            denom = dx*dx + dy*dy + dz*dz

            proj = (rx*dx + ry*dy + rz*dz)/denom
            ! proj = 0 at graded_beg, proj  = 1 at graded_end
        else if (graded_type_arr(fluid_idx) == 2) then
            ! Radial: normalized distance from center between r_beg and r_end
            dx = xi_x - graded_center_arr(1, fluid_idx)
            dy = xi_y - graded_center_arr(2, fluid_idx)
            dz = xi_z - graded_center_arr(3, fluid_idx)

            proj = (sqrt(dx*dx + dy*dy + dz*dz) - graded_r_beg_arr(fluid_idx))/(graded_r_end_arr(fluid_idx) &
                    & - graded_r_beg_arr(fluid_idx))
            ! proj  = 0 at r_beg, proj  = 1 at r_end
        else
            proj = 0._wp
        end if

    end function f_calc_projection

    !> Helper function for graded profile calculation
    function f_calc_profile(proj, fluid_idx) result(ramp)

        real(wp), intent(in) :: proj
        integer, intent(in)  :: fluid_idx
        real(wp)             :: ramp

        $:GPU_ROUTINE(parallelism='[seq]')

        if (graded_profile_arr(fluid_idx) == 1) then
            ! Linear ramp
            ramp = proj
        else if (graded_profile_arr(fluid_idx) == 2) then
            ! Sinusoidal - change to have params - should this even be included
            ramp = 0.5_wp*(1._wp - cos(pi*proj))
        else if (graded_profile_arr(fluid_idx) == 3) then
            ! Power law - change for params
            ramp = graded_k(fluid_idx)*(graded_a(fluid_idx)*proj)**(graded_n(fluid_idx))
        else
            ramp = proj
        end if

    end function f_calc_profile

    !> Evaluate graded elasticity at each cell using reference map
    impure subroutine s_grade_Ca_inv(xi_x, xi_y, xi_z, fluid_idx, Ca_inv_out)

        real(wp), intent(in)  :: xi_x, xi_y, xi_z
        integer, intent(in)   :: fluid_idx
        real(wp), intent(out) :: Ca_inv_out
        real(wp)              :: proj
        real(wp)              :: ramp
        real(wp)              :: Ca_inv_init, Ca_inv_end

        $:GPU_ROUTINE(parallelism='[seq]')

        Ca_inv_init = Ca_inv_init_arr(fluid_idx)
        Ca_inv_end = Ca_inv_end_arr(fluid_idx)

        proj = f_calc_projection(xi_x, xi_y, xi_z, fluid_idx)

        if (proj <= 0._wp) then
            Ca_inv_out = Ca_inv_init
        else if (proj >= 1._wp) then
            Ca_inv_out = Ca_inv_end
        else
            ramp = f_calc_profile(proj, fluid_idx)
            Ca_inv_out = Ca_inv_init + (Ca_inv_end - Ca_inv_init)*ramp
        end if

        Ca_inv_out = max(Ca_inv_out, 0._wp)

    end subroutine s_grade_Ca_inv

    !> Evaluate graded viscosity at each cell using reference map
    subroutine s_grade_Re(xi_x, xi_y, xi_z, fluid_idx, Re_out)

        real(wp), intent(in)  :: xi_x, xi_y, xi_z
        integer, intent(in)   :: fluid_idx
        real(wp), intent(out) :: Re_out
        real(wp)              :: proj
        real(wp)              :: ramp
        real(wp)              :: Re_init, Re_end

        $:GPU_ROUTINE(parallelism='[seq]')

        Re_init = Re_init_arr(fluid_idx)
        Re_end = Re_end_arr(fluid_idx)

        proj = f_calc_projection(xi_x, xi_y, xi_z, fluid_idx)

        if (proj <= 0._wp) then
            Re_out = Re_init
        else if (proj >= 1._wp) then
            Re_out = Re_end
        else
            ramp = f_calc_profile(proj, fluid_idx)
            Re_out = Re_init + (Re_end - Re_init)*ramp
        end if

        Re_out = max(Re_out, 0._wp)

    end subroutine s_grade_Re

    !> Deallocate graded buffer arrays allocated during module initialization
    subroutine s_finalize_graded

        @:DEALLOCATE(Ca_inv_graded_flag, Re_graded_flag, graded_type_arr, graded_profile_arr, Ca_inv_init_arr, Ca_inv_end_arr, &
                     & Re_init_arr, Re_end_arr, graded_beg_arr, graded_end_arr, graded_center_arr, graded_r_beg_arr, &
                     & graded_r_end_arr, graded_n, graded_k, graded_a)

    end subroutine s_finalize_graded

end module m_graded
