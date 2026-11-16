!>
!! @file m_graded.fpp
!! @brief Contains module m_graded

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Setup for graded (spatially-varying) elastic/viscous material properties

module m_graded

    use m_derived_types
    use m_global_parameters_common

    implicit none

    private; public :: s_initialize_graded, s_grade_G, s_grade_Re, s_finalize_graded, G_graded_flag

    ! flat arrays
    logical, allocatable, dimension(:)    :: G_graded_flag
    logical, allocatable, dimension(:)    :: Re_graded_flag
    real(wp), allocatable, dimension(:)   :: G_init_arr, G_end_arr
    real(wp), allocatable, dimension(:,:) :: Re_init_arr
    real(wp), allocatable, dimension(:,:) :: Re_end_arr
    integer, allocatable, dimension(:)    :: graded_type_arr, graded_profile_arr
    real(wp), allocatable, dimension(:)   :: graded_k, graded_n, graded_a
    real(wp), allocatable, dimension(:,:) :: graded_beg_arr, graded_end_arr, graded_center_arr
    real(wp), allocatable, dimension(:)   :: graded_r_beg_arr, graded_r_end_arr

    $:GPU_DECLARE(create='[G_graded_flag]')
    $:GPU_DECLARE(create='[Re_graded_flag]')
    $:GPU_DECLARE(create='[G_init_arr, G_end_arr]')
    $:GPU_DECLARE(create='[Re_init_arr, Re_end_arr]')
    $:GPU_DECLARE(create='[graded_type_arr, graded_profile_arr]')
    $:GPU_DECLARE(create='[graded_k, graded_n, graded_a]')
    $:GPU_DECLARE(create='[graded_beg_arr, graded_end_arr, graded_center_arr]')
    $:GPU_DECLARE(create='[graded_r_beg_arr, graded_r_end_arr]')

contains

    !> Initialize the graded module
    impure subroutine s_initialize_graded

        integer :: i

        @:ALLOCATE(G_graded_flag(1:num_fluids), Re_graded_flag(1:num_fluids), graded_type_arr(1:num_fluids), &
                   & graded_profile_arr(1:num_fluids), G_init_arr(1:num_fluids), G_end_arr(1:num_fluids), &
                   & Re_init_arr(1:2,1:num_fluids), Re_end_arr(1:2,1:num_fluids), graded_k(1:num_fluids), graded_n(1:num_fluids), &
                   & graded_a(1:num_fluids), graded_beg_arr(1:3, 1:num_fluids), graded_end_arr(1:3, 1:num_fluids), &
                   & graded_center_arr(1:3, 1:num_fluids), graded_r_beg_arr(1:num_fluids), graded_r_end_arr(1:num_fluids))

        do i = 1, num_fluids
            G_graded_flag(i) = fluid_pp(i)%graded_G
            Re_graded_flag(i) = fluid_pp(i)%graded_Re
            graded_type_arr(i) = fluid_pp(i)%graded_type
            graded_profile_arr(i) = fluid_pp(i)%graded_profile
            G_init_arr(i) = fluid_pp(i)%graded_G_init
            G_end_arr(i) = fluid_pp(i)%graded_G_end
            Re_init_arr(1, i) = fluid_pp(i)%graded_Re_init(1)
            Re_init_arr(2, i) = fluid_pp(i)%graded_Re_init(2)
            Re_end_arr(1, i) = fluid_pp(i)%graded_Re_end(1)
            Re_end_arr(2, i) = fluid_pp(i)%graded_Re_end(2)
            graded_r_beg_arr(i) = fluid_pp(i)%graded_r_beg
            graded_r_end_arr(i) = fluid_pp(i)%graded_r_end
            graded_k(i) = fluid_pp(i)%graded_pf_coeff
            graded_n(i) = fluid_pp(i)%graded_exp
            graded_a(i) = fluid_pp(i)%graded_scaling
            graded_beg_arr(1, i) = fluid_pp(i)%graded_beg_loc_x
            graded_beg_arr(2, i) = fluid_pp(i)%graded_beg_loc_y
            graded_beg_arr(3, i) = fluid_pp(i)%graded_beg_loc_z
            graded_end_arr(1, i) = fluid_pp(i)%graded_end_loc_x
            graded_end_arr(2, i) = fluid_pp(i)%graded_end_loc_y
            graded_end_arr(3, i) = fluid_pp(i)%graded_end_loc_z
            graded_center_arr(1, i) = fluid_pp(i)%graded_center_loc_x
            graded_center_arr(2, i) = fluid_pp(i)%graded_center_loc_y
            graded_center_arr(3, i) = fluid_pp(i)%graded_center_loc_z
            if (G_graded_flag(i)) then
                fluid_pp(i)%G = G_init_arr(i)
            end if
        end do

        $:GPU_UPDATE(device='[G_graded_flag, Re_graded_flag, graded_type_arr, graded_profile_arr, G_init_arr, G_end_arr, &
                     & Re_init_arr, Re_end_arr, graded_beg_arr, graded_end_arr, graded_center_arr, graded_r_beg_arr, &
                     & graded_r_end_arr, graded_n, graded_k, graded_a]')

    end subroutine s_initialize_graded

    !> Helper function for graded projection calculation
    function f_calc_projection(xi_x, xi_y, xi_z, fluid_idx) result(proj)

        real(wp), intent(in) :: xi_x, xi_y, xi_z
        integer, intent(in)  :: fluid_idx
        real(wp)             :: proj
        real(wp)             :: dx, dy, dz
        real(wp)             :: rx, ry, rz
        real(wp)             :: denom

        $:GPU_ROUTINE(parallelism='[seq]')

        if (graded_type_arr(fluid_idx) == 1) then
            ! Cartesian: project onto the segment from graded_beg_loc to graded_end_loc
            dx = graded_end_arr(1, fluid_idx) - graded_beg_arr(1, fluid_idx)
            dy = graded_end_arr(2, fluid_idx) - graded_beg_arr(2, fluid_idx)
            dz = graded_end_arr(3, fluid_idx) - graded_beg_arr(3, fluid_idx)

            rx = xi_x - graded_beg_arr(1, fluid_idx)
            ry = xi_y - graded_beg_arr(2, fluid_idx)
            rz = xi_z - graded_beg_arr(3, fluid_idx)

            denom = dx*dx + dy*dy + dz*dz

            proj = (rx*dx + ry*dy + rz*dz)/denom
            ! proj = 0 at graded_beg_loc, proj = 1 at graded_end_loc
        else if (graded_type_arr(fluid_idx) == 2) then
            ! Radial: normalized distance from graded_center_loc, between graded_r_beg and graded_r_end
            dx = xi_x - graded_center_arr(1, fluid_idx)
            dy = xi_y - graded_center_arr(2, fluid_idx)
            dz = xi_z - graded_center_arr(3, fluid_idx)

            proj = (sqrt(dx*dx + dy*dy + dz*dz) - graded_r_beg_arr(fluid_idx))/(graded_r_end_arr(fluid_idx) &
                    & - graded_r_beg_arr(fluid_idx))
            ! proj = 0 at r_beg, proj = 1 at r_end
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
            ! Sinusoidal ramp
            ramp = 0.5_wp*(1._wp - cos(pi*proj))
        else if (graded_profile_arr(fluid_idx) == 3) then
            ! Power-law ramp: graded_pf_coeff * (graded_scaling * proj)^graded_exp
            ramp = graded_k(fluid_idx)*(graded_a(fluid_idx)*proj)**(graded_n(fluid_idx))
        else
            ramp = proj
        end if

    end function f_calc_profile

    !> Evaluate graded shear modulus at a point, using the reference map as the spatial coordinate
    impure subroutine s_grade_G(xi_x, xi_y, xi_z, fluid_idx, G_out)

        real(wp), intent(in)  :: xi_x, xi_y, xi_z
        integer, intent(in)   :: fluid_idx
        real(wp), intent(out) :: G_out
        real(wp)              :: proj, ramp
        real(wp)              :: G_init, G_end
        real(wp)              :: r, dx, dy, dz

        $:GPU_ROUTINE(parallelism='[seq]')

        G_init = G_init_arr(fluid_idx)
        G_end = G_end_arr(fluid_idx)

        if (graded_profile_arr(fluid_idx) == 3) then
            ! Power law: G = G_init * (r/graded_scaling)^graded_exp
            dx = xi_x - graded_center_arr(1, fluid_idx)
            dy = xi_y - graded_center_arr(2, fluid_idx)
            dz = xi_z - graded_center_arr(3, fluid_idx)
            r = max(sqrt(dx*dx + dy*dy + dz*dz), 1.e-12_wp)
            G_out = G_init*(r/graded_a(fluid_idx))**graded_n(fluid_idx)
        else
            proj = f_calc_projection(xi_x, xi_y, xi_z, fluid_idx)

            if (proj <= 0._wp) then
                G_out = G_init
            else if (proj >= 1._wp) then
                G_out = G_end
            else
                ramp = f_calc_profile(proj, fluid_idx)
                G_out = G_init + (G_end - G_init)*ramp
            end if
        end if

        G_out = max(G_out, 0._wp)

    end subroutine s_grade_G

    !> Evaluate graded viscosity (shear or bulk, selected by j) at a point, using the reference map as the spatial coordinate
    subroutine s_grade_Re(xi_x, xi_y, xi_z, fluid_idx, j, Re_out)

        real(wp), intent(in)  :: xi_x, xi_y, xi_z
        integer, intent(in)   :: fluid_idx
        integer, intent(in)   :: j
        real(wp), intent(out) :: Re_out
        real(wp)              :: proj, ramp
        real(wp)              :: Re_init, Re_end

        $:GPU_ROUTINE(parallelism='[seq]')

        Re_init = Re_init_arr(j, fluid_idx)
        Re_end = Re_end_arr(j, fluid_idx)

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

        @:DEALLOCATE(G_graded_flag, Re_graded_flag, graded_type_arr, graded_profile_arr, G_init_arr, G_end_arr, Re_init_arr, &
                     & Re_end_arr, graded_beg_arr, graded_end_arr, graded_center_arr, graded_r_beg_arr, graded_r_end_arr, &
                     & graded_n, graded_k, graded_a)

    end subroutine s_finalize_graded

end module m_graded
