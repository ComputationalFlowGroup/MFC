!>
!! @file
!! @brief Contains module m_graded

#:include 'macros.fpp'
module m_graded

    use m_derived_types
    use m_global_parameters
    use m_bubbles
    use m_variables_conversion
    use m_helper_basic
    use m_constants

    implicit none

    private; public :: s_initialize_graded, s_graded_G !add s_graded_mu later

    ! flat arrays
    logical, allocatable, dimension(:)    :: G_graded_flag
    integer, allocatable, dimension(:)    :: graded_type_arr, graded_profile_arr
    real(wp), allocatable, dimension(:)   :: G_init_arr, G_end_arr
    real(wp), allocatable, dimension(:,:) :: graded_beg_arr, graded_end_arr, graded_center_arr
    real(wp), allocatable, dimension(:)   :: graded_r_beg_arr, graded_r_end_arr
    $:GPU_DECLARE(create='[G_graded_flag, graded_type_arr, graded_profile_arr]')
    $:GPU_DECLARE(create='[G_init_arr, G_end_arr]')
    $:GPU_DECLARE(create='[graded_beg_arr, graded_end_arr, graded_center_arr]')
    $:GPU_DECLARE(create='[graded_r_beg_arr, graded_r_end_arr]')

contains

    !> Initialize the graded module
    impure subroutine s_initialize_graded
        
        integer :: i, j
        
        @:ALLOCATE(G_graded_flag(1:num_fluids), graded_type_arr(1:num_fluids), &
                   & graded_profile_arr(1:num_fluids), G_init_arr(1:num_fluids), &
                   & G_end_arr(1:num_fluids), graded_beg_arr(1:3, 1:num_fluids), &
                   & graded_end_arr(1:3, 1:num_fluids), graded_center_arr(1:3, 1:num_fluids), &
                   & graded_r_beg_arr(1:num_fluids), graded_r_end_arr(1:num_fluids))

        do i = 1, num_fluids
            G_graded_flag(i)      = fluid_pp(i)%graded_G
            graded_type_arr(i)    = fluid_pp(i)%graded_type
            graded_profile_arr(i) = fluid_pp(i)%graded_profile
            G_init_arr(i)         = fluid_pp(i)%graded_G_init
            G_end_arr(i)          = fluid_pp(i)%graded_G_end
            graded_r_beg_arr(i)   = fluid_pp(i)%graded_r_beg
            graded_r_end_arr(i)   = fluid_pp(i)%graded_r_end
            do j = 1, 3
                graded_beg_arr(j, i)    = fluid_pp(i)%graded_beg(j)
                graded_end_arr(j, i)    = fluid_pp(i)%graded_end(j)
                graded_center_arr(j, i) = fluid_pp(i)%graded_center(j)
            end do
        end do

        $:GPU_UPDATE(device='[G_graded_flag, graded_type_arr, graded_profile_arr, &
                     & G_init_arr, G_end_arr, graded_beg_arr, graded_end_arr, &
                     & graded_center_arr, graded_r_beg_arr, graded_r_end_arr]')

    end subroutine s_initialize_graded


    !> Evaluate graded shear modulus at each cell using reference map
    subroutine s_graded_G(xi_x, xi_y, xi_z, fluid_idx, G_out)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)  :: xi_x, xi_y, xi_z
        integer,  intent(in)  :: fluid_idx
        real(wp), intent(out) :: G_out
        real(wp) :: dx, dy, dz, rx, ry, rz, denom, proj, s
        G_init = G_init_arr(fluid_idx)
        G_end = G_end_arr(fluid_idx)
          
        if (graded_type_arr(fluid_idx) == 1) then
                ! cartesian, project onto segement from graded_beg to
                ! graded_end
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
                proj = (sqrt(dx*dx + dy*dy + dz*dz) - graded_r_beg_arr(fluid_idx)) / &
                        & (graded_r_end_arr(fluid_idx) - graded_r_beg_arr(fluid_idx))
                ! proj  = 0 at r_beg, proj  = 1 at r_end
        else
            proj = 0._wp
        end if
      
        if (proj <= 0._wp) then
            G_out = G_init
        else if (proj >= 1._wp) then
            G_out = G_end
        else
            if (graded_profile_arr(fluid_idx) == 1) then
                ! Linear ramp
                ramp = proj
            else if (graded_profile_arr(fluid_idx) == 2) then
                ! Sinusoidal (smooth S-curve, 0 -> 1)
                ! change
            else if (graded_profile_arr(fluid_idx) == 3) then
                ! Power law
                ! change
            else
                ramp = proj
            end if
            G_out = G_init + (G_end - G_init)*ramp
        end if

        G_out = max(G_out, 0._wp)

        end subroutine s_graded_G

end module m_graded
