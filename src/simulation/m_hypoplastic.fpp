!>
!! @file m_hypoplastic.f90
!! @brief Contains module m_hypoplastic

#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoplastic model
module m_hypoplastic

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_helper

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_hypoplastic_module, &
 s_finalize_hypoplastic_module, &
 s_compute_hypoplastic_rhs

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), Gs)
    !$acc declare link(Gs)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), du_dx, du_dy, du_dz)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dv_dx, dv_dy, dv_dz)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), dw_dx, dw_dy, dw_dz)
    !$acc declare link(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :, :), rho_K_field, G_K_field)
    !$acc declare link(rho_K_field, G_K_field)

    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), allocatable, dimension(:, :), fd_coeff_x, fd_coeff_y, fd_coeff_z)
    !$acc declare link(fd_coeff_x,fd_coeff_y,fd_coeff_z)

#else
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    !$acc declare create(Gs)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: du_dx, du_dy, du_dz
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dv_dx, dv_dy, dv_dz
    real(kind(0d0)), allocatable, dimension(:, :, :) :: dw_dx, dw_dy, dw_dz
    !$acc declare create(du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz)

    real(kind(0d0)), allocatable, dimension(:, :, :) :: rho_K_field, G_K_field
    !$acc declare create(rho_K_field, G_K_field)

    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_x
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_y
    real(kind(0d0)), allocatable, dimension(:, :) :: fd_coeff_z
    !$acc declare create(fd_coeff_x,fd_coeff_y,fd_coeff_z)
#endif

contains
   !>   The following subroutine handles the hypoelastic evolution
        !! equation for Johnson-Cook plasticity model, which requires
        !! the Jaumann-Zaremba rate of the Kirchhoff stress deviator
    subroutine s_initialize_hypoplastic_module

        integer :: i, k, r

        @:ALLOCATE_GLOBAL(Gs(1:num_fluids))
        @:ALLOCATE_GLOBAL(rho_K_field(0:m,0:n,0:p), G_K_field(0:m,0:n,0:p))
        @:ALLOCATE_GLOBAL(du_dx(0:m,0:n,0:p))
        if (n > 0) then !2D
            @:ALLOCATE_GLOBAL(du_dy(0:m,0:n,0:p), dv_dx(0:m,0:n,0:p), dv_dy(0:m,0:n,0:p))
            if (p > 0) then !3D
                @:ALLOCATE_GLOBAL(du_dz(0:m,0:n,0:p), dv_dz(0:m,0:n,0:p))
                @:ALLOCATE_GLOBAL(dw_dx(0:m,0:n,0:p), dw_dy(0:m,0:n,0:p), dw_dz(0:m,0:n,0:p))
            end if
        end if

        do i = 1, num_fluids
            Gs(i) = fluid_pp(i)%G
        end do
        !$acc update device(Gs)

        @:ALLOCATE_GLOBAL(fd_coeff_x(-fd_number:fd_number, 0:m))
        if (n > 0) then !2D
            @:ALLOCATE_GLOBAL(fd_coeff_y(-fd_number:fd_number, 0:n))
        end if
        if (p > 0) then !3D
            @:ALLOCATE_GLOBAL(fd_coeff_z(-fd_number:fd_number, 0:p))
        end if

        ! Computing centered finite difference coefficients
        call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, &
                                                      fd_number, fd_order)
        !$acc update device(fd_coeff_x)
        if (n > 0) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_y)
        end if
        if (p > 0) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, &
                                                          fd_number, fd_order)
            !$acc update device(fd_coeff_z)
        end if

    end subroutine s_initialize_hypoplastic_module

    !>  The purpose of this procedure is to compute the source terms
        !!      that are needed for the elastic stress equations
        !!  @param idir Dimension splitting index
        !!  @param q_prim_vf Primitive variables
        !!  @param rhs_vf rhs variables
    subroutine s_compute_hypoplastic_rhs(idir, q_prim_vf, rhs_vf)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf

        real(kind(0d0)) :: rho_K, G_K
	real(kind(0d0)), dimension(num_dim**2) :: atensor, tensora

        integer :: i, k, l, q, r !< Loop variables
        integer :: ndirs  !< Number of coordinate directions
	
        ndirs = 1; if (n > 0) ndirs = 2; !if (p > 0) ndirs = 3

        ! compute velocity gradients and rho_K and G_K        
        !$acc parallel loop collapse(3) gang vector default(present)
        	do q = 0, p
                    do l = 0, n
                        do k = 0, m
			    du_dx(k, l, q) = 0d0;
                            du_dy(k, l, q) = 0d0; dv_dx(k, l, q) = 0d0; dv_dy(k, l, q) = 0d0; 
                        end do
                    end do
                end do
                !$acc end parallel loop

                !$acc parallel loop collapse(3) gang vector default(present)
                do q = 0, p
                    do l = 0, n
                        do k = 0, m
                            !$acc loop seq
                            do r = -fd_number, fd_number
				du_dx(k, l, q) = du_dx(k, l, q) &
						 + q_prim_vf(momxb)%sf(k + r, l, q)*fd_coeff_x(r, k)
                                du_dy(k, l, q) = du_dy(k, l, q) &
                                                 + q_prim_vf(momxb)%sf(k, l + r, q)*fd_coeff_y(r, l)
                                dv_dx(k, l, q) = dv_dx(k, l, q) &
                                                 + q_prim_vf(momxb + 1)%sf(k + r, l, q)*fd_coeff_x(r, k)
                                dv_dy(k, l, q) = dv_dy(k, l, q) &
                                                 + q_prim_vf(momxb + 1)%sf(k, l + r, q)*fd_coeff_y(r, l)
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop

            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
                        rho_K = 0d0; G_K = 0d0
                        do i = 1, num_fluids
                            rho_K = rho_K + q_prim_vf(i)%sf(k, l, q) !alpha_rho_K(1)
                            G_K = G_K + q_prim_vf(advxb - 1 + i)%sf(k, l, q)*Gs(i)  !alpha_K(1) * Gs(1)
                        end do
                        rho_K_field(k, l, q) = rho_K
                        G_K_field(k, l, q) = G_K

                        !TODO: take this out if not needed
                        if (G_K < verysmall) then
                            G_K_field(k, l, q) = 0
                        end if
                    end do
                end do
            end do

	    ! Compute the first additional term in rhs: -rho((SW)^T - SW)
	    ! Let atensor = SW, attensor = (SW)^T, tensora = attensor - atensor
	    do q = 0, p
		do l = 0, n
		   do k = 0, m
			atensor(1) = (1d0/4d0)*((dv_dx(k, l, q)**2 - du_dy(k, l, q)**2)
			atensor(2) = (1d0/2d0)*(du_dy(k, l, q)*du_dx(k, l, q) - &
				     dv_dx(k, l, q)*du_dx(k, l, q))
			atensor(3) = (1d0/2d0)*(dv_dx(k, l, q)*dv_dy(k, l, q) - &
				     du_dy(k, l, q)*dv_dy(k, l, q))
			atensor(4) = (1d0/4d0)*((du_dy(k, l, q)**2 - dv_dx(k, l, q)**2)
		        tensora(1) = 0d0
			tensora(2) = atensor(3) - atensor(2)
			tensora(3) = atensor(2) - atensor(3)
			tensora(4) = 0d0
		   end do
	        end do
	    end do

	    ! compute rhs source terms
            !$acc parallel loop collapse(3) gang vector default(present)
            do q = 0, p
                do l = 0, n
                    do k = 0, m
			! TODO: MISSING TERM 2 EVERYWHERE
			rhs_vf(strxb)%sf(k, l, q) = rhs_vf(strxb)%sf(k, l, q) + rho_K_field(k, l ,q)*tensora(1) 
                       
                        rhs_vf(strxb + 1)%sf(k, l, q) = rhs_vf(strxb + 1)%sf(k, l, q) + rho_K_field(k, l, q)* tensora(2)
                                                        
                        rhs_vf(strxb + 2)%sf(k, l, q) = rhs_vf(strxb + 2)%sf(k, l, q) + rho_K_field(k, l, q)*tensora(3)
                	
			rhs_vf(strxb + 3)%sf(k, l, q) = rhs_vf(strxb + 3)%sf(k, l, q) + rho_K_field(k, l, q)*tensora(4)
                    end do
                end do
            end do

          end subroutine s_compute_hypoplastic_rhs

    subroutine s_finalize_hypoplastic_module() ! --------------------

        @:DEALLOCATE_GLOBAL(Gs)
        @:DEALLOCATE_GLOBAL(rho_K_field, G_K_field)
        @:DEALLOCATE_GLOBAL(du_dx)
        @:DEALLOCATE_GLOBAL(fd_coeff_x)
        if (n > 0) then
            @:DEALLOCATE_GLOBAL(du_dy,dv_dx,dv_dy)
            @:DEALLOCATE_GLOBAL(fd_coeff_y)
            if (p > 0) then
                @:DEALLOCATE_GLOBAL(du_dz, dv_dz, dw_dx, dw_dy, dw_dz)
                @:DEALLOCATE_GLOBAL(fd_coeff_z)
            end if
        end if

    end subroutine s_finalize_hypoplastic_module

end module m_hypoplastic
