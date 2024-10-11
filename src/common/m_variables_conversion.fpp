!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'inline_conversions.fpp'
#:include 'case.fpp'

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive ones and vice versa. In
!!              addition, the module also contains the subroutines used to obtain
!!              the mixture variables and the subroutines used to compute
!!              pressure and temperature.
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    ! ==========================================================================

    implicit none

    private; 
    public :: s_initialize_variables_conversion_module, &
              s_initialize_pb, &
              s_initialize_mv, &
              s_convert_to_mixture_variables, &
              s_convert_mixture_to_mixture_variables, &
              s_convert_species_to_mixture_variables_bubbles, &
              s_convert_species_to_mixture_variables_bubbles_acc, &
              s_convert_species_to_mixture_variables, &
              s_convert_species_to_mixture_variables_acc, &
              s_convert_conservative_to_primitive_variables, &
              s_convert_primitive_to_conservative_variables, &
              s_convert_primitive_to_flux_variables, &
              s_compute_pressure, &
              s_compute_temperature, &
              s_finalize_variables_conversion_module

    !> Abstract interface to two subroutines designed for the transfer/conversion
    !! of the mixture/species variables to the mixture variables

    abstract interface ! =======================================================

        !> Structure of the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines
        !!  @param q_vf Conservative or primitive variables
        !!  @param i First-coordinate cell index
        !!  @param j First-coordinate cell index
        !!  @param k First-coordinate cell index
        !!  @param rho Density
        !!  @param gamma Specific heat ratio function
        !!  @param pi_inf Liquid stiffness function
        !!  @param qv Fluid reference energy
        subroutine s_convert_xxxxx_to_mixture_variables(q_vf, i, j, k, &
                                                        rho, gamma, pi_inf, qv, Re_K, G_K, G, jcook_K, jcook)

            ! Importing the derived type scalar_field from m_derived_types.f90
            ! and global variable sys_size, from m_global_variables.f90, as
            ! the abstract interface does not inherently have access to them
            import :: scalar_field, sys_size, num_fluids

            type(scalar_field), dimension(sys_size), intent(in) :: q_vf
            integer, intent(in) :: i, j, k
            real(kind(0d0)), intent(out), target :: rho, gamma, pi_inf, qv
            real(kind(0d0)), optional, dimension(2), intent(out) :: Re_K
            real(kind(0d0)), optional, intent(out) :: G_K
            real(kind(0d0)), optional, dimension(num_fluids), intent(in) :: G
            real(kind(0d0)), optional, dimension(10), intent(out) :: jcook_K
            real(kind(0d0)), optional, dimension(10), intent(in) :: jcook

        end subroutine s_convert_xxxxx_to_mixture_variables

    end interface ! ============================================================

    integer, public :: ixb, ixe, iyb, iye, izb, ize
    !$acc declare create(ixb, ixe, iyb, iye, izb, ize)

    ! In simulation, gammas, pi_infs, and qvs are already declared in m_global_variables
#ifndef MFC_SIMULATION
    real(kind(0d0)), allocatable, public, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    !$acc declare create(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps)
    real(kind(0d0)), allocatable, public, dimension(:) :: rho0, mg_a, mg_b, ein_cv1, ein_cv2
    !$acc declare create(rho0, mg_a, mg_b, ein_cv1, ein_cv2)
    #:for VAR in range(1,12)
      real(kind(0d0)), allocatable, public, dimension(:) :: jcook${VAR}$
    #:endfor
    !$acc declare create(jcook1,jcook2,jcook3,jcook4,jcook5,jcook6,jcook7,jcook8,jcook9,jcook10,jcook11)
#endif

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:), Gs)
    @:CRAY_DECLARE_GLOBAL(integer,         dimension(:), bubrs)
    @:CRAY_DECLARE_GLOBAL(real(kind(0d0)), dimension(:, :), Res)
    !$acc declare link(bubrs, Gs, Res)
#else
    real(kind(0d0)), allocatable, dimension(:) :: Gs
    integer, allocatable, dimension(:) :: bubrs
    real(kind(0d0)), allocatable, dimension(:, :) :: Res
    !$acc declare create(bubrs, Gs, Res)
#endif
    integer :: is1b, is2b, is3b, is1e, is2e, is3e
    !$acc declare create(is1b, is2b, is3b, is1e, is2e, is3e)

    real(kind(0d0)), allocatable, dimension(:, :, :), public :: rho_sf !< Scalar density function
    real(kind(0d0)), allocatable, dimension(:, :, :), public :: gamma_sf !< Scalar sp. heat ratio function
    real(kind(0d0)), allocatable, dimension(:, :, :), public :: pi_inf_sf !< Scalar liquid stiffness function
    real(kind(0d0)), allocatable, dimension(:, :, :), public :: qv_sf !< Scalar liquid energy reference function

    procedure(s_convert_xxxxx_to_mixture_variables), &
        pointer :: s_convert_to_mixture_variables => null() !<
    !! Pointer referencing the subroutine s_convert_mixture_to_mixture_variables
    !! or s_convert_species_to_mixture_variables, based on model equations choice

contains

    !>  This procedure conditionally calculates the appropriate pressure
        !! @param energy Energy
        !! @param alf Void Fraction
        !! @param dyn_p Dynamic Pressure
        !! @param pi_inf Liquid Stiffness
        !! @param gamma Specific Heat Ratio
        !! @param rho Density
        !! @param qv fluid reference energy
        !! @param pres Pressure to calculate
        !! @param stress Shear Stress
        !! @param mom Momentum
        !! @param G shear modulus
        !! @param alpha_K volume fraction of mixture
        !! @param alpha_rho_K conservative volume fraction of mixture
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, &
        rho, qv, pres, stress, mom, G, pref_over_gamma, rho_eref)
        !$acc routine seq

        real(kind(0d0)), intent(in) :: energy, alf
        real(kind(0d0)), intent(in) :: dyn_p
        real(kind(0d0)), intent(in) :: pi_inf, gamma, rho, qv
        real(kind(0d0)), intent(out) :: pres
        real(kind(0d0)), intent(in), optional :: stress, mom, G
        real(kind(0d0)), intent(in), optional :: pref_over_gamma, rho_eref

        real(kind(0d0)) :: E_e, pres_sg
        integer :: s !< Generic loop iterator

        ! Depending on model_eqns and bubbles, the appropriate procedure
        ! for computing pressure is targeted by the procedure pointer
        ! model_eqns = 5 corresponds to the Mie-Gruneisen EOS

        if ((model_eqns /= 4 .and. model_eqns /=5) .and. (bubbles .neqv. .true.)) then
            pres = (energy - dyn_p - pi_inf - qv)/gamma
        else if ((model_eqns /= 4 .and. model_eqns /=5) .and. bubbles) then
            pres = ((energy - dyn_p)/(1.d0 - alf) - pi_inf - qv)/gamma
        else if (model_eqns == 5) then
            pres = (energy - dyn_p + pref_over_gamma - rho_eref)/gamma
            !if (pres /= pres) then
            !    print *,energy, dyn_p, pref_over_gamma, rho_eref, gamma
            !end if

        else
            pres = (pref + pi_inf)* &
                   (energy/ &
                    (rhoref*(1.d0 - alf)) &
                    )**(1.d0/gamma + 1.d0) - pi_inf        
        end if

        if (hypoelasticity .and. present(G)) then
            ! calculate elastic contribution to Energy
            E_e = 0d0
            do s = stress_idx%beg, stress_idx%end
                if (G > 0) then
                    E_e = E_e + ((stress/rho)**2d0)/(4d0*G)
                    ! Additional terms in 2D and 3D
                    if ((s == stress_idx%beg + 1) .or. &
                        (s == stress_idx%beg + 3) .or. &
                        (s == stress_idx%beg + 4)) then
                        E_e = E_e + ((stress/rho)**2d0)/(4d0*G)
                    end if
                end if
            end do

            pres = (energy - dyn_p - pi_inf - qv - E_e)/gamma
        end if

    end subroutine s_compute_pressure

    !>  This procedure conditionally calculates the appropriate temperature of the mixture
        !! @param energy Energy
        !! @param rho Density
        !! @param pres Pressure to calculate
        !! @param alpha_K volume fraction of mixture
    subroutine s_compute_temperature(pres, Pref, gamma_inv,&
        rho, temp, alpha_K)

        !$acc routine seq
        real(kind(0d0)), intent(in) :: pres, Pref, gamma_inv, rho
        real(kind(0d0)), intent(out) :: temp
        real(kind(0d0)), dimension(num_fluids), intent(in), optional :: alpha_K
      
        ! Temporary local variables
        real(kind(0d0)) :: log_rho_mix, rho0_mix, phi, theta_E_mix
        real(kind(0d0)) :: denom, mg_a_mix, mg_b_mix, A_cv, gamma
        integer :: i !< Generic loop iterator

        ! model_eqns = 5 corresponds to the Mie-Gruneisen EOS

        if (model_eqns .eq. 5) then
            A_cv = 0d0
            theta_E_mix = 0d0
            mg_a_mix    = 0d0
            mg_b_mix    = 0d0
            rho0_mix    = 0d0
            gamma       = 0d0
            do i=1, num_fluids
               mg_b_mix = mg_b_mix  +   alpha_K(i)*mg_b(i)
               mg_a_mix = mg_a_mix  +   alpha_K(i)*mg_a(i)
               rho0_mix = rho0_mix  +   alpha_K(i)*rho0(i)
               A_cv     = A_cv + alpha_K(i)*ein_cv1(i)
               theta_E_mix = theta_E_mix + alpha_K(i)*ein_cv2(i)
               gamma    = gamma + alpha_K(i)*gammas(i)
            end do
            log_rho_mix = dlog(rho/rho0_mix)
            phi = ((rho0_mix/rho)**(-mg_a_mix))*dexp((gamma-mg_a_mix)*(1d0-(rho0_mix/rho)**mg_b_mix))
            denom = gamma_inv*(pres-pref)/(rho*phi*A_cv*theta_E_mix)+&
                1d0/(dexp(phi*theta_E_mix)-1d0)
            temp = (phi*theta_E_mix)/dlog(1d0 + 1d0/denom)
        end if
    end subroutine s_compute_temperature

    !>  This subroutine is designed for the gamma/pi_inf model
        !!      and provided a set of either conservative or primitive
        !!      variables, transfers the density, specific heat ratio
        !!      function and the liquid stiffness function from q_vf to
        !!      rho, gamma and pi_inf.
        !! @param q_vf conservative or primitive variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        !! @param rho density
        !! @param gamma  specific heat ratio function
        !! @param pi_inf liquid stiffness
        !! @param qv fluid reference energy
    subroutine s_convert_mixture_to_mixture_variables(q_vf, i, j, k, &
                                                      rho, gamma, pi_inf, qv, Re_K, G_K, G, jcook_K, jcook)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in) :: i, j, k

        real(kind(0d0)), intent(out), target :: rho
        real(kind(0d0)), intent(out), target :: gamma
        real(kind(0d0)), intent(out), target :: pi_inf
        real(kind(0d0)), intent(out), target :: qv

        real(kind(0d0)), optional, dimension(2), intent(out) :: Re_K
        real(kind(0d0)), optional, intent(out) :: G_K
        real(kind(0d0)), optional, dimension(num_fluids), intent(in) :: G
        real(kind(0d0)), optional, dimension(10), intent(out) :: jcook_K
        real(kind(0d0)), optional, dimension(10), intent(in) :: jcook

        ! Transferring the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        rho = q_vf(1)%sf(i, j, k)
        gamma = q_vf(gamma_idx)%sf(i, j, k)
        pi_inf = q_vf(pi_inf_idx)%sf(i, j, k)
        qv = 0d0 ! keep this value nill for now. For future adjustment

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(i, j, k) = rho
        gamma_sf(i, j, k) = gamma
        pi_inf_sf(i, j, k) = pi_inf
        qv_sf(i, j, k) = qv
#endif

    end subroutine s_convert_mixture_to_mixture_variables

    !>  This procedure is used alongside with the gamma/pi_inf
        !!      model to transfer the density, the specific heat ratio
        !!      function and liquid stiffness function from the vector
        !!      of conservative or primitive variables to their scalar
        !!      counterparts. Specifically designed for when subgrid bubbles
        !!      must be included.
        !! @param q_vf primitive variables
        !! @param j Cell index
        !! @param k Cell index
        !! @param l Cell index
        !! @param rho density
        !! @param gamma specific heat ratio
        !! @param pi_inf liquid stiffness
        !! @param qv fluid reference energy
    subroutine s_convert_species_to_mixture_variables_bubbles(q_vf, j, k, l, &
                                                              rho, gamma, pi_inf, qv, Re_K, G_K, G, jcook_K, jcook)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf

        integer, intent(in) :: j, k, l

        real(kind(0d0)), intent(out), target :: rho
        real(kind(0d0)), intent(out), target :: gamma
        real(kind(0d0)), intent(out), target :: pi_inf
        real(kind(0d0)), intent(out), target :: qv

        real(kind(0d0)), optional, dimension(2), intent(out) :: Re_K
        real(kind(0d0)), optional, intent(out) :: G_K
        real(kind(0d0)), optional, dimension(num_fluids), intent(in) :: G
        real(kind(0d0)), optional, dimension(10), intent(out) :: jcook_K
        real(kind(0d0)), optional, dimension(10), intent(in) :: jcook

        integer :: i, q
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        do i = 1, num_fluids
            alpha_rho_K(i) = q_vf(i)%sf(j, k, l)
            alpha_K(i) = q_vf(advxb + i - 1)%sf(j, k, l)
        end do

        if (mpp_lim) then

            do i = 1, num_fluids
                alpha_rho_K(i) = max(0d0, alpha_rho_K(i))
                alpha_K(i) = min(max(0d0, alpha_K(i)), 1d0)
            end do

            alpha_K = alpha_K/max(sum(alpha_K), 1d-16)

        end if

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively

        if (model_eqns == 4) then
            rho = q_vf(1)%sf(j, k, l)
            gamma = fluid_pp(1)%gamma    !qK_vf(gamma_idx)%sf(i,j,k)
            pi_inf = fluid_pp(1)%pi_inf   !qK_vf(pi_inf_idx)%sf(i,j,k)
            qv = fluid_pp(1)%qv
        else if ((model_eqns == 2) .and. bubbles) then
            rho = 0d0; gamma = 0d0; pi_inf = 0d0; qv = 0d0

            if (mpp_lim .and. (num_fluids > 2)) then
                do i = 1, num_fluids
                    rho = rho + q_vf(i)%sf(j, k, l)
                    gamma = gamma + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf = pi_inf + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                    qv = qv + q_vf(i)%sf(j, k, l)*fluid_pp(i)%qv
                end do
            else if (num_fluids == 2) then
                rho = q_vf(1)%sf(j, k, l)
                gamma = fluid_pp(1)%gamma
                pi_inf = fluid_pp(1)%pi_inf
                qv = fluid_pp(1)%qv
            else if (num_fluids > 2) then
                !TODO: This may need fixing for hypo + bubbles
                do i = 1, num_fluids - 1 !leave out bubble part of mixture
                    rho = rho + q_vf(i)%sf(j, k, l)
                    gamma = gamma + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf = pi_inf + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                    qv = qv + q_vf(i)%sf(j, k, l)*fluid_pp(i)%qv
                end do
                ! rho    = qK_vf(1)%sf(j,k,l)
                ! gamma_K  = fluid_pp(1)%gamma
                ! pi_inf_K = fluid_pp(1)%pi_inf
            else
                rho = q_vf(1)%sf(j, k, l)
                gamma = fluid_pp(1)%gamma
                pi_inf = fluid_pp(1)%pi_inf
                qv = fluid_pp(1)%qv
            end if
        end if

#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        if (any(Re_size > 0)) then
            if (num_fluids == 1) then ! need to consider case with num_fluids >= 2
                do i = 1, 2

                    Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0d0

                    do q = 1, Re_size(i)
                        Re_K(i) = (1 - alpha_K(Re_idx(i, q)))/fluid_pp(Re_idx(i, q))%Re(i) &
                                  + Re_K(i)
                    end do

                    Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

                end do
            end if
        end if
#endif

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(j, k, l) = rho
        gamma_sf(j, k, l) = gamma
        pi_inf_sf(j, k, l) = pi_inf
        qv_sf(j, k, l) = qv
#endif

    end subroutine s_convert_species_to_mixture_variables_bubbles

    !>  This subroutine is designed for the volume fraction model
        !!              and provided a set of either conservative or primitive
        !!              variables, computes the density, the specific heat ratio
        !!              function and the liquid stiffness function from q_vf and
        !!              stores the results into rho, gamma and pi_inf.
        !! @param q_vf primitive variables
        !! @param k Cell index
        !! @param l Cell index
        !! @param r Cell index
        !! @param rho density
        !! @param gamma specific heat ratio
        !! @param pi_inf liquid stiffness
        !! @param qv fluid reference energy
    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, r, rho, &
                                                      gamma, pi_inf, qv, Re_K, G_K, G, jcook_K, jcook)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf

        integer, intent(in) :: k, l, r

        real(kind(0d0)), intent(out), target :: rho
        real(kind(0d0)), intent(out), target :: gamma
        real(kind(0d0)), intent(out), target :: pi_inf
        real(kind(0d0)), intent(out), target :: qv

        real(kind(0d0)), optional, dimension(2), intent(out) :: Re_K
            !! Partial densities and volume fractions
        real(kind(0d0)), optional, intent(out) :: G_K
        real(kind(0d0)), optional, dimension(num_fluids), intent(in) :: G
        real(kind(0d0)), optional, dimension(10), intent(out) :: jcook_K
        real(kind(0d0)), optional, dimension(10), intent(in) :: jcook

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K !<

        integer :: i, j !< Generic loop iterator

        ! Computing the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        do i = 1, num_fluids
            alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
            alpha_K(i) = q_vf(advxb + i - 1)%sf(k, l, r)
        end do

        if (mpp_lim) then

            do i = 1, num_fluids
                alpha_rho_K(i) = max(0d0, alpha_rho_K(i))
                alpha_K(i) = min(max(0d0, alpha_K(i)), 1d0)
            end do

            alpha_K = alpha_K/max(sum(alpha_K), 1d-16)

        end if

        ! Calculating the density, the specific heat ratio function, the
        ! liquid stiffness function, and the energy reference function,
        ! respectively, from the species analogs
        rho = 0d0; gamma = 0d0; pi_inf = 0d0; qv = 0d0

        do i = 1, num_fluids
            rho = rho + alpha_rho_K(i)
            gamma = gamma + alpha_K(i)*gammas(i)
            pi_inf = pi_inf + alpha_K(i)*pi_infs(i)
            qv = qv + alpha_rho_K(i)*qvs(i)
        end do
#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        do i = 1, 2

            Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0d0

            do j = 1, Re_size(i)
                Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) &
                          + Re_K(i)
            end do

            Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

        end do
#endif

        if (present(G_K)) then
            !TODO Check our mixture rule? Replace with Cauchy numbers, make code nondimensional
            G_K = 0d0
            do i = 1, num_fluids
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0d0, G_K)
        end if

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(k, l, r) = rho
        gamma_sf(k, l, r) = gamma
        pi_inf_sf(k, l, r) = pi_inf
        qv_sf(k, l, r) = qv
#endif

    end subroutine s_convert_species_to_mixture_variables

    subroutine s_convert_species_to_mixture_variables_acc(rho_K, &
                                                          gamma_K, pi_inf_K, qv_K, &
                                                          alpha_K, alpha_rho_K, Re_K, k, l, r, &
                                                          G_K, G, jcook_K, jcook)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_convert_species_to_mixture_variables_acc
#else
        !$acc routine seq
#endif

        real(kind(0d0)), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K

        real(kind(0d0)), dimension(num_fluids), intent(inout) :: alpha_rho_K, alpha_K !<
        real(kind(0d0)), dimension(2), intent(out) :: Re_K
        !! Partial densities and volume fractions

        real(kind(0d0)), optional, intent(out) :: G_K
        real(kind(0d0)), optional, dimension(num_fluids), intent(in) :: G
        real(kind(0d0)), optional, dimension(10), intent(out) :: jcook_K
        real(kind(0d0)), optional, dimension(10), intent(in) :: jcook

        integer, intent(in) :: k, l, r

        integer :: i, j !< Generic loop iterators
        real(kind(0d0)) :: alpha_K_sum

#ifdef MFC_SIMULATION
        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        rho_K = 0d0
        gamma_K = 0d0
        pi_inf_K = 0d0
        qv_K = 0d0

        alpha_K_sum = 0d0

        if (mpp_lim) then
            do i = 1, num_fluids
                alpha_rho_K(i) = max(0d0, alpha_rho_K(i))
                alpha_K(i) = min(max(0d0, alpha_K(i)), 1d0)
                alpha_K_sum = alpha_K_sum + alpha_K(i)
            end do

            alpha_K = alpha_K/max(alpha_K_sum, sgm_eps)

        end if

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            qv_K = qv_K + alpha_rho_K(i)*qvs(i)
        end do

        if (present(G_K)) then
            G_K = 0d0
            do i = 1, num_fluids
                !TODO: change to use Gs directly here?
                !TODO: Make this changes as well for GPUs
                G_K = G_K + alpha_K(i)*G(i)
            end do
            G_K = max(0d0, G_K)
        end if

        if (any(Re_size > 0)) then

            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0d0

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res(i, j) &
                              + Re_K(i)
                end do

                Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

            end do
        end if
#endif

    end subroutine s_convert_species_to_mixture_variables_acc

    subroutine s_convert_species_to_mixture_variables_bubbles_acc(rho_K, &
                                                                  gamma_K, pi_inf_K, qv_K, &
                                                                  alpha_K, alpha_rho_K, Re_K, k, l, r)
#ifdef CRAY_ACC_WAR
        !DIR$ INLINEALWAYS s_convert_species_to_mixture_variables_bubbles_acc
#else
        !$acc routine seq
#endif

        real(kind(0d0)), intent(inout) :: rho_K, gamma_K, pi_inf_K, qv_K

        real(kind(0d0)), dimension(num_fluids), intent(in) :: alpha_K, alpha_rho_K !<
            !! Partial densities and volume fractions

        real(kind(0d0)), dimension(2), intent(out) :: Re_K
        integer, intent(in) :: k, l, r

        integer :: i, j !< Generic loop iterators

#ifdef MFC_SIMULATION
        rho_K = 0d0
        gamma_K = 0d0
        pi_inf_K = 0d0
        qv_K = 0d0

        if (mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
            do i = 1, num_fluids
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
                qv_K = qv_K + alpha_rho_K(i)*qvs(i)
            end do
        else if ((model_eqns == 2) .and. (num_fluids > 2)) then
            do i = 1, num_fluids - 1
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
                qv_K = qv_K + alpha_rho_K(i)*qvs(i)
            end do
        else
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
            qv_K = qvs(1)
        end if

        if (any(Re_size > 0)) then
            if (num_fluids == 1) then ! need to consider case with num_fluids >= 2

                do i = 1, 2
                    Re_K(i) = dflt_real

                    if (Re_size(i) > 0) Re_K(i) = 0d0

                    do j = 1, Re_size(i)
                        Re_K(i) = (1d0 - alpha_K(Re_idx(i, j)))/Res(i, j) &
                                  + Re_K(i)
                    end do

                    Re_K(i) = 1d0/max(Re_K(i), sgm_eps)

                end do
            end if
        end if
#endif

    end subroutine s_convert_species_to_mixture_variables_bubbles_acc

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_variables_conversion_module

        integer :: i, j

#ifdef MFC_PRE_PROCESS
        ixb = 0; iyb = 0; izb = 0; 
        ixe = m; iye = n; ize = p; 
#else
        ixb = -buff_size
        ixe = m - ixb

        iyb = 0; iye = 0; izb = 0; ize = 0; 
        if (n > 0) then
            iyb = -buff_size; iye = n - iyb

            if (p > 0) then
                izb = -buff_size; ize = p - izb
            end if
        end if
#endif

!$acc enter data copyin(ixb, ixe, iyb, iye, izb, ize)
!$acc enter data copyin(is1b, is1e, is2b, is2e, is3b, is3e)
!$acc update device(ixb, ixe, iyb, iye, izb, ize)

#ifdef MFC_SIMULATION
        @:ALLOCATE_GLOBAL(gammas (1:num_fluids))
        @:ALLOCATE_GLOBAL(gs_min (1:num_fluids))
        @:ALLOCATE_GLOBAL(pi_infs(1:num_fluids))
        @:ALLOCATE_GLOBAL(ps_inf (1:num_fluids))
        @:ALLOCATE_GLOBAL(cvs    (1:num_fluids))
        @:ALLOCATE_GLOBAL(qvs    (1:num_fluids))
        @:ALLOCATE_GLOBAL(qvps   (1:num_fluids))
        @:ALLOCATE_GLOBAL(Gs     (1:num_fluids))
#else
        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(gs_min (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(ps_inf (1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps   (1:num_fluids))
        @:ALLOCATE(Gs     (1:num_fluids))
#endif

        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            gs_min(i) = 1.0d0/gammas(i) + 1.0d0
            pi_infs(i) = fluid_pp(i)%pi_inf
            Gs(i) = fluid_pp(i)%G
            ps_inf(i) = pi_infs(i)/(1.0d0 + gammas(i))
            cvs(i) = fluid_pp(i)%cv
            qvs(i) = fluid_pp(i)%qv
            qvps(i) = fluid_pp(i)%qvp
        end do
!$acc update device(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs)

        if (model_eqns == 5) then 
#ifdef MFC_SIMULATION
        @:ALLOCATE_GLOBAL(rho0   (1:num_fluids))
        @:ALLOCATE_GLOBAL(mg_a   (1:num_fluids))
        @:ALLOCATE_GLOBAL(mg_b   (1:num_fluids))
        @:ALLOCATE_GLOBAL(ein_cv1(1:num_fluids))
        @:ALLOCATE_GLOBAL(ein_cv2(1:num_fluids))
        #:for VAR in range(1,12)
          @:ALLOCATE_GLOBAL(jcook${VAR}$(1:num_fluids))          
        #:endfor 
#else
        @:ALLOCATE(rho0   (1:num_fluids))
        @:ALLOCATE(mg_a   (1:num_fluids))
        @:ALLOCATE(mg_b   (1:num_fluids))
        @:ALLOCATE(ein_cv1(1:num_fluids))
        @:ALLOCATE(ein_cv2(1:num_fluids))
        #:for VAR in range(1,12)
          @:ALLOCATE(jcook${VAR}$(1:num_fluids))          
        #:endfor 
#endif
        do i = 1, num_fluids
            rho0(i) = fluid_pp(i)%rho0  
            mg_a(i) = fluid_pp(i)%mg_a
            mg_b(i) = fluid_pp(i)%mg_b
            ein_cv1(i) = fluid_pp(i)%ein_cv(1)
            ein_cv2(i) = fluid_pp(i)%ein_cv(2)
        end do
!$acc update device(rho0, mg_a, mg_b, ein_cv1, ein_cv2)
        end if
 
        if (hypoplasticity) then    
          do i = 1, num_fluids
            #:for VAR in range(1,12)
                jcook${VAR}$(i) = fluid_pp(i)%jcook(${VAR}$) 
            #:endfor
          end do
!$acc update device(jcook1,jcook2,jcook3,jcook4,jcook5,jcook6,jcook7,jcook8,jcook9,jcook10,jcook11)
        end if

#ifdef MFC_SIMULATION
        if (any(Re_size > 0)) then
            @:ALLOCATE_GLOBAL(Res(1:2, 1:maxval(Re_size)))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            !$acc update device(Res, Re_idx, Re_size)
        end if
#endif

        if (bubbles) then
#ifdef MFC_SIMULATION
            @:ALLOCATE_GLOBAL(bubrs(1:nb))
#else
            @:ALLOCATE(bubrs(1:nb))
#endif

            do i = 1, nb
                bubrs(i) = bub_idx%rs(i)
            end do
            !$acc update device(bubrs)
        end if

#ifdef MFC_POST_PROCESS
        ! Allocating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then

            ! Simulation is 3D
            if (p > 0) then

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 -buff_size:p + buff_size))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   -buff_size:p + buff_size))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    -buff_size:p + buff_size))
                allocate (qv_sf(-buff_size:m + buff_size, &
                                -buff_size:n + buff_size, &
                                -buff_size:p + buff_size))

                ! Simulation is 2D
            else

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 0:0))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    0:0))
                allocate (qv_sf(-buff_size:m + buff_size, &
                                -buff_size:n + buff_size, &
                                0:0))
            end if

            ! Simulation is 1D
        else

            allocate (rho_sf(-buff_size:m + buff_size, &
                             0:0, &
                             0:0))
            allocate (gamma_sf(-buff_size:m + buff_size, &
                               0:0, &
                               0:0))
            allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                0:0, &
                                0:0))
            allocate (qv_sf(-buff_size:m + buff_size, &
                            0:0, &
                            0:0))

        end if
#endif

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables

        else if (bubbles) then
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles

        else
            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if
    end subroutine s_initialize_variables_conversion_module

    !Initialize mv at the quadrature nodes based on the initialized moments and sigma
    subroutine s_initialize_mv(qK_cons_vf, mv)
        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf
        real(kind(0d0)), dimension(ixb:, iyb:, izb:, 1:, 1:), intent(inout) :: mv

        integer :: i, j, k, l
        real(kind(0d0)) :: mu, sig, nbub_sc

        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                    !$acc loop seq
                    do i = 1, nb
                        mu = qK_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5

                        mv(j, k, l, 1, i) = (mass_v0(i))*(mu - sig)**(3d0)/(R0(i)**(3d0))
                        mv(j, k, l, 2, i) = (mass_v0(i))*(mu - sig)**(3d0)/(R0(i)**(3d0))
                        mv(j, k, l, 3, i) = (mass_v0(i))*(mu + sig)**(3d0)/(R0(i)**(3d0))
                        mv(j, k, l, 4, i) = (mass_v0(i))*(mu + sig)**(3d0)/(R0(i)**(3d0))
                    end do

                end do
            end do
        end do

    end subroutine s_initialize_mv

    !Initialize pb at the quadrature nodes using isothermal relations (Preston model)
    subroutine s_initialize_pb(qK_cons_vf, mv, pb)
        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf
        real(kind(0d0)), dimension(ixb:, iyb:, izb:, 1:, 1:), intent(in) :: mv
        real(kind(0d0)), dimension(ixb:, iyb:, izb:, 1:, 1:), intent(inout) :: pb

        integer :: i, j, k, l
        real(kind(0d0)) :: mu, sig, nbub_sc

        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                    !$acc loop seq
                    do i = 1, nb
                        mu = qK_cons_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc
                        sig = (qK_cons_vf(bubxb + 3 + (i - 1)*nmom)%sf(j, k, l)/nbub_sc - mu**2)**0.5

                        !PRESTON (ISOTHERMAL)
                        pb(j, k, l, 1, i) = (pb0(i))*(R0(i)**(3d0))*(mass_n0(i) + mv(j, k, l, 1, i))/(mu - sig)**(3d0)/(mass_n0(i) + mass_v0(i))
                        pb(j, k, l, 2, i) = (pb0(i))*(R0(i)**(3d0))*(mass_n0(i) + mv(j, k, l, 2, i))/(mu - sig)**(3d0)/(mass_n0(i) + mass_v0(i))
                        pb(j, k, l, 3, i) = (pb0(i))*(R0(i)**(3d0))*(mass_n0(i) + mv(j, k, l, 3, i))/(mu + sig)**(3d0)/(mass_n0(i) + mass_v0(i))
                        pb(j, k, l, 4, i) = (pb0(i))*(R0(i)**(3d0))*(mass_n0(i) + mv(j, k, l, 4, i))/(mu + sig)**(3d0)/(mass_n0(i) + mass_v0(i))
                    end do
                end do
            end do
        end do

    end subroutine s_initialize_pb

    !> The following procedure handles the conversion between
        !!      the conservative variables and the primitive variables.
        !! @param qK_cons_vf Conservative variables
        !! @param qK_prim_vf Primitive variables
        !! @param gm_alphaK_vf Gradient magnitude of the volume fraction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, &
                                                             qK_prim_vf, &
                                                             gm_alphaK_vf, &
                                                             ix, iy, iz)

        type(scalar_field), dimension(sys_size), intent(in) :: qK_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: qK_prim_vf
        type(scalar_field), allocatable, optional, dimension(:), &
            intent(in) :: gm_alphaK_vf

        type(int_bounds_info), optional, intent(in) :: ix, iy, iz

        real(kind(0d0)), dimension(num_fluids) :: alpha_K, alpha_rho_K
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)) :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K
        !WENO primitive reconstruction
        real(kind(0d0)), dimension(sys_size, 3) :: P0, P2, Poly_cons, Poly_prim, omega
        !P_1 is equal to the cell average for zeta_1, zeta_2, zeta 3
        real(kind(0d0)), dimension(sys_size) :: P1, theta
        real(kind(0d0)), dimension(5,5) :: A0, A2, M0, M2, M3, M4
        real(kind(0d0)), dimension(3,5) :: zeta
        real(kind(0d0)), dimension(5) :: q_stencil, beta, zeta1, zeta2, zeta3
        real(kind(0d0)), dimension(3) :: omega_bar, density, weights, dyn_pres
        real(kind(0d0)), dimension(4) :: mu_bar
        real(kind(0d0)) :: r0, r1, r2, s0, s2, s3, s4, mu_0, theta_min, zeta_i = 0.5d0*dsqrt(3.0d0/5.0d0)
        real(kind(0d0)) :: pres_var, phi, gamma_inv, eref_gamma, pref_by_gamma
        #:if MFC_CASE_OPTIMIZATION
#ifndef MFC_SIMULATION
            real(kind(0d0)), dimension(:), allocatable :: nRtmp
#else
            real(kind(0d0)), dimension(nb) :: nRtmp
#endif
        #:else
            real(kind(0d0)), dimension(:), allocatable :: nRtmp
        #:endif

        real(kind(0d0)) :: vftmp, nR3, nbub_sc, R3tmp

        real(kind(0d0)) :: G_K

        real(kind(0d0)) :: pres, temp

        integer :: i, j, k, l, q, r !< Generic loop iterators

        real(kind(0d0)) :: ntmp

        #:if MFC_CASE_OPTIMIZATION
#ifndef MFC_SIMULATION
            if (bubbles) then
                allocate (nRtmp(nb))
            else
                allocate (nRtmp(0))
            end if
#endif
        #:else
            if (bubbles) then
                allocate (nRtmp(nb))
            else
                allocate (nRtmp(0))
            end if
        #:endif
        
        if (model_eqns == 5 .and. num_dims == 1) then 
            A0(1,1) = 3.0d0/640.0d0;A0(1,2) = -29.0d0/480.0d0;A0(1,3) = 1067.0d0/960.0d0;A0(1,4) = -29.0d0/480.0d0;A0(1,5) = 3.0d0/640.0d0
            A0(2,1) = 5.0d0/48.0d0;A0(2,2) = -34.0d0/48.0d0;A0(2,3) =0.0d0;A0(2,4) = 34.0d0/48.0d0;A0(2,5) = -5.0d0/48.0d0
            A0(3,1) = -1.0d0/16.0d0;A0(3,2) = 12.0d0/16.0d0;A0(3,3) =-22.0d0/16.0d0;A0(3,4) = 12.0d0/16.0d0;A0(3,5) = -1.0d0/16.0d0
            A0(4,1) = -1.0d0/12.0d0;A0(4,2) = -2.0d0/12.0d0;A0(4,3) = 0.0d0;A0(4,4) = -2.0d0/12.0d0;A0(4,5) = 1.0d0/12.0d0
            A0(5,1) = 1.0d0/24.0d0;A0(5,2) = -4.0d0/24.0d0;A0(5,3) = 6.0d0/24.0d0;A0(5,4) = -4.0d0/24.0d0;A0(5,5) = 1.0d0/24.0d0

              ! Assign matrix A2
            A2(1,1) = 0.0d0;A2(1,2) = -1.0d0/24.0d0;A2(1,3) = 13.0d0/12.0d0;A2(1,4) = -1.0d0/24.0d0;A2(1,5) = 0.0d0
            A2(2,1) = 0.0d0;A2(2,2) = -1.0d0/2.0d0;A2(2,3) =0.0d0;A2(2,4) = 1.0d0/2.0d0;A2(2,5) = 0.0d0
            A2(3,1) = 0.0d0;A2(3,2) = 1.0d0/2.0d0;A2(3,3) =-1.0d0;A2(3,4) = 1.0d0/2.0d0;A2(3,5) = 0.0d0
            A2(4,1) = 0.0d0;A2(4,2) = 0.0d0;A2(4,3) = 0.0d0;A2(4,4) = 0.0d0;A2(4,5) = 0.0d0
            A2(5,1) = 0.0d0;A2(5,2) = 0.0d0;A2(5,3) = 0.0d0;A2(5,4) = 0.0d0;A2(5,5) = 0.0d0

             ! Assign matrix M0 and M2
          ! Assign upper diagonal elements of M0
            M0(1,1) = 1727.0d0/1680.0d0;M0(1,2) = -51001.0d0/16800.0d0;M0(1,3) = 7547.0d0/1680.0d0
            M0(1,4) = -38947.0d0/16800.0d0;M0(1,5) = 8209.0d0/1680.0d0
            M0(2,2) = 104963.0d0/5040.0d0;M0(2,3) =-24923.0d0/1680.0d0;M0(2,4) = 89549.0d0/5040.0d0;M0(2,5) = -38947.0d0/16800.0d0
            M0(3,3) = 77051.0d0/1680.0d0;M0(3,4) = -24923.0d0/1680.0d0;M0(3,5) = 7547.0d0/1680.0d0
            M0(4,4) = 104963.0d0/5040.0d0;M0(4,5) = -51001.0d0/16800.0d0;M0(5,5) = 1727.0d0/1680.0d0

          ! Assign upper diagonal elements of M2
            M2(1,1) = 4.0d0/3.0d0;M2(1,2) = 5.0d0/3.0d0
            M2(1,3) = -13.0d0/3.0d0;M2(1,4) = 5.0d0/3.0d0;M2(1,5) = 4.0d0/3.0d0;M2(2,2) = 13.0d0/3.0d0
            M2(2,3) = -13.0d0/3.0d0;M2(2,4) = 0.0d0;M2(2,5) = 0.0d0;M2(3,3) = 13.0d0;M2(3,4) = -13.0d0/3.0d0
            M2(3,5) = -13.0d0/3.0d0;M2(4,4) = 5.0d0/3.0d0;M2(4,5) = 4.0d0/3.0d0;M2(5,5) = 0d0
            
          ! Assign upper diagonal elements of M3  
            M3(1,1) = 4.0d0/3.0d0;M3(1,2) = -19.0d0/3.0d0;M3(1,3) =11.0d0/3.0d0;M3(1,4) = 0.0d0;M3(1,5) = 0.0d0
            M3(2,2) = 25.0d0/3.0d0;M3(2,3) = -31.0d0/3.0d0;M3(2,4) = 25.0d0/3.0d0;M3(2,5) = 0.0d0
            M3(3,3) = 10.0d0/3.0d0;M3(3,4) = -31.0d0/3.0d0;M3(3,5) = 11.0d0/3.0d0
            M3(4,4) = 25.0d0/3.0d0;M3(4,5) = -19.0d0/3.0d0
            M3(5,5) = 4.0d0/3.0d0

          ! Assign upper diagonal elements of M4
            M4(1,1) = 0.0d0;M4(1,2) = 0.0d0;M4(1,3) =11.0d0/3.0d0;M4(1,4) = 0.0d0;M4(1,5) = 4.0d0/3.0d0
            M4(2,2) = -19.0d0/3.0d0;M4(2,3) = -31.0d0/3.0d0;M4(2,4) =25.0d0/3.0d0;M4(2,5) = 0.0d0
            M4(3,3) = 10.0d0/3.0d0;M4(3,4) = -31.0d0/3.0d0;M4(3,5) = 11.0d0/3.0d0
            M4(4,4) = 25.0d0/3.0d0;M4(4,5) = -19.0d0/3
            M4(5,5) = 4.0d0/3.0d0

          !Fill the lower triangular part of M0 and M2 using symmetry
            do i = 2, 5
                do j = 1, i-1
                    M0(i,j) = M0(j,i)
                    M2(i,j) = M2(j,i)
                    M3(i,j) = M3(j,i)
                    M4(i,j) = M4(j,i)
                end do
            end do
            r0 = 1100d0/1111d0; r1 = 1d0/1111d0; r2 = 10d0/1111d0; s0 = 0.97d0; s2 = 0.01d0; s3 = 0.01d0; s4 = 0.01d0 
            !calculate zeta
            zeta(1,1) = 1.0d0;zeta(1,2) = -zeta_i;zeta(1,3) = zeta_i**2d0;zeta(1,4) = -zeta_i**3d0;zeta(1,5) = zeta_i**4d0
            zeta(2,1) = 1.0d0;zeta(2,2) = 0d0; zeta(2,3) = 0d0; zeta(2,4) = 0d0; zeta(2,5) = 0d0
            zeta(3,1) = 1.0d0;zeta(3,2) = zeta_i; zeta(3,3) = zeta_i**2d0; zeta(3,4) = zeta_i**3d0; zeta(3,5) =zeta_i**4d0
            !calc weights for each gauss quadrature point
            weights(1) = 5.0d0/9.0d0; weights(2) = 8.0d0/9.0d0; weights(3) = 5.0d0/9.0d0 
        end if
        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_K, alpha_rho_K, Re_K, nRtmp, rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, R3tmp, G_K)
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe
                    dyn_pres_K = 0d0

                    !$acc loop seq
                    do i = 1, num_fluids
                        alpha_rho_K(i) = qK_cons_vf(i)%sf(j, k, l)
                        alpha_K(i) = qK_cons_vf(advxb + i - 1)%sf(j, k, l)
                    end do

                    !$acc loop seq
                    do i = 1, contxe
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                    end do

                    if (model_eqns /= 4) then
#ifdef MFC_SIMULATION
                        ! If in simulation, use acc mixture subroutines
                        if (elasticity) then
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, &
                                                                            alpha_rho_K, Re_K, j, k, l, G_K, Gs)
                        else if (bubbles) then
                            call s_convert_species_to_mixture_variables_bubbles_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                                    alpha_K, alpha_rho_K, Re_K, j, k, l)
                        else
                            call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                            alpha_K, alpha_rho_K, Re_K, j, k, l)
                        end if
#else
                        ! If pre-processing, use non acc mixture subroutines
                        if (elasticity) then
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, &
                                                                rho_K, gamma_K, pi_inf_K, qv_K, Re_K, G_K, fluid_pp(:)%G)
                        else
                            call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, &
                                                                rho_K, gamma_K, pi_inf_K, qv_K)
                        end if
#endif
                    end if

#ifdef MFC_SIMULATION
                    rho_K = max(rho_K, sgm_eps)
#endif

                    !$acc loop seq
                    do i = momxb, momxe
                        if (model_eqns /= 4) then
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                        /rho_K
                            dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf(i)%sf(j, k, l) &
                                         *qK_prim_vf(i)%sf(j, k, l)
                        else
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                        /qK_cons_vf(1)%sf(j, k, l)
                        end if
                    end do
                    ! PRESSURE CALCULATION
                    if (model_eqns/= 5) then 
                        call s_compute_pressure(qK_cons_vf(E_idx)%sf(j, k, l), &
                                                qK_cons_vf(alf_idx)%sf(j, k, l), &
                                                dyn_pres_K, pi_inf_K, gamma_K, rho_K, qv_K, pres)
                    elseif ((num_dims == 1) .and. (model_eqns == 5)) then
                           !write code for polynomial reconstruction
                           do i = 1, sys_size
                                do q = 1, 5
                                    q_stencil(q) = qK_cons_vf(i)%sf(j - 3 + q, k, l)
                                    zeta1(q) = zeta(1,q)
                                    zeta2(q) = zeta(2,q)
                                    zeta3(q) = zeta(3,q)
                                end do  
                                beta(1) = dot_product(q_stencil,matmul(M0,q_stencil))
                                beta(2) = min((qK_cons_vf(i)%sf(j+1,k,l)-qK_cons_vf(i)%sf(j,k,l))**2d0,&
                                    (qK_cons_vf(i)%sf(j,k,l)-qK_cons_vf(i)%sf(j-1,k,l))**2d0)
                                beta(3) = dot_product(q_stencil,matmul(M2,q_stencil))
                                beta(4) = dot_product(q_stencil,matmul(M3, q_stencil))
                                beta(5) = dot_product(q_stencil,matmul(M4, q_stencil))

                                !Compute the non-linear weights omega_bar_k
                                omega_bar(1) = r0/(beta(1)+verysmall)**2d0
                                omega_bar(2) = r1/(beta(2)+verysmall)**2d0
                                omega_bar(3) = r2/(beta(3)+verysmall)**2d0

                                mu_bar(1) = s0/(beta(1)+verysmall)**4d0
                                mu_bar(2) = s2/(beta(3)+verysmall)**4d0
                                mu_bar(3) = s3/(beta(4)+verysmall)**4d0
                                mu_bar(4) = s4/(beta(5)+verysmall)**4d0
                                mu_0 = mu_bar(1)/(mu_bar(1)+mu_bar(2)+mu_bar(3)+mu_bar(4))
                                
                                theta(i) = 1.0d0-(1.0d0-min(1.0d0,mu_0/s0)**4d0)**4d0
                                
                                !Calculate polynomials P0,P1,P2, theta(m)
                                P0(i,1) = dot_product(zeta1,matmul(A0,q_stencil))
                                P0(i,2) = dot_product(zeta2,matmul(A0,q_stencil))
                                P0(i,3) = dot_product(zeta3,matmul(A0,q_stencil))
                                P1(i) = qK_cons_vf(i)%sf(j, k, l)
                                P2(i,1) = dot_product(zeta1,matmul(A2,q_stencil))
                                P2(i,2) = dot_product(zeta2,matmul(A2,q_stencil))
                                P2(i,3) = dot_product(zeta3,matmul(A2,q_stencil))
                                do q = 1, 3
                                    omega(i,q) = omega_bar(q)/(omega_bar(1)+omega_bar(2)+omega_bar(3))
                                end do
                               ! do q = 1, 3
                               !     Poly_cons(i,q) = theta(i)*P0(i,q)+&
                               !     (1d0-theta(i))*(omega(i,1)*P0(i,q)+omega(i,2)*P1(i)+omega(i,3)*P2(i,q))
                               ! end do
                           end do
                           theta_min = minval(theta)
                           do i = 1, sys_size
                                do q = 1,3 
                                    Poly_cons(i,q) = theta_min*P0(i,q) +(1d0-theta_min)*P2(i,q)
                                end do
                           end do 
                                !if (j == 7 .or. j==8 .or. j==9) then 
                                !    print *,'j::',j,'mom::','P0',P0(3,1),P0(3,2),P0(3,3),'P2',P2(3,1),P2(3,2),P2(3,3)
                                !    print *,'j::',j,'dens::','P0',P0(1,1),P0(1,2),P0(1,3),'P2',P2(1,1),P2(1,2),P2(1,3)
                                !    print *,'j::',j,'dens::','P0',P0(2,1),P0(2,2),P0(2,3),'P2',P2(2,1),P2(2,2),P2(2,3)
                                !end if
                           !Step 3: Calculate the polynomials of primitive variables
                           !Calc densities
                           density(:) = 0d0
                           dyn_pres(:) = 0d0
                           pres_var = 0d0
                           !if (j == 8 .or. j == 7 .or. j==9) then 
                           !    print *, j, Poly_cons(3,1), Poly_cons(3,2), Poly_cons(3,3)
                           !end if
                           !if (j == 32 .or. j == 31) then
                           !     print *, Poly_cons(mgidxb,1),Poly_cons(mgidxb,2), Poly_cons(mgidxb,3)
                           !end if
                           do q = 1,3
                                do i = 1, contxe
                                    Poly_prim(i,q) = Poly_cons(i,q)
                                    density(q) = density(q) + Poly_cons(i,q)
                                end do
                                !Calc velocities
                                do i = momxb, momxe
                                    Poly_prim(i,q) = Poly_cons(i,q)/density(q)
                                    dyn_pres(q)= dyn_pres(q) + 0.5d0*density(q)*Poly_prim(i,q)**2d0
                                end do
                                !if (j ==7 .or. j==8 .or. j==9) then 
                                !    print *,'j::',j, 'mom::',Poly_cons(3,q),'density::',density(q)
                                !end if
                                !Calc pressure
                                !pres_var1=0d0
                               ! do i= 1, num_fluids
                               !     phi =((Poly_cons(advxb+i-1,q)*rho0(i)/Poly_cons(i,q))**(-mg_a(i)))*&
                               !         dexp((gammas(i)-mg_a(i))*(1d0-(Poly_cons(advxb+i-1,q)*rho0(i)/Poly_cons(i,q))**mg_b(i)))
                               !     pres_var1 = pres_var1 + &
                               !        Poly_cons(advxb+i-1,q)*Poly_cons(mgidxb,q)*&
                               !         pi_infs(i)*dlog(Poly_cons(i,q)/(Poly_cons(advxb+i-1,q)*rho0(i)))*&
                               !        (1d0+0.5d0*(qvs(i)-2d0)*dlog(Poly_cons(i,q)/(Poly_cons(advxb+i-1,q)*rho0(i))))-&
                               !       Poly_cons(i,q)*(0.5d0*(pi_infs(i)/rho0(i))*&
                               !         ((dlog(Poly_cons(i,q)/(Poly_cons(advxb+i-1,q)*rho0(i))))**2d0)*&
                               !         (1d0+(1d0/3d0)*(qvs(i)-2d0)*dlog(Poly_cons(i,q)/(Poly_cons(advxb+i-1,q)*rho0(i))))&
                               !         +ein_cv1(i)*(phi*ein_cv2(i)*dexp(phi*ein_cv2(i))/&
                               !         (dexp(phi*ein_cv2(i))-1d0)-dlog(dexp(phi*ein_cv2(i))-1d0)))
                                    
                                   !if (pres_var /= pres_var) then
                                   !     print *,i,Pres_var,phi,Poly_cons(advxb+i-1,q),rho0(i),Poly_cons(i,q)
                                   ! end if
                                !end do
                                gamma_inv = 0d0
                                Pref_by_gamma = 0d0
                                eref_gamma = 0d0
                                do i = 1, num_fluids
                                alpha_K(i) =Poly_cons(advxb+i-1,q)
                                alpha_rho_K(i)=Poly_cons(i,q)
                                gamma_inv = gamma_inv + &
                                alpha_K(i)/(mg_a(i)+(gammas(i)-mg_a(i))*(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i))
                                phi = (alpha_K(i)*rho0(i)/alpha_rho_K(i))**(-mg_a(i))*&
                                    dexp((gammas(i)-mg_a(i))*(1d0-(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i)))
                                Pref_by_gamma = Pref_by_gamma !+&
                                   ! alpha_K(i)*pi_infs(i)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i)))*&
                                   ! (1d0+0.5d0*(qvs(i)-2d0)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i))))/&
                                   ! (mg_a(i)+(gammas(i)-mg_a(i))*(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i))
                                eref_gamma = eref_gamma + &
                                   ! alpha_K(i)*(0.5d0*(pi_infs(i)/rho0(i))*(dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i)))**2d0)*&
                                   ! (1d0+(1d0/3d0)*(qvs(i)-2d0)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i))))+&
                                    alpha_K(i)*(ein_cv1(i)*(phi*ein_cv2(i)*dexp(phi*ein_cv2(i))/&
                                    (dexp(phi*ein_cv2(i))-1d0)-dlog(dexp(phi*ein_cv2(i))-1.d0)))
                                end do
                                !Poly_prim(E_idx,q) = (Poly_cons(E_idx,q)& 
                                !        - dyn_pres(q) +Poly_cons(mgidxb+1,q)& 
                                !    - Poly_cons(mgidxe,q))/Poly_cons(mgidxb,q)
                                Poly_prim(E_idx,q) = (Poly_cons(E_idx,q)& 
                                        - dyn_pres(q) + Pref_by_gamma &
                                    -density(q)*eref_gamma)/Poly_cons(mgidxb,q)
                                if (j == 8 .or. j==9 .or. j==7 .or. j==6) then
                                    !print *,j,alpha_rho_K(1)/(alpha_K(1)*rho0(1)),alpha_rho_K(2)/(alpha_K(2)*rho0(2))
                                    print *,j,q, eref_gamma
                                end if
                               ! print *,'Prefbygam',Pref_by_gamma,'rho_e_ref',density(q)*eref_gamma,'ene',Poly_cons(E_idx,q),'dyn_pr',dyn_pres(q)
                                do i = advxb, advxe
                                    Poly_prim(i,q) = Poly_cons(i,q)
                                end do
                                !Calc primitive variables of MG EoS
                                Poly_prim(mgidxb,q) = Poly_cons(mgidxb,q)
                                Poly_prim(mgidxb+1, q) = Poly_cons(mgidxb+1,q)/Poly_cons(mgidxb,q)
                                Poly_prim(mgidxe,q) = Poly_cons(mgidxe,q)/density(q)
                                !Calc Stress
                                if (hypoplasticity) then
                                    do i = strxb,strxe
                                        Poly_prim(i,q) = Poly_cons(i,q)/density(q)
                                    end do
                                    Poly_prim(plasidx,q) = Poly_cons(plasidx,q)/density(q)
                                end if
                           end do
                           !Step 4: Integrate to calculate cell average
                           !primitive variables
                           !print *, 'j', j
                           do i = 1, sys_size
                                    qK_prim_vf(i)%sf(j, k, l) = 0d0
                                do q = 1,3
                                    qK_prim_vf(i)%sf(j, k, l) = &
                                    qK_prim_vf(i)%sf(j, k, l) + 0.5d0*weights(q)*Poly_prim(i,q)
                                end do
                           !     print *,i,qK_prim_vf(i)%sf(j, k, l)
                           end do
                           !print *,'Ene = ', qK_cons_vf(E_idx)%sf(j, k, l)
                           pres = qK_prim_vf(E_idx)%sf(j ,k, l)
                           print * ,'pres=',j, pres
                           !Step 5: Once the code runs create function
                           !to calculate all the primitive variables
                           !Step 6: Do it for 2D
                           !Step 7: Do it for 3D
                    else
                            call s_compute_pressure(qK_cons_vf(E_idx)%sf(j, k, l), &
                                               0d0, dyn_pres_K, & 
                                               pi_inf_K, qK_cons_vf(mgidxb)%sf(j, k, l), rho_K, qv_K, &
                                               pres, 0d0, 0d0, 0d0,qK_cons_vf(mgidxb+1)%sf(j, k, l),&
                                               qK_cons_vf(mgidxe)%sf(j, k, l))
                                           
                            qK_prim_vf(mgidxb)%sf(j, k, l)   = qK_cons_vf(mgidxb)%sf(j, k, l) 
                            qK_prim_vf(mgidxb+1)%sf(j, k, l) = qK_cons_vf(mgidxb+1)%sf(j, k, l)/&
                                                               qK_cons_vf(mgidxb)%sf(j, k, l)
                            qK_prim_vf(mgidxe)%sf(j, k, l)   = qK_cons_vf(mgidxe)%sf(j, k, l)/rho_K
                           !print *, 'j',j,'rho_eref',qK_cons_vf(mgidxe)%sf(j,k,l)
                           ! if (pres /= pres)  then 
                           !     print *,'j',j,'pres',pres,'energy'
                           !     call s_mpi_abort()
                           ! end if
                            !$acc loop seq
                            do i=1, sys_size
                               if (qK_cons_vf(i)%sf(j,k,l) /= qK_cons_vf(i)%sf(j,k,l)) then
                                   print *, 'i',i,'j k',j,k,qK_cons_vf(i)%sf(j, k, l)
                                   call s_mpi_abort()
                               end if 
                            end do
                    end if
                    if (model_eqns == 5) then
#ifdef MFC_POST_PROCESS                        
                        call s_compute_temperature(qK_prim_vf(E_idx)%sf(j, k, l), &
                                                   qK_prim_vf(mgidxb+1)%sf(j, k, l), &
                                                   qK_prim_vf(mgidxb)%sf(j, k, l),&
                                                   rho_K, temp, alpha_K)
                        qK_prim_vf(plasidx+1)%sf(j, k, l) = temp
#endif
                    end if   
                       
                    qK_prim_vf(E_idx)%sf(j, k, l) = pres

                    if (bubbles) then
                        !$acc loop seq
                        do i = 1, nb
                            nRtmp(i) = qK_cons_vf(bubrs(i))%sf(j, k, l)
                        end do

                        vftmp = qK_cons_vf(alf_idx)%sf(j, k, l)

                        if (qbmm) then
                            !Get nb (constant across all R0 bins)
                            nbub_sc = qK_cons_vf(bubxb)%sf(j, k, l)

                            !Convert cons to prim
                            !$acc loop seq
                            do i = bubxb, bubxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                            !Need to keep track of nb in the primitive variable list (converted back to true value before output)
#ifdef MFC_SIMULATION
                            qK_prim_vf(bubxb)%sf(j, k, l) = qK_cons_vf(bubxb)%sf(j, k, l)
#endif

                        else
                            if (adv_n) then
                                qK_prim_vf(n_idx)%sf(j, k, l) = qK_cons_vf(n_idx)%sf(j, k, l)
                                nbub_sc = qK_prim_vf(n_idx)%sf(j, k, l)
                            else
                                call s_comp_n_from_cons(vftmp, nRtmp, nbub_sc, weight)
                            end if

                            !$acc loop seq
                            do i = bubxb, bubxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                            end do
                        end if
                    end if

                    if (elasticity) then 
                        !$acc loop seq
                        do i = strxb, strxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    if (hypoelasticity) then
                        !$acc loop seq
                        do i = strxb, strxe
                            ! subtracting elastic contribution for pressure calculation
                            if (G_K > verysmall) then !TODO: check if stable for >0
                                qK_prim_vf(E_idx)%sf(j, k, l) = qK_prim_vf(E_idx)%sf(j, k, l) - &
                                                                ((qK_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G_K))/gamma_K
                                ! extra terms in 2 and 3D
                                if ((i == strxb + 1) .or. &
                                    (i == strxb + 3) .or. &
                                    (i == strxb + 4)) then
                                    qK_prim_vf(E_idx)%sf(j, k, l) = qK_prim_vf(E_idx)%sf(j, k, l) - &
                                                                    ((qK_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G_K))/gamma_K
                                end if
                            end if
                        end do
                    end if 
                    if (hypoplasticity) then
                        !$acc loop seq
                        do i = strxb, strxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                            !if (i == 11 .or. i==12 .or. i==13 .or. i==14) then
                            !    print *,'j',j,'k',k,'i',i,qK_prim_vf(i)%sf(j,k,l)
                            !end if
                        !    if (qK_prim_vf(i)%sf(j, k, l) /= qK_prim_vf(i)%sf(j, k, l)) then
                        !         print *,'j',j,'k',k,'i',i,qK_prim_vf(i)%sf(j, k, l)
                        !         call s_mpi_abort('stress is NaN')
                        !   end if
                        end do
                        qK_prim_vf(plasidx)%sf(j, k, l) = qK_cons_vf(plasidx)%sf(j, k, l)/rho_K
                    end if

                    if (hyperelasticity) then
                        !$acc loop seq
                        do i = xibeg, xiend
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        end do
                    end if

                    !$acc loop seq
                    do i = advxb, advxe
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                    end do

                    if (.not. f_is_default(sigma)) then
                        qK_prim_vf(c_idx)%sf(j, k, l) = qK_cons_vf(c_idx)%sf(j, k, l)
                    end if
#ifdef MFC_SIMULATION
                    !if (qK_prim_vf(12)%sf(j, k, l)/= qK_prim_vf(12)%sf(j,k,l)) then
                        !print *,j,'prim' 
                    !do i=1, sys_size
                    !    print *,j,qK_prim_vf(12)%sf(j,k,l)
                    !end if
                    !end do
                    !print *,'cons'
                    !do i=1, sys_size
                    !    print *,qK_cons_vf(i)%sf(j,k,l)
                    !end do
#endif
                end do
            end do
        end do
        print *, maxval(qK_prim_vf(E_idx)%sf(ixb:ixe,0,0))
        !$acc end parallel loop
    end subroutine s_convert_conservative_to_primitive_variables ! ---------

    !>  The following procedure handles the conversion between
        !!      the primitive variables and the conservative variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Conservative variables
        !!  @param gm_alphaK_vf Gradient magnitude of the volume fractions
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                             q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function
        ! and dynamic pressure, as defined in the incompressible flow sense,
        ! respectively
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: qv
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: nbub, R3, vftmp, R3tmp
        real(kind(0d0)), dimension(nb) :: Rtmp
        real(kind(0d0)) :: G = 0d0
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)), dimension(num_fluids) :: alpha_K, alpha_rho_K
        
        !Local variables used for Mie-Gruneisen EOS
        real(kind(0d0)) :: eref_gamma, Pref_by_gamma, gamma_inv, ein_cv1_mix,&
            theta_E_mix , mg_a_mix, mg_b_mix, rho0_mix, &
            log_rho_mix, phi

        integer :: i, j, k, l !< Generic loop iterators      

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Obtaining the density, specific heat ratio function
                    ! and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                          rho, gamma, pi_inf, qv, Re_K, G, fluid_pp(:)%G)

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, contxe
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        alpha_rho_K(i) = q_prim_vf(i)%sf(j, k, l)
                    end do
                    ! Transferring the advection equation(s) variable(s)
                    do i = adv_idx%beg, adv_idx%end
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        alpha_K(i- adv_idx%beg + 1) = q_prim_vf(i)%sf(j, k, l)
                    end do
                    ! Zeroing out the dynamic pressure since it is computed
                    ! iteratively by cycling through the velocity equations
                    dyn_pres = 0d0

                    ! Computing momenta and dynamic pressure from velocity
                    do i = momxb, momxe
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)* &
                                   q_prim_vf(i)%sf(j, k, l)/2d0
                    end do

                    ! Computing the energy from the pressure
                    if ((model_eqns /= 4) .and. (model_eqns /= 5) .and. (bubbles .neqv. .true.)) then
                        ! E = Gamma*P + \rho u u /2 + \pi_inf + (\alpha\rho qv)
                        q_cons_vf(E_idx)%sf(j, k, l) = &
                            gamma*q_prim_vf(E_idx)%sf(j, k, l) + dyn_pres + pi_inf &
                            + qv
                    else if ((model_eqns /= 4) .and. (bubbles)) then
                        ! \tilde{E} = dyn_pres + (1-\alf)(\Gamma p_l + \Pi_inf)
                        q_cons_vf(E_idx)%sf(j, k, l) = dyn_pres + &
                                                       (1.d0 - q_prim_vf(alf_idx)%sf(j, k, l))* &
                                                       (gamma*q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)
                    else if (model_eqns == 5) then
                        ! Calculate the extra primitive variables
                        !   ein_cv1_mix = 0d0
                        !   theta_E_mix = 0d0
                        !   mg_a_mix = 0d0
                        !   mg_b_mix = 0d0
                        !   rho0_mix = 0d0
                        !do i=1, num_fluids
                        !   mg_b_mix = mg_b_mix+q_prim_vf(adv_idx%beg-1+i)%sf(j, k, l)*mg_b(i)
                        !   mg_a_mix = mg_a_mix+q_prim_vf(adv_idx%beg-1+i)%sf(j, k, l)*mg_a(i)
                        !   rho0_mix = rho0_mix+q_prim_vf(adv_idx%beg-1+i)%sf(j, k, l)*rho0(i)
                        !   ein_cv1_mix = ein_cv1_mix + q_prim_vf(adv_idx%beg -1 +i)%sf(j, k, l)*ein_cv1(i)
                        !   theta_E_mix = theta_E_mix + q_prim_vf(adv_idx%beg-1+i)%sf(j, k, l)*ein_cv2(i)
                        !end do
                        !   log_rho_mix = dlog(rho/rho0_mix)
                        !   phi = ((rho0_mix/rho)**(-mg_a_mix))*dexp((gamma-mg_a_mix)*(1d0-(rho0_mix/rho)**mg_b_mix))
                        !   gamma_inv = 1d0/(mg_a_mix+(gamma-mg_a_mix)*(rho0_mix/rho)**mg_b_mix) 
                        !   Pref = pi_inf*log_rho_mix*(1d0+0.5d0*(qv-2d0)*log_rho_mix)
                        !   eref = 0.5d0*pi_inf*(log_rho_mix**2d0)*(1d0+(1/3d0)*(qv-2d0)*log_rho_mix)+&
                        !   ein_cv1_mix*(phi*theta_E_mix*dexp(phi*theta_E_mix)/(dexp(phi*theta_E_mix)-1)-dlog(dexp(phi*theta_E_mix)-1))
                        Pref_by_gamma = 0d0;eref_gamma = 0d0;gamma_inv = 0d0
                        do i = 1, num_fluids
                            gamma_inv = gamma_inv + &
                            alpha_K(i)/(mg_a(i)+(gammas(i)-mg_a(i))*(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i))
                            phi = (alpha_K(i)*rho0(i)/alpha_rho_K(i))**(-mg_a(i))*&
                                dexp((gammas(i)-mg_a(i))*(1d0-(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i)))
                            Pref_by_gamma = Pref_by_gamma + alpha_K(i)*pi_infs(i)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i)))*&
                                (1d0+0.5d0*(qvs(i)-2d0)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i))))/&
                                (mg_a(i)+(gammas(i)-mg_a(i))*(alpha_K(i)*rho0(i)/alpha_rho_K(i))**mg_b(i))
                            eref_gamma = eref_gamma + &
                                alpha_K(i)*(0.5d0*(pi_infs(i)/rho0(i))*(dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i)))**2d0)*&
                                (1d0+(1d0/3d0)*(qvs(i)-2d0)*dlog(alpha_rho_K(i)/(alpha_K(i)*rho0(i))))+&
                                ein_cv1(i)*(phi*ein_cv2(i)*dexp(phi*ein_cv2(i))/&
                                (dexp(phi*ein_cv2(i))-1d0)-dlog(dexp(phi*ein_cv2(i))-1.d0)))
                        end do
                        q_prim_vf(mgidxb)%sf(j, k, l) = gamma_inv
                        q_prim_vf(mgidxb+1)%sf(j, k, l) = Pref_by_gamma/gamma_inv
                        q_prim_vf(mgidxe)%sf(j, k, l) = eref_gamma
                        ! Energy corresponding to Mie-Gruneisen EOS 
                        q_cons_vf(E_idx)%sf(j, k, l)    = rho*eref_gamma +&
                                                          gamma_inv*(q_prim_vf(E_idx)%sf(j,k,l)-Pref_by_gamma/gamma_inv) +& 
                                                          dyn_pres
                    
!                        print *,q_cons_vf(E_idx)%sf(j, k, l), &
!                            q_prim_vf(E_idx)%sf(j,k,l),&
!                        Pref_by_gamma-rho*q_prim_vf(mgidxe)%sf(j, k, l)&
!                        , Pref_by_gamma, rho*q_prim_vf(mgidxe)%sf(j,k,l)
!                        print *,'Prefbygam',Pref_by_gamma,'rho_e_ref',rho*eref_gamma,'ene',q_cons_vf(E_idx)%sf(j,k,l),'dyn_pr',dyn_pres
                        q_cons_vf(mgidxb)%sf(j, k, l)   = q_prim_vf(mgidxb)%sf(j, k, l)
                        q_cons_vf(mgidxb+1)%sf(j, k, l) = Pref_by_gamma
                        q_cons_vf(mgidxe)%sf(j, k, l)   = rho*eref_gamma
                        !print *, alpha_rho_K(1),alpha_K(1),rho0(1),&
                        !alpha_rho_K(2), alpha_K(2), rho0(2)
                    else
                        !Tait EOS, no conserved energy variable
                        q_cons_vf(E_idx)%sf(j, k, l) = 0.d0
                    end if

                    ! Computing the internal energies from the pressure and continuities
                    if (model_eqns == 3) then
                        do i = 1, num_fluids
                            ! internal energy calculation for each of the fluids
                            q_cons_vf(i + internalEnergies_idx%beg - 1)%sf(j, k, l) = &
                                q_cons_vf(i + adv_idx%beg - 1)%sf(j, k, l)* &
                                (fluid_pp(i)%gamma*q_prim_vf(E_idx)%sf(j, k, l) + &
                                 fluid_pp(i)%pi_inf) + &
                                q_cons_vf(i + cont_idx%beg - 1)%sf(j, k, l)*fluid_pp(i)%qv
                        end do
                    end if

                    if (bubbles) then
                        ! From prim: Compute nbub = (3/4pi) * \alpha / \bar{R^3}
                        do i = 1, nb
                            Rtmp(i) = q_prim_vf(bub_idx%rs(i))%sf(j, k, l)
                        end do

                        if (.not. qbmm) then
                            if (adv_n) then
                                q_cons_vf(n_idx)%sf(j, k, l) = q_prim_vf(n_idx)%sf(j, k, l)
                                nbub = q_prim_vf(n_idx)%sf(j, k, l)
                            else
                                call s_comp_n_from_prim(q_prim_vf(alf_idx)%sf(j, k, l), Rtmp, nbub, weight)
                            end if
                        else
                            !Initialize R3 averaging over R0 and R directions
                            R3tmp = 0d0
                            do i = 1, nb
                                R3tmp = R3tmp + weight(i)*0.5d0*(Rtmp(i) + sigR)**3d0
                                R3tmp = R3tmp + weight(i)*0.5d0*(Rtmp(i) - sigR)**3d0
                            end do
                            !Initialize nb
                            nbub = 3d0*q_prim_vf(alf_idx)%sf(j, k, l)/(4d0*pi*R3tmp)
                        end if

                        if (j == 0 .and. k == 0 .and. l == 0) print *, 'In convert, nbub:', nbub

                        do i = bub_idx%beg, bub_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)*nbub
                        end do
                    end if

                    if (elasticity) then 
                        ! adding the elastic contribution
                        ! Multiply \tau to \rho \tau
                        do i = strxb, strxe
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    if (hypoelasticity) then
                        do i = strxb, strxe
                            ! adding elastic contribution
                            if (G > verysmall .and. .not. hypoplasticity) then
                                q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                               (q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G)
                                ! extra terms in 2 and 3D
                                if ((i == stress_idx%beg + 1) .or. &
                                    (i == stress_idx%beg + 3) .or. &
                                    (i == stress_idx%beg + 4)) then
                                    q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                                   (q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G)
                                end if
                            end if
                        end do
                    end if

                    if (hypoplasticity) then 
                      q_cons_vf(plasidx)%sf(j, k, l) = rho*q_prim_vf(plasidx)%sf(j, k, l) 
                    end if 
  
                    ! using \rho xi as the conservative formulation stated in Kamrin et al. JFM 2022
                    if (hyperelasticity) then
                        ! Multiply \xi to \rho \xi
                        do i = xibeg, xiend
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if
                    if (.not. f_is_default(sigma)) then
                        q_cons_vf(c_idx)%sf(j, k, l) = q_prim_vf(c_idx)%sf(j, k, l)
                    end if

                end do
            end do
        end do
#else
        if (proc_rank == 0) then
            call s_mpi_abort('Conversion from primitive to '// &
                             'conservative variables not '// &
                             'implemented. Exiting ...')
        end if
#endif
    end subroutine s_convert_primitive_to_conservative_variables

    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, &
                                                     FK_vf, &
                                                     FK_src_vf, &
                                                     is1, is2, is3, s2b, s3b)

        integer, intent(in) :: s2b, s3b
        real(kind(0d0)), dimension(0:, s2b:, s3b:, 1:), intent(in) :: qK_prim_vf
        real(kind(0d0)), dimension(0:, s2b:, s3b:, 1:), intent(inout) :: FK_vf
        real(kind(0d0)), dimension(0:, s2b:, s3b:, advxb:), intent(inout) :: FK_src_vf

        type(int_bounds_info), intent(in) :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K
        real(kind(0d0)), dimension(num_fluids) :: alpha_K
        real(kind(0d0)) :: rho_K
        real(kind(0d0)), dimension(num_dims) :: vel_K
        real(kind(0d0)) :: vel_K_sum
        real(kind(0d0)) :: pres_K
        real(kind(0d0)) :: E_K
        real(kind(0d0)) :: gamma_K
        real(kind(0d0)) :: pi_inf_K
        real(kind(0d0)) :: qv_K
        real(kind(0d0)), dimension(2) :: Re_K
        real(kind(0d0)) :: G_K

        integer :: i, j, k, l !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        !$acc update device(is1b, is2b, is3b, is1e, is2e, is3e)

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity
#ifdef MFC_SIMULATION
        !$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_K, vel_K, alpha_K, Re_K)
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e

                    !$acc loop seq
                    do i = 1, contxe
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

                    !$acc loop seq
                    do i = advxb, advxe
                        alpha_K(i - E_idx) = qK_prim_vf(j, k, l, i)
                    end do
                    !$acc loop seq
                    do i = 1, num_dims
                        vel_K(i) = qK_prim_vf(j, k, l, contxe + i)
                    end do

                    vel_K_sum = 0d0
                    !$acc loop seq
                    do i = 1, num_dims
                        vel_K_sum = vel_K_sum + vel_K(i)**2d0
                    end do

                    pres_K = qK_prim_vf(j, k, l, E_idx)
                    if (elasticity) then
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                        alpha_K, alpha_rho_K, Re_K, &
                                                                        j, k, l, G_K, Gs)
                    else if (bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho_K, gamma_K, &
                                                                                pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K, j, k, l)
                    else
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, &
                                                                        alpha_K, alpha_rho_K, Re_K, j, k, l)
                    end if

                    ! Computing the energy from the pressure
                    E_K = gamma_K*pres_K + pi_inf_K &
                          + 5d-1*rho_K*vel_K_sum + qv_K

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    !$acc loop seq
                    do i = 1, contxe
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    !$acc loop seq
                    do i = 1, num_dims
                        FK_vf(j, k, l, contxe + dir_idx(i)) = &
                            rho_K*vel_K(dir_idx(1)) &
                            *vel_K(dir_idx(i)) &
                            + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, E_idx) = vel_K(dir_idx(1))*(E_K + pres_K)

                    if (riemann_solver == 1) then
                        !$acc loop seq
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = 0d0
                            FK_src_vf(j, k, l, i) = alpha_K(i - E_idx)
                        end do

                    else
                        ! Could be bubbles!
                        !$acc loop seq
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = vel_K(dir_idx(1))*alpha_K(i - E_idx)
                        end do

                        !$acc loop seq
                        do i = advxb, advxe
                            FK_src_vf(j, k, l, i) = vel_K(dir_idx(1))
                        end do

                    end if

                end do
            end do
        end do
#endif
    end subroutine s_convert_primitive_to_flux_variables

    subroutine s_finalize_variables_conversion_module() ! ------------------

        integer :: i !< Generic loop iterators

        ! Deallocating the density, the specific heat ratio function and the
        ! liquid stiffness function
#ifdef MFC_POST_PROCESS
        deallocate (rho_sf, gamma_sf, pi_inf_sf, qv_sf)
#endif

#ifdef MFC_SIMULATION
        @:DEALLOCATE_GLOBAL(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs)
        if (bubbles) then
            @:DEALLOCATE_GLOBAL(bubrs)
        end if
#else
        @:DEALLOCATE(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps, Gs)
        if (bubbles) then
            @:DEALLOCATE(bubrs)
        end if
#endif

        ! Nullifying the procedure pointer to the subroutine transferring/
        ! computing the mixture/species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module

end module m_variables_conversion
