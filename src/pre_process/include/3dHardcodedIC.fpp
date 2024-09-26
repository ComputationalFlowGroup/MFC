#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(kind(0d0)) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(kind(0d0)) :: eps

    real(kind(0d0)) :: rcoord, theta, phi, xi_sph, x_centriod, y_centriod, z_centriod
    real(kind(0d0)), dimension(num_dims) :: xi_cart
    integer :: ii

    eps = 1e-9
#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3
        rhoL = 1
        pRef = 1e5
        pInt = pRef
        h = 0.7
        lam = 0.2
        wl = 2*pi/lam
        amp = 0.025/wl

        intH = amp*(sin(2*pi*x_cc(i)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h

        alph = 5d-1*(1 + tanh((y_cc(j) - intH)/2.5e-3))

        if (alph < eps) alph = eps
        if (alph > 1 - eps) alph = 1 - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81*(1.2 - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1 - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1 - alph)*rhoL
            pInt = pref + rhoH*9.81*(1.2 - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0
        lam = 1.0
        amp = patch_icpp(patch_id)%a2
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

     case (302) ! pre_stress for hyperelasticity, bubble in material
        R0ref = 30E-6    ! equilibrium radius
        Rinit = 244.8E-6 ! initial radius
        x_centroid = patch_icpp(patch_id)%x_centroid
        y_centroid = patch_icpp(patch_id)%y_centroid
        z_centroid = patch_icpp(patch_id)%z_centroid
        rcoord = sqrt(x_cc(i)**2 + y_cc(j)**2 + z_cc(k)**2)
        theta = atan2(y_cc(j), x_cc(i))
        phi = atan2(sqrt(x_cc(i)**2 + y_cc(j)**2), z_cc(k))
        !spherical coord, assuming Rmax=1
        xi_sph = (rcoord**3 + R0ref**3 - Rinit**3)**(1d0/3d0)
        !xi_cart(1) = xi_sph*sin(phi)*cos(theta)
        !xi_cart(2) = xi_sph*sin(phi)*sin(theta)
        !xi_cart(3) = xi_sph*cos(phi)
        xi_cart(1) = (xi_sph*x_cc(i)) / &
            sqrt((x_cc(i) - x_centriod)**2 + (y_cc(j) - y_centriod)**2 + (z_cc(k) - z_centriod)**2)
        xi_cart(2) = (xi_sph*y_cc(j)) / &
            sqrt((x_cc(i) - x_centriod)**2 + (y_cc(j) - y_centriod)**2 + (z_cc(k) - z_centriod)**2)
        xi_cart(3) = (xi_sph*z_cc(k)) / &
            sqrt((x_cc(i) - x_centriod)**2 + (y_cc(j) - y_centriod)**2 + (z_cc(k) - z_centriod)**2)
 
        ! assigning the reference map to the q_prim vector field
        do ii = 1, num_dims
            q_prim_vf(ii + xibeg - 1)%sf(i, j, k) = xi_cart(ii)
        end do
!           q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(patch_id)%alpha_rho(1)
!           q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(patch_id)%alpha_rho(2)
!what about alpha_rho(3) and alpha_rho(4) [contxb + 1, contxe - 1]?
!           q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(patch_id)%pres
!           q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(patch_id)%alpha(1)
!           q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(patch_id)%alpha(2)

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
