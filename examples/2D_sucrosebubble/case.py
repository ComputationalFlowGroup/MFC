#!/usr/bin/env python3

import math
import json

ps  = 248758.567
gam = 1.4
rho = 1.
c_l = math.sqrt( 1.4*ps/rho )
vel = 0

leng = 1.
Ny = 100.
Nx = Ny*3
dx = leng/Nx

time_end = 3.E1         #5*leng/vel
cfl = 0.1

dt = cfl * dx/c_l 
Nt = int(time_end/dt)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 :  -leng/2.,
    'x_domain%end'                 :  leng/2+2*leng,
    'y_domain%beg'                 :  -leng/2.,
    'y_domain%end'                 :  leng/2.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : Nt,
    't_step_save'                  : int(Nt/20.),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 5,
    'alt_soundspeed'               : 'F',
    'hypoplasticity'               : 'T',
    'num_fluids'                   : 2,
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -6,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -6,
    'bc_y%end'                     : -6,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================
                                                                
    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 10*leng,
    'patch_icpp(1)%length_y'       : leng,
    'patch_icpp(1)%vel(1)'         : vel,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%pres'           : 8.9088E-9,
    'patch_icpp(1)%alpha_rho(1)'   : 1.E0,
    'patch_icpp(1)%alpha_rho(2)'   : 0.E+00,
    'patch_icpp(1)%alpha(1)'       : 1.E+00,
    'patch_icpp(1)%alpha(2)'       : 0.E+00,
    # ==========================================================================

    # Patch 2: Shocked state ===================================================
    #'patch_icpp(2)%geometry'       : 3,
    #'patch_icpp(2)%alter_patch(1)' : 'T',
    #'patch_icpp(2)%x_centroid'     : -3*leng/8.,
    #'patch_icpp(2)%y_centroid'     : 0.,
    #'patch_icpp(2)%length_x'       : leng/4.,
    #'patch_icpp(2)%length_y'       : leng,
    #'patch_icpp(2)%vel(1)'         : vel,
    #'patch_icpp(2)%vel(2)'         : 0.E+00,
    #'patch_icpp(2)%pres'           : ps,
    #'patch_icpp(2)%alpha_rho(1)'   : 2.4,
    #'patch_icpp(2)%alpha_rho(2)'   : 0.E+00,
    #'patch_icpp(2)%alpha(1)'       : 1.E+00,
    #'patch_icpp(2)%alpha(2)'       : 0.E+00,
    # ==========================================================================

    # Patch 3: Bubble ==========================================================
    'patch_icpp(2)%geometry'       : 2,
    'patch_icpp(2)%x_centroid'     : 0.E+00,
    'patch_icpp(2)%y_centroid'     : 0.E+00,
    'patch_icpp(2)%radius'         : leng/5,                            # was leng/5
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%pres'           : 8.9088E-9,
    'patch_icpp(2)%alpha_rho(1)'   : 0.E+00,
    'patch_icpp(2)%alpha_rho(2)'   : 9.82E-7,
    'patch_icpp(2)%alpha(1)'       : 0.E+00,
    'patch_icpp(2)%alpha(2)'       : 1.E+00,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.09E0,                           # 1.E+00/(1.4E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.9556E0,                         # isothermal bulk modulus
    'fluid_pp(2)%gamma'            : 1.4E0,                            # 1.E+00/(1.6666E+00-1.E+00),
    'fluid_pp(2)%pi_inf'           : 6.747E-6,                         # 0.0
    'fluid_pp(1)%qv'               : 3.75E0,                           # K'_theta0 for sucrose
    'fluid_pp(2)%qv'               : 2.0E0,                            #    
    'fluid_pp(1)%G'                : 0.5733E0,                         # Shear modulus
    'fluid_pp(2)%G'                : 1.0E-9,                           # Shear modulus of air taken to be a very small value
    'fluid_pp(1)%ein_cv(1)'        : 0.094647E0,                       # Can be replaced with fluid_pp(:)%cv at some point
    'fluid_pp(2)%ein_cv(1)'        : 0.0226E0,                         # 
    'fluid_pp(1)%ein_cv(2)'        : 3.775E0,                          # Can be replaced with a scalar theta_E at some point
    'fluid_pp(2)%ein_cv(2)'        : 0.335E0,
    'fluid_pp(1)%mg_a'             : 1.E0,                             #a_mg
    'fluid_pp(1)%mg_b'             : 0.E0,                             #b_mg
    'fluid_pp(2)%mg_a'             : 0.E0,                             #a_mg
    'fluid_pp(2)%mg_b'             : 1.E0,                             #b_mg
    'fluid_pp(1)%rho0'             : 1.E0,                             #Non-dimensional initial density in Birch-Murnaghan cold curve
    'fluid_pp(2)%rho0'             : 9.82E-7,
    'fluid_pp(1)%jcook(1)'         : 0.0334,                           # A, Static yield strength
    'fluid_pp(1)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
    'fluid_pp(1)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
    'fluid_pp(1)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
    'fluid_pp(1)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
    'fluid_pp(1)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient pressure
    'fluid_pp(1)%jcook(7)'         : 1.0E7,                            # Limiting strain-rate
    'fluid_pp(1)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
    'fluid_pp(1)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
    'fluid_pp(1)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limit
    'fluid_pp(2)%jcook(1)'         : 0.0334,                           # A, Static yield strength
    'fluid_pp(2)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
    'fluid_pp(2)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
    'fluid_pp(2)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
    'fluid_pp(2)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
    'fluid_pp(2)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient pressure
    'fluid_pp(2)%jcook(7)'         : 1.0E7,                            # Limiting strain-rate
    'fluid_pp(2)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
    'fluid_pp(2)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
    'fluid_pp(2)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limit
}))

# ==============================================================================
