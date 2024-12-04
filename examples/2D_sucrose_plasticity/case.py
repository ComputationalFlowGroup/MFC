#!/usr/bin/env python3
import math
import json
import numpy as np
#import scipy as sc

#Non-dimensional input file
c_l = 3.0776       #mm/us

leng = 2.0        #mm
Ny = 768
Nx = 1024
dx = leng/Nx

time_end = 0.001   #5*leng/vel
cfl = 0.1
#dt = 1.0E-12
dt = cfl * dx/c_l
#Nt = int(time_end/dt)
Nt = 20000
eps = 1E-6

#Material parameters of sucrose (dimensional)
Kt0_suc       = 14.3e9      #Pa
Kt0_prime_suc = 3.75        #-
rho_0_suc     = 1.5805E3    #kg/m^3
ein_cv1_suc   = 3279        #J/Kg-K
ein_cv2_suc   = 1125        #K
G_suc         = 8.58e9      #Pa
c_0           = 3077.6      #m/s
theta_0_suc   = 298         #K
gamma_suc     = 1.09

#Material parameters for air (dimensional)
Kt0_air       = 1.013e5     #Pa
Kt0_prime_air = 2.00        #Pa
rho_0_air     = 1.2         #kg/m^3
theta_0_air   = theta_0_suc
ein_cv1_air   = 718         #J/Kg-K
ein_cv2_air   = 100         #K
G_air         = 0           #for now
gamma_air     = 0.4         #n-1

#Initial condition
theta_0           = 310             #K
P_0               = 1.0040000406510039E-004         #GPa
compression_ratio = 1.2             #rho/rho_0 in the shocked region
rho_suc           = 1580.5          #kg/m^3
vel0              = 1.0E-6          #For seeding everything with some non-zero velocity

#Shock EoS
vel = 0.2                        #mm/us
Us = c_l+ 1.104*vel              #mm/us
rho1 = 1.5805106*Us/(Us-vel)     #density
ps = P_0 + 1.5805106*Us*vel      #GPa

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 :  -3*leng/4.,
    'x_domain%end'                 :  leng,
    'y_domain%beg'                 :  -3*leng/4.,
    'y_domain%end'                 :  3*leng/4.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    #'dt'                           : dt,
    #'t_step_start'                 : 0,
    #'t_step_stop'                  : 9000,
    #'t_step_save'                  : int(50),
    'cfl_adap_dt'                  : 'T',
    'cfl_target'                   : 0.2,
    'n_start'                      : 0,
    't_save'                       : 4.E-02,
    't_stop'                       : 4.0,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 3,             #change this to 3 for shocked state
    'model_eqns'                   : 5,
    'alt_soundspeed'               : 'F',
    'hypoplasticity'               : 'T',
    'MGEoS_model'                  : 1,
    'num_fluids'                   : 2,
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-40,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 1,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    'bc_y%beg'                     : -3,
    'bc_y%end'                     : -3,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    'E_wrt'                        :'T',
    'c_wrt'                        :'T',
    'schlieren_wrt'                :'T',
    'schlieren_alpha(1)'           : 1.00,
    'schlieren_alpha(2)'           : 0.25,
    'fd_order'                     : 2,
    # ==========================================================================

    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 10*leng,
    'patch_icpp(1)%length_y'       : 10*leng,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%vel(2)'         : 0.0,
    'patch_icpp(1)%pres'           : P_0,
    'patch_icpp(1)%alpha_rho(1)'   : (1.E+00-eps)*1580.5106E-3,
    'patch_icpp(1)%alpha_rho(2)'   : eps*0.0012,
    'patch_icpp(1)%alpha(1)'       : (1.E+00-eps),
    'patch_icpp(1)%alpha(2)'       : eps,
    # ==========================================================================

    # Patch 2: Shocked state ===================================================
    'patch_icpp(2)%geometry'       : 3,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : -2.5625*leng,
    'patch_icpp(2)%y_centroid'     : 0.,
    'patch_icpp(2)%length_x'       : 4.0*leng,
    'patch_icpp(2)%length_y'       : 10*leng,
    'patch_icpp(2)%vel(1)'         : vel,
    'patch_icpp(2)%vel(2)'         : 0.0,
    'patch_icpp(2)%pres'           : ps,
    'patch_icpp(2)%alpha_rho(1)'   : (1.E0-eps)*rho1,
    'patch_icpp(2)%alpha_rho(2)'   : eps*1.2E-3,
    'patch_icpp(2)%alpha(1)'       : 1.E+00-eps,
    'patch_icpp(2)%alpha(2)'       : eps,
    # ==========================================================================

    # Patch 3: Bubble ==========================================================
    'patch_icpp(3)%geometry'       : 2,
    'patch_icpp(3)%x_centroid'     : 0.E+00,
    'patch_icpp(3)%y_centroid'     : 0.E+00,
    'patch_icpp(3)%radius'         : leng/2.0,
    'patch_icpp(3)%alter_patch(1)' : 'T',
    #'patch_icpp(3)%smoothen'       : 'T',
    #'patch_icpp(3)%smooth_coeff'   : 1.00,
    #'patch_icpp(3)%smooth_patch_id' : 1,
    'patch_icpp(3)%vel(1)'         : 0.0,
    'patch_icpp(3)%vel(2)'         : 0.0,
    'patch_icpp(3)%pres'           : P_0,
    'patch_icpp(3)%alpha_rho(1)'   : eps*1580.5106E-3,
    'patch_icpp(3)%alpha_rho(2)'   : (1.E0-eps)*1.2E-3,
    'patch_icpp(3)%alpha(1)'       : eps,
    'patch_icpp(3)%alpha(2)'       : 1.E+00-eps,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
                    'fluid_pp(1)%gamma'            : 1.09,              # Gruneisen constant
                    'fluid_pp(1)%pi_inf'           : 0.0,               # p0
                    'fluid_pp(1)%mg_a'             : 3077.6E-3,         # c0
                    'fluid_pp(1)%mg_b'             : 1.104,             # s
                    'fluid_pp(1)%qv'               : 0.0,               # e0
                    'fluid_pp(1)%qvp'              : 1.0,               # Gruneisen exponent
                    'fluid_pp(1)%rho0'             : 1580.5E-3,         # reference density
                    'fluid_pp(1)%cv'               : 3000E-6,           # specific heat capacity
                    'fluid_pp(2)%gamma'            : 0.4,               # Gruneisen constant
                    'fluid_pp(2)%pi_inf'           : 0.0,               # p0
                    'fluid_pp(2)%mg_a'             : 0.0,               # c0
                    'fluid_pp(2)%mg_b'             : 0.0,               # s
                    'fluid_pp(2)%qv'               : 0.0,               # e0
                    'fluid_pp(2)%qvp'              : 0.0,               # Gruneisen exponent
                    'fluid_pp(2)%rho0'             : 1.2E-3,            # reference density
                    'fluid_pp(2)%cv'               : 1000E-6,           # specific heat capacity
    # Plasticity parameters =====================================================
                    'fluid_pp(1)%G'                : 8.5E0,               # shear modulus of sucrose
                    'fluid_pp(1)%jcook(1)'         : 0.5E0,               # A, Static yield strength
                    'fluid_pp(1)%jcook(2)'         : 0.5E0,               # B, Strain-Hardening coefficient
                    'fluid_pp(1)%jcook(3)'         : 0.1E0,               # n, Strain-Hardening exponent
                    'fluid_pp(1)%jcook(4)'         : 0.01E0,              # C, Strain-rate hardening #coefficient
                    'fluid_pp(1)%jcook(5)'         : 0.45E0,              # m, Thermal softening exponent
                    'fluid_pp(1)%jcook(6)'         : 459E0,               # theta_m, Melt temperature at ambient pressure
                    'fluid_pp(1)%jcook(7)'         : 10E0,                # Limiting strain-rate
                    'fluid_pp(1)%jcook(8)'         : 0.3E0,               # Parameter in Simon-Glatzel melt relation
                    'fluid_pp(1)%jcook(9)'         : 3.25E0,              # exponent in Simon-Glatzel melt relation
                    'fluid_pp(1)%jcook(10)'        : 1E-6,              # non-dimensional strain-rate limit
                    'fluid_pp(1)%jcook(11)'        : 298E0,               # Reference temperature
                    'fluid_pp(2)%G'                : 0.0E0,               # shear modulus of air
                    'fluid_pp(2)%jcook(1)'         : 0.5E0,               # A, Static yield strength
                    'fluid_pp(2)%jcook(2)'         : 0.5E0,               # B, Strain-Hardening coefficient
                    'fluid_pp(2)%jcook(3)'         : 0.1E0,               # n, Strain-Hardening exponent
                    'fluid_pp(2)%jcook(4)'         : 0.01E0,              # C, Strain-rate hardening coefficient
                    'fluid_pp(2)%jcook(5)'         : 0.45E0,              # m, Thermal softening exponent
                    'fluid_pp(2)%jcook(6)'         : 459E0,               # theta_m, Melt temperature at ambient pressure
                    'fluid_pp(2)%jcook(7)'         : 10E0,                # Limiting strain-rate
                    'fluid_pp(2)%jcook(8)'         : 0.3E0,               # Parameter in Simon-Glatzel melt relation
                    'fluid_pp(2)%jcook(9)'         : 3.25E0,              # exponent in Simon-Glatzel melt relation
                    'fluid_pp(2)%jcook(10)'        : 1E-6,              # non-dimensional strain-rate limitI
                    'fluid_pp(2)%jcook(11)'        : 0E0,               # non-dimensionalized Reference temperature
}))

# ==============================================================================
