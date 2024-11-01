#!/usr/bin/python
import math
import json

#Numerical setup
c_l     = 3077.6
Nx      = 99
cfl     = 0.1
leng    = 2.
dx      = leng/(Nx+1)
mydt    = cfl*dx/c_l
Tend    = 2E-3
Nt      = int(Tend/mydt)
#mydt   = Tend/(1.*Nt)
vel1    = 1.0
vel2    = 0.0
theta_0 = 298.0

eps = 0.0
Kt0_suc       = 14.3e9      #Pa
Kt0_prime_suc = 3.75        #-
rho_0_suc     = 1580.5      #kg/m^3
ein_cv1_suc   = 3279        #J/Kg-K
ein_cv2_suc   = 1125        #K
G_suc         = 8.58e9      #Pa
c_0           = 3077.6      #m/s
theta_0_suc   = 298         #K
gamma_suc     = 1.09

#Initial condition
theta_0           = 298             #K
P_0               = 1.0E5           #Pa
compression_ratio = 1.1             #rho/rho_0 in the shocked region
rho_0             = 1580.5          #kg/m^3

# Configuring case dictionary
print(json.dumps({
                    # Logistics ================================================
                    'run_time_info'                : 'T',
                    # ==========================================================

                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,
                    'x_domain%end'                 : 2.E+00,
                    'm'                            : Nx,
                    'n'                            : 0,
                    'p'                            : 0,
                    'dt'                           : mydt,
                    't_step_start'                 : 0,
                    't_step_stop'                  : int(Nt),
                    't_step_save'                  : int(math.ceil((Nt)/100.)),
		    # ==========================================================

                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 3,
                    'model_eqns'                   : 5,
                    'alt_soundspeed'               : 'F',
                    'num_fluids'                   : 2,
		            'mpp_lim'                      : 'F',
		            'mixture_err'                  : 'F',
		            'time_stepper'                 : 3,
                    'weno_order'                   : 5,
                    'weno_eps'                     : 1.E-16,
		            'weno_Re_flux'                 : 'F',
                    'weno_avg'                     : 'F',
                    'mapped_weno'                  : 'T',
                    'null_weights'                 : 'T',
                    'mp_weno'                      : 'T',
		            'riemann_solver'               : 2,
                    'wave_speeds'                  : 1,
                    'avg_state'                    : 2,
                    'bc_x%beg'                     : -3,
                    'bc_x%end'                     : -3,
                    # ==========================================================
                    # Hypoplasticity ================================
                     'hypoplasticity'              : 'F',
                    # ==========================================================
                     'MGEoS_model'                 : 1,
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,
                    'precision'                    : 2,
                    'prim_vars_wrt'                :'T',
		            'parallel_io'                  :'T',
		            # Patch 1 L ================================================
                    'patch_icpp(1)%geometry'       : 1,
                    'patch_icpp(1)%x_centroid'     : 1.0,
                    'patch_icpp(1)%length_x'       : leng,
                    'patch_icpp(1)%vel(1)'         : 0.0,
                    'patch_icpp(1)%pres'           : P_0,
                    'patch_icpp(1)%alpha_rho(1)'   : (1.0-eps)*1580.5,
                    'patch_icpp(1)%alpha_rho(2)'   : (eps)*1.2,
                    'patch_icpp(1)%alpha(1)'       : 1.0-eps,
                    'patch_icpp(1)%alpha(2)'       : eps,
                    # Shocked State ============================================
                    'patch_icpp(2)%geometry'       : 1,
                    'patch_icpp(2)%x_centroid'     : 0.125,
                    'patch_icpp(2)%length_x'       : 0.25,
                    'patch_icpp(2)%alter_patch(1)' : 'T',
                    'patch_icpp(2)%vel(1)'         : 59.337,
                    'patch_icpp(2)%pres'           : 303.804E6,
                    'patch_icpp(2)%alpha_rho(1)'   : (1.E0-eps)*1610,
                    'patch_icpp(2)%alpha_rho(2)'   : eps*1.5,
                    'patch_icpp(2)%alpha(1)'       : 1.0-eps,
                    'patch_icpp(2)%alpha(2)'       : eps,
                    # Patch 2 R ================================================
                    'patch_icpp(3)%geometry'       : 1,
                    'patch_icpp(3)%x_centroid'     : 1.0,
                    'patch_icpp(3)%length_x'       : 0.25,
                    'patch_icpp(3)%alter_patch(1)' : 'T',
                    'patch_icpp(3)%vel(1)'         : 0.0,
                    'patch_icpp(3)%pres'           : P_0,
                    'patch_icpp(3)%alpha_rho(1)'   : eps*1580.5,
                    'patch_icpp(3)%alpha_rho(2)'   : (1.E0-eps)*1.2,
                    'patch_icpp(3)%alpha(1)'       : eps,
                    'patch_icpp(3)%alpha(2)'       : 1.0-eps,
                    # ==========================================================
                    # Fluids Physical Parameters ===============================================
                    'fluid_pp(1)%gamma'            : 1.09,              # Gruneisen constant
                    'fluid_pp(1)%pi_inf'           : P_0,               # p0
                    'fluid_pp(1)%mg_a'             : 3077.6,            # c0
                    'fluid_pp(1)%mg_b'             : 2.71,              # s
                    'fluid_pp(1)%qv'               : 0.0,               # e0
                    'fluid_pp(1)%qvp'              : 1.0,               # Gruneisen exponent
                    'fluid_pp(1)%rho0'             : 1580.5,            # reference density
                    'fluid_pp(1)%cv'               : 3000,              # specific heat capacity
                    'fluid_pp(2)%gamma'            : 0.4,               # Gruneisen constant
                    'fluid_pp(2)%pi_inf'           : P_0,               # p0
                    'fluid_pp(2)%mg_a'             : 233,               # c0
                    'fluid_pp(2)%mg_b'             : 1.058,             # s
                    'fluid_pp(2)%qv'               : 0.0,               # e0
                    'fluid_pp(2)%qvp'              : 1E-4,               # Gruneisen exponent
                    'fluid_pp(2)%rho0'             : 1.2,               # reference density
                    'fluid_pp(2)%cv'               : 1000,              # specific heat capacity
}))
	            # ==========================================================
# ==============================================================================
