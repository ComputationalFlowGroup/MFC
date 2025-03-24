#!/usr/bin/env python3
import math
import json
import numpy as np

# import scipy as sc

# Non-dimensional input file
# ps  = 248758.567
gam = 1.4
rho = 1.0
# c_l = math.sqrt( 1.4*ps/rho )
c_l = 3.0776  # mm/us

leng = 2.0  # mm
# Ny = 768
# Nx = 1024
Ny = 50
Nx = 50
dx = leng / Nx

time_end = 0.001  # 5*leng/vel
cfl = 0.1
# dt = 1.0E-12
dt = cfl * dx / c_l
# Nt = int(time_end/dt)
Nt = 20000
eps = 1e-6

# Material parameters of sucrose (dimensional)
Kt0_suc = 14.3e9  # Pa
Kt0_prime_suc = 3.75  # -
rho_0_suc = 1.5805e3  # kg/m^3
einstein_cv1_suc = 3279  # J/Kg-K
einstein_cv2_suc = 1125  # K
G_suc = 8.58e9  # Pa
c_0 = 3077.6  # m/s
theta_0_suc = 298  # K
gamma_suc = 1.09

# Material parameters for air (dimensional)
Kt0_air = 1.013e5  # Pa
Kt0_prime_air = 2.00  # Pa
rho_0_air = 1.2  # kg/m^3
theta_0_air = theta_0_suc
einstein_cv1_air = 718  # J/Kg-K
einstein_cv2_air = 100  # K
G_air = 0  # for now
gamma_air = 0.4  # n-1

# Initial condition
theta_0 = 310  # K
P_0 = 1.0040000406510039e-004  # GPa
compression_ratio = 1.2  # rho/rho_0 in the shocked region
rho_suc = 1580.5  # kg/m^3
vel0 = 1.0e-6  # For seeding everything with some non-zero velocity

# Shock EoS
# xi = 1.0-1.0/tilde_rho
# ps = tilde_P_0 +  xi/pow(1-2.87*xi,2)       #Need to fix the shock slope everywhere in the code
# vel = vel0 +  math.sqrt((ps-tilde_P_0)*xi/rho_0_tilde)
# print('pressure=',ps*rho_0_suc*c_0*c_0*1.0E-9,'GPa')
# print('vel=',vel*c_0,'m/s')
vel = 0.2  # mm/us
Us = c_l + 1.104 * vel  # mm/us
rho1 = 1.5805106 * Us / (Us - vel)  # density
ps = P_0 + 1.5805106 * Us * vel  # GPa

# print(vel)
# print(Us)
# print(rho1)
# print(ps)
# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -2 * leng,
            "x_domain%end": 3 * leng,
            "y_domain%beg": -1.5 * leng,
            "y_domain%end": 1.5 * leng,
            "stretch_x": "T",
            "a_x": 20,
            "x_a": -1 * leng,
            "x_b": 2 * leng,
            "stretch_y": "T",
            "a_y": 20,
            "y_a": -3 * leng / 4,
            "y_b": 3 * leng / 4,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            #'dt'                           : dt,
            #'t_step_start'                 : 0,
            #'t_step_stop'                  : 9000,
            #'t_step_save'                  : int(50),
            "cfl_adap_dt": "T",
            "cfl_target": 0.3,
            "n_start": 0,
            "t_save": 4.0e-03,
            "t_stop": 4.0e-02,
            # Simulation Algorithm Parameters
            "num_patches": 3,  # change this to 3 for shocked state
            "model_eqns": 5,
            "alt_soundspeed": "F",
            "hypoplasticity": "F",
            "MGEoS_model": 1,
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-40,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # "E_wrt": "T",
            # "c_wrt": "T",
            # "schlieren_wrt": "T",
            # "schlieren_alpha(1)": 1.00,
            # "schlieren_alpha(2)": 0.25,
            "fd_order": 2,
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * leng,
            "patch_icpp(1)%length_y": 10 * leng,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": P_0,
            "patch_icpp(1)%alpha_rho(1)": (1.0e00 - eps) * 1580.5106e-3,
            "patch_icpp(1)%alpha_rho(2)": eps * 0.00116,
            "patch_icpp(1)%alpha(1)": (1.0e00 - eps),
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -2.5625 * leng,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": 4.0 * leng,
            "patch_icpp(2)%length_y": 10 * leng,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": (1.0e0 - eps) * rho1,
            "patch_icpp(2)%alpha_rho(2)": eps * 1.16e-3,
            "patch_icpp(2)%alpha(1)": 1.0e00 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Patch 3: Bubble
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": 0.0e00,
            "patch_icpp(3)%y_centroid": 0.0e00,
            "patch_icpp(3)%radius": leng / 2.0,
            "patch_icpp(3)%alter_patch(1)": "T",
            # 'patch_icpp(3)%smoothen'       : 'T',
            # 'patch_icpp(3)%smooth_coeff'   : .25,
            # 'patch_icpp(3)%smooth_patch_id' : 1,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": P_0,
            "patch_icpp(3)%alpha_rho(1)": eps * 1580.5106e-3,
            "patch_icpp(3)%alpha_rho(2)": (1.0e0 - eps) * 1.16e-3,
            "patch_icpp(3)%alpha(1)": eps,
            "patch_icpp(3)%alpha(2)": 1.0e00 - eps,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.09,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": 0.0,  # p0
            "fluid_pp(1)%mg_a": 3077.6e-3,  # c0
            "fluid_pp(1)%mg_b": 1.104,  # s
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%rho0": 1580.5e-3,  # reference density
            "fluid_pp(1)%cv": 3000e-6,  # specific heat capacity
            "fluid_pp(2)%gamma": 0.4,  # Gruneisen constant
            "fluid_pp(2)%pi_inf": 0.0,  # p0
            "fluid_pp(2)%mg_a": 0.0,  # c0
            "fluid_pp(2)%mg_b": 0.0,  # s
            "fluid_pp(2)%qv": 0.0,  # e0
            "fluid_pp(2)%qvp": 0.0,  # Gruneisen exponent
            "fluid_pp(2)%rho0": 1.16e-3,  # reference density
            "fluid_pp(2)%cv": 718e-6,  # specific heat capacity
            "fluid_pp(1)%jcook(11)": 298.0,
            "fluid_pp(2)%jcook(11)": 0.0,
        }
    )
)
