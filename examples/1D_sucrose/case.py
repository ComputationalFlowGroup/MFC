#!/usr/bin/python
import math
import json

# Numerical setup
c_l = 3077.6
Nx = 192
cfl = 0.05
leng = 1.0
dx = leng / (Nx + 1)
mydt = cfl * dx / c_l
Tend = 1e-5
Nt = int(Tend / mydt)
# mydt   = Tend/(1.*Nt)
vel0 = 0.0
vel2 = 0.0
theta_0 = 298.0


Kt0_suc = 14.3e9  # Pa
Kt0_prime_suc = 3.75  # -
rho_0_suc = 1580.5  # kg/m^3
einstein_cv1_suc = 3279  # J/Kg-K
einstein_cv2_suc = 1125  # K
G_suc = 8.58e9  # Pa
c_0 = 3077.6  # m/s
theta_0_suc = 298  # K
gamma_suc = 1.09

# Initial condition
theta_0 = 298  # K
P_0 = 1.0e5  # Pa
compression_ratio = 1.2  # rho/rho_0 in the shocked region
rho_0 = 1580.5  # kg/m^3

tilde_P0 = P_0 / (rho_0_suc * c_0 * c_0)

# phi = math.exp(gamma_suc*(1-1/tilde_rho))
# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics ================================================
            "run_time_info": "T",
            # ==========================================================
            # Computational Domain Parameters ==========================
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 100.0)),
            # ==========================================================
            # Simulation Algorithm Parameters ==========================
            "num_patches": 2,
            "model_eqns": 5,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # ==========================================================
            # Hypoplasticity ================================
            "hypoplasticity": "F",
            # ==========================================================
            "MGEoS_model": 1,
            # Formatted Database Files Structure Parameters ============
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # ==========================================================
            # Patch 1 L ================================================
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": leng,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": P_0,
            "patch_icpp(1)%alpha_rho(1)": 1580.5,
            "patch_icpp(1)%alpha(1)": 1.0,
            # ==========================================================
            # shocked state ============================================
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.25,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 59.337,
            "patch_icpp(2)%pres": 303.8043e6,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - 1.0e-7) * 1610,
            "patch_icpp(2)%alpha(1)": 1.0,
            # ==========================================================
            # Fluids Physical Parameters ===============================
            "fluid_pp(1)%gamma": 1.09,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": P_0,  # p0
            "fluid_pp(1)%mg_a": 3077.6,  # c0
            "fluid_pp(1)%mg_b": 2.87,  # s
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%rho0": 1580.5,  # reference density
            "fluid_pp(1)%cv": 3279,  # specific heat capacity
        }
    )
)
#
# ==========================================================
# ==============================================================================
