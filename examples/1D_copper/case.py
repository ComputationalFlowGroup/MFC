#!/usr/bin/python
import math
import json

# Numerical setup
c_l = 3910
Nx = 256
cfl = 0.1
leng = 1.0
dx = leng / (Nx + 1)
mydt = cfl * dx / c_l
Tend = 100.0e-6
Nt = int(Tend / mydt)
# mydt   = Tend/(1.*Nt)
vel1 = 0.0
vel2 = 0.0
theta_0 = 298.0

# Initial condition
theta_0 = 298  # K
P_0 = 1.000e5  # Pa
compression_ratio = 1.1  # rho/rho_0 in the shocked region
rho_0 = 8924.0  # kg/m^3
vel0 = 0.0
s = 1.51
tilde_rho = compression_ratio

xi = 1.0 - 1.0 / tilde_rho
ps = P_0 + c_l * c_l * rho_0 * xi / ((1.0 - s * xi) * (1.0 - s * xi))
vel = vel0 + math.sqrt((ps - P_0) * xi / rho_0)

# phi = math.exp(gamma_suc*(1-1/tilde_rho))
# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 100.0)),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 5,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "MGEoS_model": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "mapped_weno": "T",
            "null_weights": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Equilibrium state
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": leng,
            "patch_icpp(1)%vel(1)": vel0,
            "patch_icpp(1)%pres": P_0,
            "patch_icpp(1)%alpha_rho(1)": rho_0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.0625,
            "patch_icpp(2)%length_x": 0.125,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": tilde_rho * rho_0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.96e0,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": 1.0e5,  # p0
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%mg_a": 3910.0,  # c0
            "fluid_pp(1)%mg_b": 1.51,  # s
            "fluid_pp(1)%rho0": 8924.0,  # reference density
            "fluid_pp(1)%cv": 385.0,  # specific heat capacity
        }
    )
)
