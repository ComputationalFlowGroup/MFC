#!/usr/bin/python
import math
import json

# Numerical setup
Nx = 256
dx = 1.0 / (1.0 * (Nx + 1))
c0 = 4770
Tend = 85.0e-6
cfl = 0.1
mydt = cfl * dx / c0
Nt = int(Tend / mydt)

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
            "t_step_save": int(math.ceil(Nt / 50.0)),
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 5,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Turning on Hypoelasticity
            "hypoplasticity": "F",
            "MGEoS_model": 1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.2,
            "patch_icpp(1)%length_x": 0.4,
            "patch_icpp(1)%vel(1)": 543,
            "patch_icpp(1)%pres": 3e10,
            "patch_icpp(1)%alpha_rho(1)": 11042 * (1.0 - (1.0e-6)),
            "patch_icpp(1)%alpha(1)": (1.0 - (1.0e-6)),
            "patch_icpp(1)%alpha_rho(2)": 2260 * (1.0e-6),
            "patch_icpp(1)%alpha(2)": 1.0e-6,
            # Patch 2
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.2,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 9961 * (1.0 - (1.0e-6)),
            "patch_icpp(2)%alpha(1)": 1.0 - (1.0e-6),
            "patch_icpp(2)%alpha_rho(2)": 2260 * (1.0e-6),
            "patch_icpp(2)%alpha(2)": 1.0e-6,
            # Patch 3
            "patch_icpp(3)%geometry": 1,
            "patch_icpp(3)%x_centroid": 0.8,
            "patch_icpp(3)%length_x": 0.4,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%pres": 0.0,
            "patch_icpp(3)%alpha_rho(1)": 9961 * (1.0e-6),
            "patch_icpp(3)%alpha(1)": (1.0e-6),
            "patch_icpp(3)%alpha_rho(2)": 2260 * (1.0 - (1.0e-6)),
            "patch_icpp(3)%alpha(2)": 1.0 - (1.0e-6),
            # Fluids Physical Parameters for Molybdenum
            "fluid_pp(1)%gamma": 2.56,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": 0.0,  # p0
            "fluid_pp(1)%mg_a": 4770,  # c0
            "fluid_pp(1)%mg_b": 1.43,  # s
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%rho0": 9961,  # reference density
            # Fluids Physical Parameters for MORB
            "fluid_pp(2)%gamma": 1.18,  # Gruneisen constant
            "fluid_pp(2)%pi_inf": 0.0,  # p0
            "fluid_pp(2)%mg_a": 2100,  # c0
            "fluid_pp(2)%mg_b": 1.68,  # s
            "fluid_pp(2)%qv": 0.0,  # e0
            "fluid_pp(2)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(2)%rho0": 2660,  # reference density
        }
    )
)
