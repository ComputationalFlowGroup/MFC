#!/usr/bin/python
import math
import json

# Numerical setup
Nx = 640
Ny = 640
leng = 1.0
dx = leng / (1.0 * (Nx + 1))
c0 = 4770
Tend = 85e-6
cfl = 0.1
dt = cfl * dx / c0
Nt = int(Tend / dt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": leng,
            "y_domain%beg": 0.0,
            "y_domain%end": leng,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": int(Nt / 100.0),
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 5,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "T",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Turning on Hypoelasticity
            "hypoplasticity": "F",
            "MGEoS_model": 1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 4,
            "schlieren_wrt": "T",
            # Molybdenum unshocked patch
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 9961 * (1.0 - (1.0e-6)),
            "patch_icpp(1)%alpha(1)": (1.0 - (1.0e-6)),
            "patch_icpp(1)%alpha_rho(2)": 2260 * (1.0e-6),
            "patch_icpp(1)%alpha(2)": 1.0e-6,
            # MORB patch
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.55,
            "patch_icpp(2)%length_x": 0.3,
            "patch_icpp(2)%y_centroid": 0.25,
            "patch_icpp(2)%length_y": 0.5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 9961 * (1.0e-6),
            "patch_icpp(2)%alpha(1)": 1.0e-6,
            "patch_icpp(2)%alpha_rho(2)": 2260 * (1.0 - (1.0e-6)),
            "patch_icpp(2)%alpha(2)": 1.0 - (1.0e-6),
            # Patch 3
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%x_centroid": 0.15,
            "patch_icpp(3)%length_x": 0.3,
            "patch_icpp(3)%y_centroid": 0.5,
            "patch_icpp(3)%length_y": 1.0,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 543.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": 3.0e10,
            "patch_icpp(3)%alpha_rho(1)": 11042 * (1.0 - (1.0e-6)),
            "patch_icpp(3)%alpha(1)": 1.0 - (1.0e-6),
            "patch_icpp(3)%alpha_rho(2)": 2260 * (1.0e-6),
            "patch_icpp(3)%alpha(2)": 1.0e-6,
            # Fluids Physical Parameters for Molybdenum
            "fluid_pp(1)%gamma": 2.56,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": 0.0,  # p0
            "fluid_pp(1)%mg_a": 4770,  # c0
            "fluid_pp(1)%mg_b": 1.43,  # s
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%rho0": 9961,  # reference density
            "fluid_pp(1)%cv": 250,  # specific heat
            # Fluids Physical Parameters for MORB
            "fluid_pp(2)%gamma": 1.18,  # Gruneisen constant
            "fluid_pp(2)%pi_inf": 0.0,  # p0
            "fluid_pp(2)%mg_a": 2100,  # c0
            "fluid_pp(2)%mg_b": 1.68,  # s
            "fluid_pp(2)%qv": 0.0,  # e0
            "fluid_pp(2)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(2)%rho0": 2660,  # reference density
            "fluid_pp(2)%cv": 250,
        }
    )
)
