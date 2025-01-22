#!/usr/bin/python
import math
import json

# Numerical setup
Nx = 99
dx = 1.0 / (1.0 * (Nx + 1))
c0 = 5328
Tend = 50.0e-6
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
            "t_step_save": int(math.ceil(Nt / 100.0)),
            # Simulation Algorithm Parameters
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
            "mapped_weno": "F",
            "null_weights": "T",
            "mp_weno": "T",
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
            # Patch 1 L
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 2726.91,
            "patch_icpp(1)%pres": 6.816e10,
            "patch_icpp(1)%alpha_rho(1)": 4000,
            "patch_icpp(1)%alpha(1)": 1.0,
            #'patch_icpp(1)%tau_e(1)'       : 1.E-16,
            # Patch 2 R
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 2785,
            "patch_icpp(2)%alpha(1)": 1.0,
            #'patch_icpp(2)%tau_e(1)'       : 1.E-16,
            # Fluids Physical Parameters
            # Johnson-cook material parameters are taken from https://doi.org/10.1115/1.4027793
            # Parameters are for Aluminum A356
            # Simon-Glatzel parameters taken from https://doi.org/10.1016/S0925-8388(00)00736-2
            "fluid_pp(1)%gamma": 2.0,  # Gruneisen constant
            "fluid_pp(1)%pi_inf": 0.0,  # p0
            "fluid_pp(1)%mg_a": 5328.0,  # c0
            "fluid_pp(1)%mg_b": 1.338,  # s
            "fluid_pp(1)%qv": 0.0,  # e0
            "fluid_pp(1)%qvp": 1.0,  # Gruneisen exponent
            "fluid_pp(1)%rho0": 2785.0,  # reference density
            "fluid_pp(1)%cv": 903.0,  # specific heat capacity
            "fluid_pp(1)%G": 25.0e9,  # 25 GPa shear modulus
            "fluid_pp(1)%jcook(1)": 270.0e6,  # A, Static yield strength
            "fluid_pp(1)%jcook(2)": 155.0e6,  # B, Strain-Hardening coefficient
            "fluid_pp(1)%jcook(3)": 0.28,  # n, Strain-Hardening exponent
            "fluid_pp(1)%jcook(4)": 0.018,  # C, Strain-rate hardening coefficient
            "fluid_pp(1)%jcook(5)": 1.43,  # m, Thermal softening exponent
            "fluid_pp(1)%jcook(6)": 746.0,  # theta_m, Melt temperature at ambient #pressure
            "fluid_pp(1)%jcook(7)": 1.0e7,  # Limiting strain-rate
            "fluid_pp(1)%jcook(8)": 7.6e9,  # Parameter in Simon-Glatzel melt relation
            "fluid_pp(1)%jcook(9)": 1 / 0.615,  # exponent in Simon-Glatzel melt relation
            "fluid_pp(1)%jcook(10)": 1.0,  # non-dimensional strain-rate limit
            "fluid_pp(1)%jcook(11)": 298.0,  # Reference temperature
        }
    )
)
