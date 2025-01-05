#!/usr/bin/python
import math
import json

# Numerical setup
c_l = 1
Nx = 256
cfl = 0.1
leng = 1.0
dx = leng / (Nx + 1)
mydt = cfl * dx / c_l
Tend = 4.001
Nt = int(Tend / mydt)
# mydt   = Tend/(1.*Nt)
vel1 = 1.0
vel2 = 0.0
theta_0 = 298.0
eps = 1e-16

Kt0_suc = 14.3e9  # Pa
Kt0_prime_suc = 3.75  # -
rho_0_suc = 8924  # kg/m^3
einstein_cv1_suc = 3279  # J/Kg-K
einstein_cv2_suc = 1125  # K
G_suc = 8.58e9  # Pa
c_0 = 3910  # m/s
theta_0_suc = 298  # K
gamma_suc = 1.09

# Initial condition
theta_0 = 298  # K
P_0 = 1.000e5  # Pa
compression_ratio = 1.1  # rho/rho_0 in the shocked region
rho_0 = 1580.5  # kg/m^3

tilde_P0 = P_0 / (rho_0_suc * c_0 * c_0)

tilde_P_0 = P_0 / (rho_0_suc * c_0 * c_0)
tilde_rho = compression_ratio
Kt0_tilde = Kt0_suc / (rho_0_suc * c_0 * c_0)
A_tilde = einstein_cv1_suc * theta_0 / (c_0 * c_0)
theta_E_tilde = einstein_cv2_suc / theta_0
rho_0_tilde = rho_0 / rho_0_suc

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
            "num_fluids": 2,
            "mpp_lim": "T",
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
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Hypoplasticity
            "hypoplasticity": "F",
            # Mie-Gruneisen EoS model
            "MGEoS_model": 1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 L
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": leng,
            "patch_icpp(1)%vel(1)": vel1,
            "patch_icpp(1)%pres": tilde_P0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - eps),
            "patch_icpp(1)%alpha_rho(2)": eps * (1.2 / 8924),
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2 R
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": vel1,
            "patch_icpp(2)%pres": tilde_P0,
            "patch_icpp(2)%alpha_rho(1)": eps,
            "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * 1.2 / 8924,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": 1.0 - eps,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.96e0,  # 1.E+00/(1.4E+00-1.E+00),
            "fluid_pp(2)%gamma": 0.4e0,  # 1.E+00/(1.6666E+00-1.E+00),
            "fluid_pp(1)%pi_inf": tilde_P0,  # isothermal bulk modulus
            "fluid_pp(2)%pi_inf": 0.0 * 1.013e5 / (rho_0_suc * c_0 * c_0),  # 0.0
            "fluid_pp(1)%mg_a": 1.0,  # K'_theta0 for sucrose
            "fluid_pp(1)%mg_b": 1.51,  #
            "fluid_pp(2)%mg_a": 0.0,  # 237/3910,
            "fluid_pp(2)%mg_b": 1.058,
            "fluid_pp(1)%qv": tilde_P0,
            "fluid_pp(2)%qv": 0.0,
            "fluid_pp(1)%qvp": 1.0,
            "fluid_pp(2)%qvp": 0.0,
            "fluid_pp(1)%rho0": 1.0,  # Non-dimensional initial density in Birch-Murnaghan cold curve
            "fluid_pp(2)%rho0": 1.2 / 8924,
            # 'fluid_pp(2)%rho0'             : 0.429683187295/1580.5,
            # 'fluid_pp(1)%jcook(1)'         : 0.0334,                           # A, Static yield strength
            # 'fluid_pp(1)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
            # 'fluid_pp(1)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
            # 'fluid_pp(1)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
            # 'fluid_pp(1)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
            # 'fluid_pp(1)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient #pressure
            # 'fluid_pp(1)%jcook(7)'         : 1.0E7,                            # Limiting strain-rate
            # 'fluid_pp(1)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
            # 'fluid_pp(1)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
            # 'fluid_pp(1)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limit
            # 'fluid_pp(1)%jcook(11)'        : 298/theta_0,                      # Reference temperature
            # 'fluid_pp(2)%jcook(1)'         : 0.0334,                           # A, Static yield strength
            # 'fluid_pp(2)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
            # 'fluid_pp(2)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
            # 'fluid_pp(2)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
            # 'fluid_pp(2)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
            # 'fluid_pp(2)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient pressure
            # 'fluid_pp(2)%jcook(7)'         : 1.0E7,                            # Limiting strain-rate
            # 'fluid_pp(2)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
            # 'fluid_pp(2)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
            # 'fluid_pp(2)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limitI
            # 'fluid_pp(2)%jcook(11)'        : 298/theta_0,                      # non-dimensionalized Reference temperature
        }
    )
)
