#!/usr/bin/env python3
import json
import math

patmos = 101325.0
eps = 1.0e-6
# water
rho_w = 1000.0
gamma_w = 7.12
K_w = 2.2e9
B_w = K_w / gamma_w - patmos
c_w = math.sqrt(gamma_w * (patmos + B_w) / rho_w)

# graded
E0 = 1.0e9
Ca_inv = E0 / 2
rho_s = 1000.0
gamma_s = gamma_w
K_s = K_w
B_s = K_s / gamma_s - patmos

# shock
del_P_s = 200000
P_amp = del_P_s + patmos

# Characteristic scales
rho_char = rho_w
c_ref = c_w
stress_char = rho_char * c_ref**2 / gamma_w
u_char = math.sqrt(stress_char / rho_char)
L_char = 1.0
t_char = L_char / u_char

# Nondim scales
# gamma_w fine
B_w_n = B_w / stress_char
p_r_n = patmos / stress_char
p_l_n = P_amp / stress_char
rho_w_n = rho_w / rho_char
u_r_n = 0.0
c_w_n = math.sqrt(gamma_w * (p_r_n + B_w_n) / rho_w_n)
rho_s_n = rho_s / rho_char
B_s_n = B_s / stress_char
Ca_inv_nd = Ca_inv / stress_char
Ca_inv_init = Ca_inv_nd * 0.00001
Ca_inv_end = Ca_inv_nd * 2

# Numerical setup
Nx = 300
Ny = 50
Nz = 50
Lx = 1.0 * 4
dx = Lx / Nx
Ly = Ny * dx
Lz = Nz * dx
dy = Ly / Ny
dz = Lz / Nz
cfl = 0.1

# shock
g_r = gamma_w
g_l = gamma_w
p_r = p_r_n
p_l = p_l_n
B_r = B_w_n
r_r = rho_w_n
c_r = c_w_n

rho_shock_l = (((g_r + 1) / (g_r - 1)) * ((p_l + B_r) / (p_r + B_r)) + 1) / ((g_r + 1) / (g_r - 1) + (p_l + B_r) / (p_r + B_r)) * r_r
shock_l = (c_r / g_r) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) / (math.sqrt(((g_r + 1) / (2 * g_r)) * ((p_l / p_r) - 1) * (p_r / (p_r + B_r)) + 1))
s_nd = c_r * math.sqrt((g_r + 1.0) / (2.0 * g_r) * (((p_l + B_r) / (p_r + B_r)) - 1.0) + 1.0)

# Probes
n_probes = 10
probe_x = [(i + 1) * Lx / (n_probes + 1) for i in range(n_probes)]

# End time (travel time of the wave from the patch interface past the last probe)
x0 = 0.1 * Lx
x_stop = probe_x[-1] + 0.05 * Lx
t_end = 1.1 * (x_stop - x0) / s_nd
dt = 0.0005
Nt = int(math.ceil(t_end / dt))
t_step_save = max(1, Nt // 30)
case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "z_domain%beg": 0.0,
    "z_domain%end": Lz,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": Nz - 1,
    "dt": dt,
    "cfl_target": cfl,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": t_step_save,
    "num_patches": 3,
    "model_eqns": 2,  # 2 for hypo 3 for hyper
    "alt_soundspeed": "F",
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 1,  # 2 for hyper 1 for hypo
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "probe_wrt": "T",
    "fd_order": 4,
    "num_probes": n_probes,
    # Patch 1: post-shock
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.05 * Lx,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%z_centroid": Lz / 2,
    "patch_icpp(1)%length_x": 0.1 * Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": Lz,
    "patch_icpp(1)%vel(1)": shock_l,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_l_n,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_shock_l,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_s_n,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    # Patch 2: unshocked (water)
    "patch_icpp(2)%geometry": 9,
    "patch_icpp(2)%x_centroid": 0.55 * Lx,
    "patch_icpp(2)%y_centroid": Ly / 2,
    "patch_icpp(2)%z_centroid": Lz / 2,
    "patch_icpp(2)%length_x": 0.9 * Lx,
    "patch_icpp(2)%length_y": Ly,
    "patch_icpp(2)%length_z": Lz,
    "patch_icpp(2)%vel(1)": u_r_n,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_r_n,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_w_n,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_s_n,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    # Patch 3: graded
    "patch_icpp(3)%geometry": 9,
    "patch_icpp(3)%x_centroid": 0.5 * Lx,
    "patch_icpp(3)%y_centroid": Ly / 2,
    "patch_icpp(3)%z_centroid": Lz / 2,
    "patch_icpp(3)%length_x": 0.3 * Lx,
    "patch_icpp(3)%length_y": Ly,
    "patch_icpp(3)%length_z": Lz,
    "patch_icpp(3)%vel(1)": u_r_n,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%vel(3)": 0.0,
    "patch_icpp(3)%pres": p_r_n,
    "patch_icpp(3)%alpha_rho(1)": eps * rho_w_n,
    "patch_icpp(3)%alpha_rho(2)": (1.0 - eps) * rho_s_n,
    "patch_icpp(3)%alpha(1)": eps,
    "patch_icpp(3)%alpha(2)": 1.0 - eps,
    "patch_icpp(3)%alter_patch(2)": "T",
    "patch_icpp(3)%tau_e(1)": 0,
    "patch_icpp(3)%tau_e(2)": 0,
    "patch_icpp(3)%tau_e(3)": 0,
    "patch_icpp(3)%tau_e(4)": 0,
    "patch_icpp(3)%tau_e(5)": 0,
    "patch_icpp(3)%tau_e(6)": 0,
    # Fluid
    "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
    "fluid_pp(1)%pi_inf": gamma_w * B_w_n / (gamma_w - 1.0),
    "fluid_pp(1)%Ca_inv": 0,
    "fluid_pp(2)%gamma": 1.0 / (gamma_s - 1.0),
    "fluid_pp(2)%pi_inf": gamma_s * B_s_n / (gamma_s - 1.0),
    "fluid_pp(2)%graded_Ca_inv": "T",
    "fluid_pp(2)%graded_Re": "F",
    "fluid_pp(2)%graded_type": 1,
    "fluid_pp(2)%graded_profile": 1,
    "fluid_pp(2)%graded_Ca_inv_init": Ca_inv_init,
    "fluid_pp(2)%graded_Ca_inv_end": Ca_inv_end,
    "fluid_pp(2)%graded_beg_loc_x": 0.35 * Lx,
    "fluid_pp(2)%graded_beg_loc_y": Ly / 2,
    "fluid_pp(2)%graded_beg_loc_z": Lz / 2,
    "fluid_pp(2)%graded_end_loc_x": 0.65 * Lx,
    "fluid_pp(2)%graded_end_loc_y": Ly / 2,
    "fluid_pp(2)%graded_end_loc_z": Lz / 2,
    "graded": "F",
    "fluid_pp(2)%Ca_inv": Ca_inv_init,
    "hyperelasticity": "T",
    "rmt": "T",
    "hyper_model": 1,
}

for i, x in enumerate(probe_x):
    case[f"probe({i+1})%x"] = x
    case[f"probe({i+1})%y"] = Ly / 2
    case[f"probe({i+1})%z"] = Lz / 2

print(json.dumps(case))
