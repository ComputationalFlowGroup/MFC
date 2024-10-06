#!/usr/bin/env python3
import math
import json
import numpy as np
#import scipy as sc

#ps  = 248758.567
gam = 1.4
rho = 1.
#c_l = math.sqrt( 1.4*ps/rho )
c_l = 3077.6       #m/s
vel1 = 1.E0

leng = 1.
Ny = 256.
Nx = 512.
dx = leng/Nx

time_end = 2.E-1         #5*leng/vel
cfl = 0.1

dt = cfl * dx/c_l
Nt = int(time_end/dt)

#Material parameters of sucrose (dimensional)
Kt0_suc       = 14.3e9      #Pa
Kt0_prime_suc = 3.75        #-
rho_0_suc     = 1.580488803979682E3      #kg/m^3
ein_cv1_suc   = 3279        #J/Kg-K
ein_cv2_suc   = 1125        #K
G_suc         = 8.58e9      #Pa
c_0           = 3077.6      #m/s
theta_0_suc   = 298         #K
gamma_suc     = 1.09

#Material parameters for air (dimensional)
Kt0_air       = 1.013e5     #Pa
Kt0_prime_air = 2.00        #Pa
rho_0_air     = 0.429683187288245      #kg/m^3
theta_0_air   = theta_0_suc
ein_cv1_air   = 718         #J/Kg-K
ein_cv2_air   = 100         #K
G_air         = 0           #for now
gamma_air     = 0.4         #n-1

#Initial condition
theta_0           = 310             #K
P_0               = 1.013E5  #Pa
compression_ratio = 1.1             #rho/rho_0 in the shocked region
rho_suc           = 1580.5          #kg/m^3

#Calculate initial density of air at theta_0 and P_0
rho_air1 = np.linspace(1.168, 15)
P_air_values = []
for rho_air in rho_air1:
    phi_air     = np.exp(gamma_air*(1-rho_0_air/rho_air))
    P_air = Kt0_air*(rho_air/rho_0_air)*math.log(rho_air/rho_0_air)+ gamma_air*rho_0_air*ein_cv1_air*ein_cv2_air*phi_air*(1/(np.exp(phi_air*ein_cv2_air/theta_0)-1)-1/(np.exp(phi_air*ein_cv2_air/theta_0_air)-1))
    P_air_values.append(P_air)

#interp_function = interp1d(P_air_values, rho_air1, kind='linear', fill_value="extrapolate")

rho_air_at_P_0 = np.interp(P_0,P_air_values, rho_air1, left = None,right = None, period =None)
#for rho_air, P_air in zip(rho_air1, P_air_values):
#        print(f"rho_air: {rho_air:.16f}, P_air: {P_air:.16f}")
# Interpolate to find rho_air for P_air = P_0
#rho_air_at_P_0 = interp_function(P_0)
#print(f"rho_air corresponding to P_air = {P_0} is {rho_air_at_P_0}")

#Calculate initial density of suc at theta_0 and P_0
rho_suc1 = np.linspace(1577, 1600)
P_suc_values = []
for rho_suc in rho_suc1:
    phi_suc     = np.exp(gamma_suc*(1-rho_0_suc/rho_suc))
    P_suc = Kt0_suc*(rho_suc/rho_0_suc)*math.log(rho_suc/rho_0_suc)+ gamma_suc*rho_0_suc*ein_cv1_suc*ein_cv2_suc*phi_suc*(1/(np.exp(phi_suc*ein_cv2_suc/theta_0)-1)-1/(np.exp(phi_suc*ein_cv2_suc/theta_0_suc)-1))
    P_suc_values.append(P_suc)

#interp_function = interp1d(P_suc_values, rho_suc1, kind='linear', fill_value="extrapolate")

#for rho_suc, P_suc in zip(rho_suc1, P_suc_values):
#        print(f"rho_suc: {rho_suc:.16f}, P_suc: {P_suc:.16f}")
# Interpolate to find rho_suc for P_suc = P_0
rho_suc_at_P_0 = np.interp(P_0, P_suc_values,rho_suc1, left=None, right =None,period=None)
#rho_suc_at_P_0 = interp_function(P_0)
#print(f"rho_suc corresponding to P_suc = {P_0} is {rho_suc_at_P_0}")
#print(f"tilde_rho_suc_initial is {rho_suc_at_P_0/rho_0_suc}")

#Calculate bulk speed of sound at ambient (used for non-dimensionalization)
c_squared = (Kt0_suc/rho_0_suc)+ gamma_suc*P_0/rho_0_suc + math.pow(gamma_suc,2)*ein_cv1_suc*(math.pow(ein_cv2_suc/theta_0_suc,2))*math.exp(ein_cv2_suc/theta_0)/(math.pow(math.exp(ein_cv2_suc/theta_0_suc)-1,2))

c_0 = math.sqrt(c_squared)
#print(c_0)

#RH jump conditions in non-dimensional form to calculate P, u_p, U_s
tilde_P_0 = P_0/(rho_0_suc*c_0*c_0)
tilde_rho = compression_ratio
Kt0_tilde = Kt0_suc/(rho_0_suc*c_0*c_0)
A_tilde   = ein_cv1_suc*theta_0_suc/(c_0*c_0)
theta_E_tilde = ein_cv2_suc/theta_0_suc
rho_0_tilde = rho_0_suc/rho_0_suc
phi = math.exp(gamma_suc*(1-1/tilde_rho))
int_energy = 0.5*Kt0_tilde*pow(math.log(tilde_rho),2)*(1+(Kt0_prime_suc-2)*math.log(tilde_rho)/3)+A_tilde*(phi*theta_E_tilde*math.exp(phi*theta_E_tilde)/(math.exp(phi*theta_E_tilde)-1)-math.log(math.exp(phi*theta_E_tilde)-1))
int_energy0 = A_tilde*(theta_E_tilde*math.exp(theta_E_tilde)/(math.exp(theta_E_tilde)-1)-math.log(math.exp(theta_E_tilde)-1))
p_theta0 = Kt0_tilde*tilde_rho*math.log(tilde_rho)*(1+0.5*(Kt0_prime_suc-2)*math.log(tilde_rho))
ps =(-tilde_P_0*(1-1/tilde_rho+2/gamma_suc)+(2/gamma_suc)*(-p_theta0 + gamma_suc*(int_energy-int_energy0)))/(1-1/tilde_rho-2/gamma_suc)

vel = math.sqrt((ps-tilde_P_0)*(1-1/tilde_rho))


# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 :  -leng/2.,
    'x_domain%end'                 :  leng/2,
    'y_domain%beg'                 :  -leng/2.,
    'y_domain%end'                 :  leng/2.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : 0,
    'dt'                           : dt,
    't_step_start'                 : 0,
    't_step_stop'                  : Nt,
    't_step_save'                  : int(Nt/80000.),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 3,             #change this to 3 for shocked state
    'model_eqns'                   : 5,
    'alt_soundspeed'               : 'F',
    'hypoplasticity'               : 'T',
    'num_fluids'                   : 2,
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================

    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : 10*leng,
    'patch_icpp(1)%length_y'       : leng,
    'patch_icpp(1)%vel(1)'         : 1.E-06,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%pres'           : tilde_P_0,
    'patch_icpp(1)%alpha_rho(1)'   : (1.E+00-(1.E-06))*(1580.5/1.580488803979682E3),
    'patch_icpp(1)%alpha_rho(2)'   : (1.E-06)*(1.168/1.580488803979682E3),
    'patch_icpp(1)%alpha(1)'       : (1.E+00-1.E-06),
    'patch_icpp(1)%alpha(2)'       : 1.E-06,
    'patch_icpp(1)%tau_e(1)'       : 1.E-17,
    'patch_icpp(1)%tau_e(2)'       : 1.E-17,
    'patch_icpp(1)%tau_e(3)'       : 1.E-17,
    'patch_icpp(1)%tau_e(4)'       : 1.E-17,
    # ==========================================================================

    # Patch 2: Shocked state ===================================================
    'patch_icpp(2)%geometry'       : 3,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : -3*leng/8.,
    'patch_icpp(2)%y_centroid'     : 0.,
    'patch_icpp(2)%length_x'       : leng/4.,
    'patch_icpp(2)%length_y'       : leng,
    'patch_icpp(2)%vel(1)'         : vel,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%pres'           : ps,
    'patch_icpp(2)%alpha_rho(1)'   : (1.E0-(1.E-06))*tilde_rho,
    'patch_icpp(2)%alpha_rho(2)'   : (1.E-06)*(1.168/1.580488803979682E3),
    'patch_icpp(2)%alpha(1)'       : 1.E+00-(1.E-06),
    'patch_icpp(2)%alpha(2)'       : 1.E-06,
    # ==========================================================================

    # Patch 3: Bubble ==========================================================
    'patch_icpp(3)%geometry'       : 2,
    'patch_icpp(3)%x_centroid'     : 0.E+00,
    'patch_icpp(3)%y_centroid'     : 0.E+00,
    'patch_icpp(3)%radius'         : leng/5,
    'patch_icpp(3)%alter_patch(1)' : 'T',
    'patch_icpp(3)%vel(1)'         : 1.E-6,
    'patch_icpp(3)%vel(2)'         : 0.E+00,
    'patch_icpp(3)%pres'           : tilde_P_0,
    'patch_icpp(3)%alpha_rho(1)'   : 1.E-06*(1580.5/1.580488803979682E3),
    'patch_icpp(3)%alpha_rho(2)'   : (1.E0-(1.E-6))*(1.168/1.580488803979682E3),
    'patch_icpp(3)%alpha(1)'       : 1.E-06,
    'patch_icpp(3)%alpha(2)'       : 1.E+00-(1.E-06),
    'patch_icpp(3)%tau_e(1)'       : 1.E-17,
    'patch_icpp(3)%tau_e(2)'       : 1.E-17,
    'patch_icpp(3)%tau_e(3)'       : 1.E-17,
    'patch_icpp(3)%tau_e(4)'       : 1.E-17,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.09E0,                           # 1.E+00/(1.4E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : Kt0_suc/(rho_0_suc*c_0*c_0),      # isothermal bulk modulus
    'fluid_pp(2)%gamma'            : 0.4E0,                            # 1.E+00/(1.6666E+00-1.E+00),
    'fluid_pp(2)%pi_inf'           : Kt0_air/(rho_0_suc*c_0*c_0),                              # 0.0
    'fluid_pp(1)%qv'               : 3.75E0,                           # K'_theta0 for sucrose
    'fluid_pp(2)%qv'               : 2.0E0,                            # K'_theta0 for air
    'fluid_pp(1)%G'                : G_suc/(rho_0_suc*c_0*c_0),        # Shear modulus
    'fluid_pp(2)%G'                : 1.E-5,                            # Shear modulus of air taken to be a very small value
    'fluid_pp(1)%ein_cv(1)'        : A_tilde,                          # Can be replaced with fluid_pp(:)%cv at some point
    'fluid_pp(2)%ein_cv(1)'        : 0.026937087111210E0,              #
    'fluid_pp(1)%ein_cv(2)'        : theta_E_tilde,                    # Can be replaced with a scalar theta_E at some point
    'fluid_pp(2)%ein_cv(2)'        : 100E0/298E0,                      #0.335E0,
    'fluid_pp(1)%mg_a'             : 0.E0,                             #a_mg
    'fluid_pp(1)%mg_b'             : 1.E0,                             #b_mg
    'fluid_pp(2)%mg_a'             : 0.4E0,                            #a_mg
    'fluid_pp(2)%mg_b'             : 0.E0,                             #b_mg
    'fluid_pp(1)%rho0'             : 1.E0,                             #Non-dimensional initial density in Birch-Murnaghan cold curve
    'fluid_pp(2)%rho0'             : rho_0_air/(rho_0_suc),
    'fluid_pp(1)%jcook(1)'         : 0.0334,                           # A, Static yield strength
    'fluid_pp(1)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
    'fluid_pp(1)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
    'fluid_pp(1)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
    'fluid_pp(1)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
    'fluid_pp(1)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient pressure
    'fluid_pp(1)%jcook(7)'         : 3.24928515727,                            # Limiting strain-rate
    'fluid_pp(1)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
    'fluid_pp(1)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
    'fluid_pp(1)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limit
    'fluid_pp(1)%jcook(11)'        : 298/theta_0_suc,                              # Reference temperature
    'fluid_pp(2)%jcook(1)'         : 0.0334,                           # A, Static yield strength
    'fluid_pp(2)%jcook(2)'         : 0.0334,                           # B, Strain-Hardening coefficient
    'fluid_pp(2)%jcook(3)'         : 0.1,                              # n, Strain-Hardening exponent
    'fluid_pp(2)%jcook(4)'         : 0.01,                             # C, Strain-rate hardening coefficient
    'fluid_pp(2)%jcook(5)'         : 0.45,                             # m, Thermal softening exponent
    'fluid_pp(2)%jcook(6)'         : 1.5403,                           # theta_m, Melt temperature at ambient pressure
    'fluid_pp(2)%jcook(7)'         : 3.24928515727,                            # Limiting strain-rate
    'fluid_pp(2)%jcook(8)'         : 0.02,                             # Parameter in Simon-Glatzel melt relation
    'fluid_pp(2)%jcook(9)'         : 3.25,                             # exponent in Simon-Glatzel melt relation
    'fluid_pp(2)%jcook(10)'        : 3.2493E-7,                        # non-dimensional strain-rate limitI
    'fluid_pp(2)%jcook(11)'        : 298/theta_0_suc,                  # non-dimensionalized Reference temperature
}))

# ==============================================================================
