#!/usr/bin/env python3

import math
import json

pi = 3.141592653589
# material parameters


#patterson 2018
#material1 :: gas
gammag = 1.4    #unitless
Bg = 0.E0       #pascals
rhog = 1.18     #kg/m^3
c_g = 347.2     #m/s
G_g = 0         #pa

#gas Cv calculation
Ru = 8.3144598
Ra = Ru/(28.955E-3)
Cp_g = Ra*gammag/(gammag-1)
Cv_g = Cp_g/gammag

#material2 :: water
gammal = 5.5
Bl = 492.E+06
rhol = 996.0
c_l = 1540#1648.7
G_l = 1.0E4
Cv_l = 1816

#primitive vartiables
patmos = 101325. #pa

#problem specific variable
lambda_wave = 200.E-6

#define pulse
P_amp = 5E+5

#non-dim
#define characteristic density, length, time, stress material                   #make it liquid
rho_char = 1#rhol
length_char = 1#lambda_wave
c_char = 1#c_l                                                                    #should be liquid
time_char = 1#length_char/c_char
stress_char = 1#rho_char*c_char*c_char/gammal

#non-dim the properties
rhog_n  = rhog/rho_char    
c_g_n = c_g/c_char    
rhol_n = rhol/rho_char
c_l_n = c_l/c_char
Bg_n = Bg/stress_char
Bl_n = Bl/stress_char
G_g_n = G_g/stress_char
G_l_n = G_l/stress_char
patmos_n = patmos/stress_char
P_amp_n = P_amp/stress_char

#geometry
amp = 0.5
interface_amp = amp*lambda_wave
Namp = 32

dlengx = 8.*lambda_wave
dlengy = 1.*lambda_wave
dlengz = 1.*lambda_wave
Ny = int(Namp/amp)
Nx = dlengx*Ny/dlengy
Nz = dlengz*Ny/dlengy
dx = dlengx/Nx
dy = dlengy/Ny
dz = dlengz/Nz

eps = 1E-8

alphal_back = 1.0 - eps
alphag_back = 1.0 - alphal_back

alphag_lung = 1.0 - eps
alphal_lung = 1.0 - alphag_lung

# time stepping requirements
time_end = 1E-5
frames = 100
time_save = time_end/frames
cfl = 0.25

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    'sim_data'                     : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -dlengx/2.,
    'x_domain%end'                 : dlengx/2.,
    'y_domain%beg'                 : -dlengy/2.,
    'y_domain%end'                 : dlengy/2.,
    'z_domain%beg'                 : -dlengz/2.,
    'z_domain%end'                 : dlengz/2.,
    'm'                            : int(Nx),
    'n'                            : int(Ny),
    'p'                            : int(Nz),
    'stretch_x'                    : 'T',
    'a_x'                          : 1.0E+00,
    'x_a'                          : -2.5*lambda_wave,
    'x_b'                          : 2.5*lambda_wave,
    'loops_x'                      : 3,
    'cfl_adap_dt'                  : 'T',
    'cfl_target'                   : cfl,
    'n_start'                      : 0,#36,
    't_stop'                       : time_end,
    't_save'                       : time_save,
    # ==========================================================================
    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 3,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -6,
    'bc_x%end'                     : -6,
    'bc_y%beg'                     : -1,            
    'bc_y%end'                     : -1,
    'bc_z%beg'                     : -1,
    'bc_z%end'                     : -1,
    # ==========================================================================
    # Turning on Hypoelasticity ================================================
    'hypoelasticity'               : 'T',     
    #'hyperelasticity'              : 'T',
    # ==========================================================================
    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    #===========================================================================
    # Patch 1: Background ======================================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.,
    'patch_icpp(1)%y_centroid'     : 0.,
    'patch_icpp(1)%z_centroid'     : 0.,
    'patch_icpp(1)%length_x'       : dlengx,
    'patch_icpp(1)%length_y'       : dlengy,
    'patch_icpp(1)%length_z'       : dlengz,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,            
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : patmos_n,
    'patch_icpp(1)%alpha_rho(1)'   : rhol_n*alphal_back,
    'patch_icpp(1)%alpha_rho(2)'   : rhog_n*alphag_back,
    'patch_icpp(1)%alpha(1)'       : alphal_back,
    'patch_icpp(1)%alpha(2)'       : alphag_back,
    # ==========================================================================
    # Patch 2: Lung ============================================================
    'patch_icpp(2)%geometry'       : 13,
    'patch_icpp(2)%hcid'           : 302,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%x_centroid'     : 0.,
    'patch_icpp(2)%y_centroid'     : 0.,
    'patch_icpp(2)%z_centroid'     : 0.,
    'patch_icpp(2)%length_x'       : 5.,#dlengx,
    'patch_icpp(2)%length_y'       : dlengy,  
    'patch_icpp(2)%length_z'       : dlengz,                 
    'patch_icpp(2)%a2'             : interface_amp,
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%vel(3)'         : 0.E+00,
    'patch_icpp(2)%pres'           : patmos_n,
    'patch_icpp(2)%alpha_rho(1)'   : rhol_n*alphal_lung,
    'patch_icpp(2)%alpha_rho(2)'   : rhog_n*alphag_lung,
    'patch_icpp(2)%alpha(1)'       : alphal_lung,
    'patch_icpp(2)%alpha(2)'       : alphag_lung,
    # ==========================================================================
    # Patch 3: Shocked Lung ====================================================
    'patch_icpp(3)%geometry'       : 9,
    'patch_icpp(3)%alter_patch(2)' : 'T',
    'patch_icpp(3)%x_centroid'     : 7.65E-3,
    'patch_icpp(3)%y_centroid'     : 0.,
    'patch_icpp(3)%z_centroid'     : 0.,
    'patch_icpp(3)%length_x'       : 15.2E-3,#dlengx,
    'patch_icpp(3)%length_y'       : dlengy,  
    'patch_icpp(3)%length_z'       : dlengz,                 
    'patch_icpp(3)%vel(1)'         : 0.E+00,
    'patch_icpp(3)%vel(2)'         : 0.E+00,
    'patch_icpp(3)%vel(3)'         : 0.E+00,
    'patch_icpp(3)%pres'           : P_amp, 
    'patch_icpp(3)%alpha_rho(1)'   : rhol_n*alphal_back,
    'patch_icpp(3)%alpha_rho(2)'   : rhog_n*alphag_back,            
    'patch_icpp(3)%alpha(1)'       : alphal_back,
    'patch_icpp(3)%alpha(2)'       : alphag_back,
    # ==========================================================================
    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gammal-1.E+00),
    'fluid_pp(1)%pi_inf'           : gammal*Bl_n/(gammal-1.E+00),
    'fluid_pp(1)%G'                : G_l_n,
    'fluid_pp(2)%gamma'            : 1.E+00/(gammag-1.E+00),
    'fluid_pp(2)%pi_inf'           : gammag*Bg_n/(gammag-1.E+00),   
    'fluid_pp(2)%G'                : G_g_n,
 #==============================================================================
}))
# ==============================================================================
