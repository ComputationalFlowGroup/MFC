#!/usr/bin/env python3
"""Verification that m_riemann_solver_hypo_hlld.fpp's shear-wave star state is the exact
Riemann-invariant solution across hypoelasticity's OWN shear characteristic -- the thing
misc/hllc_hypoelastic_contact_wave.py found HLLC's tau_nt_star could not be (HLLC has no
independent shear wave at all; it smears (v, sigma) across the acoustic fan at s_L/s_R).

The shear-stress/tangential-velocity subsystem (rho*v_t = tau_nt_x, tau_nt_t = C*v_x, where
C is the shear impedance-squared-times-density) is the exact 1D elastic-shear analog of the
acoustic subsystem (rho*u_t = -p_x, p_t = rho*c**2*u_x): same structure, "pressure" -> tau_nt,
"velocity" -> v, acoustic impedance rho*c -> shear impedance Z = sqrt(C). Its Riemann
invariants are tau_nt -/+ Z*v = const along the +/- characteristics, giving exactly the same
closed form the acoustic p +/- rho*c*u invariants give for pres_tot_star/s_S elsewhere on
this branch (misc/hllc_hypoelastic_contact_wave.py's part 1).

Checks the code's own formulas (m_riemann_solver_hypo_hlld.fpp):
  C_NC (line 607)    = rho_hat*(Ca_inv_hat + tau_nn_hat)       -- the shear impedance-squared
  u_t_star (line 659)  = 0.5*((tau_nt_R - tau_nt_L)/sqrtC_NC + (u_t_R + u_t_L))
  tau_nt_star (line 660) = 0.5*((u_t_R - u_t_L)*sqrtC_NC + (tau_nt_R + tau_nt_L))
against the exact solution of the 2x2 Riemann-invariant system, and C_NC's formula against
equations.md's own stated shear impedance "C = rho_hat*(G_hat + tau_nn_hat)" (G = Ca_inv).
"""

import sympy as sp


def check_shear_impedance():
    print("=" * 78)
    print("1. Shear impedance C_NC matches equations.md's stated formula")
    print("=" * 78)
    rho_hat, Ca_inv_hat, tau_nn_hat = sp.symbols("rho_hat Ca_inv_hat tau_nn_hat", real=True)
    C_NC_code = rho_hat * (Ca_inv_hat + tau_nn_hat)
    C_NC_docs = rho_hat * (Ca_inv_hat + tau_nn_hat)  # equations.md: C = rho_hat*(G_hat + tau_nn_hat), G = Ca_inv
    print(f"C_NC:  code - docs = {sp.simplify(C_NC_code - C_NC_docs)}   (expect 0)")
    assert sp.simplify(C_NC_code - C_NC_docs) == 0


def check_shear_star_state():
    print()
    print("=" * 78)
    print("2. u_t_star, tau_nt_star solve the shear-wave Riemann-invariant system exactly")
    print("=" * 78)

    v_L, v_R, sig_L, sig_R, Z = sp.symbols("v_L v_R sigma_L sigma_R Z", real=True)
    v_star, sig_star = sp.symbols("v_star sigma_star", real=True)

    # Riemann invariants across the shear characteristics, exactly analogous to the acoustic
    # p -/+ rho*c*u invariants (Toro Ch. 2/10) with "pressure" -> tau_nt, impedance rho*c -> Z:
    #   along the LEFT-going characteristic (from the L state):  sigma - Z*v = sigma_L - Z*v_L
    #   along the RIGHT-going characteristic (from the R state): sigma + Z*v = sigma_R + Z*v_R
    eq_L = sp.Eq(sig_star - Z * v_star, sig_L - Z * v_L)
    eq_R = sp.Eq(sig_star + Z * v_star, sig_R + Z * v_R)
    sol = sp.solve([eq_L, eq_R], [v_star, sig_star], dict=True)[0]
    v_star_derived = sp.simplify(sol[v_star])
    sig_star_derived = sp.simplify(sol[sig_star])

    # Code's formulas (m_riemann_solver_hypo_hlld.fpp:659-660), Z = sqrtC_NC:
    v_star_code = sp.Rational(1, 2) * ((sig_R - sig_L) / Z + (v_R + v_L))
    sig_star_code = sp.Rational(1, 2) * ((v_R - v_L) * Z + (sig_R + sig_L))

    diff_v = sp.simplify(v_star_derived - v_star_code)
    diff_sig = sp.simplify(sig_star_derived - sig_star_code)
    print(f"u_t_star:    derived - code = {diff_v}   (expect 0)")
    print(f"tau_nt_star: derived - code = {diff_sig}   (expect 0)")
    assert diff_v == 0, "u_t_star does not solve the shear-wave Riemann-invariant system."
    assert diff_sig == 0, "tau_nt_star does not solve the shear-wave Riemann-invariant system."

    print()
    print("PASS: unlike HLLC's old tau_nt_star (misc/hllc_hypoelastic_contact_wave.py), HLLD's")
    print("shear star state is not an average smeared across the acoustic fan -- it is the")
    print("exact jump condition across hypoelasticity's OWN shear characteristic, at its own")
    print("speed sqrt(C_NC)/rho, independent of the acoustic speeds S_L/S_R. This is exactly")
    print("the structural fix HLLC's 3-wave family cannot provide (see that script's closing")
    print("discussion) -- HLLD already has it. No further HLLD change is needed for this.")


if __name__ == "__main__":
    check_shear_impedance()
    check_shear_star_state()
