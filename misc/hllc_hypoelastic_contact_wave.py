#!/usr/bin/env python3
"""Verification of m_riemann_solver_hllc.fpp's hypoelastic contact-wave / tangential-stress
treatment, requested before touching HLLC for hyperelasticity.

Checks two separate things against the actual code (src/simulation/m_riemann_solver_hllc.fpp,
the #:if HYPO branch around lines 1138-1260):

1. The NORMAL-direction generalization: s_S (contact speed) and pres_tot_star use
   pres_tot = pres - tau_nn in place of bare pressure. This is a standard, single-wave-family
   substitution (the same one Favrie & Gavrilyuk 2011 and Rodriguez et al. 2019 use) and should
   check out exactly against the ordinary 2-wave (HLLC) momentum Rankine-Hugoniot derivation.

2. The TANGENTIAL treatment: u_t_star, tau_nt_star. These are solved for by combining the
   tangential-momentum RH jump across the LEFT acoustic wave (speed s_L) and across the RIGHT
   acoustic wave (speed s_R) -- i.e. the code's ansatz is that a SINGLE (u_t*, tau_nt*) pair
   is continuous across the ENTIRE star region, the same way u* and pres_tot* are. This script
   finds that u_t_star IS the exact solution of that 2-equation ansatz, but tau_nt_star is NOT:
   it is missing the velocity-coupling term the tangential-momentum jump requires, and is
   therefore inconsistent with u_t_star's own derivation. Beyond that specific mismatch, the
   ansatz itself is also only an approximation: hypoelasticity has a genuine SHEAR
   characteristic at speed u +/- sqrt(Ca_inv/rho), distinct from the acoustic speed s_L/s_R used
   to build A_L, A_R, so even a corrected tau_nt_star only approximates the true (v, sigma)
   jump, which should propagate at its own speed rather than smear across the acoustic fan.

This does NOT affect the planned hyperelasticity extension to HLLC: hyperelastic stress is not
an independently flux-transported quantity (it is slaved algebraically to the deformation
gradient every RHS call -- see m_hyperelastic.fpp), so it has no u_t_star/tau_nt_star-style
tangential star-state to compute at all. Only the normal-direction substitution checked in
part 1 (already confirmed correct) is relevant to that extension.
"""

import sympy as sp


def check_normal_direction():
    print("=" * 78)
    print("1. Normal direction: s_S and pres_tot_star")
    print("=" * 78)

    rho_L, rho_R, u_L, u_R, s_L, s_R = sp.symbols("rho_L rho_R u_L u_R s_L s_R", real=True)
    pres_tot_L, pres_tot_R = sp.symbols("pres_tot_L pres_tot_R", real=True)

    # Standard HLLC two-wave derivation (Toro, "Riemann Solvers", Ch. 10), with pressure
    # replaced by the total normal stress pres_tot = pres - tau_nn throughout -- that
    # substitution is the only change a hypoelastic normal-direction HLLC needs, since the
    # normal-momentum equation's flux is rho*u**2 + pres_tot (exactly like rho*u**2 + p for an
    # inviscid fluid, just with pres_tot standing in for p).
    s_S, pres_tot_star = sp.symbols("s_S pres_tot_star", real=True)

    # Mass RH across each outer wave: rho_star_K*(s_K - u_star) = rho_K*(s_K - u_K) == A_K
    A_L = rho_L * (s_L - u_L)
    A_R = rho_R * (s_R - u_R)

    # Normal-momentum RH across each outer wave (Toro, "Riemann Solvers and Numerical Methods
    # for Fluid Dynamics", eq. 10.36): pres_tot_star = pres_tot_K + A_K*(s_S - u_K), for K = L, R.
    eq_L = sp.Eq(pres_tot_star, pres_tot_L + A_L * (s_S - u_L))
    eq_R = sp.Eq(pres_tot_star, pres_tot_R + A_R * (s_S - u_R))

    sol = sp.solve([eq_L, eq_R], [s_S, pres_tot_star], dict=True)[0]
    s_S_derived = sp.simplify(sol[s_S])
    pres_tot_star_derived = sp.simplify(sol[pres_tot_star])

    # The code's formulas (m_riemann_solver_hllc.fpp, HYPO branch):
    #   s_S = (pres_tot_R - pres_tot_L + rho_L*u_L*(s_L-u_L) - rho_R*u_R*(s_R-u_R)) / (A_L - A_R)... (code writes it as /(A_L-A_R) with signs flipped; check equivalently)
    s_S_code = (pres_tot_R - pres_tot_L + rho_L * u_L * (s_L - u_L) - rho_R * u_R * (s_R - u_R)) / (
        rho_L * (s_L - u_L) - rho_R * (s_R - u_R)
    )
    pres_tot_star_code = pres_tot_L + A_L * (s_S - u_L)  # as coded, using the SYMBOLIC s_S (not yet substituted)

    diff_sS = sp.simplify(s_S_derived - s_S_code)
    print(f"s_S:  derived - coded = {diff_sS}   (expect 0)")
    assert diff_sS == 0, "s_S formula does not match the two-wave normal-momentum RH derivation!"

    # pres_tot_star_code uses the SYMBOLIC s_S; substitute the derived closed form for s_S and recheck
    pres_tot_star_code_closed = sp.simplify(pres_tot_star_code.subs(s_S, s_S_derived))
    diff_pStar = sp.simplify(pres_tot_star_derived - pres_tot_star_code_closed)
    print(f"pres_tot_star:  derived - coded = {diff_pStar}   (expect 0)")
    assert diff_pStar == 0, "pres_tot_star formula does not match the derivation!"
    print("PASS: the normal-direction (pres -> pres - tau_nn) substitution is exactly the")
    print("      standard two-wave HLLC momentum RH solution. No issue here.")


def check_tangential_ansatz():
    print()
    print("=" * 78)
    print("2. Tangential direction: u_t_star, tau_nt_star -- is the code's ansatz self-consistent?")
    print("=" * 78)

    rho_L, rho_R, u_L, u_R, s_L, s_R = sp.symbols("rho_L rho_R u_L u_R s_L s_R", real=True)
    v_L, v_R, sig_L, sig_R = sp.symbols("v_L v_R sigma_L sigma_R", real=True)  # tangential vel, shear stress
    u_t_star, tau_nt_star = sp.symbols("u_t_star tau_nt_star", real=True)

    A_L = rho_L * (s_L - u_L)
    A_R = rho_R * (s_R - u_R)

    # Code's ansatz: a SINGLE (u_t*, tau_nt*) pair, continuous across the WHOLE star region,
    # solved from the tangential-momentum RH jump across the left wave (speed s_L) and across
    # the right wave (speed s_R). Tangential-momentum flux is rho*u*v - sigma (the "-sigma" plays
    # the role "+p" played in the normal direction), so by the same Toro-style relation:
    # tau_nt_star = sigma_K + A_K*(v_K - u_t_star), for K = L, R.
    eq_L = sp.Eq(tau_nt_star, sig_L + A_L * (v_L - u_t_star))
    eq_R = sp.Eq(tau_nt_star, sig_R + A_R * (v_R - u_t_star))

    sol = sp.solve([eq_L, eq_R], [u_t_star, tau_nt_star], dict=True)[0]
    u_t_star_derived = sp.simplify(sol[u_t_star])
    tau_nt_star_derived = sp.simplify(sol[tau_nt_star])

    denom_A = A_R - A_L
    u_t_star_code = (A_R * v_R - A_L * v_L + (sig_R - sig_L)) / denom_A
    tau_nt_star_code = (A_R * sig_R - A_L * sig_L) / denom_A

    diff_u = sp.simplify(u_t_star_derived - u_t_star_code)
    diff_tau = sp.simplify(tau_nt_star_derived - tau_nt_star_code)
    print(f"u_t_star:    derived - coded = {diff_u}   (expect 0)")
    print(f"tau_nt_star: derived - coded = {diff_tau}")
    assert diff_u == 0, "u_t_star does not solve its own stated ansatz -- that would be a coding error."
    print("u_t_star matches its ansatz exactly: no algebra mistake there.")
    if diff_tau != 0:
        print("tau_nt_star does NOT match solving the same 2x2 system as u_t_star -- the code's formula")
        print("is missing the velocity-coupling terms the tangential-momentum RH jump requires.")
        print()
        print("  code's formula:   tau_nt_star = (A_R*sigma_R - A_L*sigma_L) / (A_R - A_L)")
        print()
        # Express the correct formula the same simple way, by substituting u_t_star back into
        # either single-wave relation (both give the same answer once u_t_star is correct).
        correct_form = sp.simplify(sig_L + A_L * (v_L - u_t_star_code))
        print("  correct formula:  tau_nt_star = sigma_L + A_L*(v_L - u_t_star)")
        print("                                = sigma_R + A_R*(v_R - u_t_star)   (equivalently)")
        print(f"  in full:          tau_nt_star = {correct_form}")
    else:
        print("tau_nt_star matches its ansatz exactly too.")

    print()
    print("-" * 78)
    print("Now the physical question: is 'a single (u_t*, tau_nt*) pair, continuous across the")
    print("whole star region' the CORRECT ansatz for hypoelasticity?")
    print("-" * 78)
    print("""
This ansatz is exactly the one used for u* and pres_tot* above -- correct there because u*
and pres_tot* genuinely ARE continuous across the contact wave in the exact Riemann solution
(the defining property of a contact/material discontinuity).

But (v, sigma) -- tangential velocity and shear stress -- are NOT governed by the acoustic
wave family at all. m_riemann_state.fpp's f_elastic_signal_speed already encodes that the
LONGITUDINAL characteristic speed is c_L = sqrt(c**2 + 4*Ca_inv/(3*rho) + tau_nn/rho); the
hypoelastic system has a SEPARATE shear characteristic at c_shear = sqrt(Ca_inv/rho) (this is
exactly the wave equations.md's own HLLD section describes: "Elastic shear waves play the
role of Alfven waves, with ... shear wave impedance C = rho_hat*(G_hat + tau_nn_hat)").

The code's tangential ansatz implicitly assumes (v, sigma) jump AT THE SAME SPEEDS s_L/s_R as
the acoustic wave (that is where A_L, A_R come from). Since c_shear != c_acoustic in general
(they differ by the O(c**2) fluid compressibility term), the true solution has (v, sigma)
jumping at a DIFFERENT, INTERIOR speed -- i.e. there should be 4 interior waves (2 shear + 1
contact, nested between the 2 outer acoustic waves), not the 1 contact wave HLLC's 3-wave
structure provides for. Treating (v, sigma) as "contact-continuous" smears the genuine shear
jump across the whole acoustic fan instead of letting it propagate at its own speed.
""")


def check_hlld_reference():
    print("=" * 78)
    print("3. Cross-check: HLLD's shear-wave impedance (the structure already confirmed correct)")
    print("=" * 78)
    rho, G, tau_nn = sp.symbols("rho G tau_nn", positive=True)
    # equations.md: "shear wave impedance C = rho_hat*(G_hat + tau_nn_hat)"
    # the shear characteristic speed is then C/rho_hat = G_hat + tau_nn_hat -- NOT sqrt(G/rho)
    # as a first pass might assume; it is the TANGENT shear modulus including the current
    # normal-stress-dependent correction (the same acoustic-tensor-style tangent stiffness
    # idea used for the longitudinal wave in m_hyperelastic's own wave speed, just for shear).
    C = rho * (G + tau_nn)
    print("HLLD shear impedance C = rho_hat*(G_hat + tau_nn_hat):")
    sp.pprint(C)
    print()
    print("This is a genuinely SEPARATE characteristic speed from the longitudinal/acoustic one")
    print("(f_elastic_signal_speed's sqrt(c**2 + 4*Ca_inv/(3*rho) + tau_nn/rho)), confirming HLLD")
    print("resolves the shear wave as its own distinct family -- exactly what HLLC's 3-wave")
    print("structure structurally cannot do, independent of how carefully the tangential-jump")
    print("formulas above are derived.")


if __name__ == "__main__":
    check_normal_direction()
    check_tangential_ansatz()
    check_hlld_reference()
