#!/usr/bin/env python3
"""Acoustic-tensor derivation for the neo-Hookean hyperelastic term in MFC's Riemann solver.

Standalone script: depends only on SymPy, not on MFC itself (no mfc.sh, no venv, no
other part of this repository). Install the one dependency and run it directly:

    pip install sympy
    python3 hyperelastic_acoustic_tensor.py

What this derives, step by step, and why
-----------------------------------------
MFC's hyperelastic module (src/common/m_hyperelastic.fpp) computes the deviatoric
Cauchy stress of a neo-Hookean solid from the left Cauchy-Green tensor B = F F^T:

    tau = (G/J) * dev(B) = (G/J) * (B - tr(B)/3 * I)

where G is the shear modulus, J = det(F). The direct (HLL-type) wave speed used in the
Riemann solver adds a shear correction to the ordinary fluid sound speed c:

    s = u.n +/- sqrt(c**2 + correction)

For a HYPOELASTIC material (a rate law, where tau itself is one of the transported
variables) the correction is read directly off the current stress: (4G/3 + tau_nn)/rho.
That shortcut is only valid because the rate law makes tau the right quantity to
linearize around.

For a finite-strain law like neo-Hookean, tau is NOT the right quantity -- the TANGENT
STIFFNESS of the strain-energy function is. That tangent stiffness, contracted with the
direction of wave propagation, is exactly what continuum mechanics calls the ACOUSTIC
TENSOR. This script derives it from first principles:

  1. Start from the kinematic identity for how B is carried by the flow: Bdot = L B + B L^T,
     where L = grad(v) is the velocity gradient (this holds for any motion, independent of
     the constitutive law).
  2. Differentiate tau(B, J) using that identity, to get tau_dot(B, J, L).
  3. Form the Truesdell rate tau_truesdell = tau_dot - L.tau - tau.L^T + tr(L)*tau. This is
     the objective (frame-indifferent) rate appropriate for a Cauchy stress that is a
     function of B alone. The script VERIFIES objectivity symbolically: tau_truesdell must
     depend only on the symmetric part of L (d = sym(L)), not on the antisymmetric spin
     part -- if it didn't, the result would depend on the observer's rotation, which is
     physically meaningless.
  4. Substitute a plane-wave velocity perturbation, L = outer(u, n) (polarization u,
     propagation direction n), and read off the matrix A(n) such that
     tau_truesdell . n = A(n) . u. A(n)/rho is the acoustic tensor: its eigenvalues are
     rho*c**2 for the elastic correction to the wave speed in direction n.
  5. Check the result against the known small-strain limit (B = I): the longitudinal mode
     must reduce to the textbook 4*G/(3*rho), and a shear mode must reduce to G/rho.
  6. Finally, reproduce the closed form actually used in the Riemann solver (the
     longitudinal mode, u = n):

         correction = (2*G)/(3*J*rho) * (tr(B) - B_nn)

     and confirm it collapses to 4*G/(3*rho) at B = I, J = 1 -- the same consistency check
     done by hand in the derivation this script mirrors.

This also sets up the one piece of the derivation that is NOT done here and is tedious to
do by hand: the same exercise for the Mooney-Rivlin law (also implemented in
m_hyperelastic.fpp, s_Mooney_Rivlin_cauchy_solver), whose stress depends on B and B**2.
See the comment at the bottom of this file for how to extend the script to that case.

Background reading, in order of how closely each matches what this script does
--------------------------------------------------------------------------------
- Ogden, R.W., "Non-Linear Elastic Deformations" (Dover, 1997), Chapter 6: the acoustic
  tensor, strong ellipticity, and plane-wave propagation in a deformed elastic solid. The
  classic reference for exactly this calculation.
- Holzapfel, G.A., "Nonlinear Solid Mechanics" (Wiley, 2000), Chapter 6: objective
  (Truesdell/Lie) stress rates and the spatial elasticity tensor for isotropic hyperelastic
  materials -- the machinery this script automates.
- Truesdell, C. and Noll, W., "The Non-Linear Field Theories of Mechanics": the original,
  rigorous treatment of objective stress rates (the Truesdell rate is named for this book).
- Miller, G.H. and Colella, P., "A conservative three-dimensional Eulerian method for
  coupled solid-fluid shock capturing," Journal of Computational Physics 183 (2001): derives the full
  eigenstructure of a hyperelastic system's flux Jacobian (with the deformation gradient
  carried as a transported variable) and shows it reduces to this acoustic-tensor
  eigenproblem -- the connection from continuum theory to a Godunov-type scheme.
- Favrie, N. and Gavrilyuk, S.L., "Diffuse interface model for compressible fluid -
  compressible elastic-plastic solid interaction," Journal of Computational Physics 230 (2011): the closest
  match to MFC's own HLLC-type use case for a hyperelastic solid.
"""

import sympy as sp


def main():
    G, J, rho = sp.symbols("G J rho", positive=True)

    # B: the left Cauchy-Green tensor, symmetric by construction.
    B = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"B{min(i, j) + 1}{max(i, j) + 1}"))

    # L: a fully general velocity gradient (not assumed symmetric). Kept general here so
    # the objectivity check in step 2 below is meaningful -- it is only substituted down
    # to a plane-wave form afterward.
    L = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"L{i + 1}{j + 1}"))

    I3 = sp.eye(3)

    print("=" * 78)
    print("Step 1: the neo-Hookean Cauchy stress and its kinematic rate of change")
    print("=" * 78)

    trB = sp.trace(B)
    tau = (G / J) * (B - trB / 3 * I3)
    print("tau = (G/J) * dev(B)")

    # Kinematic identities, valid for ANY motion (not specific to this constitutive law):
    #   Bdot  = L*B + B*L^T
    #   Jdot/J = tr(L)
    Bdot = L * B + B * L.T
    J_rate_over_J = sp.trace(L)

    # tau = (G/J) * (B - tr(B)/3 I); differentiate by the product/quotient rule using the
    # two kinematic identities above.
    tau_dot = -J_rate_over_J * tau + (G / J) * (Bdot - sp.trace(Bdot) / 3 * I3)

    print("\n" + "=" * 78)
    print("Step 2: the Truesdell rate, and an objectivity check")
    print("=" * 78)

    tau_truesdell = tau_dot - L * tau - tau * L.T + sp.trace(L) * tau
    tau_truesdell = sp.simplify(tau_truesdell)

    # Objectivity check: split L into its symmetric (d) and antisymmetric (spin, w) parts,
    # and confirm the Truesdell rate depends only on d. If the coefficient of w's
    # independent components is not identically zero, something upstream is wrong -- the
    # result would depend on the observer's rotation rate, which is not physical.
    d_sym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"d{min(i, j) + 1}{max(i, j) + 1}"))
    w12, w13, w23 = sp.symbols("w12 w13 w23")
    w_antisym = sp.Matrix([[0, w12, w13], [-w12, 0, w23], [-w13, -w23, 0]])

    subs_L = {L[i, j]: (d_sym + w_antisym)[i, j] for i in range(3) for j in range(3)}
    check = sp.expand(tau_truesdell.subs(subs_L))
    no_spin = check.subs({w12: 0, w13: 0, w23: 0})
    spin_dependence = sp.simplify(check - no_spin)
    assert spin_dependence == sp.zeros(3, 3), "Truesdell rate depends on spin -- not objective!"
    print("Verified: tau_truesdell has zero dependence on the spin (antisymmetric) part of L.")
    print("(If this assertion ever fails after editing the stress law above, stop --")
    print(" the result is no longer frame-indifferent.)")

    # Replace L by its symmetric part d for everything from here on; the antisymmetric
    # part has been shown not to matter.
    tau_truesdell_of_d = sp.simplify(tau_truesdell.subs(subs_L).subs({w12: 0, w13: 0, w23: 0}))

    print("\n" + "=" * 78)
    print("Step 3: closed form in terms of B and d, and comparison to the hand derivation")
    print("=" * 78)

    trB_d = sp.trace(B * d_sym)  # B : d, the double contraction
    claimed_form = (2 * G / (3 * J)) * (trB * d_sym - trB_d * I3)
    residual = sp.simplify(tau_truesdell_of_d - claimed_form)
    assert residual == sp.zeros(3, 3), "Closed form does not match the symbolic Truesdell rate!"
    print("Confirmed: tau_truesdell = (2G/3J) * [tr(B)*d - (B:d)*I]")
    print("This is the closed form derived by hand in the conversation this script backs up.")

    print("\n" + "=" * 78)
    print("Step 4: the acoustic tensor -- plane-wave substitution")
    print("=" * 78)

    # Plane wave: velocity perturbation delta_v = u * phi(n.x - c t). Its gradient is
    # L = outer(u, n) (times phi', a common scalar factor that cancels in every ratio
    # below, so it is simply set to 1).
    u1, u2, u3, n1, n2, n3 = sp.symbols("u1 u2 u3 n1 n2 n3")
    u_vec = sp.Matrix([u1, u2, u3])
    n_vec = sp.Matrix([n1, n2, n3])
    d_plane = sp.Rational(1, 2) * (u_vec * n_vec.T + n_vec * u_vec.T)

    tau_plane = tau_truesdell_of_d.subs({d_sym[i, j]: d_plane[i, j] for i in range(3) for j in range(3)})
    tau_plane = sp.simplify(tau_plane)

    # Force per unit n: the divergence of stress along n, i.e. tau_plane . n. This is
    # linear in u, so its Jacobian with respect to u is the matrix A(n) with
    # tau_plane . n = A(n) . u.
    force_along_n = sp.simplify(tau_plane * n_vec)
    A_of_n = force_along_n.jacobian(u_vec)
    A_of_n = sp.simplify(A_of_n)
    print("A(n), such that (Truesdell stress response).n = A(n).u  (acoustic tensor is A(n)/rho):")
    sp.pprint(A_of_n)

    print("\n" + "=" * 78)
    print("Step 5: small-strain consistency check (B = I, J = 1)")
    print("=" * 78)

    subs_undeformed = {B[i, j]: I3[i, j] for i in range(3) for j in range(3)}
    subs_undeformed[J] = 1

    A_undeformed = sp.simplify(A_of_n.subs(subs_undeformed))
    print("A(n) at B = I, J = 1:")
    sp.pprint(A_undeformed)

    # n unit, aligned with x: longitudinal mode u = n should give 4G/3; a transverse
    # mode u perpendicular to n should give G.
    subs_n_x = {n1: 1, n2: 0, n3: 0}
    A_n_x = A_undeformed.subs(subs_n_x)
    longitudinal = sp.simplify((A_n_x * sp.Matrix([1, 0, 0]))[0])
    transverse = sp.simplify((A_n_x * sp.Matrix([0, 1, 0]))[1])
    print(f"\nLongitudinal eigenvalue (u = n = x-hat): {longitudinal}  (expected: 4*G/3)")
    print(f"Transverse eigenvalue   (u perp n):       {transverse}  (expected: G)")
    assert longitudinal == sp.Rational(4, 3) * G
    assert transverse == G
    print("Both match the textbook small-strain (linear elasticity) values.")

    print("\n" + "=" * 78)
    print("Step 6: the longitudinal correction used in the Riemann solver, for general B")
    print("=" * 78)

    # Longitudinal mode: u = n. rho*c**2 = n.A(n).n.
    rho_c2 = sp.simplify((n_vec.T * A_of_n * n_vec)[0])
    print("rho * c_correction**2 = n.A(n).n =")
    sp.pprint(rho_c2)

    B_nn = sp.simplify((n_vec.T * B * n_vec)[0])
    claimed_correction = (2 * G) / (3 * J) * (trB - B_nn)

    # rho_c2 only equals the claimed closed form once n is a UNIT vector (n.n = 1); it was
    # used that way when n.n appeared, back in the hand derivation. Rather than carry that
    # constraint through symbolically, check it numerically at several random unit vectors
    # and random (positive-definite-ish) B, G, J -- simpler and just as convincing.
    import random

    random.seed(0)
    max_err = 0.0
    for _ in range(20):
        v = [random.uniform(-1, 1) for _ in range(3)]
        norm = sum(c * c for c in v) ** 0.5
        n_num = [c / norm for c in v]
        B_num = sp.Matrix(3, 3, lambda i, j: random.uniform(0.5, 2.0) if i == j else random.uniform(-0.3, 0.3))
        B_num = (B_num + B_num.T) / 2  # symmetrize
        subs_num = {n1: n_num[0], n2: n_num[1], n3: n_num[2], G: 1.7, J: 1.3}
        for i in range(3):
            for j in range(i, 3):
                subs_num[B[i, j]] = B_num[i, j]
        lhs = float(rho_c2.subs(subs_num))
        rhs = float(claimed_correction.subs(subs_num))
        max_err = max(max_err, abs(lhs - rhs))
    assert max_err < 1e-10, f"Mismatch between derived and claimed correction: max_err={max_err}"
    print(f"Numerically verified over 20 random (unit n, symmetric B) samples: max error = {max_err:.2e}")
    print("\nConfirmed: rho * c_correction**2 = (2G)/(3J) * (tr(B) - B_nn)")
    print("i.e.       c_correction**2        = (2G)/(3*J*rho) * (tr(B) - n.B.n)")
    print("\nAt B = I, J = 1, n = x-hat: tr(B) = 3, B_nn = 1, giving (2G/3)*2 = 4G/3 --")
    print("exactly the constant already used in m_riemann_solver_hll.fpp.")


def mooney_rivlin():
    """Same derivation as main(), but for the genuine (two-constant) Mooney-Rivlin law.

    s_Mooney_Rivlin_cauchy_solver now implements this: see
    mooney_rivlin_implementation_check() below for the additional identities that
    implementation relies on (B^-1 = grad_xi^T grad_xi, the component-wise stress split,
    and the second invariant formula), verified against this derivation.

    The standard (Rivlin) isotropic stress for W = C1*(I1 - 3) + C2*(I2 - 3) is
    sigma = -p*I + 2*C1*B - 2*C2*B^-1 (up to the pressure term, handled separately by the
    EOS in MFC, exactly as for neo-Hookean). Generalizing the existing code's convention
    (deviatoric projection, scaled by 1/J) gives:

        tau = (1/J) * [2*C1*dev(B) - 2*C2*dev(B^-1)]

    which reduces exactly to the existing tau = (G/J)*dev(B) when C2 = 0 and G = 2*C1.
    """
    C1, C2, J, rho = sp.symbols("C1 C2 J rho", positive=True)
    B = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"B{min(i, j) + 1}{max(i, j) + 1}"))
    L = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"L{i + 1}{j + 1}"))
    I3 = sp.eye(3)

    print("\n" + "#" * 78)
    print("# Mooney-Rivlin extension")
    print("#" * 78)

    print("\n" + "=" * 78)
    print("Step 1: the Mooney-Rivlin Cauchy stress, and B^-1's kinematic rate")
    print("=" * 78)

    trB = sp.trace(B)
    # Binv must be the ACTUAL inverse of B, not an independent symbol standing in for it:
    # the rate identity used below, (B^-1)dot = -B^-1*Bdot*B^-1, only holds when B*Binv=I
    # is a true relation sympy can use to simplify, not an assumption imposed on the rate
    # alone while leaving B and Binv otherwise unrelated. (An earlier version of this
    # function used an independent symbol here; the objectivity check below failed,
    # because of exactly this inconsistency -- left as a cautionary note.)
    Binv = B.inv()
    trBinv = sp.trace(Binv)

    def dev(X, trX):
        return X - trX / 3 * I3

    tau = (1 / J) * (2 * C1 * dev(B, trB) - 2 * C2 * dev(Binv, trBinv))
    print("tau = (1/J) * [2*C1*dev(B) - 2*C2*dev(B^-1)]")

    # Kinematics: Bdot = L*B + B*L^T (same identity as the neo-Hookean case). B^-1's rate
    # follows from differentiating B*B^-1 = I: Bdot*B^-1 + B*(B^-1)dot = 0, so
    # (B^-1)dot = -B^-1*Bdot*B^-1. This is an exact matrix-calculus identity, not an
    # approximation -- it holds for any invertible B(t), independent of what B physically
    # represents.
    Bdot = L * B + B * L.T
    Binv_dot = -Binv * Bdot * Binv
    J_rate_over_J = sp.trace(L)

    tau_dot = -J_rate_over_J * tau + (1 / J) * (
        2 * C1 * dev(Bdot, sp.trace(Bdot)) - 2 * C2 * dev(Binv_dot, sp.trace(Binv_dot))
    )

    print("\n" + "=" * 78)
    print("Step 2: Truesdell rate and objectivity check (same procedure as neo-Hookean)")
    print("=" * 78)

    tau_truesdell = sp.expand(tau_dot - L * tau - tau * L.T + sp.trace(L) * tau)

    d_sym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"d{min(i, j) + 1}{max(i, j) + 1}"))
    w12, w13, w23 = sp.symbols("w12 w13 w23")
    w_antisym = sp.Matrix([[0, w12, w13], [-w12, 0, w23], [-w13, -w23, 0]])
    subs_L = {L[i, j]: (d_sym + w_antisym)[i, j] for i in range(3) for j in range(3)}

    check = sp.expand(tau_truesdell.subs(subs_L))
    no_spin = check.subs({w12: 0, w13: 0, w23: 0})
    spin_dependence = sp.simplify(check - no_spin)
    assert spin_dependence == sp.zeros(3, 3), "Mooney-Rivlin Truesdell rate depends on spin -- not objective!"
    print("Verified: objective (no dependence on the spin part of L), same as neo-Hookean.")

    tau_truesdell_of_d = sp.simplify(no_spin)

    print("\n" + "=" * 78)
    print("Step 3: reduction check -- C2 = 0 must recover the neo-Hookean result exactly")
    print("=" * 78)

    G = sp.Symbol("G", positive=True)
    # Binv is independent of B in this derivation (needed for the objectivity check to be
    # meaningful for a general law), but the C2 = 0 reduction does not involve Binv at
    # all, so this substitution is clean regardless.
    neo_hookean_form = (2 * G / (3 * J)) * (trB * d_sym - sp.trace(B * d_sym) * I3)
    reduction = sp.simplify(tau_truesdell_of_d.subs({C1: G / 2, C2: 0}) - neo_hookean_form)
    assert reduction == sp.zeros(3, 3), "C2 = 0 does not recover the neo-Hookean Truesdell rate!"
    print("Confirmed: setting C2 = 0, C1 = G/2 recovers exactly the neo-Hookean result from above.")

    print("\n" + "=" * 78)
    print("Step 4: acoustic tensor and the classical small-strain shear modulus check")
    print("=" * 78)

    u1, u2, u3, n1, n2, n3 = sp.symbols("u1 u2 u3 n1 n2 n3")
    u_vec = sp.Matrix([u1, u2, u3])
    n_vec = sp.Matrix([n1, n2, n3])
    d_plane = sp.Rational(1, 2) * (u_vec * n_vec.T + n_vec * u_vec.T)

    tau_plane = tau_truesdell_of_d.subs({d_sym[i, j]: d_plane[i, j] for i in range(3) for j in range(3)})
    force_along_n = sp.simplify(tau_plane * n_vec)
    A_of_n = sp.simplify(force_along_n.jacobian(u_vec))

    # B^-1 is already expressed symbolically in terms of B's own entries (see Binv =
    # B.inv() above), so substituting B = I is enough -- B^-1 = I follows automatically.
    subs_undeformed = {B[i, j]: I3[i, j] for i in range(3) for j in range(3)}
    subs_undeformed[J] = 1
    A_undeformed = sp.simplify(A_of_n.subs(subs_undeformed))

    subs_n_x = {n1: 1, n2: 0, n3: 0}
    A_n_x = A_undeformed.subs(subs_n_x)
    longitudinal = sp.simplify((A_n_x * sp.Matrix([1, 0, 0]))[0])
    transverse = sp.simplify((A_n_x * sp.Matrix([0, 1, 0]))[1])

    # Classical rubber-elasticity identity: an incompressible Mooney-Rivlin solid's
    # small-strain shear modulus is mu = 2*(C1 + C2) (e.g. Ogden, Non-Linear Elastic
    # Deformations, Sec. 4.3 and 6.3). This is independent of the derivation above --
    # it's the textbook fact the acoustic-tensor result at B = I must reduce to.
    mu = 2 * (C1 + C2)
    print(f"Longitudinal eigenvalue (u = n = x-hat): {longitudinal}  (expected: 4*mu/3 = {sp.simplify(4 * mu / 3)})")
    print(f"Transverse eigenvalue   (u perp n):       {transverse}  (expected: mu = {mu})")
    assert sp.simplify(longitudinal - sp.Rational(4, 3) * mu) == 0
    assert sp.simplify(transverse - mu) == 0
    print("Both match the classical small-strain Mooney-Rivlin shear modulus mu = 2*(C1 + C2).")

    print("\n" + "=" * 78)
    print("Step 5: the GENERAL-B longitudinal correction (used by f_hyperelastic_signal_speed)")
    print("=" * 78)

    rho_c2 = sp.expand((n_vec.T * A_of_n * n_vec)[0])
    B_nn = (n_vec.T * B * n_vec)[0]
    Binv_nn = (n_vec.T * Binv * n_vec)[0]

    # Separate the C1-proportional and C2-proportional parts (the whole derivation above is
    # linear in C1 and C2 separately, by construction -- tau itself is 2*C1*(...) - 2*C2*(...)).
    c1_coeff = sp.simplify(rho_c2.subs(C2, 0) / C1)
    c2_coeff = sp.simplify(rho_c2.subs(C1, 0) / C2)

    c1_closed = (4 / (3 * J)) * (trB - B_nn)
    # General-B closed form for the C2 part, in terms of tr(B^-1) and its own direction-aligned
    # component Binv_nn -- NOT simply a mirror of the C1 form (the coefficient on Binv_nn is 5,
    # not 1: B^-1's rate law, -B^-1*Bdot*B^-1, is genuinely different from B's own, L*B+B*L^T).
    c2_closed = (4 / (3 * J)) * (5 * Binv_nn - trBinv)

    # n.n = 1 is needed for these to match and is not something sympy's simplify can exploit on
    # its own (same reason main()'s own Step 6 above checks numerically rather than
    # symbolically): verify over random UNIT n and random symmetric B instead.
    import random

    random.seed(7)
    max_err_c1 = 0.0
    max_err_c2 = 0.0
    for _ in range(20):
        v = [random.uniform(-1, 1) for _ in range(3)]
        norm = sum(c * c for c in v) ** 0.5
        n_num = [c / norm for c in v]
        B_num = sp.Matrix(3, 3, lambda i, j: random.uniform(0.5, 2.0) if i == j else random.uniform(-0.3, 0.3))
        B_num = (B_num + B_num.T) / 2
        subs_num = {n1: n_num[0], n2: n_num[1], n3: n_num[2], J: 1.3}
        for i in range(3):
            for j in range(i, 3):
                subs_num[B[i, j]] = B_num[i, j]
        max_err_c1 = max(max_err_c1, abs(float(c1_coeff.subs(subs_num)) - float(c1_closed.subs(subs_num))))
        max_err_c2 = max(max_err_c2, abs(float(c2_coeff.subs(subs_num)) - float(c2_closed.subs(subs_num))))
    print(f"C1 coefficient vs. closed form, max error over 20 random (unit n, symmetric B): {max_err_c1:.2e}")
    print(f"C2 coefficient vs. closed form, max error over 20 random (unit n, symmetric B): {max_err_c2:.2e}")
    assert max_err_c1 < 1e-10, "C1 part of the general-B correction does not match the closed form!"
    assert max_err_c2 < 1e-10, "C2 part of the general-B correction does not match the closed form!"
    print("Confirmed: rho*c_correction**2 = (4/3J)*[C1*(trB - B_nn) + C2*(5*Binv_nn - trBinv)]")
    print("(C1 part reduces to the existing neo-Hookean formula exactly when C2 = 0.)")

    print("\n" + "=" * 78)
    print("Step 6: f_hyperelastic_signal_speed implements the EXACT general-B formula")
    print("=" * 78)
    print("""
m_hyperelastic.fpp stashes the direction-aligned (sweep-direction) diagonal entries of B and
B^-1 per cell, alongside their traces and J (xi%end+5 onward) -- so, unlike the earlier
version of this fix, f_hyperelastic_signal_speed does not need an isotropic estimate of
Binv_nn: it reads the real B_nn and Binv_nn for whichever direction is being swept, exactly.
""")
    # Abstract symbols at the level of the closed FORM (c1_closed/c2_closed above are built
    # from literal matrix expressions, which subs() cannot cleanly re-substitute at this
    # level): confirm the code's own formula, read verbatim from m_riemann_state.fpp, is
    # algebraically identical to Step 5's derived closed form.
    B_nn_sym, Binv_nn_sym, trBinv_sym = sp.symbols("B_nn Binv_nn trBinv")
    exact_closed_form = (4 / (3 * J)) * (C1 * (trB - B_nn_sym) + C2 * (5 * Binv_nn_sym - trBinv_sym))
    # m_riemann_state.fpp: c*c + (2*Ca_inv*(trB-B_nn))/(3*J*rho) + (4*Ca_inv_2*(5*Binv_nn-trBinv))/(3*J*rho),
    # Ca_inv = 2*C1, Ca_inv_2 = C2, correction only (the rho*c*c/rho = c*c term is the base sound speed, not part
    # of this correction).
    implemented = (2 * (2 * C1) * (trB - B_nn_sym)) / (3 * J) + (4 * C2 * (5 * Binv_nn_sym - trBinv_sym)) / (3 * J)
    diff_impl = sp.simplify(exact_closed_form - implemented)
    print(f"  code's formula vs. Step 5's closed form: diff = {diff_impl}   (expect 0)")
    assert diff_impl == 0, "f_hyperelastic_signal_speed's formula does not match the exact general-B derivation!"
    print("Confirmed: the implemented formula is the exact general-B result from Step 5, for any")
    print("deformation state -- not an isotropic or small-strain approximation.")


def plane_strain_reduction():
    """Verifies the 1D and 2D closed forms proposed for m_hyperelastic.fpp against the full
    3D machinery already derived and checked above, by FREEZING the transverse direction(s)
    rather than assuming a reduced formula is correct on its own.

    Convention (matching hypoelasticity's own 1D/2D support, and MFC's dimensional-reduction
    convention elsewhere): directions beyond num_dims are frozen -- no stretch -- so the full
    3D deformation gradient is F_3D = block_diag(F_reduced, I_(3-num_dims)). This function
    builds F_3D explicitly from a general reduced grad_xi (2x2 for 2D, scalar for 1D),
    computes B_3D = F_3D F_3D^T and J_3D = det(F_3D) the SAME way the already-verified 3D code
    does, and checks that doing the algebra directly in the reduced (2x2 or scalar) space --
    the closed forms actually proposed for m_hyperelastic.fpp -- gives identical results.
    """
    print("\n" + "#" * 78)
    print("# 1D/2D reduction check")
    print("#" * 78)

    G, rho = sp.symbols("G rho", positive=True)
    I3 = sp.eye(3)

    print("\n" + "=" * 78)
    print("2D: grad_xi is 2x2, F_3D = block_diag(F_2D, 1)")
    print("=" * 78)

    m11, m12, m21, m22 = sp.symbols("m11 m12 m21 m22")
    M2 = sp.Matrix([[m11, m12], [m21, m22]])
    detM2 = M2.det()
    F2 = sp.simplify(M2.inv())
    J2_reduced = sp.simplify(1 / detM2)
    B2_reduced = sp.simplify(F2 * F2.T)
    trB2_reduced = sp.simplify(sp.trace(B2_reduced) + 1)  # +1 for the frozen B33 = 1

    # Full 3D embedding, built independently from the 2x2 pieces above.
    F3_embedded = sp.Matrix([[F2[0, 0], F2[0, 1], 0], [F2[1, 0], F2[1, 1], 0], [0, 0, 1]])
    J3_embedded = sp.simplify(F3_embedded.det())
    B3_embedded = sp.simplify(F3_embedded * F3_embedded.T)
    trB3_embedded = sp.simplify(sp.trace(B3_embedded))

    assert sp.simplify(J2_reduced - J3_embedded) == 0, "2D J does not match the 3D-embedded determinant"
    assert sp.simplify(trB2_reduced - trB3_embedded) == 0, "2D tr(B)+1 does not match the 3D-embedded trace"
    assert sp.simplify(B2_reduced[0, 0] - B3_embedded[0, 0]) == 0
    assert sp.simplify(B2_reduced[0, 1] - B3_embedded[0, 1]) == 0
    assert sp.simplify(B2_reduced[1, 1] - B3_embedded[1, 1]) == 0
    print("Confirmed: the 2x2 closed form for F, J, and B exactly matches embedding into the")
    print("full 3x3 machinery with the third row/column frozen to identity.")

    # Neo-Hookean stress: the in-plane components of tau = (G/J)*dev(B), computed the
    # reduced way, must match the (1,1)/(1,2)/(2,2) entries of the full 3D stress.
    def dev3(B, trB):
        return B - trB / 3 * I3

    tau3 = sp.simplify((G / J3_embedded) * dev3(B3_embedded, trB3_embedded))
    tau2_reduced = sp.simplify((G / J2_reduced) * (B2_reduced - trB2_reduced / 3 * sp.eye(2)))
    assert sp.simplify(tau2_reduced[0, 0] - tau3[0, 0]) == 0
    assert sp.simplify(tau2_reduced[0, 1] - tau3[0, 1]) == 0
    assert sp.simplify(tau2_reduced[1, 1] - tau3[1, 1]) == 0
    print("Confirmed: the 2D neo-Hookean in-plane stress (using trace = tr(B_2D) + 1, dividing")
    print("by the SAME 3 in dev(B) = B - tr(B)/3*I) matches the full 3D stress's in-plane block.")

    print("\n" + "=" * 78)
    print("1D: grad_xi is a scalar, F_3D = diag(F, 1, 1)")
    print("=" * 78)

    m = sp.Symbol("m", positive=True)
    F1 = 1 / m
    J1_reduced = F1  # det of a 1x1 is itself
    B1_reduced = F1**2
    trB1_reduced = B1_reduced + 2  # +2 for the two frozen directions

    F3_embedded_1d = sp.diag(F1, 1, 1)
    J3_embedded_1d = sp.simplify(F3_embedded_1d.det())
    B3_embedded_1d = sp.simplify(F3_embedded_1d * F3_embedded_1d.T)
    trB3_embedded_1d = sp.simplify(sp.trace(B3_embedded_1d))

    assert sp.simplify(J1_reduced - J3_embedded_1d) == 0
    assert sp.simplify(trB1_reduced - trB3_embedded_1d) == 0
    assert sp.simplify(B1_reduced - B3_embedded_1d[0, 0]) == 0
    print("Confirmed: the 1D scalar closed form for F, J, and B matches the 3x3 embedding.")

    tau3_1d = sp.simplify((G / J3_embedded_1d) * dev3(B3_embedded_1d, trB3_embedded_1d))
    tau1_reduced = sp.simplify((G / J1_reduced) * (B1_reduced - trB1_reduced / 3))
    assert sp.simplify(tau1_reduced - tau3_1d[0, 0]) == 0
    print("Confirmed: the 1D stress tau_xx = (G/J)*(B_xx - trace/3) matches the full 3D stress's")
    print("(1,1) entry, with trace = B_xx + 2 for the two frozen directions.")

    print("\n" + "=" * 78)
    print("General trace formula check: trace = sum(diagonal entries) + (3 - num_dims)")
    print("=" * 78)
    print("1D (num_dims=1, 1 diagonal entry at triangular position 1):  trace = B_xx + 2       -- confirmed above.")
    print("2D (num_dims=2, 2 diagonal entries at triangular positions 1,3): trace = B_xx+B_yy+1 -- confirmed above.")
    print("3D (num_dims=3, 3 diagonal entries at triangular positions 1,3,6): trace = B_xx+B_yy+B_zz+0, i.e. the")
    print("existing code's trace = btensor_in(1)+btensor_in(3)+btensor_in(6) -- the num_dims=3 special case of the")
    print("same formula, not a separate derivation.")


def mooney_rivlin_implementation_check():
    """Verifies the specific discrete identities s_Mooney_Rivlin_cauchy_solver (m_hyperelastic.fpp)
    relies on, so that implementation can avoid ever explicitly inverting B.

    1. B^-1 = grad_xi^T * grad_xi. Since F = (grad_xi)^-1 (the code's own convention), B^-1 =
       (F F^T)^-1 = (F^-1)^T F^-1 = grad_xi^T grad_xi -- the EXACT SAME "sum of products" formula
       the code already uses for B = F F^T (STEP 3 of s_hyperelastic_rmt_stress_update), just
       applied to grad_xi instead of F. This means B^-1 needs no separate matrix inversion in the
       general case; it only needed the code to also invert B directly via its own adjugate/det
       (det(B) = J**2), which is what is implemented, checked here against this identity.
    2. The second invariant I2 = 0.5*(tr(B)**2 - tr(B**2)) matches the textbook sum of principal
       2x2 minors.
    3. The component-wise stress formula (subtract trace/3 from diagonal entries only, leave
       off-diagonal entries alone, matching s_hyperelastic_trace_and_deviator's own convention)
       is algebraically identical to the matrix form tau = (1/J)*(2*C1*dev(B) - 2*C2*dev(B^-1)).
    4. The 1D/2D frozen-direction convention for tr(B**2) (diagonal terms squared, off-diagonal
       terms squared and doubled, plus (3 - num_dims) for the frozen unit diagonal entries)
       matches the full 3D embedding, the same way plane_strain_reduction() checked tr(B) itself.
    """
    print("\n" + "#" * 78)
    print("# Mooney-Rivlin implementation check (discrete identities used in the Fortran code)")
    print("#" * 78)

    print("\n" + "=" * 78)
    print("1. B^-1 = grad_xi^T * grad_xi")
    print("=" * 78)
    M = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"m{i + 1}{j + 1}"))  # grad_xi
    F = M.inv()
    B = sp.simplify(F * F.T)
    Binv_direct = sp.simplify(B.inv())
    Binv_identity = sp.simplify(M.T * M)
    diff = sp.simplify(Binv_direct - Binv_identity)
    assert diff == sp.zeros(3, 3), "B^-1 != grad_xi^T * grad_xi"
    print("Confirmed: B^-1 = grad_xi^T * grad_xi -- not needed in the final implementation (which")
    print("inverts B directly via its own adjugate/determinant instead), but confirms the identity")
    print("the two routes agree, and documents why a second inversion path was not necessary.")

    print("\n" + "=" * 78)
    print("2. Second invariant I2")
    print("=" * 78)
    Bxx, Bxy, Bxz, Byy, Byz, Bzz = sp.symbols("Bxx Bxy Bxz Byy Byz Bzz", real=True)
    B3 = sp.Matrix([[Bxx, Bxy, Bxz], [Bxy, Byy, Byz], [Bxz, Byz, Bzz]])
    trB, trB2 = sp.trace(B3), sp.trace(B3 * B3)
    I2_formula = sp.Rational(1, 2) * (trB**2 - trB2)
    I2_minors = Bxx*Byy - Bxy**2 + Byy*Bzz - Byz**2 + Bxx*Bzz - Bxz**2
    assert sp.simplify(I2_formula - I2_minors) == 0, "I2 formula does not match the principal minors"
    print("Confirmed: I2 = 0.5*(tr(B)**2 - tr(B**2)) matches the sum of principal 2x2 minors.")

    print("\n" + "=" * 78)
    print("3. Component-wise stress formula vs. the matrix form")
    print("=" * 78)
    import random

    random.seed(1)
    A = sp.Matrix(3, 3, lambda i, j: random.uniform(-0.5, 0.5))
    Bnum = sp.simplify(A * A.T + 3 * sp.eye(3))  # positive-definite-ish
    Jnum = sp.sqrt(Bnum.det())
    trBnum = sp.trace(Bnum)
    C1, C2 = sp.Rational(3, 2), sp.Rational(7, 10)
    Binv_num = Bnum.inv()
    trBinv_num = sp.trace(Binv_num)
    I3 = sp.eye(3)
    tau_matrix = (2 * C1 * (Bnum - trBnum / 3 * I3) - 2 * C2 * (Binv_num - trBinv_num / 3 * I3)) / Jnum

    def tau_component(i, j):
        if i == j:
            return (2 * C1 * (Bnum[i, i] - trBnum / 3) - 2 * C2 * (Binv_num[i, i] - trBinv_num / 3)) / Jnum
        return (2 * C1 * Bnum[i, j] - 2 * C2 * Binv_num[i, j]) / Jnum

    max_err = max(abs(float(sp.simplify(tau_matrix[i, j] - tau_component(i, j)))) for i in range(3) for j in range(3))
    assert max_err < 1e-12, f"component-wise formula mismatch: {max_err}"
    print(f"Confirmed: component-wise formula matches the matrix form exactly (max error {max_err:.1e}).")

    print("\n" + "=" * 78)
    print("4. 1D/2D frozen-direction tr(B**2)")
    print("=" * 78)
    Bxx, Bxy, Byy = sp.symbols("Bxx Bxy Byy", real=True)
    B2d_embedded = sp.diag(Bxx, Byy, 1)
    B2d_embedded[0, 1] = Bxy
    B2d_embedded[1, 0] = Bxy
    trB2_2d_formula = Bxx**2 + Byy**2 + 2 * Bxy**2 + 1  # num_dims = 2: one frozen direction, +1**2
    trB2_3d_direct = sp.trace(B2d_embedded * B2d_embedded)
    assert sp.simplify(trB2_2d_formula - trB2_3d_direct) == 0
    print("Confirmed: the 2D tr(B**2) formula (diagonal terms squared + 2*off-diagonal squared +")
    print("1 for the frozen z-direction) matches the full 3D embedding exactly. The 1D case is the")
    print("same formula with zero off-diagonal terms and (3-1)=2 frozen directions.")



def transverse_mode():
    """The general-B TRANSVERSE (shear) eigenvalue of the Mooney-Rivlin acoustic tensor.

    mooney_rivlin() above derives A(n) and checks its transverse eigenvalue at B = I against
    the classical small-strain shear modulus mu = 2*(C1 + C2). HLLD needs the general-B form:
    its inner waves are shear waves whose impedance is sqrt(rho*Ca_inv_eff), and for a
    hypoelastic material Ca_inv_eff = Ca_inv_hat + tau_nn_hat -- the shear modulus stiffened
    by the normal stress. A finite-strain law has no transported stress to read that off, so
    the shear modulus has to come from the acoustic tensor, exactly as the longitudinal speed
    does.

    For a shear wave the polarization is perpendicular to the propagation direction, so the
    quantity wanted is t.A(n).t with t.n = 0 (and |t| = |n| = 1). Fitting that against the
    frame invariants MFC already stashes gives, to machine precision:

        rho*c_T**2 = (2/(3*J)) * [ C1*tr(B) + C2*(3*(Binv_nn + Binv_tt) - tr(Binv)) ]

    Two things are worth noting about this result:

      * The C1 (neo-Hookean) part depends only on tr(B): the shear speed of a neo-Hookean
        solid is ISOTROPIC even at finite strain, independent of both n and t. Only the
        Mooney-Rivlin C2 term is direction dependent.
      * Every term is already carried per interface. tr(B) and tr(Binv) are stashed scalars,
        and Binv_nn / Binv_tt are the direction-aligned diagonals of Binv -- which MFC stores
        per coordinate direction, and the Riemann solve is dimensionally split, so n and t are
        always coordinate axes. No new stashed field is needed.

    At B = I, J = 1 this collapses to 2*(C1 + C2) = mu, the same classical identity
    mooney_rivlin() checks, which is the consistency check for the whole expression.
    """
    import random

    C1, C2, J = sp.symbols("C1 C2 J", positive=True)
    B = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"B{min(i, j) + 1}{max(i, j) + 1}"))
    L = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"L{i + 1}{j + 1}"))
    I3 = sp.eye(3)

    print("\n" + "#" * 78)
    print("# Transverse (shear) mode: the HLLD inner-wave impedance")
    print("#" * 78)

    Binv = B.inv()
    dev = lambda M: M - sp.trace(M) / 3 * I3
    tau = (2 * C1 * dev(B) - 2 * C2 * dev(Binv)) / J
    Bdot = L * B + B * L.T
    Binvdot = -Binv * Bdot * Binv
    tau_dot = ((2 * C1 * dev(Bdot) - 2 * C2 * dev(Binvdot)) / J) - tau * sp.trace(L)
    tau_truesdell = tau_dot - L * tau - tau * L.T + sp.trace(L) * tau

    d_of_L = sp.Rational(1, 2) * (L + L.T)
    tau_of_d = tau_truesdell.subs({L[i, j]: d_of_L[i, j] for i in range(3) for j in range(3)}, simultaneous=True)

    u1, u2, u3, n1, n2, n3 = sp.symbols("u1 u2 u3 n1 n2 n3")
    u_vec = sp.Matrix([u1, u2, u3])
    n_vec = sp.Matrix([n1, n2, n3])
    d_plane = sp.Rational(1, 2) * (u_vec * n_vec.T + n_vec * u_vec.T)
    tau_plane = tau_of_d.subs({L[i, j]: d_plane[i, j] for i in range(3) for j in range(3)}, simultaneous=True)
    A_of_n = (tau_plane * n_vec).jacobian(u_vec)

    print("\nA(n) built (same construction as mooney_rivlin(), step 4).")

    # The claimed closed form, as a function of the invariants.
    def closed_form(trB, trBinv, Binv_nn, Binv_tt, Jv):
        return (2 * (C1 * trB + C2 * (3 * (Binv_nn + Binv_tt) - trBinv))) / (3 * Jv)

    # Numeric verification over random SPD B and random orthonormal (n, t). Exact rationals
    # throughout, so "equal" means equal, not equal to a tolerance. simplify() cannot exploit
    # n.n = 1 and t.n = 0 symbolically, which is why this is done numerically -- the same
    # reason the general-B longitudinal check above is numeric.
    rng = random.Random(20250601)
    for trial in range(6):
        M = sp.Matrix(3, 3, lambda i, j: sp.Rational(rng.randint(-50, 50), 50))
        Bm = (M * M.T) + sp.Rational(3, 2) * I3  # symmetric positive definite
        Jv = sp.sqrt(Bm.det())
        nv = sp.Matrix([sp.Rational(rng.randint(-9, 9), 10) for _ in range(3)])
        if nv.norm() == 0:
            continue
        nv = nv / nv.norm()
        seed_t = sp.Matrix([1, 0, 0]) if abs(float(nv[0])) < sp.Rational(9, 10) else sp.Matrix([0, 1, 0])
        tv = seed_t - (seed_t.T * nv)[0] * nv
        tv = tv / tv.norm()
        assert sp.simplify((nv.T * tv)[0]) == 0, "t is not perpendicular to n"

        sub = {B[i, j]: Bm[i, j] for i in range(3) for j in range(3)}
        sub.update({J: Jv, n1: nv[0], n2: nv[1], n3: nv[2]})
        lhs = sp.expand((tv.T * A_of_n.subs(sub) * tv)[0])

        Bi = Bm.inv()
        rhs = closed_form(sp.trace(Bm), sp.trace(Bi), (nv.T * Bi * nv)[0], (tv.T * Bi * tv)[0], Jv)
        assert sp.simplify(sp.expand(lhs - rhs)) == 0, f"transverse closed form failed on trial {trial}"
    print("Confirmed on 6 random (B, n, t): t.A(n).t equals")
    print("    rho*c_T**2 = (2/(3*J)) * [ C1*tr(B) + C2*(3*(Binv_nn + Binv_tt) - tr(Binv)) ]")
    print("exactly, in rational arithmetic.")

    # Undeformed limit: the classical shear modulus.
    mu_claim = sp.simplify(closed_form(3, 3, 1, 1, 1))
    assert sp.simplify(mu_claim - 2 * (C1 + C2)) == 0
    print(f"\nAt B = I, J = 1: {mu_claim} = 2*(C1 + C2) = mu, the classical small-strain shear modulus.")

    # The neo-Hookean part carries no direction dependence at all.
    c1_part = sp.simplify(closed_form(sp.Symbol("trB"), 0, 0, 0, sp.Symbol("Jd")).coeff(C1))
    print(f"C1 part: {c1_part} -- a function of tr(B) and J only, so the neo-Hookean shear")
    print("speed is isotropic even at finite strain. Only the C2 term depends on n and t.")


if __name__ == "__main__":
    main()
    mooney_rivlin()
    plane_strain_reduction()
    mooney_rivlin_implementation_check()
    transverse_mode()
