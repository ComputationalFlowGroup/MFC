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


if __name__ == "__main__":
    main()


# --------------------------------------------------------------------------------------
# Extending this to Mooney-Rivlin (src/common/m_hyperelastic.fpp,
# s_Mooney_Rivlin_cauchy_solver):
#
# Mooney-Rivlin's stress depends on the SECOND invariant of B as well as the first, so
# tau is a function of B and B**2 (there are several equivalent ways to write this
# depending on convention -- check against the exact form in
# s_Mooney_Rivlin_cauchy_solver before changing anything below). The only changes needed
# are:
#   1. Replace the `tau = (G/J) * (B - trB/3*I3)` line with the Mooney-Rivlin stress.
#   2. Bdot is unchanged (B's kinematics do not depend on the constitutive law), but you
#      will also need d(B**2)/dt = Bdot*B + B*Bdot (product rule), since tau now depends
#      on B**2.
# Steps 2-6 (the Truesdell rate, objectivity check, plane-wave substitution, and
# small-strain check) carry over unchanged -- they are generic to any stress law
# expressed as a function of B and J.
# --------------------------------------------------------------------------------------
