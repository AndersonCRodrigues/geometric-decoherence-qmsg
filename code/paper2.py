#!/usr/bin/env python3
"""
Geometric Response and State-Resolved Decoherence:
Microscopic Born-Markov Derivation for a Commuting Qubit-Bath Coupling

Publication-grade verification + figure generation
Author: Anderson Costa Rodrigues

v2 additions (closing gaps found in a blind multi-referee review round):
  [A1c]  symbolic check that A(theta)^2 = I identically, at any theta and
         any temperature -- the actual reason the imaginary (Lamb-shift)
         part of the bath correlation function is dynamically trivial
         here, independent of the high-temperature approximation.
  [A2x]  exact spectral computation of sup_X ||L(X)||_1/||X||_1 via direct
         linear algebra on the 3 Hermitian traceless basis operators, with
         no optimizer. This is the primary numerical evidence for
         Theorem 1 (angle-independence); the pre-existing Nelder-Mead
         global search (A2) is kept only as a heuristic, independent
         confirmation, not as the main argument.
  [A3b]  Gaussian-approximation error reported BOTH relative to the local
         exact value (as before) and relative to the peak Gamma=L^2. The
         two are different numbers (~4.95% vs ~2.86%) and the manuscript
         previously conflated them; both are now printed explicitly and
         unambiguously labeled.
  [A4]   worked numerical example of Obstruction 1 (Section X.B of the
         manuscript): coherent-precession contamination of the Frobenius
         witness at omega_0/gamma = 0, 1, 10, 50.
  [A5]   exact vs. leading-order integration of purity loss along a fixed
         trajectory (Section X.A, Eq. 23), by directly solving the Bloch
         ODE and comparing against the naive extrapolation of the t=0
         rate. Quantifies the regime of validity (Delta P << 1) claimed
         in the revised manuscript, including the point at which the
         leading-order formula predicts an unphysical purity below 1/2.
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize

# ================================================================
# GLOBAL CONSTANTS
# ================================================================

print("=" * 80)
print("GEOMETRIC DECOHERENCE MODEL — PUBLICATION VERIFICATION")
print("=" * 80)

# ================================================================
# A1 — BATH CORRELATION & MARKOV COEFFICIENT
# ================================================================

eta, beta, wc, s = sp.symbols("eta beta wc s", positive=True)

C = (2 * eta / beta) * wc / (1 + wc**2 * s**2)

I_C = sp.integrate(C, (s, 0, sp.oo))
gamma_sym = 2 * I_C
norm_L_sym = 2 * gamma_sym

print("\n[A1] Analytical results")
print("∫C(s)ds =", I_C)
print("gamma =", gamma_sym)
print("||L|| =", norm_L_sym)

# Numerical check
eta_v = beta_v = wc_v = 1.0
C_num = lambda t: (2 * eta_v / beta_v) * wc_v / (1 + (wc_v * t) ** 2)

I_num, _ = quad(C_num, 0, np.inf)
gamma_num = 2 * I_num
norm_num = 2 * gamma_num

print("\n[A1] Numerical check")
print("||L|| =", norm_num)
print(f"Expected = {4*np.pi:.8f}")

tau_B = 1.0 / wc_v
print(f"\ntau_B = 1/omega_c = {tau_B:.4f}  (these units)")
print(f"gamma = {gamma_num:.4f}  (these units)")
print(
    "NOTE: the Markov validity condition is a separation of TIME SCALES,\n"
    "      t >> tau_B  (equivalently omega_c*t >> 1), quantified in [A1b]\n"
    "      below. It is NOT the same condition as gamma*t >> 1: gamma and\n"
    "      omega_c are independent parameters of the model, and the two\n"
    "      inequalities coincide only if gamma ~ omega_c, which is not\n"
    "      assumed anywhere. (This distinction was mislabeled in an\n"
    "      earlier manuscript draft and is corrected here.)"
)

# ----------------------------------------------------------------
# A1b — Markov tail-integral error bound (supports Appendix C,
# Step 4: separation of time scales tau_B << tau_S)
# ----------------------------------------------------------------

print("\n[A1b] Markov extension error: tail integral int_t^inf C(s) ds")
t_sym = sp.Symbol("t", positive=True)
tail_sym = sp.integrate(C, (s, t_sym, sp.oo))
tail_sym_simpl = sp.simplify(tail_sym)
print("Symbolic tail (int_t^inf C(s) ds) =", tail_sym_simpl)

for t_val in [1, 5, 10, 50]:
    tail_num, _ = quad(C_num, t_val, np.inf)
    rel_err = tail_num / I_num
    print(
        f"  t = {t_val:>3}: tail = {tail_num:.6e}  "
        f"(relative to total Markov integral: {rel_err:.4%})"
    )
print(
    "Tail vanishes as t -> infinity, confirming the Markov "
    "extension is controlled for tau_B = wc^-1 << tau_S."
)

# ----------------------------------------------------------------
# A1c — A(theta)^2 = I identically: the actual (temperature-independent)
# reason the imaginary part of C(t) never affects the dynamics.
# ----------------------------------------------------------------

print("\n[A1c] Symbolic check: A(theta)^2 = I for all theta, all temperature")
sx_sym = sp.Matrix([[0, 1], [1, 0]])
sz_sym = sp.Matrix([[1, 0], [0, -1]])
theta_sym = sp.Symbol("theta", real=True)
A_sym = sp.cos(theta_sym) * sx_sym + sp.sin(theta_sym) * sz_sym
A2_sym = sp.simplify(A_sym * A_sym)
print("A(theta)^2 =")
sp.pprint(A2_sym)
print(
    "A^2 = I identically (any theta). Any Lamb-shift term is therefore\n"
    "proportional to A^2 = I, hence commutes with everything and is\n"
    "dynamically trivial -- at ANY temperature, not only in the high-T\n"
    "limit used below to get the closed-form real part of C(t). The\n"
    "high-T limit is used only to obtain the specific Lorentzian shape\n"
    "of the real part (Eq. eq:Ct in the manuscript); it is not the\n"
    "reason the imaginary part drops out."
)

# ================================================================
# A2 — ANGLE-INDEPENDENCE OF THE GENERATOR NORM
# ================================================================

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def A(theta):
    return np.cos(theta) * sx + np.sin(theta) * sz


def apply_L(X, theta):
    Aop = A(theta)
    return gamma_num * (Aop @ X @ Aop - X)


def trace_norm(M):
    return np.sum(np.abs(np.linalg.eigvalsh(M)))


def bloch_X(r):
    return r[0] * sx + r[1] * sy + r[2] * sz


print("\n[A2x] EXACT spectral verification of Theorem 1 (primary evidence)")
print(
    "sup_X ||L(X)||_1/||X||_1 is, by Lemma 1, exactly the operator norm of\n"
    "the real 3x3 matrix (R_theta - I) acting on Bloch vectors. Evaluating\n"
    "on {sigma_x,sigma_y,sigma_z} is NOT generally enough to find this max\n"
    "(quadratic-form ratios are generally maximized on the eigenbasis of\n"
    "the operator, not on an arbitrary spanning set). It IS enough here\n"
    "for a model-specific reason: n_theta=(cos,0,sin) has zero y-component\n"
    "for every theta, so sigma_y is always orthogonal to n_theta, hence\n"
    "(Lemma 2) always lies in the maximal (-2) eigenspace of R_theta-I --\n"
    "sigma_y alone already attains the exact supremum for every theta;\n"
    "sigma_x, sigma_z are included only for redundancy. This does NOT\n"
    "generalize to a coupling family where n_theta leaves the xz-plane."
)
thetas_exact = np.linspace(0, np.pi, 9)
exact_sup = []
for th in thetas_exact:
    Aop = A(th)
    ratios = []
    for X in (sx, sy, sz):
        LX = gamma_num * (Aop @ X @ Aop - X)
        ratios.append(trace_norm(LX) / trace_norm(X))
    exact_sup.append(max(ratios))
exact_sup = np.array(exact_sup)
print("theta values:      ", np.round(thetas_exact, 4))
print("exact sup (3 evals):", np.round(exact_sup, 10))
print(
    f"All equal to 4*pi = {4*np.pi:.10f}? ",
    np.allclose(exact_sup, 4 * np.pi, atol=1e-10),
)
print(
    "This is the number reported as Theorem 1 in the manuscript: an exact\n"
    "algebraic identity, evaluated here on 3 fixed operators per angle\n"
    "(9 x 3 = 27 linear-algebra evaluations total), not a search."
)

print("\n[A2] Global numerical search (SECONDARY, heuristic confirmation only)")
print(
    "Kept for independent cross-validation of [A2x] by a method that does\n"
    "not presuppose Lemma 1 / Theorem 1; NOT the primary evidence."
)
thetas_test = np.linspace(0, np.pi, 9)
norms = []
for th in thetas_test:
    best = -1
    for _ in range(20):
        r0 = np.random.randn(3)
        r0 /= np.linalg.norm(r0)
        res = minimize(
            lambda r: -trace_norm(apply_L(bloch_X(r), th)) / trace_norm(bloch_X(r)),
            r0,
            method="Nelder-Mead",
            options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 10000},
        )
        val = -res.fun
        if val > best:
            best = val
    norms.append(best)

norms = np.array(norms)
print("Max ratio for each theta:", np.round(norms, 6))
print("All equal?", np.allclose(norms, norm_num, atol=1e-8))

# ================================================================
# A3 — GEOMETRIC RESPONSE FUNCTION Γ(θ)
# ================================================================

θ = sp.Symbol("θ")
ex, ez = sp.symbols("ex ez", real=True)

Gamma = (sp.cos(θ) * ex + sp.sin(θ) * ez) ** 2
G1 = sp.diff(Gamma, θ)
G2 = sp.diff(G1, θ)
G3 = sp.diff(G2, θ)

# Explicit critical point: solve dΓ/dθ = 0
critical_solutions = sp.solve(G1, θ)
print("\n[A3] Critical points from dΓ/dθ=0:", critical_solutions)
# SymPy may return atan(ez/ex) etc.

# Define μ from the correct branch
L_s = sp.sqrt(ex**2 + ez**2)
μ = sp.atan2(ez, ex)

# Verify that μ satisfies dΓ/dθ = 0 and is a maximum
check_G1 = sp.simplify(G1.subs(θ, μ))
check_G2 = sp.simplify(G2.subs(θ, μ))
print("Check dΓ/dθ(μ) =", check_G1)
print("Check d²Γ/dθ²(μ) =", check_G2)

# Geometric substitution (equivalent to evaluation at μ)
subs_mu = {sp.sin(θ): ez / L_s, sp.cos(θ): ex / L_s}
Gamma_mu = sp.simplify(Gamma.subs(subs_mu))
G2_mu = sp.simplify(G2.subs(subs_mu))
G3_mu = sp.simplify(G3.subs(subs_mu))

sigma_sq = sp.simplify(-Gamma_mu / G2_mu)
sigma_geo_sym = sp.simplify(sp.sqrt(sigma_sq))

print("Γ''(μ) =", G2_mu)
print("Γ'''(μ) =", G3_mu)
print("σ_geo =", sigma_geo_sym)

sigma_geo_val = float(sp.N(sigma_geo_sym))
print(f"σ_geo numerical = {sigma_geo_val:.8f}")

# ================================================================
# NUMERICAL RESPONSE FUNCTION
# ================================================================


def Gamma_func(theta, ex_val, ez_val):
    return (np.cos(theta) * ex_val + np.sin(theta) * ez_val) ** 2


alpha = np.radians(50)
ex_val = np.sin(2 * alpha)
ez_val = np.cos(2 * alpha)

theta_vals = np.linspace(0, np.pi, 2000)
Gamma_vals = Gamma_func(theta_vals, ex_val, ez_val)

mu_val = np.arctan2(ez_val, ex_val)
L_val = np.sqrt(ex_val**2 + ez_val**2)

Gamma_gauss = L_val**2 * np.exp(-((theta_vals - mu_val) ** 2) / (2 * sigma_geo_val**2))
dP = 2 * gamma_num * (Gamma_vals - 1)
dP_gauss = 2 * gamma_num * (Gamma_gauss - 1)

idx_1sigma = np.argmin(np.abs(theta_vals - (mu_val + sigma_geo_val)))
exact_1sigma = Gamma_vals[idx_1sigma]
gauss_1sigma = Gamma_gauss[idx_1sigma]
abs_error_1sigma = abs(exact_1sigma - gauss_1sigma)
peak_val = L_val**2

# Two DIFFERENT quantities -- do not conflate them (this conflation was
# the source of a real error in an earlier manuscript draft, caught in
# peer review: it reported "~4.9% of the peak value" when the number
# actually computed here is relative to the LOCAL value).
rel_error_1sigma_local = abs_error_1sigma / exact_1sigma * 100
rel_error_1sigma_peak = abs_error_1sigma / peak_val * 100

print(f"\nAbsolute error of Gaussian at 1σ: {abs_error_1sigma:.6f}  (peak L^2 = {peak_val:.6f})")
print(f"Relative error vs LOCAL exact value at 1σ: {rel_error_1sigma_local:.2f}%")
print(f"Relative error vs PEAK value (L^2):        {rel_error_1sigma_peak:.2f}%")
print(
    "These are two different reference points and give two different\n"
    "percentages; the manuscript now reports both explicitly and labels\n"
    "which is which."
)
# kept for backward compatibility with the rest of this script / older refs
rel_error_1sigma = rel_error_1sigma_local

# ================================================================
# A4 — OBSTRUCTION 1: COHERENT-PRECESSION CONTAMINATION OF W_bio
# (worked example supporting Section X.B / Appendix E of the manuscript)
# ================================================================

print("\n" + "=" * 80)
print("[A4] Obstruction 1: coherent contamination of the Frobenius witness")
print("=" * 80)
print(
    "theta=0 (A=sigma_x, dephasing preserves x), plus H_S=(omega_0/2)sigma_z\n"
    "added on TOP of the commuting model (this add-on breaks [H_S,A]=0 by\n"
    "construction; it is used here only as an illustrative counterexample,\n"
    "not as part of the commuting model itself -- see manuscript Appendix E).\n"
    "State: equatorial, r=(0,1,0). Bloch eqs: rdot_x=-w0*r_y,\n"
    "rdot_y=w0*r_x-2*gamma*r_y, rdot_z=-2*gamma*r_z."
)
gamma_A4 = 1.0
r_eq = np.array([0.0, 1.0, 0.0])
print(f"\n{'omega_0/gamma':>14} | {'||rhodot||_F / gamma':>22} | {'dP/dt':>10}")
for w0_over_g in [0, 1, 10, 50]:
    w0 = w0_over_g * gamma_A4
    rdot = np.array(
        [-w0 * r_eq[1], w0 * r_eq[0] - 2 * gamma_A4 * r_eq[1], -2 * gamma_A4 * r_eq[2]]
    )
    frob = np.linalg.norm(rdot) / np.sqrt(2)
    dPdt = np.dot(r_eq, rdot)
    print(f"{w0_over_g:14.1f} | {frob/gamma_A4:22.4f} | {dPdt:10.4f}")
print(
    "dP/dt is exactly -2*gamma for every omega_0 (purity is blind to\n"
    "unitary precession); ||rhodot||_F grows without bound. This is the\n"
    "numerical content behind Obstruction 1 in the manuscript."
)

# ================================================================
# A5 — EXACT vs. LEADING-ORDER TRAJECTORY INTEGRATION
# (quantifies the Delta P << 1 regime of validity of Eq. 23)
# ================================================================

from scipy.integrate import solve_ivp

print("\n" + "=" * 80)
print("[A5] Exact vs. leading-order (Eq. 23) purity loss along a trajectory")
print("=" * 80)
print(
    "Fixed theta=0 (n_hat=(1,0,0)), purely dissipative (no coherent term),\n"
    "worst-case initial state r=(0,1,0) (maximal decay rate, Gamma=0\n"
    "throughout since the n_hat-component of r never changes at fixed\n"
    "theta). Exact: solve dr/dt = gamma*(R-I)r. Leading-order/naive:\n"
    "extrapolate the t=0 rate linearly, Delta P_naive = 2*gamma*T*(Gamma-1),\n"
    "i.e. exactly Eq. 23 evaluated as if the state stayed pure throughout."
)

gamma_A5 = 1.0
r0_A5 = np.array([0.0, 1.0, 0.0])
Rmat_A5 = np.diag([0.0, -2.0, -2.0])  # R_theta - I at theta=0, in the sx,sy,sz basis


def rhs_A5(t, r):
    return gamma_A5 * Rmat_A5 @ r


P0_A5 = (1 + np.dot(r0_A5, r0_A5)) / 2
Gamma0_A5 = r0_A5[0] ** 2

T_list = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
print(f"\n{'T (1/gamma)':>12} | {'dP exact':>10} | {'dP naive (Eq.23)':>17} | {'ratio':>7} | naive physical?")
traj_rows = []
for T in T_list:
    sol = solve_ivp(rhs_A5, [0, T], r0_A5, dense_output=True, rtol=1e-11, atol=1e-13)
    rT = sol.y[:, -1]
    PT = (1 + np.dot(rT, rT)) / 2
    dP_exact = PT - P0_A5
    dP_naive = 2 * gamma_A5 * T * (Gamma0_A5 - 1)
    physical = "yes" if -0.5 <= dP_naive <= 0 else "NO (< -0.5, impossible for a qubit)"
    print(f"{T:12.2f} | {dP_exact:10.4f} | {dP_naive:17.4f} | {dP_naive/dP_exact:7.3f} | {physical}")
    traj_rows.append((T, dP_exact, dP_naive))
print(
    "\nEq. 23 is accurate only while |Delta P| << 1/2 (the full range available\n"
    "to a qubit); for this worst-case trajectory it becomes unphysical\n"
    "(predicts purity below the maximally-mixed floor of 1/2) beyond\n"
    "T ~ 0.3/gamma. The manuscript now states this regime explicitly rather\n"
    "than presenting Eq. 23 as valid for arbitrary T."
)

# Figure 4: exact vs naive trajectory purity loss
T_fine = np.linspace(0.001, 1.2, 400)
dP_exact_fine = []
for T in T_fine:
    sol = solve_ivp(rhs_A5, [0, T], r0_A5, rtol=1e-10, atol=1e-12)
    rT = sol.y[:, -1]
    dP_exact_fine.append((1 + np.dot(rT, rT)) / 2 - P0_A5)
dP_exact_fine = np.array(dP_exact_fine)
dP_naive_fine = 2 * gamma_A5 * T_fine * (Gamma0_A5 - 1)

plt.figure(figsize=(7, 5))
plt.plot(T_fine, dP_exact_fine, label="Exact (ODE)", linewidth=2)
plt.plot(T_fine, dP_naive_fine, "--", label="Leading-order (Eq. 23)", linewidth=2)
plt.axhline(-0.5, color="red", linestyle=":", linewidth=1, label="unphysical floor (ΔP=-1/2)")
plt.xlabel("T (units of 1/γ)")
plt.ylabel("ΔP(T)")
plt.title("Exact vs. leading-order purity loss (worst-case trajectory)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("figure_trajectory.png", dpi=200)
plt.savefig("figure_trajectory.pdf", bbox_inches="tight")
print("\nFigure 4 saved: figure_trajectory.png/pdf")

# ================================================================
# FIGURES
# ================================================================

# Figure 1: Gamma
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(np.degrees(theta_vals), Gamma_vals, label="Exact Γ(θ)", linewidth=2)
plt.plot(
    np.degrees(theta_vals),
    Gamma_gauss,
    "--",
    label="Local Gaussian approx",
    linewidth=2,
)
plt.axvline(
    np.degrees(mu_val),
    linestyle=":",
    color="gray",
    label=f"μ = {np.degrees(mu_val):.1f}°",
)
plt.xlabel("θ (degrees)")
plt.ylabel("Γ(θ)")
plt.title("Geometric Response Function")
plt.legend()
plt.grid()

plt.subplot(1, 2, 2)
mask = np.abs(theta_vals - mu_val) <= 1.5 * sigma_geo_val
plt.plot(np.degrees(theta_vals[mask]), Gamma_vals[mask], label="Exact", linewidth=2)
plt.plot(
    np.degrees(theta_vals[mask]),
    Gamma_gauss[mask],
    "--",
    label="Gaussian approx",
    linewidth=2,
)
plt.axvline(np.degrees(mu_val), linestyle=":", color="gray")
plt.axvline(
    np.degrees(mu_val + sigma_geo_val),
    linestyle="--",
    color="purple",
    alpha=0.7,
    label="±σ_geo",
)
plt.axvline(
    np.degrees(mu_val - sigma_geo_val), linestyle="--", color="purple", alpha=0.7
)
plt.xlabel("θ (degrees)")
plt.ylabel("Γ(θ)")
plt.title("Zoom around maximum")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("figure_gamma.png", dpi=200)
plt.savefig("figure_gamma.pdf", bbox_inches="tight")
print("\nFigure 1 saved: figure_gamma.png/pdf")

# Figure 2: Purity decay
# NOTE: not included in the main text (paper2.tex). dP/dt is a fixed
# rescaling of Gamma(theta) [Eq. (eq:pureloss): dP/dt = 2*gamma*(Gamma-1)],
# so it carries no information beyond Fig. 1 (figure_gamma) and was
# dropped from the manuscript to avoid a redundant panel. Kept here as
# supplementary output for readers who want the purity-rate view directly.
plt.figure(figsize=(10, 5))
plt.plot(np.degrees(theta_vals), dP, label="Exact purity loss", linewidth=2)
plt.plot(np.degrees(theta_vals), dP_gauss, "--", label="Gaussian approx", linewidth=2)
plt.axhline(0, color="black", linewidth=0.8)
plt.axvline(
    np.degrees(mu_val),
    linestyle=":",
    color="gray",
    label=f"μ = {np.degrees(mu_val):.1f}°",
)
plt.xlabel("θ (degrees)")
plt.ylabel("dP/dt|₀")
plt.title("Initial Purity Decay Rate")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("figure_purity.png", dpi=200)
plt.savefig("figure_purity.pdf", bbox_inches="tight")
print("Figure 2 saved: figure_purity.png/pdf")

# Figure 3: Error
error_abs = np.abs(Gamma_vals - Gamma_gauss)
error_rel = error_abs / (Gamma_vals + 1e-15)
plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.semilogy(np.degrees(theta_vals), error_abs)
plt.axvline(np.degrees(mu_val), linestyle=":", color="gray")
plt.axvline(np.degrees(mu_val + sigma_geo_val), linestyle="--", color="purple")
plt.axvline(np.degrees(mu_val - sigma_geo_val), linestyle="--", color="purple")
plt.xlabel("θ (degrees)")
plt.ylabel("Absolute error")
plt.title("Gaussian Approximation Error")
plt.grid()

plt.subplot(1, 2, 2)
plt.semilogy(np.degrees(theta_vals), error_rel)
plt.axvline(np.degrees(mu_val), linestyle=":", color="gray")
plt.axvline(
    np.degrees(mu_val + sigma_geo_val), linestyle="--", color="purple", label="±σ_geo"
)
plt.axvline(np.degrees(mu_val - sigma_geo_val), linestyle="--", color="purple")
plt.axhline(0.05, color="black", linestyle=":", linewidth=0.8, label="5%")
plt.xlabel("θ (degrees)")
plt.ylabel("Relative error")
plt.title("Relative Error (Gaussian vs Exact)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("figure_error.png", dpi=200)
plt.savefig("figure_error.pdf", bbox_inches="tight")
print("Figure 3 saved: figure_error.png/pdf")

# ================================================================
# FINAL SUMMARY
# ================================================================

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)
print(f"gamma = {gamma_num:.8f}")
print(f"tau_B = {tau_B:.8f}")
print(f"||L|| = {norm_num:.8f}  (exact spectral check: {np.allclose(exact_sup, 4*np.pi, atol=1e-10)})")
print(f"σ_geo = {sigma_geo_val:.8f}")
print(f"Gaussian error at 1σ, vs local value = {rel_error_1sigma_local:.2f}%")
print(f"Gaussian error at 1σ, vs peak value  = {rel_error_1sigma_peak:.2f}%")
print(f"A(theta)^2 = I identically: confirmed symbolically (all theta, all T)")
print(f"Obstruction-1 example: dP/dt = -2*gamma independent of omega_0; ||rhodot||_F unbounded")
print(f"Eq.23 leading-order trajectory formula: physical only for T ~< 0.3/gamma (worst case)")
