# Geometric Response and State-Resolved Decoherence

[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#)
<!-- Anderson: substitua o badge acima pelo badge real do Zenodo assim que o DOI existir,
     no mesmo formato do badge usado no repositório do Paper 1. -->

Microscopic Born–Markov derivation for a commuting qubit–bath coupling.

This repository contains the complete derivation, symbolic/numerical
verification, and figure-generation pipeline for **Paper 2** of the
biographical-irreducibility program (companion to
[Path-Dependent Decoherence in Open Quantum Systems](https://doi.org/10.5281/zenodo.20896384)).

Within an exactly solvable model — a qubit coupled to a thermal
bosonic bath through an interaction operator that commutes with the
system Hamiltonian — this work derives the Born–Markov dissipative
generator in closed form, without a secular approximation, and
establishes:

1. a **closed-form proof** (not a numerical search) that the induced
   trace norm of the generator is independent of the coupling angle;
2. an **exact geometric response**, $\Gamma(\theta) = L^2\cos^2(\theta-\mu)$,
   for the state-resolved purity-loss rate, with a curvature-derived
   angular scale $\sigma_{\rm geo} = 1/\sqrt2$ that is a parameter-free
   invariant of the coupling family;
3. an explicit, quantified account of the **algebraic (not
   exponential) decay** of the Markov truncation error for this bath;
4. a direct test of whether this microscopic result reduces to the
   phenomenological trajectory witness $W_{\rm bio}$ of the companion
   paper — with the answer proven **negative**, and the specific
   obstruction identified rather than left implicit.

> **Scope note.** The commuting condition $[H_S(\theta),A(\theta)]=0$
> is a significant restriction: it reduces the model to pure dephasing
> in a fixed basis, and this is stated explicitly in the manuscript's
> introduction, not left for the reader to infer. This paper does
> **not** establish that $\gamma \propto W_{\rm bio}$; it proves the
> opposite — that the two are not proportional in general — and
> identifies exactly why (Section X.B, Proposition). If you are
> reading this after updating the companion repository's own scope
> note, make sure the two notes agree with each other and with this
> one; they currently do not (see below).

---

## A note on consistency with the companion repository

The [companion repository](https://github.com/AndersonCRodrigues/path-dependent-decoherence-qmsg)'s
README currently states that this paper "shows that
γ ≈ α_eff · W_bio emerges as a leading-order approximation." That is
not what this paper shows. Before both papers are public in their
current form, the companion repository's scope note should be
rewritten to match the actual result below — a proven
non-proportionality, with an explicitly narrower open question about
whether any restricted or reformulated relation can be established —
rather than the leading-order-approximation claim it currently makes.

---

## The Paper

[**View the manuscript (PDF)**](paper/paper2.pdf) · [LaTeX source](paper/paper2.tex) · [Bibliography](paper/references.bib)

---

## Key Results

All values below are reproducible by running the script in `/code`
(deterministic; no random seed needed except for the secondary,
non-essential Nelder–Mead cross-check).

| Quantity | Value |
|---|---|
| Dissipative rate $\gamma$ | $2\pi\eta/\beta$ |
| Induced trace norm $\lVert\mathcal L^{(D)}_\theta\rVert$ | $4\pi\eta k_BT$, exact for every $\theta$ |
| Angular scale $\sigma_{\rm geo}$ | $1/\sqrt2 \approx 40.5^\circ$ |
| Gaussian approx. error at $1\sigma_{\rm geo}$, vs. local value | $4.95\%$ |
| Gaussian approx. error at $1\sigma_{\rm geo}$, vs. peak value | $2.86\%$ |
| Markov tail at $t=10\,\tau_B$ | $6.3\%$ uncorrected |
| Markov tail at $t=50\,\tau_B$ | $1.3\%$ uncorrected |
| Obstruction 1: $\lVert\dot\rho\rVert_F/\gamma$ at $\omega_0/\gamma=0,1,10,50$ | $1.41,\,1.58,\,7.21,\,35.4$ |
| Obstruction 1: $dP/dt$ across the same range | exactly $-2\gamma$, unchanged |
| Eq. (23) leading-order trajectory formula, valid for | $\gamma T \lesssim 0.2$–$0.3$ (worst case), unphysical beyond |

---

## Reproducing All Results

```bash
git clone https://github.com/AndersonCRodrigues/geometric-decoherence-qmsg.git
cd geometric-decoherence-qmsg

pip install -r requirements.txt

python code/paper2.py
```

Figures are written to the working directory; compare console output
against `code/run_final.log`. The script performs every symbolic
derivation and numerical check cited in the manuscript, including the
exact spectral verification of Theorem 1 (no optimizer), the symbolic
confirmation that $A(\theta)^2=I$ for all $\theta$ (all temperatures),
and the exact-vs-leading-order trajectory comparison behind Figure 4.

---

## Repository Structure

```
/paper
  paper2.tex          LaTeX source (revtex4-2)
  references.bib      Bibliography (22 entries, all cited, all verified)
  paper2.pdf           Compiled manuscript
/code
  paper2.py            Full derivation + verification + figure generation
  requirements.txt     numpy, scipy, sympy, matplotlib (only what is used)
  run_final.log        Console output from the run above
/figures
  figure_gamma.pdf      Fig. 1 — exact vs. Gaussian response
  figure_error.pdf      Fig. 2 — approximation error, local vs. peak
  figure_trajectory.pdf Fig. 4 — exact vs. leading-order trajectory integration
/docs
  AUDIT_TRAIL.md        What was checked, what was fixed, and why
```

---

## Validation Summary

| Check | Result | Status |
|---|---|---|
| Angle-independence, exact spectral evaluation (3 basis operators) | matches $4\pi\eta k_BT$ to $10^{-10}$ | ✓ |
| Angle-independence, heuristic Nelder–Mead cross-check | matches to $10^{-8}$ (secondary, non-essential) | ✓ |
| $A(\theta)^2=I$, symbolic, all $\theta$ | confirmed identically | ✓ |
| $\sigma_{\rm geo}=1/\sqrt2$, symbolic Taylor matching | confirmed identically | ✓ |
| Rotating-frame identity behind the adiabatic condition (App. G) | confirmed identically, SymPy | ✓ |
| Exact vs. leading-order trajectory integration (Fig. 4) | naive formula unphysical beyond $\gamma T\approx0.3$ | ✓ (bound quantified, not hidden) |
| Non-proportionality of $\Delta P$ and $W_{\rm bio}$ | proven, two independent obstructions | ✓ |
| BibTeX compilation against `apsrev4-2.bst` | 0 errors, 22/22 references cited | ✓ |

---

## Citation

```bibtex
@article{rodrigues2026geometric,
  title   = {Geometric Response and State-Resolved Decoherence:
             Microscopic Born-Markov Derivation for a Commuting
             Qubit-Bath Coupling},
  author  = {Rodrigues, Anderson Costa},
  year    = {2026},
  note    = {DOI to be assigned on deposit / acceptance.
             Code: https://github.com/AndersonCRodrigues/geometric-decoherence-qmsg}
}
```
<!-- Anderson: substitua o bloco acima pela entrada final assim que o DOI existir
     (Zenodo, no mesmo padrão do Paper 1, e/ou a versão aceita pela revista). -->

---

## License

MIT — see `LICENSE`. (Matches the license used in the companion
repository; change here if you'd rather use CC-BY for the paper text
specifically and keep MIT only for the code.)
