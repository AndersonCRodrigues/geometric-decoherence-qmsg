# Audit trail

This paper went through multiple independent verification passes before
submission, documented here for transparency rather than left implicit.

## Round 1 — internal audit
Full read of the manuscript against its own reproducibility script.
Findings included: an unjustified reliance on numerical search where a
closed-form proof was available (Theorem 1); an unquantified Markov
truncation error already computed by the script but never reported in
the text; a missing mathematical bridge to the companion paper's
phenomenological witness; an incorrectly typed self-citation causing
BibTeX compilation errors; four orphaned bibliography entries.

## Round 2 — blind review, 6 independent models
Prompted with no knowledge of Round 1's findings. Cross-checked every
claim against the manuscript source and the script before accepting or
rejecting it. Confirmed and fixed: an incomplete derivation of the
central norm theorem's key algebraic step; an unstated adiabatic
condition; a mischaracterized "open problem" that was in fact a proven
negative result; an overstated decoherence-free-subspace analogy.
Rejected as unfounded after direct verification: at least one specific
claim in each of two reports that did not match the actual source when
checked.

## Round 3 — blind review, second batch (5 independent models, different prompt)
Found: a genuine sign error in the trajectory purity-loss equation,
confirmed by direct inspection to contradict the manuscript's own
script, figure, and an earlier equation in the same text; a
mathematically invalid general justification for a numerical
cross-check that happened to give the right answer for a model-specific
reason; a PDF text-extraction (ligature) defect; missing connections to
the geometric-phase and filter-function literature. Two specific
numerical claims in the reports were checked against the source and
found to be fabricated; both are documented as rejected rather than
silently ignored.

## Round 4 — derivation on request
A reviewer noted that the adiabatic condition had been justified by
analogy rather than derived. Re-derived it via an exact change of
frame (an appendix dedicated to it), verified symbolically at every step before writing
it up; explicitly stated what remains at the level of analogy
(the two secondary conditions) rather than claiming the derivation is
now complete in every respect.

## A related, incorrect derivation not adopted
A proposed extension analyzing robustness to a transverse Hamiltonian
perturbation was checked and found to rest on an incorrect frequency
scale (it should scale with the dominant coupling frequency, not the
perturbation strength, when the perturbation is small); the extension
was not included in the manuscript for this reason, independent of
the separate decision not to pursue that direction at all in this
version (see main text discussion of scope).

## What this means for a reader
Every numerical value in the manuscript is reproducible by running
`code/paper2.py`; console output should be compared against
`code/run_final.log`. Where the manuscript reports a limitation or an
open problem, it does so because a specific check failed or a specific
gap was identified, not as a hedge.
