# M3R-NOV-02 — Public v1.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23063122.svg)](https://doi.org/10.5281/zenodo.23063122)

## An Exact Symmetry-Locked Physical-Time Counterconstruction to First-Order Angular Turnover with Matched-Path and Fixed-Scale Transfer-Threshold Classifications

**Author:** Alexanja Senke  
**Resource type:** Public technical note / independent mathematical research  
**Public v1.0 DOI:** [10.5281/zenodo.23063122](https://doi.org/10.5281/zenodo.23063122)  
**Prior restricted v0.3 DOI:** [10.5281/zenodo.22673450](https://doi.org/10.5281/zenodo.22673450)

This package contains the canonical public v1.0 manuscript, its LaTeX source, the exact public claim boundary, verification-status notes, and reproducible symbolic replay material.

The result is **family-specific and local in physical time**. It distinguishes:

1. a calorically matched observation path $\tau(t)=T-t$, and
2. physical-time differentiation at fixed positive heat scale.

For the named explicit periodic Navier–Stokes family, the two protocols produce different exact initial transfer classifications. The matched-path threshold has a unique certified interior minimum near `5.788127221357525...`, whereas the fixed-scale threshold is strictly decreasing with limiting infimum `3`.

## Claim boundary

This release does **not** claim:

- global Navier–Stokes regularity;
- finite-time singularity formation;
- a Millennium Prize solution;
- a universal Reynolds threshold;
- a heat-scale-uniform persistence interval;
- worldwide novelty or priority;
- independent line-by-line certification of the complete manuscript.

See `verification/CLAIM_BOUNDARY.md`.

## Verification

The public package includes two independent symbolic replay scripts inherited from the frozen review state:

- `verification/replay/independent_symbolic_verification.py`
- `verification/replay/r2_theorem_line_audit.py`

Both were replayed successfully on 30 September 2026. See `verification/VERIFICATION_STATUS.md`.

## Repository/package structure

- `paper/main.pdf` — canonical public v1.0 manuscript.
- `paper/main.tex` — LaTeX source.
- `verification/CLAIM_BOUNDARY.md` — exact public claim boundary.
- `verification/VERIFICATION_STATUS.md` — verification state and hostile-review targets.
- `verification/replay/` — reproducible symbolic replay scripts and frozen result summaries.
- `CITATION.cff` — machine-readable citation metadata.
- `ZENODO_METADATA.md` — public Zenodo metadata.
- `RELEASE_NOTES_v1.0.md` — changes from restricted v0.3.
- `CONTRIBUTING.md` — defect-reporting guidance.
- `LICENSE.txt` — rights statement.
- `SHA256_MANIFEST.txt` — integrity hashes.

## Provenance

The earlier v0.3 record remains the restricted historical review state at DOI `10.5281/zenodo.22673450`.

The public v1.0 release uses DOI `10.5281/zenodo.23063122`.

Private correspondence, internal handovers, protected parallel research branches, and superseded exploratory artifacts are intentionally excluded from this public package.
