# Verification Status — M3R-NOV-02 Public v1.0

**DOI:** `10.5281/zenodo.23063122`

## Current status

- Internal algebraic reconstruction: **PASS**
- Independent dual-protocol symbolic replay: **PASS**
- Theorem-line replay: **PASS**
- Exact Sturm root count: **PASS**
- Rational bracket/enclosure checks: **PASS**
- External line-by-line review: **NOT OBTAINED**
- Peer review: **NOT OBTAINED**
- Formal proof-assistant verification of this manuscript: **NOT PERFORMED**
- Worldwide novelty / priority: **NOT CLAIMED**

## Replayed on 30 September 2026

The two public verification scripts terminate with:

`ALL MATCHED AND FIXED-SCALE SYMBOLIC CHECKS PASS`

and

`ALL R2 THEOREM-LINE ALGEBRAIC CHECKS PASS`

The replay certifies the exact rational brackets

`0.63651056784173490941 < q* < 0.63651056784173490942`

and

`5.78812722135752519109224751394 < R_min < 5.78812722135752519155412626161`.

## Suggested hostile checks

A useful external review should attempt to identify the earliest failure in one of the following:

1. divergence freedom or symmetry equivariance of the explicit initial datum;
2. symmetry forcing of the gradient or deviatoric covariance tensor shape;
3. the fixed-$\tau$ quantifier structure for the nonzero physical-time alignment interval;
4. the matched-path chain rule and sign of $\dot\tau=-1$;
5. the exact matched threshold and degree-10 derivative polynomial;
6. the exact Sturm count and rational enclosure;
7. the fixed-scale derivative and its distinct threshold;
8. the proof that the fixed-scale threshold is strictly decreasing with endpoint limits $+\infty$ and $3$;
9. any hidden overclaim across the matched-path / fixed-scale boundary;
10. known prior art containing essentially the same family-specific theorem.

A precise earliest failing identity or implication is more useful than a general assessment.
