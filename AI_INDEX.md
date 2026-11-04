# AI index: divergent isolated Kalman covariances (v0.2.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) and the
paper [isokf_manuscript.pdf](isokf_manuscript.pdf) (LaTeX source [manuscript.tex](manuscript.tex), bibliography [refs.bib](refs.bib)).
Version DOI 10.5281/zenodo.23941896 (0.2.0; corrects 0.1.0, 10.5281/zenodo.23935675); concept DOI 10.5281/zenodo.23935674. Target: the convergence conjecture for the
isolated Kalman filter (IsoKF) of Jung, Luft and Weiss, *Autonomous Robots* 49(1):7, 2025, DOI 10.1007/s10514-025-10191-x
(Conjecture 1 in the manuscript, quoted verbatim).

## Notation

Linear time-invariant system x_{k+1} = A x_k + w_k, block-diagonal A and Q, node states (scalar in the counterexamples).
An observation with participant set U, matrix H_U and noise R updates the participant block jointly (S = H_U B H_U^T + R,
C = B - B H_U^T S^{-1} H_U B with B = P[U,U]); every cross-covariance between a participant i and a nonparticipant j is
rescaled by the marginal ratio D_i = C_ii B_ii^{-1} (equations (2.1)-(2.3); published eqs. (47)-(50), Table 2).
Triangle schedule: H rows (1,1,0), (0,1,1), (1,0,1) processed in that order each epoch, R = I. m_k = min_i P_ii after epoch k.

## Claims and where they are proved

| ID | Claim | Location | Evidence type |
|---|---|---|---|
| C1 | Recurrence (2.1)-(2.3) is the published IsoKF update; buffered implementation equals the dense array for LTI models | Section 2, Table 2 | transcription check and argument |
| C2 | Schur-stable block-diagonal A: reported marginals uniformly bounded for every schedule | Proposition 3.1 | analytic proof |
| C3 | Triangle, A = 3I, Q = R = P_0 = I: completely observable (one epoch) and controllable, all local blocks positive definite, m_k >= (5/4)(9/5)^k - 1/4 -> infinity, while the centralised KF keeps every marginal <= 3/4 | Theorem 4.1, Lemmas 4.2-4.3 | analytic proof + exact check of 4 epochs |
| C4 | Any connected non-bipartite pair-sum graph with A = aI diverges when a^2 >= Perron root of the schedule matrix (2+sqrt5 for the triangle); extra attenuation or deletion of cross-covariances does not restore boundedness | Theorem 5.1, Remark 5.2, Corollary 5.3 | analytic proof |
| C5 | For a^2 in (1, 1.38] large priors diverge at the open-loop rate a^{2k} once in a certified region | Theorem 6.2 | computer-assisted (certificate C1, interval arithmetic) |
| C6 | Entry from P_0 = 10^12 I into that region at a^2 in {1+10^-8, 1+10^-6, 1.001, 1.01, 21/20, 1.1, 1.2, 1.38} | Corollary 6.3 | computer-assisted (certificate C2) |
| C7 | At a^2 = 21/20, P_0 = I converges to a bounded fixed point while 10^12 I diverges: no globally attracting steady covariance | Theorem 6.4, Corollary 6.5 | computer-assisted (certificate C3) + proof |
| C8 | Actual estimation error of the IsoKF estimate stays bounded in the divergent cases tested (gross pessimism, up to ~10^190 at a^2 = 9) | Finding 7.1 | numerical (80 digits), not a theorem |
| C9 | Divergence threshold from P_0 = I near a^2 = 1.1225; deletion variant threshold close to 2+sqrt5 | Finding 7.2 | numerical (double precision) |
| C10 | Every local block stays positive definite for every schedule iff participant groups overlap in at most one node (|U cap V| <= 1) | Theorem 8.1 | analytic proof + exact witness |
| C11 | Anchored relative (difference) measurements with a directly observed anchor give bounded marginals | Proposition 9.1 | analytic proof + exact identity check |
| C12 | Not claimed: any failure of the ikf_lib software (inspected, not run); anything for a^2 in (1.38, 2+sqrt5); uniqueness for stable A; boundedness of the actual error as a theorem | Sections 9-11, README Boundaries | statement of scope |

## Verification artefacts

- Exact rational (standard library): [scripts/verify_pairwise_divergence.py](scripts/verify_pairwise_divergence.py) (C3, C4 schedule matrix), [scripts/verify_exact.py](scripts/verify_exact.py) (C10 witnesses, Appendix A), [scripts/verify_anchored_identity.py](scripts/verify_anchored_identity.py) (C11), [scripts/fingerprints.py](scripts/fingerprints.py).
- Interval certificates (mpmath 1.3.0 `iv`, outward rounding): [scripts/cert_x1_tail.py](scripts/cert_x1_tail.py) (C1), [scripts/cert_x1_entry.py](scripts/cert_x1_entry.py) (C2), [scripts/cert_x2_bounded.py](scripts/cert_x2_bounded.py) (C3); runner [scripts/run_ext.sh](scripts/run_ext.sh).
- Numerical: [scripts/x3_true_error.py](scripts/x3_true_error.py) (C8), [scripts/x4_thresholds.py](scripts/x4_thresholds.py) (C9).
- Figures: [scripts/plot_growth.py](scripts/plot_growth.py), [scripts/plot_coexistence.py](scripts/plot_coexistence.py) -> [figures/](figures/).
- Outputs in normal and `python -O` mode: [results/normal/](results/normal/), [results/O/](results/O/). Claim-by-claim ledger: [VERIFICATION.md](VERIFICATION.md).
- Every script requires `--out NEWDIR`, uses explicit checks (no `assert`) and exits non-zero on failure, also under `-O`. CI: [.github/workflows/replay.yml](.github/workflows/replay.yml).

## Sources, review and provenance

- [refs.bib](refs.bib), generated by [tools/make_refs.py](tools/make_refs.py) from Crossref and arXiv records in [tools/records/](tools/records/).
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal review objections and resolutions (no scores). [RESEARCH_GATES.json](RESEARCH_GATES.json): contribution map, sources read, limitation triage.
- Integrity: [MANIFEST.sha256](MANIFEST.sha256). Licences: [LICENSE.md](LICENSE.md) (text, figures, data CC0-1.0; code MIT). Citation: [CITATION.cff](CITATION.cff); archive metadata: [.zenodo.json](.zenodo.json).

## Limits

Unrefereed candidate: internal AI-assisted review only; no human peer review, independent reproduction or formal
verification. Certificates assume mpmath interval arithmetic is correctly rounded. Numerical findings are evidence, not
proofs. The counterexamples use pair-sum observation channels, not plain differences. Priority holds only within a bounded
search of 1-4 November 2026. Luft et al. (2018, IJRR), read in full for 0.2.0, introduce the correction (their eq. 26), its
scaled variant (eq. 28) and prove two-robot blocks stay positive definite (Proposition 1, credited at Theorem 8.1); they make
no covariance boundedness or convergence claim and give no divergence example. Version 0.2.0 changes only this comparison
(manuscript Sections 10-11); no theorem, proof, certificate or computation changed.
