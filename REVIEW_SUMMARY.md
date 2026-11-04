# Summary of the internal review and the revision

The 2 November 2026 version of this paper went through one internal review round on 4 November 2026. Five reviewer roles took part: an editor, a methodology reviewer, a domain reviewer, a cross-disciplinary reviewer and a devil's advocate. A separate adversarial check and an independent replay of all scripts were also run. All reviews were AI-assisted. This was **not** external peer review.

This file lists the reviewers' objections and how the revision answers each one. It contains no scores or grades.

## Overall verdict on correctness

No reviewer found an error in a stated result. All load-bearing algebra was independently recomputed in exact arithmetic: Lemma 1, identity (11), the bound table, the Perron data, the participant-family witnesses, the singularity witnesses and the anchored identity. The domain reviewer checked the recurrence against the published equations and confirmed that it is a faithful transcription.

## Objections and resolutions

| # | Objection | Resolution in this version |
|---|---|---|
| 1 | The refutation attacks a paraphrase. The conjecture and the update equations must be quoted, and the reading of the correction factor justified. | Conjecture 1 is quoted verbatim, together with the dissertation's extra clause. Table 2 sets the published eqs. (47)–(50) beside ours. The factor is the marginal ratio, applied with the marginal prior. The journal's own remark admitting unstable dynamics is quoted. An unconfirmed section locator was removed. |
| 2 | A = 3I is strongly unstable, and the proved threshold a² ≥ 4.24 is far from what simulations show (about 1.12). | New Section 6 gives computer-assisted proofs. Large priors diverge at the open-loop rate for every a² in (1, 1.38] once they reach a certified region. Entry into that region is certified down to a² = 1 + 10⁻⁸. Section 7 reports the numerical threshold table. |
| 3 | Divergence needs instability. Say so, and say what remains open for stable dynamics. | New Proposition 3.1: Schur-stable dynamics always give bounded marginals. The uniqueness question for stable A is stated as open. |
| 4 | The "unique" part of the conjecture was not addressed. | Theorem 6.4 and Corollary 6.5: at a² = 21/20 one prior provably converges and another provably diverges, so there is no globally attracting steady covariance. Figure 1 shows both. |
| 5 | Is the estimate itself harmed, or only the reported covariance? | Finding 7.1 computes the exact error covariance of the IsoKF estimate in 80-digit arithmetic. It stays bounded while the reported covariance diverges. The abstract and conclusions say plainly that the failure is gross pessimism, not estimate divergence, in these tests. They also note that the estimate is still worse than the centralised one. |
| 6 | Nothing IsoKF-specific: deleting cross-covariances also diverges. | Finding 7.2: the deletion variant stays bounded up to a threshold close to 2+√5 (finite-horizon bias unresolved), while the IsoKF rule diverges from about 1.12. The rule itself is less robust than deletion. |
| 7 | Pair-sum measurements are not the robotics model. | Section 9 explains that unanchored difference graphs are unobservable and that anchored differences are proved bounded. The risk lies in sum-type or mixed-sign channels. A physical example is given (delay biases in two-way ranging), stated without claiming that any deployed system diverges. |
| 8 | Missing related work: Schmidt/consider filter, covariance intersection, split CI, Roumeliotis–Bekey, double counting. | New Section 10, with references taken from Crossref records. |
| 9 | The buffered implementation might differ from the dense recurrence. | Section 2 gives a one-paragraph argument that the two are identical for linear time-invariant models. Section 11 adds why the library's positive-semidefinite projection would not mask the divergence. |
| 10 | Forward citations of the 2025 article were not checked. The 2018 journal article was not read in full. | Two related citing works were checked by title and abstract (D-GVIO, De Carli–Dimarogonas). The limitation for the 2018 article is kept. |
| 11 | Minor points: unit pair sums only, "interval" wording, order dependence, strictness after prediction, closedness in the limit argument, rank-deficient H_U, k ≥ D+1, deletion inside the attenuation corollary, labelling of the singularity section, the u_k device, figure caption, abstract length, title, author line. | All addressed. Remark 5.2 gives the rescaling, the open condition and all six orders. The singularity witnesses and the attenuation repair moved to appendices, and the witnesses are labelled as outside the algorithm's domain. The input device was dropped. The author line is now Anonymous. |
| 12 | Scripts: one verifier imported a sibling module and failed under `python -I`. `assert` statements vanish under `-O`. | Helpers inlined. All asserts replaced with explicit checks that exit with code 1. Every script needs a new `--out` directory. Outputs were regenerated in normal and `-O` mode and are identical. |

| 13 | Verification round (certificate audit, proof audit, re-review): the coexistence proof named a displayed box that is not invariant; some certificate boundaries were not rounded outward; the limit array was said to annihilate observed rows globally; the a = 1 remark and some numerical wording overstated the evidence. | Contraction argument now runs on the exact certificate box. Rounding fixed in C1–C3 and all certificates re-run (still pass). Statements narrowed: local annihilation only, strict hypothesis for (d), no global claim at a = 1, "no globally attracting steady covariance", divergent-case wording, deletion threshold "close to" 2+√5. Failures now emit a JSON record. |

## What remains open

- Certified entry for every a² in (1, 1.38].
- Any proof for a² in (1.38, 2+√5).
- A proof that the actual estimation error stays bounded.
- Uniqueness of the steady state for stable dynamics.
- Behaviour of the compiled software.
