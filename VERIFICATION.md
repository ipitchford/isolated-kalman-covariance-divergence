# Verification ledger

The authoritative statement is `manuscript.tex` and its build `isokf_manuscript.pdf` (4 November 2026). Numbering below refers to that PDF.

Evidence types:
- **Proof**: a written analytic argument. No proof is machine-verified.
- **Exact**: rational arithmetic (Python `fractions`), standard library only.
- **Certificate**: outward-rounded interval arithmetic (mpmath 1.3.0 `iv`). It checks the finite hypotheses of a written argument.
- **Numerical**: high-precision or double-precision evaluation, not certified.

Every script needs `--out NEWDIR` and refuses to overwrite. Checks are explicit `if`/exit tests, not `assert`, so they also run under `python -O`, and any failure exits with code 1. Each exact script records its own SHA-256 and the number of checks passed. Outputs are in `results/normal/` and `results/O/`. All JSON and CSV result files are byte-identical between the two modes. (The `-O` figure runs used `--no-figure`.) After the verification-round rounding fixes, C1–C3 were re-run and still pass, with the same entry epochs and bounds. Any unexpected exception is caught and reported as a JSON failure record with exit code 1.

## Claims, proofs and checks

| Claim | Where proved | Check (script → output) | Recorded result |
|---|---|---|---|
| Recurrence (2.1)–(2.3) equals published eqs. (47)–(50) | Sec. 2, Table 2 | Verbatim source excerpts from the internal domain review | Factor is the marginal ratio, not a square root |
| Prop. 3.1: Schur-stable A ⇒ bounded marginals | Proof (Sec. 3) | — | Proof only |
| Thm. 4.1: m_k ≥ (5/4)(9/5)^k − 1/4; central KF ≤ 3/4 | Lemmas 4.2–4.3 | `verify_pairwise_divergence.py` → `verify_pairwise_divergence/` | 4 exact epochs; cone, identities (4.4)–(4.5) and bound checked at every update (240 checks) |
| Comparison matrix M, ρ = 2+√5 | Sec. 5 | same script; `fingerprints.py` | char. poly λ(λ²−4λ−1) |
| Thm. 5.1, Rem. 5.2, Cor. 5.3 | Proof | — | Proof only |
| Lemma 6.1 (deficit form, cycle identity) | Proof | used by all certificates | — |
| Certificate C1: invariant boxes, a² ∈ [1, 1.38] | Sec. 6.2 | `cert_x1_tail.py` → `x1_tail/` | θ upper endpoint rounded up from 10⁻⁶; 190 intervals, all_ok; max L_m ≤ 3.1717; L₁⁻ ≥ 2.0104 |
| Thm. 6.2 (mild-instability divergence) | Proof + C1 | as above | — |
| Cor. 6.3 (entry from 10¹²I at 8 values of a²) | C2 | `cert_x1_entry.py` → `x1_entry/` | exact strict test min p > π_a (rational); all_ok; k₀ = 29, 29, 29, 29, 31, 32, 35, 43 |
| Thm. 6.4 (P₀ = I converges at a² = 21/20) | Proof + C3 on the exact box B_cert (endpoints recorded in binary in `x2_bounded.json`); the displayed box only reports bounds | `cert_x2_bounded.py` → `x2/` | entry epoch 80; ‖D(F⁴)‖∞ ≤ 0.5955 (interval row sums, upper endpoints); all_ok |
| Cor. 6.5 (no unique attracting steady state) | Thm. 6.2–6.4 + classical Riccati theory | — | — |
| Finding 7.1 (actual error bounded) | Numerical (80 digits, exact linear recursion) | `x3_true_error.py` → `x3/` | table in Sec. 7; cross-checked by the independent dense Decimal recursion in `plot_coexistence.py` (agree to all printed digits) |
| Finding 7.2 (thresholds) | Numerical (double precision) | `x4_thresholds.py` → `x4/` | IsoKF from I: a² ∈ [1.12251, 1.12256]; deletion ≈ 4.234 |
| Thm. 8.1 necessity witness, (zero innovation) | Proof + exact | `verify_exact.py` → `verify_exact/` | det P₁ = −7/162; hᵀP₁h = −1/3 |
| Prop. 9.1 identity | Proof + exact | `verify_anchored_identity.py` → `verify_anchored_identity/` | exact equality with non-commuting 2×2 blocks |
| App. A witnesses | Exact | `verify_exact.py` | det P_* = −257375/270336; interval 155/132 < R₃ < 2265/1804 |
| App. B criterion | Proof | — | Proof only |
| Fig. 1 | Illustration | `plot_coexistence.py` → `plot_coexistence/` | P₀=I ends at certified fixed point; 2706 checks |
| Fig. 2 | Illustration | `plot_growth.py` → `plot_growth/` | bound and central ≤ 3/4 checked at each epoch |

## Fault-injection check

A copy of `verify_exact.py` with one expected determinant changed (−7/162 → 7/162) was run under `python -O`. It printed `{"status": "FAILED", "failed_check": "verify_exact:111"}` and exited with code 1. This confirms that the checks are not removed by `-O`.

## Build

`pdflatex` + `bibtex` (TinyTeX): 18 pages, no undefined references or citations, no overfull boxes. Every page was inspected visually after rendering.

## Limits of this ledger

- The fidelity row relies on the domain reviewer's verbatim excerpts, which were read from the open-access article on 4 November 2026. A direct automated re-fetch during revision was blocked by the publisher site, so the excerpts were not re-read a second time.
- Certificates C1–C3 assume that mpmath's interval arithmetic is correctly rounded.
