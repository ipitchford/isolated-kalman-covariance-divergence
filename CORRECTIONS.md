# Corrections

## 11 November 2026: version 0.3.0-candidate (doi:10.5281/zenodo.24142193)

Numerical finding 7.1 (Section 7) said that for a² = 9, P₀ = I the actual error variances converge to the limit (17/4, 13/4, 13/4), that the asymptotic error-transition matrix has spectral radius ≈ 0, and that the estimate is worse than the centralised one by a factor of about five to six. Replaying `scripts/x3_true_error.py` (80 digits) at consecutive epochs shows that the variances do not converge. After a short transient they repeat with period two: (17/4, 13/4, 13/4) after even epochs, which is the value recorded at k = 200, and (3/4, 21/4, 3/4) after odd epochs. The largest actual error variance therefore alternates between 4.25 and 5.25. The one-epoch transition matrix has spectral radius alternating between 0 and 1.5; the two-epoch product has spectral radius ≈ 0. Against the centralised filter (0.678 per node) the factor is about five to eight. The other tabulated cases are constant between consecutive epochs.

The actual error remains bounded in these tests, and the conclusion that the failure is pessimism of the reported covariance is unchanged. No theorem, proof, certificate, script or recorded output changed. Versions 0.2.0 (doi:10.5281/zenodo.23941896, tag `v0.2.0-candidate`) and 0.1.0 (doi:10.5281/zenodo.23935675, tag `v0.1.0-candidate`) remain archived unchanged.

## 4 November 2026: version 0.2.0-candidate

Full-text comparison with Luft et al. (2018, IJRR); see the README version history.
