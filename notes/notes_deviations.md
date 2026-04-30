# Code Deviations from Lab Methods and Differences from Paper Figures

---

## 1. Non-Lab Functions Used

### 1.1 `sandwich::NeweyWest()` — All LP Regression Files

**Files affected:** `Fig_2.R`, `Fig_3.R`, `Fig_4.R`, `Fig_5.R`, `Fig_6.R`, `Fig_7.R`, `Fig_10.R`, `Fig_11.R`, `Fig_13.R`, `Table_1.R`

The course labs do not cover Newey-West HAC standard errors. However, local projections on time-series data require heteroskedasticity- and autocorrelation-consistent (HAC) standard errors to produce valid inference. `NeweyWest()` from the `sandwich` package is the standard R implementation and cannot be replaced by any function introduced in the labs.

Note: the default behavior of `sandwich::NeweyWest()` applies a small-sample degrees-of-freedom correction (`adjust = TRUE`). EViews, which the original paper likely used for estimation, does not apply this correction by default. This causes a small difference in the width of confidence bands between our output and the paper.

### 1.2 `glmnet` (Elastic Net VAR) — Connectedness Computation Files

**Files affected:** `Compute_Ct.R`, `Compute_Ct_macro.R`

The Diebold-Yilmaz (2012) connectedness index requires estimating a high-dimensional VAR over 51 state-level housing return series. The paper follows Lee (2025) in using an Elastic Net penalty to regularize the VAR coefficients. This is implemented via `glmnet` with `alpha = 0.5`. There is no equivalent method in the course labs, and no substitution is possible without fundamentally changing the paper's methodology.

---

## 2. Differences Between Generated Figures and the Original Paper

### 2.1 Confidence Interval Width — NOT an issue ✅

**Earlier versions of this document incorrectly flagged this.** The paper explicitly states in the caption of every LP figure (Figures 2, 3, 4, 5, 6, 7, 10, 11, 13): *"68 percent confidence intervals."* Our code constructs bands as `lower = est - se` and `upper = est + se`, which is exactly ±1 standard error = 68% CI. This is correct and intentional.

### 2.2 P-values in Table 1 — Two compounding sources

**Root cause 1: Newey-West bandwidth formula not specified in paper.**

The paper states only: *"we apply the correction method proposed in Newey and West (1987) to our standard errors."* No bandwidth is specified. Our code uses:
```r
bw_fix <- floor(0.75 * T_smpl^(1/3))
```
For the baseline sample of T = 108 quarters, this gives bandwidth = 3. EViews's default automatic Newey-West bandwidth formula is `floor(4*(T/100)^(2/9))`, which gives bandwidth = 4. A different bandwidth directly changes the HAC covariance matrix, changing standard errors and therefore p-values across all LP regressions and the cumulative tests in Table 1.

**Root cause 2: `NeweyWest(adjust = TRUE)` vs. EViews default.**

R's `sandwich::NeweyWest()` applies a small-sample degrees-of-freedom correction (`adjust = TRUE`) by default, scaling the covariance matrix by n/(n-k). EViews does not apply this correction. This makes our standard errors slightly larger than the paper's at every horizon, compounding the bandwidth difference.

Both sources push in the same direction — our standard errors are slightly larger than the paper's — which explains why some p-values in Table 1 are larger in our replication than in the paper while cumulative IRF point estimates match closely.

**Partially fixable:** Setting `adjust = FALSE` and implementing EViews's bandwidth formula would bring p-values closer to the paper's values. The group chose not to do this because `adjust = TRUE` with bandwidth `floor(0.75*T^(1/3))` is the standard academic R implementation and is defensible on its own terms.

### 2.3 IRF Shape and Magnitude Differences (Figures 3, 4, 6, 7, 13)

**Root cause: cross-software Elastic Net differences → different Ct series → different regime dummies → different IRFs.**

This is a chain reaction with four links:

**Link 1 — Solver differences:** The paper was produced using MATLAB (the replication package is available at OPENICPSR E216361V3). R's `glmnet` and MATLAB's Elastic Net solver differ in their coordinate descent algorithms, convergence tolerances, and numerical precision. Even with identical data, these produce slightly different VAR coefficients at each rolling window.

**Link 2 — Random seed:** Our `connectedness_helpers.R` sets `set.seed(2)` before each `cv.glmnet` call to ensure reproducibility within R. The paper does not mention a random seed. The paper's MATLAB or EViews implementation likely used a different fold assignment for 10-fold cross-validation, producing different λ values and therefore different regularized coefficients.

**Link 3 — HP filter propagation:** A different Ct series produces a slightly different HP trend. Some quarters that were classified as "above trend" (high connectedness) in the paper may be classified as "below trend" (low connectedness) in our replication, or vice versa, changing the composition of the p1/n1 dummy variables.

**Link 4 — IRF propagation:** All LP coefficient estimates in Equations (2) and (3) are conditional on these regime dummies. Even a small misclassification of a few quarters changes the effective subsamples for High/Low estimation, shifting all IRF point estimates and bands.

The overall qualitative finding — monetary policy has stronger effects in high-connectedness states — is preserved across all replicated figures. Quantitative differences are most visible in Figures 4 and 7 (Employment and Real Wage, which have lower signal-to-noise ratios) and Figure 13 (which uses the macro-augmented Ct from a 57-variable VAR, compounding the glmnet differences further). House Price and Residential Investment responses tend to match more closely because their signals are large enough to dominate any small misclassification.

**Not fixable** without exactly replicating MATLAB's Elastic Net solver behavior and the original cross-validation fold assignments.

### 2.4 Figures 6 and 7 — Romer-Romer Shock Vintage

**Root cause: shock series vintage not specified in paper.**

The paper states: *"we use the extended series for the Romer and Romer (2004) monetary shock."* Our replication uses `macro$RRSHOCK` from `All_data.xls`. The Romer-Romer shock series has been updated and extended multiple times since 2004. Without access to the original replication package (OPENICPSR E216361V3), we cannot confirm which vintage is in `All_data.xls` matches the one the authors used. If the series differ even slightly, all Figures 6 and 7 IRFs will differ independently of the Elastic Net issue.

**Not verifiable** without the original replication data.

### 2.5 Quarterly Aggregation of Monthly Ct — Method Not Specified in Paper

**Root cause: paper is ambiguous.**

Section 2.3 states: *"We then convert our monthly measures to quarterly ones."* No aggregation method is given. Our code (in both `Fig_1.R` and `Fig_13.R`) uses the mean of the 3 months within each quarter:
```r
ct_q <- sapply(1:T_q, function(j) mean(ct_raw[(3*(j-1)+1):(3*j)]))
```
This choice is internally consistent across all scripts. However, if the paper's authors used end-of-quarter values instead, the quarterly Ct series would differ slightly, which propagates into the HP filter cycle and therefore into the regime dummies, shifting IRF estimates for all LP figures.

**Not fixable** without clarification from the authors.

### 2.6 Visual Style Differences (All Figures)

The paper's figures are produced in MATLAB, which uses serif axis labels, specific marker sizes, and a white background with light gray gridlines. Our figures use `ggplot2`'s `theme_bw`, which produces a similar but not identical appearance (different fonts, slightly different gridline weight, different legend placement). The scientific content is equivalent; only the aesthetics differ.

---

## Summary Table

| Issue | Affected | Cause | In original notes? | Fixable? |
|---|---|---|---|---|
| 68% CI width | — | Not an issue — paper explicitly states 68% CI | ❌ Incorrect | N/A |
| P-values in Table 1 | Table 1, all LP | NW bandwidth mismatch (bw=3 vs EViews bw=4) + adjust=TRUE | ❌ Missing | Partially |
| P-values in Table 1 (secondary) | Table 1, all LP | `NeweyWest(adjust=TRUE)` vs EViews no correction | ✅ | Partially |
| IRF shape/magnitude | Fig 3, 4, 6, 7, 13 | glmnet vs MATLAB solver + set.seed(2) → different Ct → different dummies | Partial | No |
| RR shock vintage | Fig 6, 7 | Paper uses "extended series" — vintage unverifiable | ❌ Missing | No |
| Quarterly Ct aggregation | All LP figures | Paper says "convert to quarterly" without specifying mean vs end-of-quarter | ❌ Missing | No |
| Visual style | All | ggplot2 vs MATLAB aesthetics | ✅ | Partially |