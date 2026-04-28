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

### 2.1 IRF Shape Differences (Figures 2–7, 10, 13)

**Root cause: different Elastic Net implementation across software.**

The paper was produced using MATLAB (figures show MATLAB-style formatting) with the Elastic Net VAR likely estimated in MATLAB or EViews. R's `glmnet` package and MATLAB's Elastic Net solver differ in their coordinate descent algorithms, convergence tolerances, and cross-validation routines. Even with identical data, these differences produce a slightly different connectedness index (Ct), which in turn shifts the HP-filter cycle and changes the classification of quarters into High/Low connectedness regimes (the p1/n1 dummies).

Because the entire identification of state-dependent responses depends on these dummies, small changes in regime classification propagate into all IRF estimates. The overall patterns (e.g., High connectedness producing larger GDP responses) are preserved, but the exact point estimates and the ordering of High vs. Low curves at specific horizons can differ from the paper, particularly for Employment and Real Wage (Figures 4 and 7).

House Price responses tend to match the paper more closely because the signal is large enough to dominate any small misclassification.

### 2.2 Confidence Band Width (All LP Figures)

**Root cause: `NeweyWest(adjust = TRUE)` vs. EViews default.**

The `sandwich` package applies a small-sample correction to the Newey-West covariance matrix by default. EViews does not. This makes our confidence bands slightly wider than those in the paper, most visibly at short horizons where the sample used in each regression is smaller.

### 2.3 Figure 12 — X-Axis Range

**Root cause: `Ct_monthly.csv` spans 1981–2016, but the paper's Figure 12 only displays 1981–2007.**

`Compute_Ct.R` computes a rolling connectedness index over the full data range available in `EI_DATA.csv` (1981–2016). The current `Fig_12.R` plots the entire series (`limits = range(df$date)`), making the x-axis roughly 10 years longer than in the paper. The fix is to clip the x-axis upper limit to December 2007.

### 2.4 Visual Style

The paper's figures are produced in MATLAB, which uses serif axis labels, specific marker sizes, and a white-background style with light gray gridlines. Our figures use `ggplot2`'s `theme_bw`, which produces a similar but not identical appearance (different fonts, slightly different gridline weight, different legend placement). The scientific content is equivalent; only the aesthetics differ.

---

## Summary Table

| Issue | Affected Files | Cause | Fixable? |
|---|---|---|---|
| `NeweyWest()` not in labs | All LP figures + Table_1 | No lab equivalent for HAC SE | N/A — necessary |
| `glmnet` not in labs | `Compute_Ct.R`, `Compute_Ct_macro.R` | No lab equivalent for high-dim VAR | N/A — necessary |
| IRF shape differences | Fig_2 – Fig_7, Fig_10, Fig_13 | R vs. MATLAB Elastic Net → different Ct → different dummies | Only by replicating MATLAB solver exactly |
| Confidence band width | All LP figures | `NeweyWest(adjust=TRUE)` vs. EViews | Yes: use `adjust=FALSE` |
| Fig_12 x-axis too long | `Fig_12.R` | Ct series extends to 2016; no date clipping | Yes: set `xlim` to 2007-12-31 |
| Visual style | All figures | ggplot2 vs. MATLAB aesthetics | Partially, with theme customization |
