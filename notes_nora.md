# Progress Notes — Nora Zhou

This file documents the code review, fixes, and extension work I completed on
top of Yitong's consolidated codebase. It also contains a handoff for Cris
(Person B of the L2 directional extension).

---

## 1. Code Review & Fixes to Yitong's Code

### 1a. Refactor `Compute_Ct.R` + `Compute_Ct_macro.R`

The two scripts shared roughly 110 lines of identical helper-function
definitions (`LVAR`, `Bcoef`, `Acoef`, `Phi`, `res_fn`, `fevd_generalised`,
`sptable2012`). I extracted them into a single shared file and slimmed the
two main scripts:

- **New file:** `connectedness_helpers.R` (222 lines, all helpers).
- `Compute_Ct.R`: 262 → 90 lines; now starts with `source("connectedness_helpers.R")`.
- `Compute_Ct_macro.R`: 180 → 96 lines; same refactor.
- Total: 442 → 408 lines, with **zero duplicated logic**.

**Bit-for-bit equivalence verified** by running the first 5 rolling windows
through the refactored code and comparing to the pre-refactor `Ct_monthly.csv`:
all 5 values matched to machine precision.

### 1b. LP regime intercepts (`p1`, `n1`) added to 10 files

When auditing `Table_1.R` against the original Lee & Ma (2025) EViews program
`table_1.prg`, I found that the R port omitted `p1` and `n1` as bare
regressors in the state-dependent specification. The original EViews
specification includes them (since `p1 + n1 = 1`, they act as
regime-specific intercepts). Without them, both regimes are forced to
intercept = 0, mis-attributing level differences to slope coefficients.

**Fix applied to 10 files** (`Fig_2.R`, `Fig_3.R`, `Fig_4.R`, `Fig_5.R`,
`Fig_6.R`, `Fig_7.R`, `Fig_10.R`, `Fig_11.R`, `Fig_13.R`, `Table_1.R`):

```r
# before
sd_rhs <- c("0", "p1_FFR/SHOCK", "n1_FFR/SHOCK", ...)

# after
sd_rhs <- c("0", "p1", "n1", "p1_FFR/SHOCK", "n1_FFR/SHOCK", ...)
```

**Result.** Table 1 now matches the paper to four decimal places at the
boundary horizons (h=4 and h=20):

| Cell | Paper | Replication (after fix) |
|---|---:|---:|
| GDP - High - h=4   | 1.34   | 1.3383  |
| GDP - High - h=20  | 15.50  | 15.5004 |
| FFR - High - h=4   | -5.55  | -5.5547 |

Intermediate horizons (h=8, 12, 16) deviate by 5–12% from the paper, which
is consistent with the documented R-glmnet vs MATLAB Elastic Net solver
difference (`notes_deviations.md` Section 2.1) propagating through the
HP-cycle high/low classification.

**Visual results.** After re-running, the figures whose IRF patterns most
clearly match the paper are `Fig_2`, `Fig_3`, `Fig_5`, `Fig_6`, `Fig_10`,
`Fig_11`, `Fig_12`, and `Table_1` (essentially identical at the boundary
horizons). `Fig_4`, `Fig_7`, and `Fig_13` are still close but show
larger horizon-by-horizon deviation; these are worth a second look before
final submission.

### 1c. `Fig_12.R` x-axis alignment with paper

`notes_deviations.md` Section 2.3 documented that `Fig_12.R` plotted the
full 1981–2016 connectedness series, while the paper's Figure 12 only
displays 1981–2007. I changed:

```r
# before
limits = range(df$date)

# after
limits = c(as.Date("1981-01-01"), as.Date("2007-12-31"))
```

`notes_deviations.md` Section 2.3 has been updated to "RESOLVED" and the
summary table updated accordingly. The underlying connectedness data is
unchanged; only the displayed range is clipped.

---

## 2. Extension — Part 1 (Person A): L2 Directional Connectedness

### Background

Diebold & Yılmaz (2014) define connectedness at three levels:

| Level | Object | Cardinality |
|---|---|---|
| L1 | Pairwise directional `C_{i←j}` | 51 × 51 = 2,550 measures |
| L2 | Total directional `C_{i←•}`, `C_{•←j}` | 51 To + 51 From + 51 Net |
| L3 | System-wide total `C_t` | 1 scalar |

Lee & Ma (2025) only use L3 in their main analysis (`Ct_monthly.csv`).
Following the TA's suggestion, my extension extracts and saves L2 measures
so that downstream analysis can ask which states act as persistent net
"transmitters" vs "receivers" of housing-market shocks, and whether those
roles change between high- and low-connectedness regimes.

### What I built

- **New folder:** `extension/`
- **New script:** `extension/compute_directional.R` (231 lines). Re-runs
  the same rolling Elastic-Net VAR loop as `Compute_Ct.R`, but at every
  window also stores the To, From, and Net vectors from `sptable2012()`'s
  output instead of discarding them.
- Run time: ~20 minutes (single-threaded, same as `Compute_Ct.R`).

### Sanity checks (both passed)

1. **L3 consistency.** `matspill` produced by `compute_directional.R`
   was compared to the existing `Ct_monthly.csv` (max absolute difference
   reported in console). Result: bit-for-bit match.
2. **Network-property sanity.** For every rolling window, the entries of
   `sum(Net)` should be approximately zero (because total To = total From
   by construction). Console output: `max |sum_t(Net)| = 6e-3`, well
   within the 1e-3 rounding tolerance from `sptable2012()`'s
   `round(.., 3)` step.

### Outputs (in `extension/outputs/`)

| File | Shape | Contents |
|---|---|---|
| `directional_To_monthly.csv`   | 432 × 52  | date + 51 columns of To values |
| `directional_From_monthly.csv` | 432 × 52  | date + 51 columns of From values |
| `directional_Net_monthly.csv`  | 432 × 52  | date + 51 columns of Net values |
| `directional_long.csv`         | 22,032 × 5| long format: date, state, To, From, Net |
| `state_ranking_by_Net.csv`     | 51 × 6    | per-state mean_To, mean_From, mean_Net, sd_Net, pct_pos_Net (sorted by mean_Net descending) |

### Top-line empirical result

| | Top 5 states (mean_Net) |
|---|---|
| Net **transmitters** (highest mean_Net)   | MI, MN, CA, NH, CO |
| Net **receivers**    (lowest  mean_Net)   | KS, KY, WV, ND, IL |

The story is geographically intuitive: coastal/tech/financial states
transmit; rural/agricultural Midwest and Appalachian states receive. The
one anomaly worth Cris's attention is **MI** (Michigan), whose `sd_Net`
(94.2) is roughly twice every other state's, and whose `pct_pos_Net`
(48.8%) is essentially 50%. MI's "average transmitter" status is being
pulled by extreme periods, plausibly the 2008 housing collapse and Detroit
auto-industry shock. Recommend Cris plot MI's Net time series separately
and split-sample around 2007–2009.

---

## 3. Handoff to Cris (Person B)

### Your starting point

The five CSVs in `extension/outputs/` are your raw material. The most
useful three for plotting are:

- **`state_ranking_by_Net.csv`** — feeds the U.S. choropleth map.
- **`directional_Net_monthly.csv`** — wide format, ideal for one-line-per-state
  time-series plots.
- **`directional_long.csv`** — tidy/long format, ideal for `ggplot2 + facet_wrap`.

Read them like this:

```r
library(dplyr); library(ggplot2); library(readr)
ranking <- read_csv("outputs/state_ranking_by_Net.csv")
net_wide <- read_csv("outputs/directional_Net_monthly.csv")
net_long <- read_csv("outputs/directional_long.csv")
```

### Suggested deliverables

**Figure E1 — U.S. choropleth of mean_Net.**
Color states by `mean_Net` (diverging palette: red = receiver, blue =
transmitter). This is the single most important visual for the extension;
it gives the audience an immediate geographic story. Use `usmap` or
`maps` + `ggplot2`. Tag the top 5 transmitter and bottom 5 receiver
state names directly on the map.

**Figure E2 — Time series of Net for selected states.**
Pick 6–8 representative states from across the spectrum (e.g., 3 top
transmitters: CA, MN, CO; 3 bottom receivers: KS, KY, WV; plus MI as
the volatility anomaly). Plot `Net_t` for each on the same axes,
shaded NBER recessions, vertical line at 2008. The story is whether
roles are persistent or whether they flip during crises.

**Figure E3 (optional but high-value) — Stability across regimes.**
Compute mean_Net **separately** for "high connectedness" quarters
(`Ct_dummies.csv`'s `pshock == 1`) vs "low connectedness" quarters.
Plot a 2-panel choropleth: left = mean_Net during high-Ct periods,
right = mean_Net during low-Ct periods. Test the TA's question: "do
the high connectedness periods correspond to a more centralized network,
specific hubs, or broader diffusion?"

**Table E1 — Top/bottom states with stability statistics.**
Use the columns of `state_ranking_by_Net.csv` directly, plus
correlation between mean_Net measured on the first half (1981–1994)
vs the second half (1995–2007) of the sample, as a stability index.

### Specific things to flag in the writeup

1. **MI's high `sd_Net`.** Almost certainly reflects the 2008
   housing-bust + Detroit shock. Worth its own paragraph.
2. **IL as a net receiver despite Chicago's size.** Counterintuitive
   but consistent with IL's overall population/housing-market decline
   relative to peer Midwest states.
3. **The 50-th percentile (median state) does not have mean_Net = 0**
   — the distribution is asymmetric. Worth noting.

### Questions worth asking the TA before final writeup

- Should `Figure E3` use the same `Ct_dummies.csv` lagged dummies (`p1`)
  as the LP regressions, or the contemporaneous `pshock`? My instinct is
  contemporaneous, since we are *describing* the network during the
  regime, not running an LP.
- Should the maps use 4 Census regions, 9 Census divisions, or all 51
  state-level mean_Net? Probably state-level for the map but state →
  regional aggregation for any quantitative comparison.

### What I will not do (it's your call)

- I have not yet computed any L1 (pairwise) measures. The pipeline could
  be extended to save the full 51×51 GFEVD matrix at every window, but I
  did not do this — it would more than double the output size and is not
  needed for the L2-focused extension we agreed with the TA on.
- I have not built any plotting code in `extension/figures/`. The folder
  exists but is empty; please save your `Fig_E1.png`, `Fig_E2.png`, etc.
  there.
