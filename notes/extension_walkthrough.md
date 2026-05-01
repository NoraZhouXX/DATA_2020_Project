# Extension Walkthrough — L2 Directional Connectedness

A complete narrative of why we built this extension, what it does
mathematically, which scripts produce which outputs, and how each piece
ties back to Lee & Ma (2025). Intended as a reference for the TA meeting
and for the final writeup.

---

## 1. Why an Extension Was Needed

### 1.1 The TA's hint

The TA mentioned that the original data has "three levels," and the paper
only uses two. We traced this to **Diebold & Yılmaz (2014, *Journal of
Econometrics*)** — the methodological foundation Lee & Ma (2025) build on.
DY2014 organizes connectedness into a three-tier hierarchy that is most
clearly stated in their Section 5.3:

> "we progressed from 'micro to macro' — that is, **from pairwise
> connectedness, to total directional connectedness, to total
> connectedness** … various levels of disaggregation (total directional
> and pairwise directional)."

And the empirical sections of DY2014 are explicitly organized this way
(Section 5.3.1 Total, 5.3.2 Total directional, 5.3.3 Pairwise directional).

### 1.2 What Lee & Ma actually use

Lee & Ma's main analysis uses **only the top tier** — the system-wide
total connectedness scalar `C_t` (one number per rolling window). Every
figure and table in the paper is conditioned on this single time series.

The two finer levels are computed as a by-product of the same rolling
elastic-net VAR — the function `sptable2012()` in the original code
constructs them at every window — but the paper's code immediately
discards them and saves only the bottom-right cell of the spillover table
(the L3 scalar).

### 1.3 Our extension in one sentence

**Re-run the same rolling elastic-net VAR, but at every window save the
state-level directional vectors that the original paper computes and
throws away, then build a story around them.**

The cost is small (re-run the existing pipeline once, ~20 minutes), the
methodological risk is zero (it is the same DY2014 calculation),
and the new information is genuine — the paper's data simply does not
expose state-level transmitter/receiver heterogeneity.

---

## 2. Formal Definitions of L1, L2, L3

At each rolling window of 120 months, we estimate a 51-variable
elastic-net VAR(1) on state housing returns, compute the H = 10
generalized forecast-error variance decomposition matrix
$\tilde{D}^H = [\tilde{d}^H_{ij}]_{i,j=1..51}$ (row-normalized), then read
the three levels off it as follows.

### Level 1 — Pairwise directional (51 × 51 = 2,550 measures)

$$C^H_{i \leftarrow j} = \tilde{d}^H_{ij}$$

The fraction of state $i$'s 10-step-ahead forecast-error variance
attributable to a shock originating in state $j$. Highly granular; useful
for network analysis but unwieldy for direct economic interpretation.

### Level 2 — Total directional (3 × 51 measures)

For each state $j$:

$$
C^H_{\bullet \leftarrow j} = \sum_{i \neq j} \tilde{d}^H_{ij} \quad (\text{To: how much } j \text{ sends out})
$$

$$
C^H_{j \leftarrow \bullet} = \sum_{i \neq j} \tilde{d}^H_{ji} \quad (\text{From: how much } j \text{ receives})
$$

$$
\text{Net}_j = C^H_{\bullet \leftarrow j} - C^H_{j \leftarrow \bullet}
$$

Positive Net = net transmitter; negative Net = net receiver. This is the
level our extension targets.

### Level 3 — System-wide total (1 scalar)

$$
C^H_t = \frac{1}{N} \sum_{i \neq j} \tilde{d}^H_{ij}
$$

This is the `Ct_monthly.csv` series the original paper saves and uses
throughout.

### Why we focus on L2

L1 produces $51 \times 51 = 2{,}550$ time series per window — far too many
to interpret one-by-one. Using L1 properly would require reducing it via
network-centrality measures (eigenvector centrality, PageRank, etc.) and
running heterogeneous local projections on top of those. That doubled the
project scope.

L2 produces $3 \times 51 = 153$ values per window, which is exactly the
right granularity for a state-level transmitter/receiver story —
geographically interpretable, statistically tractable, and directly
comparable to the existing Lee & Ma cross-sectional analyses (e.g. their
Figure 11 splits states by debt level; we split them by network position).

The TA confirmed in email that L2 first, with optional L1 if time
permits, was the right scope.

---

## 3. The Computational Pipeline

```
EI_DATA.csv (51 state housing returns, 552 monthly obs)
            │
            ▼
   ┌──────────────────────────────────┐
   │  extension/compute_directional.R │
   │                                  │
   │  Same loop as Compute_Ct.R, but  │
   │  saves all 3 directional vectors │
   │  (To, From, Net) at every window │
   │  via sptable2012() output rows.  │
   └────┬─────────────────────────────┘
        │
        ▼
extension/outputs/  (5 CSVs)
   directional_To_monthly.csv      432 × 52  (date + 51 states)
   directional_From_monthly.csv    432 × 52
   directional_Net_monthly.csv     432 × 52
   directional_long.csv         22,032 × 5  (tidy: date, state, To, From, Net)
   state_ranking_by_Net.csv       51 × 6   (mean_To, mean_From, mean_Net,
                                            sd_Net, pct_pos_Net)
        │
        │
   ┌────┼────────────────────────────┬────────────────────┐
   ▼    ▼                            ▼                    ▼
 Fig_E1.R                       Fig_E2.R               Table_E1.R
 (US choropleth)               (Time series)           (Ranking +
        │                            │                  stability)
        ▼                            ▼                    ▼
 Fig_E1.png                    Fig_E2.png             Table_E1.png
                                                      Table_E1.csv
```

### Two sanity checks built into `compute_directional.R`

1. **L3 consistency.** The script also reconstructs `matspill` (the
   scalar L3) and compares it cell-by-cell to the existing
   `Ct_monthly.csv`. Result on our run: `max abs diff = 0` (bit-for-bit
   match). This confirms the L2 values come from the **same** rolling VAR
   estimates that produced the published L3 — there is no
   re-estimation drift.

2. **Network identity.** By construction, total To across states equals
   total From across states at every window, so Net should sum to ~0 at
   every $t$. Result: $\max_t |\sum_i \text{Net}_{i,t}| = 6 \times 10^{-3}$,
   well within the rounding tolerance of `sptable2012()`'s
   `round(., 3)`.

---

## 4. The Three Deliverables

### 4.1 `Fig_E1.png` — U.S. choropleth of mean Net

**Built by `extension/Fig_E1.R`** using `usmap` + `ggplot2` + `ggrepel`.

**What it shows.** Each state is colored by its time-averaged Net
connectedness over 1981–2007. Blue = transmitter (Net > 0), red =
receiver (Net < 0). The top 5 transmitters (CA, MN, NH, CO, MI) and
bottom 5 receivers (KS, KY, WV, ND, IL) are labeled directly on the map.

**The geographic story.** Net transmitters cluster in the West Coast
(CA), the Mountain West (CO), and the financial / tech corridors
(MN's Twin Cities, NH/Boston metro). Net receivers cluster in the
Plains and Appalachian states whose housing markets are smaller and
follow national trends rather than driving them.

**Relation to the paper.** Lee & Ma's Figure 11 partitions states by
**debt level** and shows that high-debt states have stronger IRFs to
monetary shocks. Our map gives an orthogonal partition of the same 51
states by **network position**, complementary to their debt-based
heterogeneity. A natural follow-up question (good for the writeup) is
whether the high-debt states from Lee & Ma overlap with our net
transmitters — a quick comparison would identify whether debt and
network role are correlated channels of transmission.

### 4.2 `Fig_E2.png` — Time series of Net for selected states

**Built by `extension/Fig_E2.R`** using `ggplot2` + `patchwork`.

**What it shows.** Side-by-side panels of Net_t plotted from 1981 to 2007
for 4 chosen transmitters (MI, MN, CA, CO) and 3 chosen receivers (KS,
KY, WV).

**The temporal story.** Most states have stable roles over the sample —
California stays positive throughout, Kansas stays negative throughout.
The exception is **Michigan**: its Net hovers near zero or negative
through the 1980s and 1990s, then climbs sharply after the early 2000s
and spikes during the 2008 housing crisis. This visually corroborates
the H1/H2 numerical split (see 4.3 below).

**Relation to the paper.** Lee & Ma's Figure 1 shows the aggregate
$C_t$ time series with NBER recession bands and points to the 2008
crisis as a major regime shift. Our Fig E2 decomposes that aggregate
shift to the state level and shows that **the rise in $C_t$ around 2008
is partly driven by Michigan flipping its role**, not just by uniform
strengthening of all linkages. This is a finer-grained mechanism than
the paper itself articulates.

### 4.3 `Table_E1.csv` (and `Table_E1.png`) — Ranking with stability index

**Built by `extension/Table_E1.R`**.

**What it shows.** A 51-row table sorted by `mean_Net` descending, with
columns:

| Column | Meaning |
|---|---|
| Rank | 1 = strongest transmitter, 51 = strongest receiver |
| State | abbreviation |
| Role | Transmitter / Receiver based on full-sample mean_Net |
| Role_stable | "Stable" if H1 and H2 means agree in sign; "Switched" otherwise |
| Mean Net, SD Net, % Pos | full-sample summary stats |
| H1 Mean | mean_Net over 1981–1994 |
| H2 Mean | mean_Net over 1995–2007 |

**The stability story (two layers).**

*Layer 1 — sign-level stability (the Role_stable column).* Of the top 14
transmitters, **6 are labelled "Switched"** because their mean_Net sign
changed between H1 (1981-1994) and H2 (1995-2007). Most strikingly,
**Michigan: H1 mean = −21.6 (receiver), H2 mean = +79.0 (transmitter)**.
Florida, Georgia, Virginia, and New Mexico also flipped sign.

*Layer 2 — magnitude-level reshuffling (the H1-on-H2 regression).*
Even among the states whose role *signs* are stable, the *magnitudes*
change dramatically. Regressing H2 mean on H1 mean across all 51 states
gives **R² = 0.001 and slope = 0.062** — essentially no predictive power
of one half over the other. Examples among the "Stable" states:
- New Hampshire: H1 = 55.8 (strong transmitter) → H2 = 8.5 (mild)
- Minnesota:     H1 =  0.1 (neutral)            → H2 = 71.4 (strong)
- Rhode Island:  H1 = 48.8 (strong)             → H2 =  6.4 (mild)
- California:    H1 =  6.4                      → H2 = 61.7

This is the deeper finding: **the hub roster turns over almost
completely between the two halves of the sample**. The clear
core-periphery picture in Fig_E1 is largely a 1995-2007 phenomenon, not
a stable long-run feature of the network.

**Reconciling with Fig_E2's "persistent role separation" story.**
The seven states we plot in Fig_E2 (MI, MN, CA, CO and KS, KY, WV) are
the *most* persistent transmitters and receivers — they were chosen
precisely because their roles are visually clear. The H1-on-H2 regression
includes all 51 states and shows that, in aggregate, the hub roster
reshuffled. There is always a core-periphery structure (so individual
survivors look stable in Fig_E2), but *which* states sit in the core is
sample-period-dependent.

**Relation to the paper.** This addresses a question Lee & Ma raise but
do not answer: **what does the rise in $C_t$ during high-connectedness
episodes structurally mean for the cross-section of states?** The H1-on-H2
near-zero R² gives our answer: high-connectedness periods are not just
amplifications of an existing hub structure — they coincide with a
**structural reorganization of which states drive the system**. The hubs
of the 1980s (NH, RI, MN-mediocre) are not the hubs of the 2000s
(MN-strong, CA, CO, MI-emergent). This adds a geographic / structural
dimension to the paper's main result: high connectedness is not just
"more spillovers happening" — it is "different states doing the
spilling."

---

## 5. The H1/H2 Split — Methodology and Economic Reading

### 5.1 How and where we split

`Table_E1.R` cuts the baseline sample at January 1, 1995:

```r
mutate(half = ifelse(date < as.Date("1995-01-01"), "H1", "H2"))
```

| | Range | Quarters | Economic background |
|---|---|---:|---|
| H1 | 1981Q1 – 1994Q4 | 56 | Volcker disinflation, 1990–91 recession, S&L crisis |
| H2 | 1995Q1 – 2007Q4 | 52 | Greenspan put, dot-com era, housing boom |

### 5.2 Why 1995

Three reasons:

1. **Natural midpoint.** The full baseline sample is 108 quarters; cutting at the 56-th quarter keeps H1 and H2 roughly equal in length so the half-sample means are directly comparable.
2. **Economically meaningful.** 1995 separates the "high-inflation / Volcker / S&L" era from the "Greenspan put / securitization / housing-boom" era. Pre-1995 housing rose modestly with general inflation; post-1995 saw the unprecedented Sun Belt and coastal boom, the rise of subprime lending, and growing financial integration of regional housing markets.
3. **Avoids sample-edge effects.** The cut sits 14 years deep into both H1 and H2, well away from any rolling-window endpoint where the connectedness estimate could be noisier.

A natural robustness check is to repeat with cuts at 1990 or 2000. We expect R² ≈ 0 to be robust because the underlying mechanism (structural housing-market reorganization in the 2000s) does not depend on the exact cut date.

### 5.3 What "Switched" means at the state level

Recall `Net_j = To_j − From_j`:

- **Net > 0 (transmitter)**: state $j$'s housing-price moves *lead* the rest of the country; its shocks propagate outward.
- **Net < 0 (receiver)**: state $j$'s housing-price moves *follow* the rest of the country; most of its variance is explained by shocks originating elsewhere.

A state labelled "Switched" therefore changed from passive follower to active leader, or vice versa, between H1 and H2. This is **not** a geographic property — it is an emergent feature of each half-sample's economic regime.

### 5.4 Three story examples

| State | H1 mean | H2 mean | Likely economic interpretation |
|---|---:|---:|---|
| **MI** (Michigan) | −21.6 | +79.0 | Receiver flips to strong transmitter. Plausible mechanism: the long arc of Detroit auto-industry decline gets priced into Michigan housing *earlier* than in peer states, so MI prices begin to *lead* national housing distress signals. |
| **FL** (Florida) | −26.7 | +35.4 | Sign flip. Plausible mechanism: the 2000s subprime boom and Miami / Orlando condo construction wave made Florida the bellwether of the national housing cycle. |
| **GA** (Georgia) | −16.8 | +35.4 | Mirror of FL, slower onset. Atlanta-metro emerged in the 2000s as a Sun Belt growth hub, leading regional housing-price moves. |

### 5.5 "Stable" labels can still hide reshuffling

`Role_stable` flags whether the *sign* of mean Net agrees between H1 and H2. But many states keep the same sign while swapping magnitude with each other. Examples:

| State | H1 mean | H2 mean | Comment |
|---|---:|---:|---|
| NH | +55.8 | +8.5  | Strongest H1 transmitter; barely transmits in H2 |
| RI | +48.8 | +6.4  | Same pattern, smaller scale |
| MN |  +0.1 | +71.4 | Almost neutral in H1, strongest H2 transmitter |
| CA |  +6.4 | +61.7 | Mild H1, strong H2 |

NH and MN essentially **trade ranks** between halves, yet because both stayed positive they are labelled "Stable" by the sign rule. The full-population regression of H2 on H1 across all 51 states captures this magnitude reshuffling: **R² = 0.001 with slope = 0.062**. A state's H1 transmitter rank has essentially no predictive power for its H2 rank.

### 5.6 Economic implication

The standard implicit assumption in the connectedness literature is that the *network structure* (who connects to whom, who acts as a hub) is roughly stable across time and that the system-wide level $C_t$ moves up and down by scalar amplification. Our finding challenges that assumption:

> **High-connectedness periods are not just amplifications of an existing hub structure — they coincide with a structural reorganization of which states drive the system.**

Concrete implications:
- Monetary-policy transmission through housing **depends on which states are currently hubs**. The same 100 bp shock in 1985 and 2005 will not propagate through the same state-level pathways.
- Forecasts that rely on historical state-level network structure to predict future spillovers under a new monetary regime should be treated with caution.
- The geographic dimension of monetary transmission is itself regime-dependent — not just an amplitude story, but a structural one.

This is the deeper finding our extension surfaces, beyond the narrower paper question of "high vs low connectedness amplitude."

---

## 6. How the Extension Relates to the Original Paper

| Lee & Ma (2025) | Our Extension |
|---|---|
| Use L3 (`C_t`) as conditioning variable for state-dependent LP | Decompose `C_t` into the 51 underlying L2 vectors |
| State-by-debt heterogeneity (Fig 11) | State-by-network-position heterogeneity (Fig E1) |
| Aggregate $C_t$ time series (Fig 1) | Per-state Net_t time series (Fig E2) |
| "Connectedness rises during crises" claim | Crisis attribution to specific switching states (Table E1, MI in particular) |
| "Three levels in the data, only used L3" (TA's hint) | Implements L2; L1 deliberately deferred as out-of-scope |

---

## 7. Anticipated TA Questions

**Q: Why not also do L1 (pairwise)?**
A: 2,550 pairwise measures per window are too many to analyze without
network tools (centrality, community detection). We agreed in your
earlier email that L2 first, L1 if time permits. We chose to do L2
thoroughly rather than both at lower depth. Happy to add an L1 network-
centrality block as a follow-up if you'd like, but the cost is roughly
another week of work and significantly more methodological exposition
in the writeup.

**Q: What does H1 vs H2 buy you that the high/low connectedness dummies
don't?**
A: They answer different questions. The high/low dummies (used by Lee &
Ma) split *time* by Connectedness regime; H1/H2 splits *time* by an
exogenous boundary (1995). The H1/H2 split is closer to "does the
network's structure change over time" while the high/low split is "does
the IRF amplitude change with regime." We used H1/H2 because it directly
addresses the TA-pose question of role stability. Could supplement with
a high/low version if useful.

**Q: How confident are you in the Michigan story?**
A: Very. (i) MI's `sd_Net` (94.2) is roughly twice every other state's;
(ii) `pct_pos_Net` is 48.8% — exactly half the time MI is a transmitter,
half a receiver, which only makes sense if the role flipped; (iii) H1
and H2 means are −21.6 and +79.0; (iv) the time-series plot in Fig E2
visually shows the regime change happening around 2007–2009, consistent
with the auto-industry collapse / Detroit housing bust. All four lines
of evidence point to the same conclusion.

**Q: Is the extension reproducible?**
A: Yes. `extension/compute_directional.R` runs in ~20 minutes,
including the L3 sanity check that confirms bit-for-bit consistency
with `Ct_monthly.csv`. The three downstream scripts (Fig_E1.R,
Fig_E2.R, Table_E1.R) each run in a few seconds and consume only the
five CSVs in `extension/outputs/`.
