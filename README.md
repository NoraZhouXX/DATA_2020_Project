# Housing Market Connectedness and Transmission of Monetary Policy

This repository contains a replication and extension study of **Lee and Ma (2025), “Housing market connectedness and transmission of monetary policy.”** The project reproduces the paper’s main connectedness measures and impulse-response results, then extends the analysis by examining **state-level directional connectedness** in the U.S. housing market network. :contentReference[oaicite:0]{index=0}

## Project Goals

This project has two main objectives:

1. **Replication**
   - Reproduce the paper’s housing-market connectedness measure and main empirical results.
   - Recreate the connectedness index and selected figures/tables from the published analysis. The repository includes scripts such as `Compute_Ct.R`, `Compute_Ct_macro.R`, `Fig_1.R` through `Fig_13.R`, and `Table_1.R` for this purpose. :contentReference[oaicite:1]{index=1}

2. **Extension**
   - Go beyond the paper’s aggregate connectedness measure by studying **state-level directional spillovers**.
   - Identify which states act as net transmitters or receivers of housing-market shocks using the extension pipeline in `extension/compute_directional.R`, along with `Fig_E1.R`, `Fig_E2.R`, and `Table_E1.R`. :contentReference[oaicite:2]{index=2}

## Repository Structure

```text
Housing-market-connectedness-and-transmission-of-monetary-policy/
├── data/               # Input datasets, paper, connectedness series, shadow rate, etc.
├── scripts/            # Main replication scripts
├── tables/             # Generated replication tables
├── output/             # Replication outputs
├── extension/          # Extension scripts, outputs, and figures
├── notes/              # Notes and project documentation
├── legacy/             # Older or archived code
└── Lab_export/         # Exported materials
