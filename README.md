# Echo parakeet — Multistate CMR Survival Analysis

An R script that fits a multistate (Arnason-Schwarz) Cormack-Jolly-Seber (CJS) model to capture-mark-recapture data using **RMark** (Laake and Rexstad 2022; Laake 2013) and **Program MARK** (White and Burnham 1999), and outputs age- and time-dependent survival, detection and transition estimates, a survival figure, and model summary.

It runs **locally on your PC** — (no timeouts or memory limits). You need **R** and **Program MARK** installed (see below).

------------------------------------------------------------------------

## What the script does

Given a CSV of encounter histories, the script:

1.  Builds a multistate capture history from the yearly state columns (pre-breeder / breeder / post-breeder / not seen).
2.  Fits an Arnason-Schwarz multistate CJS model in MARK, with:
    - **Survival (Φ)** varying by age class (juvenile 0-1 & adult 2+), state (breeders and non-breeders) and year;
    - **Detection (p)** varying by age and state over time;
    - **Transitions (Ψ)** between states, varying by age and population growth phases;
    - Biologically impossible age × state combinations for this species fixed to zero (Gellé et al. 2026).
3.  Writes a model summary and CSVs of the survival, detection and transition estimates.
4.  Produces a three-panel figure of juvenile vs adult survival (breeders and pre-breeders) over time, with a dashed threshold line level below which the population is expected to decline over the long term. Derived from a Population Viability Analysis (the point where the stochastic growth rate r \< 0). Sustained survival above the line supports a stable or growing population.
5.  Produces a figure of recapture rates over time (index of field team performance).

------------------------------------------------------------------------

## Prerequisites

### 1. R (and RStudio)

- **R** (≥ 4.4 recommended): <https://cran.r-project.org/>
- **RStudio** (recommended): <https://posit.co/download/rstudio-desktop/>

### 2. Program MARK

RMark does not fit models itself — it communicates with **Program MARK**, which must be installed on your computer.

- Download MARK: <http://www.phidot.org/software/mark/>
- Run the installer (Windows). On Windows the installer normally places `mark.exe` where RMark can find it automatically.
- macOS / Linux: you must install the MARK command-line executable and make it available on your `PATH`. See the phidot download page for platform notes.

### 3. R packages

Run once in R:

``` r
install.packages(c("RMark", "dplyr", "ggplot2"))
```

------------------------------------------------------------------------

## Getting the script

**Option A — clone from R (recommended):** In RStudio → **File** → **New Project...** → **Version Control** → **Git** → Paste repo URL:<https://github.com/adriengelle/Echo-CMR.git> → **Create Project**

This clones the repo and opens it as an RStudio project (if any future updates, only "Pull" required in the 'Git' tab)

**Option B — download:** click **Code → Download ZIP** on the GitHub page and unzip it.

------------------------------------------------------------------------

## Data format

Put your CSV (**only one file!**) in the `data/` folder. It must have:

- **One column per survey year**, with the year as the column header (e.g. `1994`, `1995`, … `2022`).
- **State codes** in those columns: `PB` (pre-breeder), `B` (breeder), `PO` (post-breeder), and `0` (not seen).
- A **`chrt`** column giving each individual's cohort (ringing) year.
- A **`subpop`** column with the subpopulation (`BO` = Bel Ombre, `GG` = Gorges).
- Individual 'ID' are unimportant as long as it is one row per individual history
- The data file name does not matter
- For more information about how to pivot re-sightings data into encounter histories, visit my other GitHub [repository](https://github.com/adriengelle/Formatting-data-for-CMR) or this [Shiny app](https://adrien-gelle.shinyapps.io/cmr-pivot/)

Example (see also `examples/`):

| id  | subpop | chrt | 1994 | 1995 | 1996 | […] |
|-----|--------|------|------|------|------|-----|
| 1   | GG     | 1994 | PB   | PB   | B    | […] |
| 2   | BO     | 1995 | 0    | PB   | 0    | […] |

> The full echo parakeet re-sighting dataset is not included in this repository. Supply your own CSV in the format above.

------------------------------------------------------------------------

## Running the analysis

1.  Open the project in RStudio (double-click the `.Rproj`, or use Option A above) so the working directory is the repository root.
2.  Open `cmr_multistate_analysis.R`.
3.  Source the whole script: **Ctrl + Shift + S** (or click **Source**).

Fitting calls MARK and can take from under a minute to several minutes depending on your data.

------------------------------------------------------------------------

## Outputs

All written to the `output/` folder:

| File                       | Contents                                       |
|---------------------------|---------------------------------------------|
| `survival_plot.png`        | Juvenile & adult survival over time (300 dpi)  |
| `model_summary.txt`        | Model summary, AICc, deviance, parameter count |
| `survival_estimates.csv`   | Survival (Φ) estimates with SE and 95% CI      |
| `detection_estimates.csv`  | Detection (p) estimates with SE and 95% CI     |
| `transition_estimates.csv` | Transition (Ψ) estimates with SE and 95% CI    |

------------------------------------------------------------------------

## Troubleshooting

**"MARK not found" / RMark cannot locate the executable.** MARK is not installed, or `mark.exe` is not where RMark looks. Reinstall MARK, or point RMark at it, e.g.:

``` r
MarkPath <- "C:/Program Files (x86)/MARK"   # folder containing mark.exe
```

set before calling the script. If your installer produced `mark64.exe` rather than `mark.exe`, copy or rename it to `mark.exe`, or set `MarkPath` to its folder.

------------------------------------------------------------------------

## References

- Gellé, A. et al., 2026. Demographic responses to population recovery illustrated by 30 years of monitoring a once critically endangered parrot. Journal of Applied Ecology, 63 (8), e70523. [10.1111/1365-2664.70523](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2664.70523).
- Laake, J. and Rexstad, E., 2022. RMark – an alternative approach to building linear models in MARK.
- Laake, J.L. (Jeffrey L., 2013. RMark : an R Interface for analysis of capture-recapture data with MARK [online]. Available at: <https://repository.library.noaa.gov/view/noaa/4372> [Accessed 10 September 2024].
- White, G.C. and Burnham, K.P., 1999. Program MARK: survival estimation from populations of marked animals. Bird Study, 46 (sup1), S120–S139. [10.1080/00063659909477239](https://www.tandfonline.com/doi/abs/10.1080/00063659909477239).
