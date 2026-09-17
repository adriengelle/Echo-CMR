# Echo Parakeet — Multi-state CMR Survival Analysis

An R script that fits a multi-state (Arnason-Schwarz) Cormack-Jolly-Seber model
to capture-mark-recapture data using **RMark** and **Program MARK**, and outputs
age- and state-dependent survival, detection and transition estimates plus a
survival figure.

It runs **locally on your own machine** — there is no web app and no server, so
there are no timeouts or memory limits. You need R and Program MARK installed
(see below).

---

## What the script does

Given a CSV of encounter histories, the script:

1. Builds a multi-state capture history from the yearly state columns
   (pre-breeder / breeder / post-breeder / not seen).
2. Fits an Arnason-Schwarz multi-state CJS model in MARK, with:
   - **Survival (Φ)** varying by age class (juvenile 0–1 vs adult 2+) and year;
   - **Detection (p)** varying by age and state over time;
   - **Transitions (Ψ)** between states, varying by age and population growth phase;
   - biologically impossible age × state combinations fixed to zero.
3. Writes a model summary and CSVs of the survival, detection and transition
   estimates.
4. Produces a two-panel figure of juvenile vs adult survival over time, with a
   dashed decline threshold drawn on the juvenile panel.

---

## Prerequisites

### 1. R (and RStudio)

- **R** (≥ 4.2 recommended): https://cran.r-project.org/
- **RStudio Desktop** (recommended): https://posit.co/download/rstudio-desktop/

### 2. Program MARK

RMark does not fit models itself — it drives **Program MARK**, which must be
installed on your computer.

- Download MARK: http://www.phidot.org/software/mark/downloads/
- Run the installer (Windows). On Windows the installer normally places
  `mark.exe` where RMark can find it automatically.
- macOS / Linux: you must install the MARK command-line executable and make it
  available on your `PATH`. See the phidot download page for platform notes.

### 3. R packages

Run once in R:

```r
install.packages(c("RMark", "dplyr", "ggplot2"))
```

---

## Getting the script

**Option A — clone from R (recommended):**

```r
install.packages("usethis")   # if needed
usethis::create_from_github("https://github.com/<USERNAME>/<REPO>.git",
                            destdir = "~/projects")
```

This clones the repo and opens it as an RStudio project.

**Option B — clone with git:**

```bash
git clone https://github.com/<USERNAME>/<REPO>.git
```

**Option C — download:** click **Code → Download ZIP** on the GitHub page and
unzip it.

---

## Data format

Put your CSV in the `data/` folder. It must have:

- **One column per survey year**, with the year as the column header
  (e.g. `1994`, `1995`, … `2022`).
- **State codes** in those columns: `PB` (pre-breeder), `B` (breeder),
  `PO` (post-breeder), and `0` (not seen).
- A **`chrt`** column giving each individual's cohort (ringing) year.
- A **`subpop`** column with the subpopulation (`BO` = Bel Ombre, `GG` = Gorges).

Example (abridged):

| id | subpop | chrt | 1994 | 1995 | 1996 | … |
|----|--------|------|------|------|------|---|
| 1  | GG     | 1994 | PB   | B    | 0    | … |
| 2  | BO     | 1995 | 0    | PB   | PO   | … |

> The real echo parakeet dataset is not included in this repository. Supply your
> own CSV in the format above, or add a small synthetic example to `data/`.

---

## Running the analysis

1. Open the project in RStudio (double-click the `.Rproj`, or use Option A above)
   so the working directory is the repository root.
2. Open `cmr_multistate_analysis.R` and edit the **CONFIG** block at the top:
   - `data_file` — path to your CSV;
   - `start_year` / `end_year` — the first and last year columns to use;
   - `last_cohort_year` — drop cohorts ringed after this year;
   - `threshold` — the juvenile decline threshold for the figure.
3. Source the whole script: **Ctrl + Shift + S** (or click **Source**).

Fitting calls MARK and can take from under a minute to several minutes depending
on your data.

---

## Outputs

All written to the `output/` folder:

| File | Contents |
|------|----------|
| `survival_plot.png` | Juvenile vs adult survival over time (300 dpi) |
| `model_summary.txt` | Model summary, AICc, deviance, parameter count |
| `survival_estimates.csv` | Real survival (Φ) estimates with SE and 95% CI |
| `detection_estimates.csv` | Real detection (p) estimates |
| `transition_estimates.csv` | Real transition (Ψ) estimates |

---

## Troubleshooting

**"MARK not found" / RMark cannot locate the executable.** MARK is not installed,
or `mark.exe` is not where RMark looks. Reinstall MARK, or point RMark at it, e.g.:

```r
MarkPath <- "C:/Program Files (x86)/MARK"   # folder containing mark.exe
```

set before calling the script. If your installer produced `mark64.exe` rather
than `mark.exe`, copy or rename it to `mark.exe`, or set `MarkPath` to its folder.

**"CSV is missing year columns".** `start_year`/`end_year` in the CONFIG block
don't match the year headers in your file. Check the exact column names.

**The model runs but survival looks flat or fails to converge.** Check the state
codes and the cohort/year settings; a very sparse dataset may not support the
fully time-varying structure.

---

## Contact

Adrien Gelle — adriengelle@gmail.com
