# Reproducibility archive: Which update does a clinical risk model need? Localization of misfit in logistic risk models

Ebrahim Khaled Ebrahim, Osama Abd El-Aziz Hussein and Ahmed El-Kotory, Department of Statistics, Faculty of Business,
Alexandria University, Egypt.

## Version 2.1.2: what the paper uses

```
option_a/R/localize_update.R        localize_external_update(): the default, updated-model reference with Holm's correction;
                                    options cluster and tau (centre-effect reference) and main (see below)
option_a/R/localize_size.R          localize_size(): size of each part, chance level and noise interval
option_a/R/localize_groups.R        group bases, frozen-model reference and the in-sample check
option_a/R/localize_plot.R          plot() for a verdict
option_a/run_sim3.R                 main simulation (twelve designs, eight settings); sim3_rows/, one row per data set
option_a/run_davis.R                Brier-score procedure and closed testing with a flexible step, same data sets; sim2d_rows/, rows_davis/
option_a/run_rivals_update.R        focused comparison; rows_update/; run_vergouwe.R, vergouwe_rows/ for the rivals
option_a/analyze_sim2.py, analyze_davis.py   summaries (SIM2_*.csv, RIVALS_DAVIS.csv)
option_a/apps_update.R, size_apps.R          SUPPORT and burn applications
option_a/burn_iecv2.R               burn internal-external cross-validation with all refits repeated in the simulations
option_a/theory/check_*.R           stress tests, clustering, machine-learning predictor, size checks (CHECK_*.csv)
option_a/worked_example.R           the worked example of the paper, with the functions of the R package ebrahim.gof; output in worked_example.txt
option_a/paper_sim/                 scripts that write the tables and figures
```
Set `main = FALSE` for risks pooled by cross-validation: they are not fully frozen, and the main-effects member responds to the
differences between folds. Scripts carry the local path of the study in `HERE`; change it to the archive folder before running.
`run_sim2.R` is kept because `run_sim3.R` was generated from it; the data sets of the old cells are identical.

## Earlier material (frozen-model reference)

This archive regenerates every table and figure of the paper from deposited per-data-set results, and holds the
code that produced those results.

## Contents

```
option_a/R/localize_groups.R      the procedure: localize_external() and localize_insample()
option_a/R/localize_plot.R        plot() for a verdict: misfit compass and lattice of closed tests (class gof_localize)
option_a/run_sim_groups.R          main simulation (Table 2); results in option_a/rows/, summary SIM_GROUPS.csv
option_a/analyze_sim_groups.py     builds SIM_GROUPS.csv from option_a/rows/
option_a/run_purity_prob.R         probability-scale purity check; option_a/rows_prob/, PURITY_PROB.csv
option_a/support_v2.R              SUPPORT application (in-sample, de-aliased, random split, transported, neural network)
option_a/support_split_first.R     SUPPORT split-first update-and-retest (lower rows of Table 3)
option_a/support_localize.R        first SUPPORT run (superseded by support_v2.R; kept for the record)
option_a/support_timing.R          timing on the SUPPORT data
option_a/burn_localize.R, burn_revise.R   burn application (in-sample, revisions, 20 cross-centre facility splits)
option_a/burn_null.R                     the same splits with outcomes simulated from a correct model (estimation-error baseline)
option_a/compare_detect.R          single tests on the main simulation's data sets; analyze_compare.py -> COMPARE.csv
option_a/analyze_multiplicity.py   closure vs naive, Bonferroni and Holm naming, from option_a/rows/ -> MULTIPLICITY.csv
option_a/run_sim_small.R           small samples (n = 200, 300); analyze_sim_small.py -> SIM_SMALL.csv
option_a/sim_robust_power.R        the robust option on the main data sets -> ROBUST_POWER.csv
option_a/paper/make_support_figure.R   Web Figure B.1
option_a/theory/check_*.R          numerical checks of the theory (Web Appendix A) and the additional studies (Web Appendix B)
option_a/theory/make_numerics.py   tables of Web Appendix A
option_a/paper/make_tables.py      Table 1 of the paper (and the table versions of Figures 1 and 3)
option_a/paper/make_table_figures.R   Figures 1 and 3
option_a/paper/make_compass_figure.R  Figure 4 (from support_closure_full.rds, written by option_a/support_closure_full.R)
option_a/paper/make_appendix_b.py  tables of Web Appendix B
option_a/paper/make_figure.R       Figure 2
localize_robust/                   corrupted-records check (run_groups_corruption.R, rows_groups/, CORRUPT_GROUPS.csv;
                                   the robust option: run_groups_robust.R, rows_groups_robust/, CORRUPT_ROBUST.csv)
data/support2.csv                  SUPPORT2 (UCI Machine Learning Repository, CC BY 4.0)
data/battery/                      stored data-set seeds and fingerprints for the corrupted-records check
MANIFEST.sha256                    sha256 of every file
```

## How to run

R 4.4 or later with the packages `splines`, `parallel` and `nnet` (all shipped with R), `aplore3` (burn data),
`ebrahim.gof` (2.8.0 or later) and `givitiR` (the comparison with single tests); Python 3 with `pandas` and `numpy`
for the table builders. The procedure is also in the R package `ebrahim.gof` as `localize.external()` and
`localize.gof()` (from version 2.10.0).

Every R script finds its files from the archive root, taken from the environment variable `LOCALIZE_ARCHIVE_ROOT`
or, if it is unset, the working directory:

```
cd <archive root>
export LOCALIZE_ARCHIVE_ROOT="$PWD"        # Windows: set LOCALIZE_ARCHIVE_ROOT=%CD%
Rscript option_a/smoke_localize.R          # a one-minute check of the procedure
Rscript option_a/support_v2.R              # the SUPPORT application (about 30 minutes)
```

The table builders run from their own folders, for example `cd option_a/paper && python make_tables.py`.

**Simulations.** Each simulation script writes one file per chunk of 10 or 20 data sets and skips files that exist, so
an interrupted run resumes. Seeds are fixed per cell and data set inside each script. To rerun a study from scratch,
delete its result folder first. The main runs used 4 to 12 worker processes; the results do not depend on the number
of workers.

## Using the procedure on your own data

```r
source("option_a/R/localize_groups.R")
# external validation of frozen predictions p (any model) for outcomes y and covariates X
localize_external(y, p, X, M = 999)
# in-sample checking of a fitted logistic glm
localize_insample(fit, X, B = 199)
```
Both return the parts named, the p-value of every intersection test and of every member test.
