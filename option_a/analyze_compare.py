"""analyze_compare.py -- detection: the procedure against single tests on the same data sets (n = 1,000).
Joins rows_compare/ (rival p-values, compare_detect.R) with rows/ (the procedure, run_sim_groups.R) by setting, family,
size and data set, checks that the event counts agree (same data), and reports the rejection rate at level .05 of
each test and of the procedure's global test (the intersection of all groups; it rejects exactly when some group can
be named). Writes COMPARE.csv."""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
rd = lambda pat: pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", ""]) for f in glob.glob(os.path.join(HERE, pat))])
c = rd("rows_compare/*.csv"); r = rd("rows/*.csv")
r = r[r.n == 1000]
key = ["setting", "fam", "C", "n", "rep"]
c["C"] = c["C"].round(2); r["C"] = r["C"].round(2)
m = c.merge(r[key + ["events", "global_p"]], on=key, suffixes=("", "_proc"))
assert (m.events == m.events_proc).all(), "the data sets differ"
tests = {"procedure (any part)": "global_p", "Hosmer-Lemeshow": "p_HL", "Spiegelhalter": "p_Spiegelhalter",
         "calibration belt": "p_belt", "le Cessie-van Houwelingen": "p_leCessie"}
rows = []
for (s, fam, C), x in m.groupby(["setting", "fam", "C"], sort=False):
    row = dict(setting=s, fam=fam, C=C, reps=len(x))
    for lab, col in tests.items():
        v = x[col].dropna()
        row[lab] = np.mean(v <= 0.05) if len(v) else np.nan
    rows.append(row)
R = pd.DataFrame(rows)
order = ["null", "intercept", "slope", "link", "ushape", "thresh", "inter"]
R["o"] = R.fam.map({f: i for i, f in enumerate(order)}); R = R.sort_values(["setting", "o"]).drop(columns="o")
R.to_csv(os.path.join(HERE, "COMPARE.csv"), index=False)
pd.set_option("display.width", 220)
print(R.round(3).to_string(index=False))
print("matched data sets:", len(m))
