# analyze_davis.py -- merges the two further rivals (sim2d_rows/, rows_davis/) with the main simulation rows and reports
# false updates under no misfit, right action, over-updating and risk error. Writes SIM2_DAVIS.csv, RIVALS_DAVIS.csv.
import glob, os
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
rd = lambda p: pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in glob.glob(os.path.join(H, p, "*.csv"))], ignore_index=True)
key = ["design", "prev", "events", "dep", "rep"]
d = rd("sim3_rows").merge(rd("sim2d_rows"), on=key, validate="one_to_one")
print("merged", len(d))
ORDER = ["null", "intercept", "slope", "transport", "link", "effect", "inter", "ushape"]
need = dict(null=0, intercept=1, slope=2, transport=2, link=3, effect=4, inter=4, ushape=4)
rank = dict(none=0, intercept=1, recalibration=2, flexible=3, revision_cov=4, revision_linear=4)
RULES = ["loc", "vergouwe", "ctp_flex", "davis", "aic"]
def err(col, metric="err", remap=None):
    c = d[col].replace(remap) if remap else d[col]
    o = np.empty(len(d))
    for k in c.unique():
        m = (c == k).values
        o[m] = d.loc[m, metric + "_" + k].values
    return o
for r in RULES:
    d["e_" + r] = err(r); d["l_" + r] = err(r, "ll"); d["k_" + r] = d[r].map(rank)
    d["ec_" + r] = err(r, "err", {"revision_linear": "revision_cov"})
d["nr"] = d.dep.map(need)
rows = []
for dep in ORDER:
    g = d[d.dep == dep]
    row = dict(dep=dep)
    for r in RULES:
        row["right_" + r] = (g["k_" + r] == g.nr).mean(); row["over_" + r] = (g["k_" + r] > g.nr).mean()
        row["err_" + r] = g["e_" + r].mean(); row["errcov_" + r] = g["ec_" + r].mean(); row["ll_" + r] = g["l_" + r].mean()
    rows.append(row)
t = pd.DataFrame(rows)
t.to_csv(os.path.join(H, "SIM2_DAVIS.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
for p in ["right_", "over_", "err_", "errcov_"]:
    c = [p + r for r in RULES]
    x = t[["dep"] + c].copy()
    print("\n", p, "\n", x.round(3).to_string(index=False), "\n mean", t[c].mean().round(3).to_dict())
n0 = d[d.dep == "null"]
cells = n0.groupby(["design", "prev", "events"]).apply(lambda g: pd.Series({"false_" + r: (g["k_" + r] > 0).mean() for r in RULES})).reset_index()
cells.to_csv(os.path.join(H, "SIM2_DAVIS_CELLS.csv"), index=False)
print("\nfalse update under no misfit by cell, max:", {r: round(n0.groupby(["design", "prev", "events"])["k_" + r].apply(lambda s: (s > 0).mean()).max(), 3) for r in RULES},
      "overall", {r: round((n0["k_" + r] > 0).mean(), 3) for r in RULES})
print("paired risk error, rule minus ours (pp, SE):")
for dep in ORDER:
    m = d.dep == dep
    print(dep.ljust(10), "  ".join(f"{r} {(d.loc[m, 'e_' + r] - d.loc[m, 'e_loc']).mean():+.2f} ({(d.loc[m, 'e_' + r] - d.loc[m, 'e_loc']).std() / np.sqrt(m.sum()):.2f})" for r in ["davis", "ctp_flex"]))
# focused
f = rd("rows_davis")
rows = []
for (fam, C), g in f.groupby(["fam", "C"]):
    n = {"null": 0, "intercept": 1, "slope": 2, "link": 3, "ushape": 4, "thresh": 4, "inter": 4}[fam]
    rows.append(dict(fam=fam, C=C, reps=len(g), right_davis=(g.davis.map(rank) == n).mean(), right_ctp_flex=(g.ctp_flex.map(rank) == n).mean(),
                     over_davis=(g.davis.map(rank) > n).mean(), over_ctp_flex=(g.ctp_flex.map(rank) > n).mean()))
ff = pd.DataFrame(rows); ff.to_csv(os.path.join(H, "RIVALS_DAVIS.csv"), index=False); print(ff.round(3).to_string(index=False))
