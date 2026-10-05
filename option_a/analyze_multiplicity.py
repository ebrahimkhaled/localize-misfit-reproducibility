"""analyze_multiplicity.py -- closed testing against the simpler ways of naming groups, on the SAME data sets of the main
simulation (rows/, which store each group's own p-value). Per data set and rule, the groups named:
  closure     the procedure (stored)
  naive       every group whose own p-value is <= .05
  Bonferroni  every group whose own p-value is <= .05 / (number of groups)
  Holm        Holm's step-down over the groups' own p-values
Reported per cell: the rate of naming any group with no part (pure families and the null; nominal .05), the right-action
rate (highest group named = the rung the departure needs), and the rate of naming the departure's own group.
Writes MULTIPLICITY.csv."""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", ""]) for f in glob.glob(os.path.join(HERE, "rows", "*.csv"))])
RUNG = ["INTERCEPT", "SLOPE", "LINK", "COV"]
NEED = {"null": None, "intercept": "INTERCEPT", "slope": "SLOPE", "slope_pure": "SLOPE", "link": "LINK", "link_pure": "LINK",
        "ushape": "COV", "thresh": "COV", "inter": "COV", "cov_pure": "COV"}
PURE = {"intercept": {"INTERCEPT"}, "slope_pure": {"SLOPE"}, "link_pure": {"LINK"}, "cov_pure": {"COV"}, "null": set()}
for g in RUNG:
    c = "named_" + g
    if c in d:
        d[c] = d[c].astype(str).str.upper().eq("TRUE")


def rules(row, groups):
    p = {g: row["p_" + g] for g in groups}
    out = {"closure": {g for g in groups if row["named_" + g]},
           "naive": {g for g in groups if p[g] <= .05},
           "Bonferroni": {g for g in groups if p[g] <= .05 / len(groups)}}
    holm, k = set(), len(groups)
    for i, (g, pv) in enumerate(sorted(p.items(), key=lambda kv: kv[1])):
        if pv <= .05 / (k - i):
            holm.add(g)
        else:
            break
    out["Holm"] = holm
    return out


rows = []
for (s, fam, C, n), x in d.groupby(["setting", "fam", "C", "n"], sort=False):
    groups = RUNG if s == "external" else ["LINK", "COV"]
    need = NEED[fam]
    acc = {r: {"wrong": [], "right": [], "own": []} for r in ["closure", "naive", "Bonferroni", "Holm"]}
    for _, row in x.iterrows():
        for r, named in rules(row, groups).items():
            hi = [g for g in RUNG if g in named]
            hi = hi[-1] if hi else None
            acc[r]["right"].append(hi is None if need is None else hi == need)
            if need is not None:
                acc[r]["own"].append(need in named)
            if fam in PURE:
                acc[r]["wrong"].append(bool(named - PURE[fam]))
    for r, a in acc.items():
        rows.append(dict(setting=s, fam=fam, C=C, n=n, rule=r, reps=len(x), wrong=np.mean(a["wrong"]) if a["wrong"] else np.nan,
                         right=np.mean(a["right"]), own=np.mean(a["own"]) if a["own"] else np.nan))
M = pd.DataFrame(rows)
M.to_csv(os.path.join(HERE, "MULTIPLICITY.csv"), index=False)
pd.set_option("display.width", 220)
P = M[M.n == 1000].pivot_table(index=["setting", "fam", "C"], columns="rule", values=["wrong", "right"])
print(P.round(3).to_string())
