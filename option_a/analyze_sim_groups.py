"""analyze_sim_groups.py -- operating characteristics of the orthogonal-group localization (run_sim_groups.R).
Per cell: the rate each group is named; any group named; the RIGHT-ACTION rate (the highest rung named equals the
rung the departure needs: INTERCEPT < SLOPE < LINK < COV); and, for the PURE departures, the rate at which a group
whose part of the misfit is zero is named anyway (strong familywise error; nominal .05). Writes SIM_GROUPS.csv."""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in glob.glob(os.path.join(HERE, "rows", "*.csv"))])
RUNG = ["INTERCEPT", "SLOPE", "LINK", "COV"]
NEED = {"null": None, "intercept": "INTERCEPT", "slope": "SLOPE", "slope_pure": "SLOPE", "link": "LINK", "link_pure": "LINK",
        "ushape": "COV", "thresh": "COV", "inter": "COV", "cov_pure": "COV"}
PURE = {"intercept": {"INTERCEPT"}, "slope_pure": {"SLOPE"}, "link_pure": {"LINK"}, "cov_pure": {"COV"}, "null": set()}
for g in RUNG:
    c = "named_" + g
    d[c] = d[c].astype(str).str.upper().eq("TRUE") if d[c].dtype == object else d[c].astype(bool)


def highest(row):
    named = [g for g in RUNG if row["named_" + g]]
    return named[-1] if named else None


d["highest"] = d.apply(highest, axis=1)
rows = []
for (s, fam, C, n), g in d.groupby(["setting", "fam", "C", "n"], sort=False):
    need = NEED[fam]
    r = dict(setting=s, fam=fam, C=C, n=n, reps=len(g), any_named=(g.highest.notna()).mean())
    for k in RUNG:
        if s == "insample" and k in ("INTERCEPT", "SLOPE"):
            continue
        r[k] = g["named_" + k].mean()
    r["right_action"] = (g.highest.isna() if need is None else (g.highest == need)).mean()
    if fam in PURE:
        wrong = [k for k in RUNG if k not in PURE[fam] and not (s == "insample" and k in ("INTERCEPT", "SLOPE"))]
        r["wrong_group_named"] = g[["named_" + k for k in wrong]].any(axis=1).mean()
    r["secs_median"] = g.secs.median()
    rows.append(r)
R = pd.DataFrame(rows)
order = ["null", "intercept", "slope", "slope_pure", "link", "link_pure", "ushape", "thresh", "inter", "cov_pure"]
R["o"] = R.fam.map({f: i for i, f in enumerate(order)})
R = R.sort_values(["setting", "o", "C", "n"]).drop(columns="o")
R.to_csv(os.path.join(HERE, "SIM_GROUPS.csv"), index=False)
pd.set_option("display.width", 220)
print(R.round(3).to_string(index=False))
