# analyze_sim2.py -- error rates, right action and consequences of the main simulation (sim2_rows/).
# Outputs SIM2_ERROR.csv, SIM2_RIGHT.csv, SIM2_CONSEQ.csv and prints the summaries.
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
files = glob.glob(os.path.join(HERE, "sim3_rows", "*.csv"))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in files], ignore_index=True)
print("rows", len(d), "files", len(files))
print("rows per dep:", d.groupby("dep").size().to_dict())
print("median secs_loc", round(d.secs_loc.median(), 1), "max", round(d.secs_loc.max(), 1))

G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
for g in G:
    d["named_" + g] = d["named_" + g].astype(str).str.upper().eq("TRUE")
cell = ["design", "prev", "events"]
order = ["null", "intercept", "slope", "transport", "link", "effect", "inter", "ushape"]
NEED = dict(null="none", intercept="intercept", slope="recalibration", transport="recalibration",
            link="flexible", effect="revision_cov", inter="revision_cov", ushape="revision_cov")

# 1. family-wise error under the null, per cell
n0 = d[d.dep == "null"].copy()
n0["any"] = n0[["named_" + g for g in G]].any(axis=1)
n0["ver_any"] = n0.vergouwe.ne("none")
n0["aic_any"] = n0.aic.ne("none")
e = n0.groupby(cell).agg(reps=("any", "size"), fwer=("any", "mean"), ver=("ver_any", "mean"), aic=("aic_any", "mean")).reset_index()
e.to_csv(os.path.join(HERE, "SIM2_ERROR.csv"), index=False)
print("\nFWER under null by cell (SE about .01):\n", e.round(3).to_string(index=False))
print("overall", round(n0["any"].mean(), 4), "ver", round(n0.ver_any.mean(), 4), "aic", round(n0.aic_any.mean(), 4))

# 2. false verdicts for a MORE demanding part than the departure needs
rank = dict(none=0, intercept=1, recalibration=2, flexible=3, revision_cov=4, revision_linear=4)
d["need_rank"] = d.dep.map(NEED).map(rank)
d["loc_rank"] = d["loc"].map(rank)
d["ver_rank"] = d.vergouwe.map(rank)
d["aic_rank"] = d.aic.map(rank)
r = []
for dep in order:
    g = d[d.dep == dep]
    nr = rank[NEED[dep]]
    r.append(dict(dep=dep, need=NEED[dep],
                  loc_right=(g.loc_rank == nr).mean(), ver_right=(g.ver_rank == nr).mean(), aic_right=(g.aic_rank == nr).mean(),
                  loc_over=(g.loc_rank > nr).mean(), ver_over=(g.ver_rank > nr).mean(), aic_over=(g.aic_rank > nr).mean(),
                  loc_under=(g.loc_rank < nr).mean(), ver_under=(g.ver_rank < nr).mean(), aic_under=(g.aic_rank < nr).mean()))
r = pd.DataFrame(r)
print("\nRight action / too much / too little, all cells pooled:\n", r.round(3).to_string(index=False))
rr = []
for (dep, des, ev), g in d.groupby(["dep", "design", "events"]):
    nr = rank[NEED[dep]]
    rr.append(dict(dep=dep, design=des, events=ev, loc=(g.loc_rank == nr).mean(), ver=(g.ver_rank == nr).mean(), aic=(g.aic_rank == nr).mean()))
rr = pd.DataFrame(rr)
rr.to_csv(os.path.join(HERE, "SIM2_RIGHT.csv"), index=False)
print("\nRight action by design and events:\n", rr.round(3).to_string(index=False))

# 3. consequences: the error of the update each policy prescribes
def chosen(col, metric, remap=None):
    c = d[col].replace(remap) if remap else d[col]
    out = np.empty(len(d))
    for k in c.unique():
        m = (c == k).values
        out[m] = d.loc[m, metric + "_" + k].values
    return out

UPD = ["none", "intercept", "recalibration", "flexible", "revision_linear", "revision_cov"]
rows = []
for metric in ["err", "ll"]:
    P = {"ours": chosen("loc", metric),
         "vergouwe": chosen("vergouwe", metric),
         "vergouwe_covrev": chosen("vergouwe", metric, {"revision_linear": "revision_cov"}),
         "aic": chosen("aic", metric),
         "aic_covrev": chosen("aic", metric, {"revision_linear": "revision_cov"})}
    for u in UPD:
        P["always_" + u] = d[metric + "_" + u].values
    best = d[[metric + "_" + u for u in UPD]].min(axis=1).values
    df = pd.DataFrame(P)
    df["dep"] = d.dep.values
    df["design"] = d.design.values
    df["best"] = best
    for dep in order + ["ALL"]:
        for des in ["indep5", "mixed10", "both"]:
            g = df
            if dep != "ALL":
                g = g[g.dep == dep]
            if des != "both":
                g = g[g.design == des]
            row = dict(metric=metric, dep=dep, design=des)
            for k in P:
                row[k] = g[k].mean()
            row["best"] = g["best"].mean()
            rows.append(row)
c = pd.DataFrame(rows)
c.to_csv(os.path.join(HERE, "SIM2_CONSEQ.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
show = ["dep", "ours", "vergouwe", "vergouwe_covrev", "aic", "aic_covrev", "always_none", "always_recalibration", "always_flexible", "always_revision_linear", "always_revision_cov", "best"]
print("\nMean abs risk error (pp), both designs:\n", c[(c.metric == "err") & (c.design == "both")][show].round(2).to_string(index=False))
print("\nmixed10 only:\n", c[(c.metric == "err") & (c.design == "mixed10")][show].round(2).to_string(index=False))
print("\nExpected log-loss x 100, both designs:\n", c[(c.metric == "ll") & (c.design == "both")][show].assign(**{k: lambda x, k=k: 100 * x[k] for k in show[1:]}).round(2).to_string(index=False))

# paired difference of ours against each rival, by departure, with its SE
pe = {k: chosen(*v) for k, v in dict(ours=("loc", "err"), ver=("vergouwe", "err"), aic=("aic", "err")).items()}
pe["vcov"] = chosen("vergouwe", "err", {"revision_linear": "revision_cov"})
pe["acov"] = chosen("aic", "err", {"revision_linear": "revision_cov"})
pe["recal"] = d.err_recalibration.values
pe["revc"] = d.err_revision_cov.values
print("\nPaired difference in risk error, rival minus ours (pp), SE in brackets; positive = ours better:")
for dep in order:
    m = (d.dep == dep).values
    s = [dep.ljust(10)]
    for k in ["ver", "vcov", "aic", "acov", "recal", "revc"]:
        x = pe[k][m] - pe["ours"][m]
        s.append(f"{k} {x.mean():+.2f} ({x.std() / np.sqrt(len(x)):.2f})")
    print("  ".join(s))
