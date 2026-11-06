"""check_figures.py -- every plotted number of Figures 2 to 5 recomputed from the stored per-data-set results, separately from the analysis
scripts, and compared with the numbers printed in the figure files. Writes FIGURE_CHECK.md."""
import glob, os, re, subprocess, sys
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); OA = os.path.dirname(H)
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, min(1, (m - i) * p[k])); adj[k] = run
    return adj
def pdftext(f):
    return subprocess.run(["pdftotext", "-enc", "UTF-8", "-layout", os.path.join(H, f), "-"], capture_output=True, stdin=subprocess.DEVNULL).stdout.decode("utf-8", "replace")

# ---- Figure 2
rd = lambda d: pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in glob.glob(os.path.join(OA, d, "*.csv"))], ignore_index=True)
a = rd("sim3_rows"); b = rd("sim2d_rows")
key = ["design", "prev", "events", "dep", "rep"]
say("rows", len(a), len(b), "duplicated keys", a.duplicated(key).sum(), b.duplicated(key).sum())
d = a.merge(b, on=key, validate="one_to_one")
say("merged", len(d), "per setting", d.groupby("dep").size().to_dict(), "designs", d[["design", "prev", "events"]].drop_duplicates().shape[0])
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
for g in G:
    d["n_" + g] = d["named_" + g].astype(str).str.upper().eq("TRUE")
P = d[["p_" + g for g in G]].values
A = np.apply_along_axis(holm, 1, P)
say("named flags equal Holm on the stored p-values: mismatches", int(((A <= .05) != d[["n_" + g for g in G]].values).sum()))
top = np.where(d.n_COV, "revision_cov", np.where(d.n_LINK, "flexible", np.where(d.n_SLOPE, "recalibration", np.where(d.n_INTERCEPT, "intercept", "none"))))
say("stored action equals the most demanding part named: mismatches", int((top != d["loc"].values).sum()))
rank = {"none": 0, "intercept": 1, "recalibration": 2, "flexible": 3, "revision_cov": 4, "revision_linear": 4}
need = {"null": 0, "intercept": 1, "slope": 2, "transport": 2, "link": 3, "effect": 4, "inter": 4, "ushape": 4}
R = ["loc", "vergouwe", "davis", "aic"]
ORDER = ["null", "intercept", "slope", "transport", "link", "effect", "inter", "ushape"]
d["need"] = d.dep.map(need)
T = pd.DataFrame({r: [(d[d.dep == k][r].map(rank) == need[k]).mean() for k in ORDER] for r in R}, index=ORDER)
F = pd.read_csv(os.path.join(OA, "SIM2_DAVIS.csv"), keep_default_na=False).set_index("dep").loc[ORDER]
say("panel B: largest difference from the plotted file", float(np.abs(T.values - F[["right_" + r for r in R]].astype(float).values).max()))
say(T.round(2).to_string())
n0 = d[d.dep == "null"]
C = n0.groupby(["design", "prev", "events"]).apply(lambda g: pd.Series({r: (g[r].map(rank) > 0).mean() for r in R}))
Fc = pd.read_csv(os.path.join(OA, "SIM2_DAVIS_CELLS.csv")).set_index(["design", "prev", "events"]).sort_index()
say("panel A: largest difference", float(np.abs(C.sort_index().values - Fc[["false_" + r for r in R]].values).max()), "| range by rule", {r: (round(C[r].min(), 3), round(C[r].max(), 3)) for r in R})
def err(r):
    o = np.empty(len(d))
    for k in d[r].unique():
        m = (d[r] == k).values; o[m] = d.loc[m, "err_" + k].values
    return o
E = pd.DataFrame({r: err(r) for r in R}); E["dep"] = d.dep.values
E["recal"] = d.err_recalibration.values; E["rev"] = d.err_revision_cov.values
E["best"] = d[["err_" + u for u in ["none", "intercept", "recalibration", "flexible", "revision_linear", "revision_cov"]]].min(axis=1).values
M = E.groupby("dep").mean().loc[ORDER]
Cq = pd.read_csv(os.path.join(OA, "SIM2_CONSEQ.csv"), keep_default_na=False); Cq = Cq[(Cq.metric == "err") & (Cq.design == "both")].set_index("dep").loc[ORDER]
say("panel C: largest difference, rules", float(np.abs(M[R].values - F[["err_" + r for r in R]].astype(float).values).max()),
    "| fixed rules and best", float(np.abs(M[["recal", "rev", "best"]].values - Cq[["always_recalibration", "always_revision_cov", "best"]].astype(float).values).max()))
say("mean over the eight settings", M[R + ["recal", "rev"]].mean().round(2).to_dict())
say("is any rule's mean below the bar?", bool((M[R + ["recal", "rev"]].min(axis=1) < M.best - 1e-9).any()))
tx = pdftext("fig_sim2.pdf")
nums = re.findall(r"(?<![\d.])(?:1\.00|\.\d\d)(?![\d])", tx)
want = [("1.00" if v >= .995 else f"{v:.2f}"[1:]) for v in T.values.ravel()]
say("panel B numbers printed in fig_sim2.pdf:", len(nums), "| all 32 present in order:", sorted(nums) == sorted(want) or sorted(want))

# ---- Figure 3
ap = pd.read_csv(os.path.join(OA, "apps_update.csv"), keep_default_na=False); sz = pd.read_csv(os.path.join(OA, "size_apps.csv"), keep_default_na=False)
bad = 0
for _, r in ap.iterrows():
    p = [float(r[g]) for g in G]; adj = holm(p)
    nm = [g for g, q in zip(G, adj) if q <= .05]
    if (" + ".join(nm) or "none") != r.named or np.abs(adj - [float(r["adj_" + g]) for g in G]).max() > 1e-9:
        bad += 1; say("MISMATCH", r.data, r.model, r.setting, r.named, nm)
say("Figure 3 and 4 source: rows", len(ap), "verdict or adjusted p differs from Holm on the stored p:", bad)
MAP = [("random split", "random split"), ("transported", "transported, all chronic"), ("transported, half A", "transported, half A"), ("half B, intercept update", "half B after intercept update"),
       ("half B, logistic recalibration", "half B after logistic recalibration"), ("half B, flexible recalibration", "half B after flexible recalibration")]
cells = []
for sa, ss in MAP:
    for m in ["M0", "M2"]:
        for g in G:
            q = sz[(sz.data == "SUPPORT") & (sz.model == m) & (sz.setting == ss) & (sz.part == g)]; assert len(q) == 1
            cells.append(f"{max(float(q.excess.iloc[0]), 0):.1f}")
tv = re.findall(r"\d+\.\d", pdftext("fig_verdicts.pdf"))
say("Figure 3: 48 cells; printed", len(tv), "same multiset:", sorted(tv) == sorted(cells))
for m, st in [("M0", "half B, model re-estimated"), ("M2", "half B, model re-estimated")]:
    say("Web figure label", m, st, ap[(ap.model == m) & (ap.setting == st)].named.iloc[0])

# ---- Figure 4
for st in ["transported, half A", "half B, intercept update"]:
    r = ap[(ap.data == "SUPPORT") & (ap.model == "M0") & (ap.setting == st)].iloc[0]; say("Figure 4", st, [float(r[g]) for g in G], r.named, "n", r.n)

# ---- Figure 5
bp = pd.read_csv(os.path.join(OA, "burn_iecv_baseline_p.csv")); it = pd.read_csv(os.path.join(OA, "burn_iecv.csv"), keep_default_na=False).set_index("model")
for m in ["B0", "B3"]:
    g = bp[bp.model == m]; Pm = g[["p_" + k for k in G]].values; Am = np.apply_along_axis(holm, 1, Pm)
    say("Figure 5", m, "simulated sets", len(g), "familywise rate", round((Am <= .05).any(axis=1).mean(), 4), "stored flags agree", bool(((Am <= .05) == g[G].values).all()))
    for j, k in enumerate(G):
        say("  ", k, "centre-effect p", it.loc[m, k], "| plain p", it.loc[m, "plain_" + k], "| share of simulated sets at or below plain p", round((Pm[:, j] <= float(it.loc[m, "plain_" + k])).mean(), 3))
open(os.path.join(H, "FIGURE_CHECK.md"), "w", encoding="utf-8", newline="\n").write("# Figure check\n\n```\n" + "\n".join(out) + "\n```\n")
