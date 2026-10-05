"""make_appendix_b.py -- Web Appendix B tables from the stored per-data-set results:
tab_appendix_b.tex (multiplier calibration, de-aliasing, corrupted records) and tab_plan.tex (SUPPORT planning).
Also writes ROBUST.csv, DEALIAS.csv and CORRUPT.csv next to the results."""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.path.join(HERE, "..", "theory")
LR = os.path.join(HERE, "..", "..", "localize_robust")
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
OWN = {"null": set(), "link_b": {"LINK"}, "cov_b": {"COV"}, "cov": {"COV"}, "link": {"LINK"},
       "transport": {"INTERCEPT", "SLOPE", "LINK"}}
B = lambda s: s.astype(str).str.upper().eq("TRUE")
read = lambda pat: pd.concat([pd.read_csv(f, keep_default_na=False) for f in glob.glob(pat) if "SMOKE" not in f])
out = []

# ---- multiplier calibration
d = read(os.path.join(TH, "rows_robust", "*.csv"))
for g in G:
    d["named_" + g] = B(d["named_" + g])
rows = []
for (f, n, cal), x in d.groupby(["fam", "n", "calibration"]):
    own = [g for g in G if g in OWN[f]]; oth = [g for g in G if g not in OWN[f]]
    w = x[["named_" + g for g in oth]].any(axis=1).mean()
    rows.append(dict(fam=f, n=n, calibration=cal, reps=len(x), own=x[["named_" + g for g in own]].all(axis=1).mean() if own else np.nan,
                     wrong=w, se=np.sqrt(w * (1 - w) / len(x))))
R = pd.DataFrame(rows); R.to_csv(os.path.join(TH, "ROBUST.csv"), index=False)
LAB = {"null": "none", "link_b": "link (probability scale)", "cov_b": "covariate (probability scale)",
       "transport": "gross miscalibration$^a$"}
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{Monte Carlo versus multiplier calibration, external validation, departures of fixed size that are pure "
        r"on the probability scale. Own: every group with a part named; wrong: any other group named (nominal bound .05). "
        r"Same data sets for both calibrations; 400 data sets per cell (200 at $n=16{,}000$).}\label{tab:robust}",
        r"\begin{tabular}{lrrrrr}\toprule & & \multicolumn{2}{c}{Monte Carlo} & \multicolumn{2}{c}{multiplier}\\",
        r"\cmidrule(lr){3-4}\cmidrule(lr){5-6} departure & $n$ & own & wrong & own & wrong\\ \midrule"]
for f in ["null", "link_b", "cov_b", "transport"]:
    for n in sorted(R[R.fam == f].n.unique()):
        c = []
        for cal in ["montecarlo", "multiplier"]:
            y = R[(R.fam == f) & (R.n == n) & (R.calibration == cal)]
            c += ["--" if y.empty or np.isnan(y.own.iloc[0]) else f"{y.own.iloc[0]:.2f}", "--" if y.empty else f"{y.wrong.iloc[0]:.3f}"]
        out.append(f"{LAB[f]} & {n:,} & ".replace(",", "{,}") + " & ".join(c) + r"\\")
out += [r"\bottomrule\end{tabular}",
        r"\par\footnotesize{$^a$$\mathrm{logit}\,\pi=-2+0.8\eta_0$: about 11\% observed against 40\% predicted. "
        r"\textsc{intercept}, \textsc{slope} and \textsc{link} have parts on the probability scale, so only \textsc{cov} counts "
        r"as wrong; own: all three named.}", r"\end{table}", ""]

# ---- de-aliasing
d = read(os.path.join(TH, "rows_dealias", "*.csv"))
for g in ["LINK", "COV"]:
    d["named_" + g] = B(d["named_" + g])
d["dealias"] = B(d["dealias"])
D = d.groupby(["fam", "dealias"]).agg(reps=("rep", "size"), LINK=("named_LINK", "mean"), COV=("named_COV", "mean")).reset_index()
D.to_csv(os.path.join(TH, "DEALIAS.csv"), index=False)
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{In-sample checking with a skewed and a binary covariate ($\E(X\mid\eta_0)$ not affine), $n=1{,}000$, $B=199$: "
        r"rate at which each group is named, default and de-aliased link group (Proposition~\ref{prop:dealias}); "
        + f"{int(D.reps.min())}" + r" data sets per cell.}\label{tab:dealias}",
        r"\begin{tabular}{lrrrr}\toprule & \multicolumn{2}{c}{default} & \multicolumn{2}{c}{de-aliased}\\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5} departure & LINK & COV & LINK & COV\\ \midrule"]
for f, lab in [("null", "none"), ("cov", "covariate structure"), ("link", "link shape")]:
    c = []
    for da in [False, True]:
        y = D[(D.fam == f) & (D.dealias == da)]
        c += [f"{y.LINK.iloc[0]:.3f}", f"{y.COV.iloc[0]:.3f}"] if len(y) else ["--", "--"]
    out.append(lab + " & " + " & ".join(c) + r"\\")
out += [r"\bottomrule\end{tabular}\end{table}", ""]

# ---- corrupted records
d = read(os.path.join(LR, "rows_groups", "*.csv"))
for g in ["named_LINK", "named_COV", "dealias"]:
    d[g] = B(d[g])
d["any"] = d.named_LINK | d.named_COV
d["block"] = d["block"].astype(str)
C = d.groupby(["block", "file", "dealias"]).agg(reps=("rep", "size"), any=("any", "mean"), LINK=("named_LINK", "mean"),
                                               COV=("named_COV", "mean")).reset_index()
C.to_csv(os.path.join(LR, "CORRUPT_GROUPS.csv"), index=False)
CL = [("9b", "clean", "$n=500$, clean"), ("9b", "x4", "$n=500$, $x\\times4$ in 3 records"), ("9b", "x8", "$n=500$, $x\\times8$ in 3 records"),
      ("9", "logit_clean_n1000", "$n=1{,}000$, clean"), ("9", "logit_C1_r010_n1000", "$n=1{,}000$, $x\\times4$ in 10 records")]
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{Corrupted records, in-sample, correct model: rate at which any group, LINK or COV is named (all are "
        r"false); default and de-aliased; 100 data sets per cell.}\label{tab:corrupt}",
        r"\begin{tabular}{lrrrrrr}\toprule & \multicolumn{3}{c}{default} & \multicolumn{3}{c}{de-aliased}\\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7} data & any & LINK & COV & any & LINK & COV\\ \midrule"]
for b, f, lab in CL:
    c = []
    for da in [False, True]:
        y = C[(C.block == b) & (C.file == f) & (C.dealias == da)]
        c += [f"{y['any'].iloc[0]:.2f}", f"{y.LINK.iloc[0]:.2f}", f"{y.COV.iloc[0]:.2f}"] if len(y) else ["--"] * 3
    out.append(lab + " & " + " & ".join(c) + r"\\")
out += [r"\bottomrule\end{tabular}\end{table}", ""]
open(os.path.join(HERE, "tab_appendix_b.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(out))

# ---- burn application (burn_localize.csv, burn_revise.csv, burn_revise_splits.csv)
bl = pd.read_csv(os.path.join(HERE, "..", "burn_localize.csv"), keep_default_na=False, na_values=["NA", ""])
br = pd.read_csv(os.path.join(HERE, "..", "burn_revise.csv"), keep_default_na=False, na_values=["NA", ""])
bs = pd.read_csv(os.path.join(HERE, "..", "burn_revise_splits.csv"), keep_default_na=False, na_values=["NA", ""])
for g in ["named_INTERCEPT", "named_SLOPE", "named_LINK", "named_COV", "any"]:
    bs[g] = B(bs[g])
fp = lambda x: "--" if pd.isna(x) else (f"{float(x):.3f}" if float(x) >= .001 else "$<$0.001")
nm = lambda s: "none" if s == "none" else " + ".join(r"\textsc{" + g.lower() + "}" for g in str(s).split(" + "))
out_b = [r"\begin{table}[ht]\centering\small",
         r"\caption{Burn injury mortality (1,000 patients, 40 facilities, 150 deaths). In-sample checking ($B=199$) of the "
         r"main-effects model B0 and of three revisions, and cross-centre external validation over 20 random splits of the "
         r"facilities into development and validation halves ($M=999$): parts named and single-group $p$-values in-sample; "
         r"share of the 20 splits naming each part externally.}\label{tab:burn}",
         r"\begin{tabular}{llrrr}\toprule model & in-sample: parts named & AIC & L & C\\ \midrule"]
LABM = {"B0": "B0 main effects", "B1": "B1: spline in age", "B2": "B2: B1 + age $\\times$ inhalation", "B3": "B3: B2 + spline in burn area"}
for _, r in br.iterrows():
    out_b.append(f"{LABM[r.model]} & {nm(r.named)} & {r.aic:.1f} & {fp(r.LINK)} & {fp(r.COV)}\\\\")
da = bl[bl.analysis == "in-sample, de-aliased"].iloc[0]
out_b.append(f"B0, de-aliased link group & {nm(da.named)} & & {fp(da.LINK)} & {fp(da.COV)}\\\\")
out_b += [r"\midrule & \multicolumn{4}{l}{cross-centre validation: share of 20 facility splits naming}\\",
          r"model & I & S & L & C\\ \midrule"]
for m in ["B0", "B3"]:
    x = bs[bs.model == m]
    out_b.append(f"{LABM[m]} & " + " & ".join(f"{x[c].mean():.2f}" for c in ["named_INTERCEPT", "named_SLOPE", "named_LINK", "named_COV"]) + r"\\")
bn = pd.read_csv(os.path.join(HERE, "..", "burn_null.csv"))
out_b += [r"\midrule & \multicolumn{4}{l}{same splits, outcomes simulated from the model fitted to all patients}\\",
          r"& \multicolumn{4}{l}{(model right by construction; mean over 25 data sets)}\\ \midrule"]
for m in ["B0", "B3"]:
    x = bn[bn.model == m]
    out_b.append(f"{LABM[m]} & " + " & ".join(f"{B(x[c]).mean():.2f}" for c in ["named_INTERCEPT", "named_SLOPE", "named_LINK", "named_COV"]) + r"\\")
out_b += [r"\bottomrule\end{tabular}\end{table}", ""]
open(os.path.join(HERE, "tab_appendix_b.tex"), "a", encoding="utf-8", newline="\n").write("\n".join(out_b))

# ---- detection against single tests (COMPARE.csv), naming rules (MULTIPLICITY.csv), small samples (SIM_SMALL.csv)
OA = os.path.join(HERE, "..")
FL = {"null": "none", "intercept": "intercept shift", "slope": "slope 0.6", "link": "curved link", "ushape": "U-shape",
      "thresh": "threshold", "inter": "interaction"}
cp = pd.read_csv(os.path.join(OA, "COMPARE.csv"), keep_default_na=False, na_values=["NA", ""])
f3 = lambda v: "--" if pd.isna(v) else f"{v:.2f}"
out_c = [r"\begin{table}[ht]\centering\small",
         r"\caption{Detection at $n=1{,}000$, the same data sets as Table~2 of the paper: rejection rate at level 0.05 of the "
         r"procedure's global test (the intersection of all groups) and of single tests; under no misfit the rate is the level. A part is named in marginally fewer data sets than the global test rejects (Table 2 of the paper). The le Cessie test was run in-sample only. "
         r"Externally the Hosmer--Lemeshow test uses deciles of the frozen predictions with 10 degrees of freedom; the "
         r"calibration belt is that of GiViTI. 300 data sets per cell (1,000 under no misfit).}\label{tab:compare}",
         r"\begin{tabular}{llrrrrr}\toprule setting & departure & procedure & Hosmer--Lemeshow & Spiegelhalter & belt & le Cessie\\ \midrule"]
for _, r in cp.iterrows():
    out_c.append(f"{'external' if r.setting == 'external' else 'in-sample'} & {FL[r.fam]} & {f3(r['procedure (any part)'])} & "
                 f"{f3(r['Hosmer-Lemeshow'])} & {f3(r['Spiegelhalter'])} & {f3(r['calibration belt'])} & {f3(r['le Cessie-van Houwelingen'])}\\\\")
out_c += [r"\bottomrule\end{tabular}\end{table}", ""]
mu = pd.read_csv(os.path.join(OA, "MULTIPLICITY.csv"), keep_default_na=False, na_values=["NA", ""])
mn = mu[mu.fam != "null"].groupby(["setting", "n", "rule"])[["right", "own"]].mean().reset_index()
m0 = mu[mu.fam == "null"].set_index(["setting", "n", "rule"])["wrong"]
out_c += [r"\begin{table}[ht]\centering\small",
          r"\caption{Ways of naming the parts, on the same data sets as the main simulation: the closed procedure, naming "
          r"every group whose own test has $p\le0.05$ (naive), and Bonferroni and Holm corrections over the groups' own "
          r"tests. Error: rate of naming any part under no misfit; right action and own part: averages over the departure "
          r"cells of the main simulation.}\label{tab:multiplicity}",
          r"\begin{tabular}{llrrrrrrrrrrrr}\toprule & & \multicolumn{4}{c}{error under no misfit} & \multicolumn{4}{c}{right action} & \multicolumn{4}{c}{own part named}\\",
          r"\cmidrule(lr){3-6}\cmidrule(lr){7-10}\cmidrule(lr){11-14} setting & $n$ & clos. & naive & Bonf. & Holm & clos. & naive & Bonf. & Holm & clos. & naive & Bonf. & Holm\\ \midrule"]
for (s, n), x in mn.groupby(["setting", "n"]):
    g = x.set_index("rule")
    cells = [f"{m0[(s, n, r)]:.3f}" for r in ["closure", "naive", "Bonferroni", "Holm"]] + \
            [f"{g.loc[r, 'right']:.3f}" for r in ["closure", "naive", "Bonferroni", "Holm"]] + \
            [f"{g.loc[r, 'own']:.3f}" for r in ["closure", "naive", "Bonferroni", "Holm"]]
    out_c.append(f"{'external' if s == 'external' else 'in-sample'} & {n:,} & ".replace(",", "{,}") + " & ".join(cells) + r"\\")
out_c += [r"\bottomrule\end{tabular}\end{table}", ""]
sm = pd.read_csv(os.path.join(OA, "SIM_SMALL.csv"), keep_default_na=False, na_values=["NA", ""])
out_c += [r"\begin{table}[ht]\centering\small",
          r"\caption{Small samples, $n=200$ and $300$ (design of the main simulation): rate of naming each part and the "
          r"right-action rate. 300 data sets per cell (1,000 under no misfit).}\label{tab:small}",
          r"\begin{tabular}{llrrrrrr}\toprule setting & departure & $n$ & I & S & L & C & right\\ \midrule"]
for _, r in sm.iterrows():
    out_c.append(f"{'external' if r.setting == 'external' else 'in-sample'} & {FL[r.fam]} & {int(r.n)} & " +
                 " & ".join(f3(r.get(g, np.nan)) for g in ["INTERCEPT", "SLOPE", "LINK", "COV"]) + f" & {f3(r.right_action)}\\\\")
out_c += [r"\bottomrule\end{tabular}\end{table}", ""]
open(os.path.join(HERE, "tab_appendix_b.tex"), "a", encoding="utf-8", newline="\n").write("\n".join(out_c))

# ---- the robust option (CORRUPT_ROBUST.csv: corrupted records; ROBUST_POWER.csv: main departures, same data sets)
cr = pd.read_csv(os.path.join(LR, "CORRUPT_ROBUST.csv"), keep_default_na=False, na_values=["NA", ""])
cr["block"] = cr["block"].astype(str); cr["robust"] = B(cr["robust"])
rp = pd.read_csv(os.path.join(OA, "ROBUST_POWER.csv"), keep_default_na=False, na_values=["NA", ""])
out_r = [r"\begin{table}[ht]\centering\small",
         r"\caption{The robust option (bases on normal scores) against the default. Top: corrupted records, correct model, "
         r"in-sample, $B=199$: rate of naming any part (all false), 100 data sets per cell. Bottom: the departures of the main "
         r"simulation at $n=1{,}000$, the same data sets: rate of naming some part and right-action rate.}\label{tab:robustopt}",
         r"\begin{tabular}{lrrrr}\toprule corrupted records & \multicolumn{2}{c}{default} & \multicolumn{2}{c}{robust}\\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5} & any part & \textsc{link} & any part & \textsc{link}\\ \midrule"]
for b, f, lab in CL:
    a = cr[(cr.block == b) & (cr.file == f) & ~cr.robust].iloc[0]; z = cr[(cr.block == b) & (cr.file == f) & cr.robust].iloc[0]
    out_r.append(f"{lab} & {a['any']:.2f} & {a.named_LINK:.2f} & {z['any']:.2f} & {z.named_LINK:.2f}\\\\")
out_r += [r"\midrule main departures & \multicolumn{2}{c}{default} & \multicolumn{2}{c}{robust}\\",
          r"\cmidrule(lr){2-3}\cmidrule(lr){4-5} & some part & right action & some part & right action\\ \midrule"]
for _, r in rp.iterrows():
    out_r.append(f"{'external' if r.setting == 'external' else 'in-sample'}, {FL[r.fam]} & {r.any_def:.3f} & {r.right_def:.3f} & "
                 f"{r.any_rob:.3f} & {r.right_rob:.3f}\\\\")
out_r += [r"\bottomrule\end{tabular}\end{table}", ""]
open(os.path.join(HERE, "tab_appendix_b.tex"), "a", encoding="utf-8", newline="\n").write("\n".join(out_r))

# ---- SUPPORT planning
v = pd.read_csv(os.path.join(HERE, "..", "support_v2.csv"), keep_default_na=False, na_values=["NA"])
v = v[v.detect_delta.notna() & ~v.analysis.str.contains("multiplier")]
P = [r"\begin{table}[ht]\centering\small",
     r"\caption{SUPPORT, external analyses: smallest mean miscalibration $\delta$ and logit-scale slope error $\epsilon$ that the "
     r"single-member tests detect with power .80 at level .05 (Corollary~\ref{cor:power}), computed from the covariates and "
     r"frozen predictions only, with the observed and mean predicted risk.}\label{tab:plan}",
     r"\begin{tabular}{llrrrrr}\toprule setting & model & $n$ & $\delta$ & $\epsilon$ & observed & predicted\\ \midrule"]
for _, r in v.iterrows():
    P.append(f"{r.analysis.replace(', montecarlo', '')} & {r.model} & {int(r.n):,} & {r.detect_delta:.3f} & {r.detect_slope:.3f} & "
             f"{r.observed:.3f} & {r.predicted:.3f}\\\\".replace(",", "{,}", 1) if False else
             f"{r.analysis.replace(', montecarlo', '')} & {r.model} & {int(r.n)} & {r.detect_delta:.3f} & {r.detect_slope:.3f} & "
             f"{r.observed:.3f} & {r.predicted:.3f}\\\\")
P += [r"\bottomrule\end{tabular}\end{table}"]
open(os.path.join(HERE, "tab_plan.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(P) + "\n")
print(R.round(3).to_string(index=False)); print(D.round(3).to_string(index=False)); print(C.round(3).to_string(index=False))
