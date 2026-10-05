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
