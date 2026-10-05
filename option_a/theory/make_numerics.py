"""make_numerics.py -- writes the tables of Section A.9 (numerics_tables.tex) from the result files, so no number in the
appendix is typed by hand: CHECK_T1.csv (Lemma A.1 / Theorem A.1d), SPLINE_APPROX.csv (Remark A.1), DECOMP.csv
(Proposition A.1) and rows_local/ (Theorems A.3-A.4 under local and fixed departures). Also writes LOCAL.csv."""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
OWN = {"null": set(), "slope_pure": {"SLOPE"}, "link_pure": {"LINK"}, "cov_pure": {"COV"}, "combo": {"LINK", "COV"},
       "slope_mod": {"SLOPE"}, "link_mod": {"LINK"}, "cov_mod": {"COV"}}
FAM = {"null": "none", "slope_pure": "slope only", "link_pure": "link only", "cov_pure": "covariate only",
       "combo": "link + covariate",
       "slope_mod": "slope only", "link_mod": "link only", "cov_mod": "covariate only"}
out = []

# ---- Table A.1: what the in-sample fit erases
t1 = pd.read_csv(os.path.join(HERE, "CHECK_T1.csv"))
lab = {"logit": "logistic (null)", "curved": "curved link", "inter": "interaction"}
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{In-sample identities after a logit and a probit fit (200 data sets of size 1,000 per truth; largest "
        r"absolute value over data sets). Sums are $\sum_i(Y_i-\hat p_i)$ and $\sum_i(Y_i-\hat p_i)\hat\eta_i$; $a$ and $b-1$ are "
        r"the intercept and slope error of the logistic recalibration of $\hat p$ on the same data.}\label{tab:t1}",
        r"\begin{tabular}{lrrrrrrrr}\toprule",
        r"& \multicolumn{4}{c}{logit fit} & \multicolumn{4}{c}{probit fit}\\ \cmidrule(lr){2-5}\cmidrule(lr){6-9}",
        r"truth & sum & sum$\cdot\hat\eta$ & $a$ & $b-1$ & sum & sum$\cdot\hat\eta$ & $a$ & $b-1$\\ \midrule"]
for tr in ["logit", "curved", "inter"]:
    d = t1[t1.truth == tr].set_index("stat")["max_abs"]
    cells = [d[f"{m}_{s}"] for m in ("logit", "probit") for s in ("sum", "sum_eta", "recal_a", "recal_b")]
    out.append(lab[tr] + " & " + " & ".join(("%.1e" % v) if v < 1e-3 else ("%.3f" % v if v < 1 else "%.2f" % v)
                                             for v in cells) + r"\\")
out += [r"\bottomrule\end{tabular}\end{table}", ""]

# ---- Table A.2: spline approximation of functions of eta
ap = pd.read_csv(os.path.join(HERE, "SPLINE_APPROX.csv"))
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{Remark~\ref{rem:impl}: relative size $\|(\Pi_{\mathcal G}-\Pi_{\mathcal G_5})b\|_w/\|b-\Pi_{\mathcal G_5}b\|_w$ of the part "
        r"of a COV basis column that the five-degree-of-freedom spline leaves in $\mathcal G$ ($\mathcal G$ approximated by a "
        r"40-degree-of-freedom spline; simulation design, $n=10^5$).}\label{tab:approx}",
        r"\begin{tabular}{lrrr}\toprule COV member & columns & median & largest\\ \midrule"]
for _, r in ap.iterrows():
    out.append(f"{r.member} & {r.columns} & {r.median_rel:.3f} & {r.max_rel:.3f}\\\\")
out += [r"\bottomrule\end{tabular}\end{table}", ""]

# ---- Table A.3: fixed departures, groups named outside the departure's own piece
dc = pd.read_csv(os.path.join(HERE, "DECOMP.csv"))
own_dc = {"slope_p": "SLOPE", "link_p": "LINK", "cov_p": "COV", "slope_pure": "SLOPE", "link_pure": "LINK", "cov_pure": "COV"}
dc = dc[dc.apply(lambda r: r.group != own_dc[r.fam] and (r.named >= .10 or r.member_power >= .10), axis=1)]
out += [r"\begin{table}[ht]\centering\small",
        r"\caption{Fixed-size departures of the simulations (external validation, $n=1{,}000$, the same 100 data sets per "
        r"cell that the procedure analyzed): groups other than the departure's own that were named in at least 10\% of "
        r"data sets or carry a member with power at least .10. Noncentrality: the largest member noncentrality of the "
        r"group computed from the true risks; power: that member's power alone at level .05; named: rate at which the "
        r"procedure named the group.}\label{tab:decomp}",
        r"\begin{tabular}{llrlrrr}\toprule scale & departure & size & group & noncentrality & power & named\\ \midrule"]
for _, r in dc.iterrows():
    out.append(f"{r.study} & {r.fam.replace('_', ' ')} & {r.C:g} & {r.group} & {r.mean_ncp:.2f} & {r.member_power:.2f} & {r.named:.2f}\\\\")
out += [r"\bottomrule\end{tabular}\end{table}", ""]

# ---- Table A.4: local versus fixed departures
files = glob.glob(os.path.join(HERE, "rows_local", "*.csv"))
if files:
    d = pd.concat([pd.read_csv(f, keep_default_na=False) for f in files if "SMOKE" not in f])
    for g in G:
        d["named_" + g] = d["named_" + g].astype(str).str.upper().eq("TRUE")
    rows = []
    for (s, f, tr, n), x in d.groupby(["setting", "fam", "track", "n"]):
        groups = G if s == "external" else ["LINK", "COV"]
        own = [g for g in groups if g in OWN[f]]
        oth = [g for g in groups if g not in OWN[f]]
        r = dict(setting=s, fam=f, track=tr, n=n, reps=len(x),
                 own_all=x[["named_" + g for g in own]].all(axis=1).mean() if own else np.nan,
                 wrong=x[["named_" + g for g in oth]].any(axis=1).mean())
        r["se_wrong"] = np.sqrt(r["wrong"] * (1 - r["wrong"]) / len(x))
        rows.append(r)
    L = pd.DataFrame(rows)
    L.to_csv(os.path.join(HERE, "LOCAL.csv"), index=False)
    for s, ns, fams in (("external", [500, 2000, 8000], ["null", "slope_pure", "link_pure", "cov_pure", "combo"]),
                        ("insample", [500, 2000, 4000], ["null", "link_pure", "cov_pure", "combo"]),
                        ("moderate", [1000, 4000, 16000, 32000], ["slope_mod", "link_mod", "cov_mod"])):
        Ls = L[(L.setting == ("external" if s == "moderate" else s)) & L.fam.isin(fams)]
        out += [r"\begin{table}[ht]\centering\small",
                r"\caption{" + {"external": "External validation, large departures", "insample": "In-sample checking",
                              "moderate": "External validation, moderate departures followed to $n=32{,}000$"}[s] +
                r": rate at which every group with a part of the departure is named (own) and at which any other "
                r"group is named (wrong; nominal bound .05; Monte Carlo standard error at most "
                + f"{Ls.se_wrong.max():.3f}" + r"). Local track: size $C\sqrt{1000/n}$; fixed track: size $C$. "
                + {"external": "400", "insample": "200", "moderate": "400 (200 at $n=32{,}000$)"}[s] + r" data sets per cell.}\label{tab:local" + s[:2] + "}",
                r"\begin{tabular}{ll" + "rr" * len(ns) + r"}\toprule",
                " & & " + " & ".join(r"\multicolumn{2}{c}{$n=" + f"{n:,}".replace(",", "{,}") + "$}" for n in ns) + r"\\",
                " ".join(r"\cmidrule(lr){%d-%d}" % (3 + 2 * i, 4 + 2 * i) for i in range(len(ns))),
                "departure & track & " + " & ".join(["own & wrong"] * len(ns)) + r"\\ \midrule"]
        for f in fams:
            for tr in ["local", "fixed"]:
                x = Ls[(Ls.fam == f) & (Ls.track == tr)]
                if x.empty:
                    continue
                cells = []
                for n in ns:
                    y = x[x.n == n]
                    if y.empty:
                        cells += ["--", "--"]
                    else:
                        y = y.iloc[0]
                        cells += ["--" if np.isnan(y.own_all) else f"{y.own_all:.2f}", f"{y.wrong:.3f}"]
                out.append(f"{FAM[f]} & {tr} & " + " & ".join(cells) + r"\\")
        out += [r"\bottomrule\end{tabular}\end{table}", ""]
    print(L.round(3).to_string(index=False))

open(os.path.join(HERE, "numerics_tables.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(out))
print("written numerics_tables.tex")
