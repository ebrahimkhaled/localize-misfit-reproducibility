"""make_sim2_outputs.py -- the tables of the main simulation, written from ../SIM2_*.csv (no typed numbers):
tab_sim2.tex (main text: risk error for new patients under each rule) and tab_sim2_web.tex (Web Tables: error rate by
cell, right action by design and events, and the rules with the same covariate revision given to every procedure)."""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OA = os.path.join(HERE, "..")
c = pd.read_csv(os.path.join(OA, "SIM2_CONSEQ.csv"), keep_default_na=False)
e = pd.read_csv(os.path.join(OA, "SIM2_ERROR.csv"))
r = pd.read_csv(os.path.join(OA, "SIM2_RIGHT.csv"), keep_default_na=False)
ORDER = ["null", "intercept", "slope", "transport", "link", "effect", "inter", "ushape"]
LAB = dict(null="no misfit", intercept="intercept $+0.4$", slope="slope $0.7$", transport="transport",
           link="curved link", effect="changed effect", inter="interaction", ushape="U-shape")
num = ["ours", "vergouwe", "vergouwe_covrev", "aic", "aic_covrev", "always_none", "always_recalibration", "always_flexible",
       "always_revision_linear", "always_revision_cov", "best"]
for k in num:
    c[k] = c[k].astype(float)
dv = pd.read_csv(os.path.join(OA, "SIM2_DAVIS.csv"), keep_default_na=False).set_index("dep")
for m_, col in [("err", "err_"), ("ll", "ll_")]:
    for k, src in [("davis", "davis"), ("ctp_flex", "ctp_flex")]:
        c.loc[(c.metric == m_) & (c.design == "both") & (c.dep != "ALL"), k] = c[(c.metric == m_) & (c.design == "both") & (c.dep != "ALL")].dep.map(dv[col + src].astype(float))


def block(metric, design, cols, scale=1.0, nd=2):
    g = c[(c.metric == metric) & (c.design == design) & (c.dep != "ALL")].set_index("dep").loc[ORDER]
    lines = []
    for d in ORDER:
        vals = [g.loc[d, k] * scale for k in cols]
        lo = min(vals[:-1])
        s = [("\\textbf{%.*f}" % (nd, v)) if abs(v - lo) < 1e-9 else "%.*f" % (nd, v) for v in vals[:-1]] + ["%.*f" % (nd, vals[-1])]
        lines.append(LAB[d] + " & " + " & ".join(s) + " \\\\")
    m = [g[k].mean() * scale for k in cols]
    lo = min(m[:-1])
    s = [("\\textbf{%.*f}" % (nd, v)) if abs(v - lo) < 1e-9 else "%.*f" % (nd, v) for v in m[:-1]] + ["%.*f" % (nd, m[-1])]
    lines.append("\\midrule\nmean of the eight rows & " + " & ".join(s) + " \\\\")
    return "\n".join(lines)


main = r"""\begin{table}[t]
\centering
\caption{Consequences for new patients: mean absolute error of the predicted risk against the true risk (percentage points; 20,000 new patients per data set), after each rule chooses and fits an update on the validation data. Mean over the twelve designs; 6,000 data sets under no misfit and 2,400 for each departure. Smallest error of the six rules in bold. Best possible: the better of the six candidate updates in each data set, known only to an oracle.}
\label{t:sim2}
\small
\begin{tabular}{lrrrrrrr}
\toprule
Setting & Localization & Closed testing & Brier score & AIC & Always recalibrate & Always revise & Best possible \\
\midrule
""" + block("err", "both", ["ours", "vergouwe", "davis", "aic", "always_recalibration", "always_revision_cov", "best"]) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(os.path.join(HERE, "tab_sim2.tex"), "w", encoding="utf-8", newline="\n").write(main)

web = []
web.append(r"""\begin{table}[ht]
\centering
\caption{Main simulation: familywise error rate under no misfit, by design (500 data sets per row; standard error about 0.01). Localization: the share naming any part. Closed testing and AIC: the share choosing any update.}
\label{tab:sim2err}
\small
\begin{tabular}{lrrrrr}
\toprule
Covariates & Prevalence & Events & Localization & Closed testing & AIC \\
\midrule
""")
for _, q in e.iterrows():
    web.append("%s & %.2f & %d & %.3f & %.3f & %.3f \\\\\n" % (q.design, q.prev, q.events, q.fwer, q.ver, q.aic))
web.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
web.append(r"""\begin{table}[ht]
\centering
\caption{Main simulation: right-action rate by covariate design and expected events, pooled over prevalence. For closed testing and AIC, a linear revision counts as the right action for the two covariate departures.}
\label{tab:sim2right}
\small
\begin{tabular}{llrrrr}
\toprule
Setting & Covariates & Events & Localization & Closed testing & AIC \\
\midrule
""")
r["o"] = r.dep.map({d: i for i, d in enumerate(ORDER)})
for _, q in r.sort_values(["o", "design", "events"]).iterrows():
    web.append("%s & %s & %d & %.2f & %.2f & %.2f \\\\\n" % (LAB[q.dep], q.design, q.events, q["loc"], q.ver, q.aic))
web.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
web.append(r"""\begin{table}[ht]
\centering
\caption{Main simulation: risk error (percentage points) when closed testing and AIC use the same ridge covariate revision as localization in place of their linear revision, and for the other fixed rules. Smallest error of the rules in bold.}
\label{tab:sim2fair}
\small
\begin{tabular}{lrrrrrrr}
\toprule
Setting & Localization & Closed testing & AIC & Never update & Always flexible & Always linear revision & Best possible \\
\midrule
""" + block("err", "both", ["ours", "vergouwe_covrev", "aic_covrev", "always_none", "always_flexible", "always_revision_linear", "best"]) + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[ht]
\centering
\caption{Main simulation: expected log-loss on new patients ($\times100$), same rules as Table~2 of the paper. Smallest of the five rules in bold.}
\label{tab:sim2ll}
\small
\begin{tabular}{lrrrrrr}
\toprule
Setting & Localization & Closed testing & AIC & Always recalibrate & Always revise & Best possible \\
\midrule
""" + block("ll", "both", ["ours", "vergouwe", "aic", "always_recalibration", "always_revision_cov", "best"], scale=100) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
ex = [r"""\\begin{table}[ht]
\\centering
\\caption{Main simulation: the Brier-score procedure of Davis et al.\\ and closed testing with a flexible step: right-action rate, share prescribing too much, risk error (percentage points), and risk error with the covariate revision in place of the linear revision.}
\\label{tab:sim2extra}
\\small
\\begin{tabular}{lrrrrrrrr}
\\toprule
 & \\multicolumn{4}{c}{Brier score} & \\multicolumn{4}{c}{Closed testing, flexible step} \\\\
Setting & right & too much & error & error, cov & right & too much & error & error, cov \\\\
\\midrule
"""]
for d_ in ORDER:
    q = dv.loc[d_].astype(float)
    ex.append("%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n" % (LAB[d_], q.right_davis, q.over_davis, q.err_davis, q.errcov_davis, q.right_ctp_flex, q.over_ctp_flex, q.err_ctp_flex, q.errcov_ctp_flex))
ex.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
ex[0] = ex[0].replace("\\\\", "\\")
web += ex
open(os.path.join(HERE, "tab_sim2_web.tex"), "w", encoding="utf-8", newline="\n").write("".join(web))
print(main)
