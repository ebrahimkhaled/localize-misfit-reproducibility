"""make_sim4_table.py -- Web Table of the simulation with estimated frozen models, written from SIM4_ERROR.csv and SIM4_COV.csv."""
import io, os
import pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.dirname(H)
e = pd.read_csv(os.path.join(A, "SIM4_ERROR.csv"), keep_default_na=False); c = pd.read_csv(os.path.join(A, "SIM4_COV.csv"), keep_default_na=False)
m = e.merge(c, on=["dep", "ndev", "events"], how="left")
NAME = dict(null="Correct form", transport="Transport", link="Curved link", inter="Interaction")
lines = [r"\begin{table}[p]\centering\footnotesize\setlength{\tabcolsep}{4pt}",
         r"\caption{A frozen model with estimated coefficients: a main-effects model fitted on a development sample ($n_d$ patients) from the baseline population, then validated at 20\% prevalence with 100, 200 or 400 expected events; both covariate designs, 200 data sets each, 9,600 in all. \textsc{cov}: share of data sets naming it. Gain: when \textsc{cov} was named, the error of the update the verdict would otherwise prescribe minus the error of the covariate revision, with its standard error (positive: revision helped). The remaining columns give the mean absolute error of the risks of new patients, in percentage points: the frozen model, following the verdict, closed testing with the covariate revision as its last step, always recalibrating, always revising the covariates, and always refitting every coefficient.}",
         r"\label{tab:estdev}", r"\begin{tabular}{llrrlrrrrrr}", r"\toprule",
         r"Setting & $n_d$ & Events & \textsc{cov} & Gain (SE) & Frozen & Verdict & Closed & Recal. & Revise & Refit \\", r"\midrule"]
last = None
for _, r in m.iterrows():
    if last is not None and r.dep != last:
        lines.append(r"\addlinespace")
    last = r.dep
    rate = r.named / 400
    lines.append(f"{NAME[r.dep]} & {int(r.ndev):,} & {int(r.events)} & {rate:.3f} & ${r.gain:+.2f}$ ({r.se:.2f}) & {r.err_none:.2f} & {r.e_ours:.2f} & {r.e_vercov:.2f} & {r.err_recalibration:.2f} & {r.err_revision_cov:.2f} & {r.err_revision_linear:.2f} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
io.open(os.path.join(H, "tab_sim4.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(lines))
print(len(m), "rows")
