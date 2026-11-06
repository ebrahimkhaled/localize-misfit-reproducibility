"""build_sim_tex.py -- assembles localize_sim.tex from its sections, adds the Web Tables of the updated-model reference
(from ../theory/CHECK_UPDATE.csv, ../theory/CHECK_CLUSTER.csv, ../burn_update_*.csv) to tab_sim2_web.tex, and adds the
archive entry to references.bib. Run make_sim2_outputs.py first."""
import io, os
import pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
OA = os.path.join(H, "..")
rd = lambda f: io.open(os.path.join(H, f), encoding="utf-8").read()
t = rd("sim_front.tex") + rd("sim_simdisc.tex").replace("%FOCUSED%", rd("sim_focused.tex")).replace("%APPS%", rd("sim_apps.tex"))
ml = pd.read_csv(os.path.join(OA, "theory", "CHECK_ML.csv")).set_index("setting")
pc = lambda v: "%.0f\\%%" % (100 * v)
t = (t.replace("%%ML0%%", "%.3f" % ml.loc["none", "any"]).replace("%%MLL%%", pc(ml.loc["transport", "LINK"])).replace("%%MLC%%", pc(ml.loc["transport", "COV"]))
      .replace("%%MLI%%", pc(ml.loc["cov", "COV"])))
io.open(os.path.join(H, "localize_sim.tex"), "w", encoding="utf-8", newline="\n").write(t)

u = pd.read_csv(os.path.join(OA, "theory", "CHECK_UPDATE.csv"), keep_default_na=False)
L = [r"\begin{table}[ht]", r"\centering",
     r"\caption{The two references on the same data sets ($n=1{,}000$, 200 data sets per row, Holm's correction): share naming each part and right-action rate. Frozen: every outcome drawn from the frozen model. Updated: each part tested against the model updated for the simpler parts. Transport: $\mathrm{logit}\,\pi=-2+0.8\eta_0$, which a logistic recalibration repairs exactly.}",
     r"\label{tab:updatecheck}", r"\small", r"\begin{tabular}{llrrrrrrrrrr}", r"\toprule",
     r" & & \multicolumn{5}{c}{Frozen} & \multicolumn{5}{c}{Updated} \\",
     r"Departure & Needs & int & slope & link & cov & right & int & slope & link & cov & right \\", r"\midrule"]
for _, q in u.iterrows():
    v = [q.mc_INTERCEPT, q.mc_SLOPE, q.mc_LINK, q.mc_COV, q.right_mc, q.up_INTERCEPT, q.up_SLOPE, q.up_LINK, q.up_COV, q.right_up]
    L.append(" & ".join([q.fam, q.needs.lower()] + ["%.2f" % float(x) for x in v]) + r" \\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
c = pd.read_csv(os.path.join(OA, "theory", "CHECK_CLUSTER.csv")).sort_values(["centres", "tau"])
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{A correct frozen model with a centre random effect of standard deviation $\tau$ on the logit scale ($n=1{,}000$, 400 data sets per row, standard error about 0.01). Share naming each part with the plain reference, and familywise error rate (any) of the plain and the cluster-aware reference, with $\tau$ estimated on each data set.}",
      r"\label{tab:cluster}", r"\small", r"\begin{tabular}{rrrrrrrrr}", r"\toprule",
      r"Centres & $\tau$ & Mean $\hat\tau$ & int & slope & link & cov & Any, plain & Any, cluster-aware \\", r"\midrule"]
for _, q in c.iterrows():
    L.append("%d & %.1f & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\" % (q.centres, q.tau, q.tau_hat, q.plain_INTERCEPT, q.plain_SLOPE, q.plain_LINK, q.plain_COV, q.plain_any, q.clust_any))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
s = pd.read_csv(os.path.join(OA, "burn_update_splits.csv"))
n = pd.read_csv(os.path.join(OA, "burn_update_null.csv"))
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{Burn data, half of the facilities used for development (20 random facility splits; about 496 patients and 74 deaths per validation set): share of splits naming each part. Baseline: the same splits applied to 25 data sets simulated from each model fitted to all patients, so that its form is right and only the estimation error of the development half remains.}",
      r"\label{tab:burnupdate}", r"\small", r"\begin{tabular}{llrrrr}", r"\toprule", r"Model & & int & slope & link & cov \\", r"\midrule"]
for m in ["B0", "B3"]:
    L.append(m + " & observed & " + " & ".join("%.2f" % s[s.model == m]["named_" + g].astype(str).str.upper().eq("TRUE").mean() for g in G) + r" \\")
    L.append(" & baseline & " + " & ".join("%.2f" % n[n.model == m]["named_" + g].astype(str).str.upper().eq("TRUE").mean() for g in G) + r" \\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
z = pd.read_csv(os.path.join(OA, "theory", "CHECK_SIZE_COVER.csv"))
z["o"] = z.part.map({g: i for i, g in enumerate(G)})
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{The size of each part under no misfit (1,000 data sets per row, 200 chance draws, five normal covariates): mean size above the chance level, in percentage points, with its standard error, and the share of 95\% noise intervals that contain zero.}",
      r"\label{tab:sizecover}", r"\small", r"\begin{tabular}{rlrrr}", r"\toprule", r"$n$ & Part & Mean size & Standard error & Contains zero \\", r"\midrule"]
for _, q in z.sort_values(["n", "o"]).iterrows():
    L.append("%s & %s & %.3f & %.3f & %.3f \\\\" % (format(int(q.n), ","), q.part.lower(), q.excess, q.se, q.cover))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
h = pd.read_csv(os.path.join(OA, "SIM_GROUPS_HOLM.csv"), keep_default_na=False, na_values=[""])
h = h[(h.setting == "insample") & (h.fam.isin(["null", "link", "ushape", "thresh", "inter"]))]
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{Checks on the development data with Holm's correction (five normal covariates, parametric bootstrap, $B=199$): share of data sets naming each part. 300 data sets per row, 1,000 under no misfit.}",
      r"\label{tab:insampleholm}", r"\small", r"\begin{tabular}{lrrrrr}", r"\toprule", r"Departure & $n$ & link & cov & any & right action \\", r"\midrule"]
for _, q in h.sort_values(["n", "fam"]).iterrows():
    L.append("%s & %d & %.2f & %.2f & %.2f & %.2f \\\\" % (q.fam, q.n, float(q.LINK), float(q.COV), float(q.any_named), float(q.right_action)))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{A frozen neural network (five hidden units, trained on 3,000 patients), $n=1{,}000$: share of data sets naming each part. No misfit: outcomes drawn from the network's risks. Transport: logit risk $-1+0.8\eta_0$. Omitted interaction: $\eta_0+0.8x_1x_3$.}",
      r"\label{tab:ml}", r"\small", r"\begin{tabular}{lrrrrrr}", r"\toprule", r"Setting & Data sets & int & slope & link & cov & any \\", r"\midrule"]
for s_, n_ in [("none", "no misfit"), ("transport", "transport"), ("cov", "omitted interaction")]:
    q = ml.loc[s_]
    L.append("%s & %d & %.3f & %.3f & %.3f & %.3f & %.3f \\\\" % (n_, q.reps, q.INTERCEPT, q.SLOPE, q.LINK, q.COV, q["any"]))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
a_ = pd.read_csv(os.path.join(OA, "theory", "CHECK_AUDIT.csv"), keep_default_na=False)
NAME = dict(coef="Changed effect", ladder="Exact step of the ladder", outside="Curve outside the spline family", few="Few events, no misfit", dev="Estimated development model")
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{Stress tests of the updated-model reference: share of data sets naming each part. Changed effect: tags give the covariate, the change and $n$; closed testing: share choosing revision. Exact steps of the ladder ($n=1{,}000$): flat, no relation with the score; nearflat, slope 0.1; steep spline, a spline with large coefficients, for which the flexible recalibration is enough; extreme transport, $-3+0.3\eta_0$. Curves outside the spline family: tag gives the shape and $n$. Few events: $n$ and prevalence. Estimated development model: design and development size, $n=1{,}000$ for validation.}",
      r"\label{tab:audit}", r"\small", r"\begin{tabular}{llrrrrrr}", r"\toprule", r"Study & Setting & int & slope & link & cov & any & closed testing \\", r"\midrule"]
for e_ in ["coef", "ladder", "outside", "few", "dev"]:
    for _, q in a_[a_.exp == e_].iterrows():
        ct = "" if q.ctp_revision in ("", "NA") else "%.2f" % float(q.ctp_revision)
        L.append("%s & %s & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\" % (NAME[e_], str(q.tag).replace("_", " "), float(q.INTERCEPT), float(q.SLOPE), float(q.LINK), float(q.COV), float(q["any"]), ct))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
st = pd.read_csv(os.path.join(OA, "theory", "CHECK_SIZE_TRUE.csv"))
L += [r"\begin{table}[ht]", r"\centering",
      r"\caption{The \textsc{cov} size against its true value, omitted interaction of strength $c$ (100 data sets per row, 60 chance draws): true mean absolute difference between the risks and the best spline of the score, mean size above chance, and mean noise interval, in percentage points.}",
      r"\label{tab:sizetrue}", r"\small", r"\begin{tabular}{rrrrrr}", r"\toprule", r"$c$ & $n$ & True & Size above chance & Lower & Upper \\", r"\midrule"]
for _, q in st.iterrows():
    L.append("%.1f & %s & %.2f & %.2f & %.2f & %.2f \\\\" % (q.c, format(int(q.n), ","), q.true, q.excess, q.lo, q.hi))
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
w = rd("tab_sim2_web.tex")
io.open(os.path.join(H, "tab_web_update.tex"), "w", encoding="utf-8", newline="\n").write((w + "\n" + "\n".join(L)).replace(r"\small", r"\footnotesize"))

b = rd("references.bib")
if "localizearchive2026" not in b:
    b = b.rstrip() + "\n\n" + """@misc{localizearchive2026,
  author = {Ebrahim, Ebrahim Khaled},
  title = {Localization of misfit in logistic risk models: {R} functions, simulation code, per-data-set results and scripts},
  year = {2026},
  howpublished = {Zenodo},
  note = {Dataset and software archive},
  doi = {10.5281/zenodo.23192535}
}
"""
    io.open(os.path.join(H, "references.bib"), "w", encoding="utf-8", newline="\n").write(b)
print("built", len(t.splitlines()), "lines")
