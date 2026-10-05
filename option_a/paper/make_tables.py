"""make_tables.py -- writes the three tables of the main text: tab_decision.tex (definitions), tab_sim.tex (from
../SIM_GROUPS.csv) and tab_support.tex (from ../support_localize.csv and ../support_v2.csv). No number is typed by hand."""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(HERE, "..")
SC = lambda s: r"\textsc{" + s.lower() + "}"

# ---- Table 1: the decision table
dec = [
    ("none", "no misfit detected", "", "keep the model; the verdict is limited by the sample size (Corollary~\\sref{cor:power})"),
    (SC("INTERCEPT"), "overall risk too high or too low", "different prevalence, case mix or period", "update the intercept"),
    (SC("SLOPE"), "risks too extreme or too modest", "overfitting at development; a different spread of risk",
     "logistic recalibration $a+b\\eta_0$; shrink future models"),
    (SC("LINK"), "curved relation between score and risk", "asymmetric tails; a transformation missing from the score",
     "flexible recalibration $f(\\eta_0)$ or another link"),
    (SC("COV"), "misfit among patients who share a score", "omitted interaction, curved covariate, an effect that differs in the new population",
     "revise the model; recalibration cannot repair it"),
]
L = [r"\begin{table}", r"\caption{Reading a verdict: act on the highest part named. In-sample only \textsc{link} and \textsc{cov} "
     r"can be named (Section~\ref{s:parts}).}", r"\label{t:decision}", r"\begin{center}\footnotesize",
     r"\begin{tabular}{p{1.5cm}p{3.6cm}p{4.2cm}p{4.6cm}}", r"\Hline",
     r"Highest part named & What it means & Typical cause & Update \\ \hline"]
L += [" & ".join(r) + r" \\" for r in dec]
L += [r"\hline", r"\end{tabular}", r"\end{center}", r"\end{table}"]
open(os.path.join(HERE, "tab_decision.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")

# ---- Table 2: simulation, n = 1000
s = pd.read_csv(os.path.join(UP, "SIM_GROUPS.csv"), keep_default_na=False, na_values=[""])
s = s[s.n == 1000]
FAM = [("null", 0.0, "none", "none"), ("intercept", 0.4, "intercept shift $+0.4$", SC("INTERCEPT")),
       ("slope", 0.6, "slope $0.6$", SC("SLOPE")), ("link", 0.25, "curved link $0.25(\\eta_0^2-1)$", SC("LINK")),
       ("ushape", 0.8, "U-shape $0.8(x_1^2-1)$", SC("COV")), ("thresh", 4.0, "threshold $4\\max(x_2-1,0)$", SC("COV")),
       ("inter", 1.0, "interaction $x_1x_3$", SC("COV"))]
f2 = lambda v: "--" if pd.isna(v) else f"{v:.2f}"
L = [r"\begin{table}", r"\caption{Simulation, $n=1{,}000$: proportion of data sets in which each part was named, and the "
     r"right-action rate (highest part named $=$ the part the departure needs; under no misfit, nothing named). External "
     r"validation: Monte Carlo reference, $M=999$; in-sample: parametric bootstrap, $B=199$. 300 data sets per cell, 1,000 "
     r"under no misfit.}", r"\label{t:sim}", r"\begin{center}\footnotesize", r"\begin{tabular}{llrrrrrrrr}", r"\Hline",
     r" & & \multicolumn{5}{c}{External validation} & \multicolumn{3}{c}{In-sample} \\",
     r"Departure & Needs & I & S & L & C & Right & L & C & Right \\ \hline"]
for fam, C, lab, need in FAM:
    e = s[(s.setting == "external") & (s.fam == fam) & ((s.C - C).abs() < 1e-9)]
    i = s[(s.setting == "insample") & (s.fam == fam) & ((s.C - C).abs() < 1e-9)]
    ec = [f2(e[g].iloc[0]) if len(e) else "--" for g in ("INTERCEPT", "SLOPE", "LINK", "COV")] + [f2(e.right_action.iloc[0]) if len(e) else "--"]
    ic = [f2(i[g].iloc[0]) if len(i) else "--" for g in ("LINK", "COV")] + [f2(i.right_action.iloc[0]) if len(i) else "--"]
    L.append(" & ".join([lab, need] + ec + ic) + r" \\")
L += [r"\hline", r"\end{tabular}", r"\end{center}",
      r"\vspace{10pt}\par\noindent\footnotesize{I, S, L, C: \textsc{intercept}, \textsc{slope}, \textsc{link}, \textsc{cov}. Dashes: the part cannot be "
      r"seen in-sample, so the departure was run externally only.}", r"\end{table}"]
open(os.path.join(HERE, "tab_sim.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")

# ---- Table 3: SUPPORT
v2p = os.path.join(UP, "support_v2.csv")
if os.path.exists(v2p):
    v = pd.read_csv(v2p, keep_default_na=False, na_values=["NA"])
    ## the split-first update-and-retest rows (support_split_first.R): verdict on half A, updates fitted on A, retests on B
    sf = pd.read_csv(os.path.join(UP, "support_split_first.csv"), keep_default_na=False, na_values=["NA"])
    sf["analysis"] = ["split: verdict on A" if s == "verdict on A" else "split: " + r for s, r in zip(sf.stage, sf.rung)]
    sf["prescribed"] = sf.stage.str.contains("prescribed")
    v = pd.concat([v, sf], ignore_index=True)
    keep = ["in-sample", "in-sample, de-aliased", "external random, montecarlo",
            "external transported, montecarlo", "split: verdict on A",
            "split: intercept update", "split: logistic recalibration", "split: flexible recalibration", "split: model revision"]
    LAB = {"in-sample": "in-sample, all patients", "in-sample, de-aliased": "\\quad de-aliased link group",
           "external random, montecarlo": "random split",
           "external transported, montecarlo": "transported to chronic diagnoses",
           "split: verdict on A": "\\quad verdict on half A",
           "split: intercept update": "\\quad half B after intercept update$^a$",
           "split: logistic recalibration": "\\quad half B after logistic recalibration$^a$",
           "split: flexible recalibration": "\\quad half B after flexible recalibration$^a$",
           "split: model revision": "\\quad half B after re-estimation$^a$"}
    fp = lambda x: "--" if x == "" or pd.isna(x) else (f"{float(x):.3f}" if float(x) >= .001 else "$<$0.001")
    L = [r"\begin{table}", r"\caption{SUPPORT: parts named ($\alpha=0.05$) and the $p$-value of each part's own test, for the "
         r"linear model M0 and the spline model M2, in-sample ($n=8{,}873$), on a random half after development on the "
         r"other half ($n=4{,}437$), and transported from the acute to the chronic diagnoses ($n=4{,}130$), with the split-first "
         r"update-and-retest: verdict on a random half A of the chronic patients, updates fitted on A, retests on the other half B. "
         r"Monte Carlo (external) or bootstrap (in-sample) reference throughout. A part can have $p\le0.05$ and not be named: "
         r"naming also needs every larger set of parts to reject.}", r"\label{t:support}",
         r"\begin{center}\footnotesize", r"\begin{tabular}{llp{3.2cm}rrrr}", r"\Hline",
         r"Setting & Model & Parts named & I & S & L & C \\ \hline"]
    for a in keep:
        for m in ["M0", "M2"]:
            r = v[(v.analysis == a) & (v.model == m)]
            if not len(r):
                continue
            r = r.iloc[0]
            named = " + ".join(SC(x) for x in str(r.named).split(" + ")) if r.named != "none" else "none"
            if r.get("prescribed", False) is True:
                named += "$^b$"
            L.append(" & ".join([LAB[a] if m == "M0" else "", m, named] +
                                [fp(r.get(g, "")) for g in ("INTERCEPT", "SLOPE", "LINK", "COV")]) + r" \\")
    nn = v[(v.model == "neural network") & ~v.analysis.str.contains("multiplier")]
    for _, r in nn.iterrows():
        named = " + ".join(SC(x) for x in str(r.named).split(" + ")) if r.named != "none" else "none"
        lab = "neural network, " + ("random split" if "random" in r.analysis else "transported") + \
              (", multiplier" if "multiplier" in r.analysis else "")
        L.append(" & ".join([lab, "NN", named] + [fp(r.get(g, "")) for g in ("INTERCEPT", "SLOPE", "LINK", "COV")]) + r" \\")
    L += [r"\hline", r"\end{tabular}", r"\end{center}",
          r"\vspace{10pt}\par\noindent\footnotesize{$^a$Tested on the 2,065 chronic patients of half B. $^b$The update that the verdict on A prescribes. I, S, L, C: $p$-values of the "
          r"single-group tests of \textsc{intercept}, \textsc{slope}, \textsc{link} and \textsc{cov}. NN: a neural network "
          r"with five hidden units fitted to the same development data.}", r"\end{table}"]
    open(os.path.join(HERE, "tab_support.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("tables written; support_v2:", os.path.exists(v2p))
