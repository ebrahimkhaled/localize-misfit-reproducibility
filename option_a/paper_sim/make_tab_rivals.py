"""make_tab_rivals.py -- the focused comparison of the SiM version: how often each procedure takes the update the departure
needs, on the same external-validation data sets (n = 1000): localization with the updated-model reference and Holm's
correction (../rows_update, run_rivals_update.R), the closed testing procedure of Vergouwe et al. (2017), and AIC and BIC
(../VERGOUWE_COMPARE.csv). No number is typed by hand. Writes tab_rivals.tex and prints the named rates."""
import glob, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OA = os.path.join(HERE, "..")
c = pd.read_csv(os.path.join(OA, "VERGOUWE_COMPARE.csv"), keep_default_na=False, na_values=["NA", ""])
u = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in glob.glob(os.path.join(OA, "rows_update", "*.csv"))])
dv = pd.read_csv(os.path.join(OA, "RIVALS_DAVIS.csv"), keep_default_na=False)
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
for g in G:
    u["named_" + g] = u["named_" + g].astype(str).str.upper().eq("TRUE")
u["top"] = u[["named_" + g for g in G]].apply(lambda r: max([i + 1 for i, v in enumerate(r) if v], default=0), axis=1)
ROWS = [("null", 0.0, "none", "keep the model", 0), ("intercept", 0.4, "intercept shift $+0.4$", "intercept update", 1),
        ("slope", 0.6, "slope $0.6$", "logistic recalibration", 2), ("link", 0.25, "curved link $0.25(\\eta_0^2-1)$", "flexible recalibration", 3),
        ("ushape", 0.8, "U-shape $0.8(x_1^2-1)$", "revision", 4), ("thresh", 4.0, "threshold $4\\max(x_2-1,0)$", "revision", 4),
        ("inter", 1.0, "interaction $x_1x_3$", "revision", 4)]
f2 = lambda v: f"{v:.2f}"
L = [r"\begin{table}[t]", r"\centering",
     r"\caption{Focused comparison on the same data sets ($n=1{,}000$, five normal covariates): the share of data sets in "
     r"which each procedure chose the update the departure needs. Closed testing has no flexible step, so it cannot take the "
     r"right action for a curved link; closed, flexible is our extension with that step. BIC: as AIC, with the penalty of the Bayesian information criterion. 300 data sets "
     r"per row, 1,000 under no misfit.}",
     r"\label{t:rivals}", r"\small", r"\begin{tabular}{llrrrrrr}", r"\toprule",
     r"Departure & Update needed & Localization & Closed testing & Closed, flexible & Brier score & AIC & BIC \\", r"\midrule"]
for n in (1000, 500):
    print("n", n)
    for fam, C, lab, need, rung in ROWS:
        g = u[(u.fam == fam) & ((u.C - C).abs() < 1e-9) & (u.n == n)]
        r = c[(c.fam == fam) & ((c.C - C).abs() < 1e-9) & (c.n == n)].iloc[0]
        right = (g.top == rung).mean()
        print(f"  {fam:10s} reps {len(g)} right {right:.3f} over {(g.top > rung).mean():.3f} any {(g.top > 0).mean():.3f} named",
              {x: round(g['named_' + x].mean(), 3) for x in G}, "holm-mc", round(r.right_holm, 3), "ctp", round(r.right_vergouwe, 3))
        if n == 1000:
            dd = dv[(dv.fam == fam) & ((dv.C - C).abs() < 1e-9)].iloc[0]
            L.append(" & ".join([lab, need, f2(right), f2(r.right_vergouwe), f2(dd.right_ctp_flex), f2(dd.right_davis), f2(r.right_aic), f2(r.right_bic)]) + r" \\")
L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(os.path.join(HERE, "tab_rivals.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("\n".join(L[10:-3]))
