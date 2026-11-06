"""make_figs_py.py -- the figures of the paper in one style (matplotlib), from the stored result files (no typed numbers):
fig_concept: the four parts on the ladder of updates (illustrative curves); fig_sim2: the main simulation;
fig_verdicts: SUPPORT verdicts and sizes; fig_burn: burn part sizes with noise intervals. Each as .pdf and .png."""
import os
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

H = os.path.dirname(os.path.abspath(__file__)); OA = os.path.join(H, "..")
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8983", "#e6e5e0"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
RAMP = ["#dbe9f9", "#a9cbf1", "#6ea6e6", "#2a78d6", "#1c5aa6"]   # one hue, light to dark
mpl.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.edgecolor": MUTED, "axes.linewidth": .6, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": .6, "ytick.major.width": .6, "xtick.major.size": 2.5,
                     "ytick.major.size": 0, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.dpi": 200})
save = lambda f, n: (f.savefig(os.path.join(H, n + ".pdf"), bbox_inches="tight", pad_inches=.03),
                     f.savefig(os.path.join(H, n + ".png"), bbox_inches="tight", pad_inches=.03), plt.close(f))
lg = lambda p: np.log(p / (1 - p)); ex = lambda e: 1 / (1 + np.exp(-e))


def concept():
    f = plt.figure(figsize=(7.2, 2.95))
    p = np.concatenate([np.geomspace(1e-5, .01, 80), np.linspace(.01, .99, 300)[1:-1], 1 - np.geomspace(.01, 1e-5, 80)]); e = lg(p)
    ## link: the misfit g = (true - predicted) / (p(1 - p)) of equation (1) is an even function of the score with its weighted mean
    ## removed, so its intercept and slope parts are zero (predicted risks spread evenly over 0 to 1)
    q = np.linspace(.0005, .9995, 2000); sh = lambda v: 1 / np.cosh(v / 1.5); w = q * (1 - q)
    gl = lambda v: -1.1 * (sh(v) - np.sum(w * sh(lg(q))) / np.sum(w))
    panels = [("Intercept", "overall risk too high", [(ex(e - .9), BLUE, "-")], "Update the intercept"),
              ("Slope", "risks too extreme", [(ex(.55 * e), BLUE, "-")], "Logistic recalibration"),
              ("Link", "curved relation with risk", [(p + p * (1 - p) * gl(e), BLUE, "-")], "Flexible recalibration"),
              ("Covariates", "equal predicted risk, unequal true risk", [(p + p * (1 - p), ORANGE, "-"), (p - p * (1 - p), AQUA, "-"), ((2 * p - p ** 2 + p ** 2) / 2, BLUE, "-")], "Revise the model")]
    L, W, G = .055, .205, .036
    for i, (nm, sub, curves, up) in enumerate(panels):
        ax = f.add_axes([L + i * (W + G), .36, W, .47])
        ax.plot([0, 1], [0, 1], color=MUTED, lw=.8, ls=(0, (3, 2)), zorder=1)
        for y, c, ls in curves:
            ax.plot(p, y, color=c, lw=2, ls=ls, solid_capstyle="round", zorder=3)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal", "box")
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1]); ax.tick_params(labelsize=7.5, length=2)
        ax.set_xlabel("Predicted risk", labelpad=-5, fontsize=7.5)
        if i == 0: ax.set_ylabel("Observed risk", labelpad=-4, fontsize=7.5)
        ax.text(0.0, 1.20, nm, transform=ax.transAxes, fontweight="bold", fontsize=9.5, color=INK)
        ax.text(0.0, 1.07, sub, transform=ax.transAxes, fontsize=7.4, color=INK2)
        if i == 3:
            ax.text(.05, .84, "group A", color=INK2, fontsize=7); ax.text(.60, .13, "group B", color=INK2, fontsize=7)
            ax.annotate("average", xy=(.72, .72), xytext=(.76, .47), fontsize=7, color=INK2, arrowprops=dict(arrowstyle="-", lw=.5, color=MUTED))
    ax = f.add_axes([0, 0, 1, .22]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    for i, (nm, sub, c, up) in enumerate(panels):
        x0 = L + i * (W + G) - .012; x1 = x0 + W + G - .004; y0 = .12; y1 = .50 + .0 * i
        ax.add_patch(Polygon([(x0, y0), (x1 - .018, y0), (x1, (y0 + y1 + .34) / 2 - .17), (x1 - .018, y1 + .34 - .34), (x0, y1)], closed=True, fc=RAMP[i + 1], ec="white", lw=1.5))
        ax.text((x0 + x1) / 2 - .006, (y0 + y1) / 2, up, ha="center", va="center", fontsize=8, fontweight="bold", color="white" if i >= 2 else INK)
    ax.text(.5, .80, "the ladder of updates, from the simplest to the most demanding", ha="center", fontsize=7.6, color=INK2, style="italic")
    save(f, "fig_concept")


def sim():
    t = pd.read_csv(os.path.join(OA, "SIM2_DAVIS.csv"), keep_default_na=False)
    c = pd.read_csv(os.path.join(OA, "SIM2_CONSEQ.csv"), keep_default_na=False)
    c = c[(c.metric == "err") & (c.design == "both")].set_index("dep")
    cell = pd.read_csv(os.path.join(OA, "SIM2_DAVIS_CELLS.csv"))
    ORDER = ["null", "intercept", "slope", "transport", "link", "effect", "inter", "ushape"]
    LAB = ["No misfit", "Intercept +0.4", "Slope 0.7", "Transport", "Curved link", "Changed effect", "Interaction", "U-shape"]
    t = t.set_index("dep").loc[ORDER]
    RULES = [("loc", "Localization", BLUE, "o"), ("vergouwe", "Closed testing", ORANGE, "^"), ("davis", "Brier-score procedure", AQUA, "s"), ("aic", "AIC", YELLOW, "D")]
    f = plt.figure(figsize=(7.2, 3.6))
    gs = f.add_gridspec(1, 3, width_ratios=[.72, 1.02, 1.12], wspace=.13, left=.09, right=.995, top=.87, bottom=.20)
    # A: false updates in the twelve designs
    a = f.add_subplot(gs[0]); rng = np.random.default_rng(1)
    a.axvspan(0, .05, color=RAMP[0], alpha=.55, lw=0); a.axvline(.05, color=BLUE, lw=.8)
    for i, (k, nm, col, mk) in enumerate(RULES):
        v = cell["false_" + k].values
        a.scatter(v, i + rng.uniform(-.16, .16, len(v)), s=17, color=col, marker=mk, ec="white", lw=.4, zorder=3)
    a.set_yticks(range(4)); a.set_yticklabels(["Localization", "Closed\ntesting", "Brier\nscore", "AIC"], fontsize=7.6); a.invert_yaxis()
    a.set_xlim(0, max(.4, cell[[c_ for c_ in cell if c_.startswith("false_")]].values.max() * 1.05)); a.set_xlabel("Correct model wrongly updated")
    a.xaxis.grid(True, color=GRID, lw=.6); a.set_axisbelow(True); a.spines["left"].set_visible(False)
    a.set_title("A  Error in 12 designs", loc="left", fontsize=9, fontweight="bold", color=INK)
    # B: right-update tile
    b = f.add_subplot(gs[1]); M = np.array([[t.loc[d, "right_" + k] for k, *_ in RULES] for d in ORDER])
    cmap = mpl.colors.LinearSegmentedColormap.from_list("b", ["#f4f8fd"] + RAMP)
    b.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            bold = M[i, j] >= M[i].max() - 1e-9
            b.text(j, i, f"{M[i, j]:.2f}".replace("0.", ".") if M[i, j] < 1 else "1.00", ha="center", va="center", fontsize=7.8,
                   color="white" if M[i, j] > .6 else INK, fontweight="bold" if bold else "normal")
    b.set_xticks(range(4)); b.set_xticklabels(["Local.", "Closed\ntesting", "Brier\nscore", "AIC"], fontsize=7.4); b.set_yticks(range(8)); b.set_yticklabels(LAB, fontsize=7.8)
    b.set_xticks(np.arange(-.5, 4), minor=True); b.set_yticks(np.arange(-.5, 8), minor=True); b.grid(which="minor", color="white", lw=2); b.tick_params(which="both", length=0)
    for s in b.spines.values(): s.set_visible(False)
    b.set_title("B  Share choosing the right update", loc="left", fontsize=9, fontweight="bold", color=INK)
    # C: risk error
    d = f.add_subplot(gs[2]); y = np.arange(8)
    fixed = [("always_recalibration", "Always recalibrate", MUTED, "|"), ("always_revision_cov", "Always revise", INK2, "x")]
    best = c.loc[ORDER, "best"].astype(float).values
    d.barh(y, best, color=RAMP[0], height=.72, lw=0, zorder=1, label="Best possible")
    for k, nm, col, mk in fixed:
        d.scatter(c.loc[ORDER, k].astype(float), y, marker=mk, s=26, color=col, lw=1.1, zorder=3, label=nm)
    for k, nm, col, mk in RULES:
        d.scatter(t["err_" + k].astype(float), y, s=20, color=col, marker=mk, ec="white", lw=.4, zorder=4 if k == "loc" else 3)
    d.set_ylim(7.5, -.5); d.set_yticks(y); d.set_yticklabels([]); d.tick_params(axis="y", length=0); d.set_xlim(0, 6)
    d.xaxis.grid(True, color=GRID, lw=.6); d.set_axisbelow(True); d.spines["left"].set_visible(False)
    d.set_xlabel("Error in predicted risk (percentage points)")
    d.set_title("C  Risk error for new patients", loc="left", fontsize=9, fontweight="bold", color=INK)
    hs = [plt.Line2D([], [], marker=mk, ls="", color=col, mec="white", mew=.4, ms=5.5, label=nm) for k, nm, col, mk in RULES]
    hs += [plt.Line2D([], [], marker=mk, ls="", color=col, mew=1.1, ms=5.5, label=nm) for k, nm, col, mk in fixed]
    hs += [mpl.patches.Patch(color=RAMP[0], label="Best possible")]
    f.legend(handles=hs, loc="lower center", ncol=7, fontsize=7.3, handletextpad=.1, columnspacing=.9, bbox_to_anchor=(.53, -.015))
    save(f, "fig_sim2")


def verdicts():
    a = pd.read_csv(os.path.join(OA, "apps_update.csv"), keep_default_na=False); s = pd.read_csv(os.path.join(OA, "size_apps.csv"), keep_default_na=False)
    G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
    MAP = [("random split", "random split", "Similar patients"), ("transported", "transported, all chronic", "Transported"),
           ("transported, half A", "transported, half A", "Transported, half A"), ("half B, intercept update", "half B after intercept update", "Half B: intercept update"),
           ("half B, logistic recalibration", "half B after logistic recalibration", "Half B: logistic recalibration"),
           ("half B, flexible recalibration", "half B after flexible recalibration", "Half B: flexible recalibration")]
    f, axs = plt.subplots(1, 2, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=.07, left=.215, right=.995, top=.84, bottom=.03))
    for ax, m in zip(axs, ["M0", "M2"]):
        ax.set_xlim(0, 4); ax.set_ylim(6, 0); ax.axis("off")
        ax.text(2, -.72, "Model " + m, ha="center", fontweight="bold", fontsize=9.5, color=INK)
        for j, g in enumerate(G):
            ax.text(j + .5, -.2, g.lower(), ha="center", fontsize=8, style="italic", color=INK2)
        for i, (sa, ss, lab) in enumerate(MAP):
            r = a[(a.data == "SUPPORT") & (a.model == m) & (a.setting == sa)].iloc[0]; nm = r.named.split(" + ")
            if m == "M0": ax.text(-.12, i + .5, lab, ha="right", va="center", fontsize=8, color=INK)
            for j, g in enumerate(G):
                q = s[(s.data == "SUPPORT") & (s.model == m) & (s.setting == ss) & (s.part == g)].iloc[0]; v = max(float(q.excess), 0); on = g in nm
                ax.add_patch(FancyBboxPatch((j + .05, i + .07), .9, .86, boxstyle="round,pad=0,rounding_size=.1", fc=(ORANGE if g == "COV" else BLUE) if on else "#f1f0ec", ec="none", mutation_aspect=1))
                ax.text(j + .5, i + .52, f"{v:.1f}", ha="center", va="center", fontsize=8.4, color="white" if on else INK2, fontweight="bold" if on else "normal")
    save(f, "fig_verdicts")


def burn():
    bp = pd.read_csv(os.path.join(OA, "burn_iecv_baseline_p.csv")); it = pd.read_csv(os.path.join(OA, "burn_iecv.csv"), keep_default_na=False).set_index("model")
    G = ["INTERCEPT", "SLOPE", "LINK", "COV"]; rng = np.random.default_rng(3)
    f, axs = plt.subplots(1, 2, figsize=(6.4, 2.0), sharey=True, sharex=True, gridspec_kw=dict(wspace=.06, left=.1, right=.99, top=.88, bottom=.21))
    for ax, m, col, mk, nm in zip(axs, ["B0", "B3"], [ORANGE, BLUE], ["o", "s"], ["B0, main effects", "B3, revised"]):
        g = bp[bp.model == m]
        ax.axvline(-np.log10(.05), color=MUTED, lw=.8, ls=(0, (3, 2)))
        for i, k in enumerate(G):
            v = -np.log10(g["p_" + k].values)
            ax.scatter(v + rng.uniform(-.015, .015, len(v)), i + rng.uniform(-.2, .2, len(v)), s=4, color="#c9c8c2", lw=0, zorder=2)
            ob = -np.log10(float(it.loc[m, "plain_" + k]))   # the reference of the grey points: no centre effect
            ax.scatter(ob, i, s=46, color=col, marker=mk, ec="white", lw=.7, zorder=4)
        ax.set_yticks(range(4)); ax.set_yticklabels([x.lower() for x in G], style="italic"); ax.set_ylim(3.5, -.5); ax.set_xlim(0, 3.1)
        ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["1", "0.1", "0.01", "0.001"])
        ax.xaxis.grid(True, color=GRID, lw=.6); ax.set_axisbelow(True); ax.spines["left"].set_visible(False)
        ax.set_title(nm, loc="left", fontsize=9, fontweight="bold", color=INK); ax.set_xlabel("$p$-value of the part")
    save(f, "fig_burn")


concept(); verdicts(); burn()
if os.path.exists(os.path.join(OA, "SIM2_DAVIS.csv")):
    sim()
print("figures written")
