"""analyze_vergouwe.py -- head-to-head on the same external-validation data sets: our procedure (closure, and Holm over the
four groups' own tests) against Vergouwe et al.'s closed testing procedure and AIC/BIC choice over the update ladder.
The 'right action' is the update the departure needs: none for no misfit; intercept update; logistic recalibration for a
slope error; flexible recalibration for a curved link; revision for misfit among patients who share a score. Vergouwe's
ladder has no flexible recalibration, so for the link departures it cannot take the right action (reported as such).
Writes VERGOUWE_COMPARE.csv."""
import glob, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
B = lambda s: s.astype(str).str.upper().eq("TRUE")
NEED = {"null": "none", "intercept": "intercept", "slope": "recalibration", "slope_pure": "recalibration",
        "link": "flexible", "link_pure": "flexible", "ushape": "revision", "thresh": "revision", "inter": "revision",
        "cov_pure": "revision"}
PART2ACT = {"INTERCEPT": "intercept", "SLOPE": "recalibration", "LINK": "flexible", "COV": "revision"}
ORDER = ["INTERCEPT", "SLOPE", "LINK", "COV"]

ours = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", ""]) for f in glob.glob(os.path.join(HERE, "rows", "external_*.csv"))], ignore_index=True)
for g in ORDER:
    ours["named_" + g] = B(ours["named_" + g])


def top_action(row, named):
    hits = [g for g in ORDER if named(row, g)]
    return PART2ACT[hits[-1]] if hits else "none"


ours["act_closure"] = ours.apply(lambda r: top_action(r, lambda rr, g: rr["named_" + g]), axis=1)


def holm(row, alpha=0.05):
    p = sorted((row["p_" + g], g) for g in ORDER)
    named = set()
    for k, (pv, g) in enumerate(p):
        if pv <= alpha / (len(p) - k):
            named.add(g)
        else:
            break
    return top_action(row, lambda rr, g: g in named)


ours["act_holm"] = ours.apply(holm, axis=1)
ver = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", ""]) for f in glob.glob(os.path.join(HERE, "vergouwe_rows", "*.csv"))], ignore_index=True)
ver = ver.rename(columns={"vergouwe": "act_vergouwe", "aic": "act_aic", "bic": "act_bic"})
ver["act_vergouwe"] = ver["act_vergouwe"].replace({"original": "none"})
ver["act_aic"] = ver["act_aic"].replace({"original": "none"})
ver["act_bic"] = ver["act_bic"].replace({"original": "none"})
m = ours.merge(ver, on=["fam", "C", "n", "rep"], suffixes=("", "_v"))
assert len(m) == len(ours) == len(ver), (len(m), len(ours), len(ver))
assert (m["events"] == m["events_v"]).all()

rows = []
for (fam, C, n), d in m.groupby(["fam", "C", "n"]):
    need = NEED[fam]
    r = {"fam": fam, "C": C, "n": n, "reps": len(d), "needs": need}
    for meth in ["closure", "holm", "vergouwe", "aic", "bic"]:
        a = d["act_" + meth]
        r["right_" + meth] = (a == need).mean()
        r["any_update_" + meth] = (a != "none").mean()
        r["revise_" + meth] = (a == "revision").mean()
    rows.append(r)
out = pd.DataFrame(rows).sort_values(["n", "fam", "C"])
out.to_csv(os.path.join(HERE, "VERGOUWE_COMPARE.csv"), index=False)
pd.set_option("display.width", 220)
print(out[out.n == 1000][["fam", "C", "needs", "right_closure", "right_holm", "right_vergouwe", "right_aic", "right_bic",
                          "any_update_closure", "any_update_vergouwe"]].round(3).to_string(index=False))
print()
print(out[out.n == 500][["fam", "C", "needs", "right_closure", "right_holm", "right_vergouwe", "right_aic", "right_bic"]].round(3).to_string(index=False))
