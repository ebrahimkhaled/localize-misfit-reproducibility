# analyze_sim4.py -- the frozen model with estimated coefficients (sim4_rows/). Output SIM4_RATES.csv, SIM4_ERROR.csv, SIM4_COV.csv.
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
files = glob.glob(os.path.join(HERE, "sim4_rows", "*.csv"))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in files], ignore_index=True)
print("rows", len(d), "files", len(files))
G = ["INTERCEPT", "SLOPE", "LINK", "COV"]
for g in G:
    d["named_" + g] = d["named_" + g].astype(str).str.upper().eq("TRUE")
UPD = ["none", "intercept", "recalibration", "flexible", "revision_linear", "revision_cov"]
key = ["dep", "ndev", "events"]
ordr = {"null": 0, "transport": 1, "link": 2, "inter": 3}


def chosen(col, metric="err", remap=None):
    c = d[col].replace(remap) if remap else d[col]
    out = np.empty(len(d))
    for k in c.unique():
        m = (c == k).values
        out[m] = d.loc[m, metric + "_" + k].values
    return out


# the error of the update the verdict prescribes if a COV verdict is ignored
def below_cov():
    top = np.where(d.named_LINK, "flexible", np.where(d.named_SLOPE, "recalibration", np.where(d.named_INTERCEPT, "intercept", "none")))
    return np.select([top == u for u in UPD[:4]], [d["err_" + u].values for u in UPD[:4]])


d["e_ours"] = chosen("loc")
d["e_ours_nocov"] = below_cov()
d["e_ver"] = chosen("vergouwe")
d["e_vercov"] = chosen("vergouwe", remap={"revision_linear": "revision_cov"})
d["e_aic"] = chosen("aic")
d["e_aiccov"] = chosen("aic", remap={"revision_linear": "revision_cov"})

r = d.groupby(key + ["design"]).agg(reps=("rep", "size"), I=("named_INTERCEPT", "mean"), S=("named_SLOPE", "mean"), L=("named_LINK", "mean"), C=("named_COV", "mean")).reset_index()
r.to_csv(os.path.join(HERE, "SIM4_RATES.csv"), index=False)
r2 = d.groupby(key).agg(reps=("rep", "size"), I=("named_INTERCEPT", "mean"), S=("named_SLOPE", "mean"), L=("named_LINK", "mean"), C=("named_COV", "mean")).reset_index()
r2 = r2.sort_values(by=key, key=lambda s: s.map(ordr) if s.name == "dep" else s)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 200)
print("\nShare naming each part, both designs:\n", r2.round(3).to_string(index=False))

cols = ["e_ours", "e_ours_nocov", "e_ver", "e_vercov", "e_aic", "e_aiccov"] + ["err_" + u for u in UPD]
e = d.groupby(key)[cols].mean().reset_index()
e = e.sort_values(by=key, key=lambda s: s.map(ordr) if s.name == "dep" else s)
e.to_csv(os.path.join(HERE, "SIM4_ERROR.csv"), index=False)
print("\nMean abs risk error (pp), both designs:\n", e.round(2).to_string(index=False))

# does a COV verdict help new patients? among data sets where COV was named: revision against what the rule would give without it
c = d[d.named_COV].copy()
c["gain"] = c.e_ours_nocov - c.err_revision_cov
o = c.groupby(key).agg(named=("gain", "size"), gain=("gain", "mean"), se=("gain", lambda x: x.std() / np.sqrt(len(x))), helped=("gain", lambda x: (x > 0).mean())).reset_index()
o = o.sort_values(by=key, key=lambda s: s.map(ordr) if s.name == "dep" else s)
o.to_csv(os.path.join(HERE, "SIM4_COV.csv"), index=False)
print("\nCOV named: error without the revision minus with it (pp; positive = revision helped):\n", o.round(3).to_string(index=False))

print("\nPaired, always revise minus ours, and always flexible minus ours (positive = ours better):")
for k, g in d.groupby(key):
    a = g.err_revision_cov - g.e_ours; b = g.err_flexible - g.e_ours; v = g.e_vercov - g.e_ours
    print(k, f"rev {a.mean():+.2f} ({a.std()/np.sqrt(len(a)):.2f})  flex {b.mean():+.2f} ({b.std()/np.sqrt(len(b)):.2f})  ver {v.mean():+.2f} ({v.std()/np.sqrt(len(v)):.2f})")
