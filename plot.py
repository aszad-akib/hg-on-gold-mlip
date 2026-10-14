#!/usr/bin/env python
"""Plot Hg adsorption energies from results_mace.csv and results_chgnet.csv"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for m in ["mace", "chgnet"]:
    rows += list(csv.DictReader(open(f"results/results_{m}.csv")))
for r in rows:
    r["E"] = float(r["E_ads_kJmol"]); r["nb"] = int(r["n_Au_neighbours"])

labels = [f"{r['facet']}\n{r['site']}" for r in rows if r["model"] == "mace"]
x = range(len(labels))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [2.2, 1]})
for i, (m, c, name) in enumerate([("mace", "tab:blue", "MACE-MP-0"), ("chgnet", "tab:orange", "CHGNet")]):
    E = [r["E"] for r in rows if r["model"] == m]
    ax1.bar([j + (i - 0.5) * 0.4 for j in x], E, 0.4, color=c, label=name)
    nb = [r["nb"] for r in rows if r["model"] == m]
    ax2.scatter(nb, E, color=c, s=50, label=name)
for ax in (ax1, ax2):
    ax.axhline(-98, color="k", ls="--", lw=1, label="Hg/Au experiment (-98)")
    ax.axhline(-52, color="r", ls=":", lw=1, label="Cn/Au experiment (-52)")
    ax.set_ylabel("E_ads (kJ/mol)")
ax1.set_xticks(list(x)); ax1.set_xticklabels(labels, fontsize=8)
ax1.set_title("Hg adsorption on gold: site and facet")
ax1.legend(fontsize=8, loc="lower left")
ax2.set_xlabel("number of Au neighbours of Hg"); ax2.set_xticks([1, 2, 3, 4, 5])
ax2.set_title("Binding vs coordination")
plt.tight_layout(); plt.savefig("figures/Eads_facets.png", dpi=150)
print("saved Eads_facets.png")

print(f"\n{'facet':8s} {'best site':12s} {'MACE':>8s} {'CHGNet':>8s}  (kJ/mol)")
for f in ["Au(111)", "Au(100)", "Au(110)"]:
    best = {m: min((r for r in rows if r["model"] == m and r["facet"] == f), key=lambda r: r["E"]) for m in ["mace", "chgnet"]}
    print(f"{f:8s} {best['mace']['site']:12s} {best['mace']['E']:8.1f} {best['chgnet']['E']:8.1f}")
print(f"{'experiment':21s} {-98:8.1f}")
