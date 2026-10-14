#!/usr/bin/env python
"""Production: Hg on Au(111), Au(100), Au(110), all sites.  usage: python production.py mace|chgnet"""
import sys, csv, time, warnings
import numpy as np
from ase import Atoms
from ase.build import bulk, fcc111, fcc100, fcc110, add_adsorbate
from ase.constraints import FixAtoms, FixCartesian
from ase.eos import EquationOfState
from ase.optimize import BFGS
warnings.filterwarnings("ignore")
model = sys.argv[1] if len(sys.argv) > 1 else "mace"
if model == "mace":
    from mace.calculators import mace_mp
    calc = mace_mp(model="medium", default_dtype="float64", device="cpu")
else:
    from chgnet.model.dynamics import CHGNetCalculator
    calc = CHGNetCalculator(use_device="cpu")
t0 = time.time()
vols, ens = [], []
for a in np.linspace(3.95, 4.35, 9):
    au = bulk("Au", "fcc", a=a, cubic=True); au.calc = calc
    vols.append(au.get_volume()); ens.append(au.get_potential_energy())
a0 = EquationOfState(vols, ens).fit()[0] ** (1 / 3)
hg = Atoms("Hg", cell=[18, 18, 18], pbc=True); hg.center(); hg.calc = calc
e_hg = hg.get_potential_energy()
print(f"[{model}] a0 = {a0:.3f} A")
FACETS = {"Au(111)": (fcc111, 4, ["ontop", "bridge", "fcc", "hcp"]),
          "Au(100)": (fcc100, 4, ["ontop", "bridge", "hollow"]),
          "Au(110)": (fcc110, 6, ["ontop", "shortbridge", "longbridge", "hollow"])}
rows = []
for facet, (build, n, sites) in FACETS.items():
    slab = build("Au", size=(3, 3, n), a=a0, vacuum=12.0)
    slab.set_constraint(FixAtoms(mask=[t > n // 2 for t in slab.get_tags()]))
    slab.calc = calc; BFGS(slab, logfile=None).run(fmax=0.02, steps=400)
    e_slab = slab.get_potential_energy()
    for site in sites:
        s = slab.copy(); add_adsorbate(s, "Hg", 2.6, site)
        s.set_constraint(s.constraints + [FixCartesian(len(s) - 1, mask=(1, 1, 0))])
        s.calc = calc; BFGS(s, logfile=None).run(fmax=0.02, steps=500)
        e = s.get_potential_energy() - e_slab - e_hg
        d = s.get_distances(len(s) - 1, range(len(s) - 1), mic=True)
        nb = int((d < d.min() + 0.25).sum())
        rows.append(dict(model=model, facet=facet, site=site, E_ads_eV=round(e, 4),
                         E_ads_kJmol=round(e * 96.485, 1), n_Au_neighbours=nb))
        print(f"{facet:8s} {site:12s} {e:8.3f} eV {e * 96.485:7.1f} kJ/mol  {nb} Au neighbours", flush=True)
        s.write(f"Hg_{facet[3:6]}_{site}_{model}.xyz")
with open(f"results_{model}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"saved results_{model}.csv   ({time.time() - t0:.0f} s)")
