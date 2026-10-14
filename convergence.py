#!/usr/bin/env python
"""Convergence and strain tests for Hg on the Au(111) fcc site.
usage: python convergence.py mace|chgnet      -> conv_<model>.csv"""
import sys, csv, time, warnings
import numpy as np
from ase import Atoms
from ase.build import bulk, fcc111, add_adsorbate
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

def e_ads(layers=4, size=3, vacuum=12.0, a=None):
    slab = fcc111("Au", size=(size, size, layers), a=a or a0, vacuum=vacuum)
    slab.set_constraint(FixAtoms(mask=[t > layers // 2 for t in slab.get_tags()]))
    slab.calc = calc; BFGS(slab, logfile=None).run(fmax=0.02, steps=400)
    e_slab = slab.get_potential_energy()
    s = slab.copy(); add_adsorbate(s, "Hg", 2.6, "fcc")
    s.set_constraint(s.constraints + [FixCartesian(len(s) - 1, mask=(1, 1, 0))])
    s.calc = calc; BFGS(s, logfile=None).run(fmax=0.02, steps=500)
    return s.get_potential_energy() - e_slab - e_hg

tests = [("layers", n, dict(layers=n)) for n in [3, 4, 5, 6]] + \
        [("cell", m, dict(size=m)) for m in [2, 3, 4]] + \
        [("vacuum", v, dict(vacuum=v)) for v in [8.0, 12.0, 16.0]] + \
        [("lattice", "own", {}), ("lattice", 4.08, dict(a=4.08))]
rows = []
for test, value, kw in tests:
    e = e_ads(**kw)
    rows.append(dict(model=model, test=test, value=value, E_ads_eV=round(e, 4), E_ads_kJmol=round(e * 96.485, 1)))
    print(f"{test:8s} {str(value):6s} {e:8.3f} eV {e * 96.485:7.1f} kJ/mol", flush=True)
with open(f"conv_{model}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"saved conv_{model}.csv   ({time.time() - t0:.0f} s)")
