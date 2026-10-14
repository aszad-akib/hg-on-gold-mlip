#!/usr/bin/env python
"""Hg adsorption on Au(111): 4 sites, own lattice constant, site-locked Hg.
Modified from miroi/open-collection Hg_on_gold baseline.  usage: python hg_sites.py mace|chgnet"""
import sys, warnings
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
# 1. this model's own lattice constant
vols, ens = [], []
for a in np.linspace(3.95, 4.35, 9):
    au = bulk("Au", "fcc", a=a, cubic=True); au.calc = calc
    vols.append(au.get_volume()); ens.append(au.get_potential_energy())
a0 = EquationOfState(vols, ens).fit()[0] ** (1 / 3)
print(f"[{model}] lattice constant a0 = {a0:.3f} A  (experiment 4.078 A)")
# 2. relaxed Au(111) slab, 3x3x4, bottom 2 layers frozen
slab = fcc111("Au", size=(3, 3, 4), a=a0, vacuum=12.0)
slab.set_constraint(FixAtoms(mask=[t > 2 for t in slab.get_tags()]))
slab.calc = calc; BFGS(slab, logfile=None).run(fmax=0.02)
e_slab = slab.get_potential_energy()
# 3. isolated Hg atom
hg = Atoms("Hg", cell=[18, 18, 18], pbc=True); hg.center(); hg.calc = calc
e_hg = hg.get_potential_energy()
# 4. Hg on each site (held laterally on the site)
print(f"{'site':8s} {'E_ads(eV)':>10s} {'kJ/mol':>8s} {'height(A)':>10s}")
for site in ["ontop", "bridge", "fcc", "hcp"]:
    s = slab.copy(); add_adsorbate(s, "Hg", 2.6, site)
    s.set_constraint(s.constraints + [FixCartesian(len(s) - 1, mask=(1, 1, 0))])
    s.calc = calc; BFGS(s, logfile=None).run(fmax=0.02, steps=500)
    e = s.get_potential_energy() - e_slab - e_hg
    h = s.positions[-1, 2] - slab.positions[:, 2].max()
    print(f"{site:8s} {e:10.3f} {e * 96.485:8.1f} {h:10.2f}")
    s.write(f"Hg_{site}_{model}.xyz")
print("experiment (thermochromatography): about -98 kJ/mol")
