import warnings
warnings.filterwarnings('ignore')
import hoomd
import gsd
import gsd.hoomd
import unyt as u
import mbuild as mb
import numpy as np
import os
import flowermd
from flowermd.base import Simulation, Molecule
from flowermd.library import FF_from_file
from flowermd.base.system import Pack
import time

def espaloma_mol(file_path):
     mol = mb.load(file_path)
     for p in mol.particles():
           p.name = f"_{p.name}"
     return mol


system_file = "/bsuhome/jacobbieri/repos/pl-validation/mol2/polymers/PCPDTFBT_C5_BO_250mer.mol2"
ff_filepath = "/bsuhome/jacobbieri/repos/pl-validation/xml/PCPDTFBT_C5_BO.xml"
gsd_path = "test-C5-monomer"
espmol = espaloma_mol(system_file)
molecule = Molecule(num_mols=1, compound=espmol)

molff = FF_from_file(ff_filepath)
start_time = time.time()
system = Pack(molecules=molecule,density=0.000001 * u.g/u.cm**3)
print("Time:", time.time() - start_time, "s to pack")
system.apply_forcefield(r_cut=2.5, force_field=molff, auto_scale=True,remove_charges=True, remove_hydrogens=True)
system.hoomd_snapshot
sim = Simulation.from_system(system=system, gsd_write_freq=5000, log_write_freq=5000, gsd_file_name=gsd_path+".gsd",log_file_name=gsd_path+".txt")
sim.run_NVT(n_steps=5e5, kT=5.0, tau_kt=0.05)
sim.flush_writers()
