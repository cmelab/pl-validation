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


system_file = "/home/jacobbieri/repos/persistence-length/mol2/10_mers/PCPDTFBT_C5_BO_10mer.mol2"
ff_filepath = "/home/jacobbieri/repos/persistence-length/xml/PCPDTFBT_C5_BO.xml"
gsd_path = "PCPDTFBT-C5-BO-10_mer"
espmol = espaloma_mol(system_file)
molecule = Molecule(num_mols=1, compound=espmol)

molff = FF_from_file(ff_filepath)
start_time = time.time()
system = Pack(molecules=molecule,density=0.01 * u.g/u.cm**3)
print("Time:", time.time() - start_time, "s to pack")
system.apply_forcefield(r_cut=2.5, force_field=molff, auto_scale=True,remove_charges=True, remove_hydrogens=True)
system.hoomd_snapshot
sim = Simulation.from_system(system=system, gsd_write_freq=10000, log_write_freq=10000, gsd_file_name=gsd_path+".gsd",log_file_name=gsd_path+".txt")
sim.run_NVT(n_steps=5e7, kT=5.0, tau_kt=1)
sim.flush_writers()
