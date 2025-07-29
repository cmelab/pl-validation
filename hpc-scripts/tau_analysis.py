from l_p2 import persistence_length
import matplotlib.pyplot as plt
import numpy as np
import scipy
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
import glob
import os
import pandas as pd
from pl_com import persistence_length1
from polymer_dictionary import polymer_dictionary

# Do not change
ref_length = 0.3563594872561357

# Change based on experiment
num_monomers = 10
variable_being_tested = "C1_tau_sweep"

key_list = sorted(list(polymer_dictionary.keys()))
path = os.getcwd()
molecule_list = sorted(glob.glob(path+"/gsd_files/tau_sweep/"+"*.gsd"))

p_lens = []
exp_fits = []
for i in range(len(molecule_list)):
    h = persistence_length(filepath=molecule_list[i],
                       atom_index=polymer_dictionary.get(key_list[1]).get("Sulfur index"),
                       monomer_count=num_monomers,
                       start=50)
    l_p = h[0]
    # l_b = h[1]
    # x_values = h[2]
    # C_n = h[3]
    exp_fit = h[4]
    # decorr = h[-2]
    i += 1
    p_lens.append(l_p*ref_length*10)
    exp_fits.append(exp_fit)


p_lens1 = []
exp_fits1 = []
for i in range(len(molecule_list)):
    h = persistence_length(filepath=molecule_list[i],
                       monomer_count=num_monomers,
                       start=50)
    l_p = h[0]
    # l_b = h[1]
    # x_values = h[2]
    # C_n = h[3]
    exp_fit = h[4]
    # decorr = h[-2]
    i += 1
    p_lens1.append(l_p*ref_length*10)
    exp_fits1.append(exp_fit)

kts = [0.5,2,3,4,5]
df = pd.DataFrame({"Polymer": kts, "Atom index pl": p_lens, "COM pl": p_lens1})
df.to_csv(variable_being_tested)
