from grits import utils
from MDAnalysis.analysis import polymer
import MDAnalysis as mda
import hoomd
import numpy as np

def pl(gsd_file, start=0, stop=None, stride=1):
    persistence_lengths = []
    frames = []
    u = mda.Universe(gsd_file)
    with gsd.hoomd.open(gsd_file) as traj:
            for snap in traj[start:stop:stride]: # Looping through each frame
                com = utils.get_com(snap.particles.position,snap.particles.mass)
                frames.append(com)
                print(type(com))
                persistence_lengths.append(polymer.PersistenceLength([com]))
    return persistence_lengths, frames
