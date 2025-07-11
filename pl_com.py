import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
import gsd
from grits import utils

def autocorr1D(array):
    """Takes in a linear np array, performs autocorrelation
    function and returns normalized array with half the length
    of the input.

    Parameters
    ----------
    data : numpy.typing.Arraylike, required
        1-D series of data to perform autocorrelation on.

    Returns
    -------
    1D np.array

    """
    ft = np.fft.rfft(array - np.average(array))
    acorr = np.fft.irfft(ft * np.conjugate(ft)) / (len(array) * np.var(array))
    return acorr[0 : len(acorr) // 2]  # noqa: E203

def get_decorr(acorr):
    """
    Returns the decorrelation time of the autocorrelation, a 1D np array
    """
    return np.argmin(acorr > 0)

def get_monomers(filepath, frame_num, monomer_count=1):
    """
    Returns a list of GSD frames of each monomer in a polymer. Use this if you need GSD Frame Particle data for only SOME atoms.
    """
    with gsd.hoomd.open(filepath, "r") as traj:
        frame = traj[frame_num]
    
    n_atoms_total = frame.particles.N
    atoms_per_monomer = n_atoms_total // monomer_count
    
    monomer_frames = []
    
    for i in range(monomer_count):
        start = i * atoms_per_monomer # First atom in a monomer
        end = (i + 1) * atoms_per_monomer # Last atom in the monomer
        new_indicies = np.arange(start, end)
        
        # Add all properties of the particles into the new frame
        new_frame = gsd.hoomd.Frame()
        new_frame.particles.N = len(new_indicies)
        new_frame.particles.position = frame.particles.position[new_indicies]
        new_frame.particles.typeid = frame.particles.typeid[new_indicies]
        new_frame.particles.types = frame.particles.types
        new_frame.particles.velocity = frame.particles.velocity[new_indicies]
        new_frame.particles.mass = frame.particles.mass[new_indicies]
        new_frame.particles.charge = frame.particles.charge[new_indicies]
        new_frame.particles.body = frame.particles.body[new_indicies]
    
        monomer_frames.append(new_frame)
    return monomer_frames


def persistence_length1(filepath, monomer_count, start=0, stop=None, interval=1):
    """
    filepath needs to be a format in which you can
    create an mdanalysis universe from, we mostly use gsd files
    """
    u = mda.Universe(topology=filepath)

    n_monomers = monomer_count
    total_atoms = len(u.atoms)
    atoms_per_monomer = total_atoms // n_monomers

    """create bonds list"""
    autocorrelation = []
    bond_len = []
    for t in u.trajectory[start:stop:interval]:
        particle_positions = []
        bonds = []
        unit_bonds = []
        bond_lengths = []
        angles = []
        coms = []
        mon_particles = get_monomers(filepath, t.frame, monomer_count)

        for frame in mon_particles:
            particles = frame.particles
            pos,mass = utils.get_heavy_atoms(particles)
            coms.append(utils.get_com(pos,mass))
        for com in coms: # Looping through atoms in a monomer
            particle_positions.append(com)
        for i in range(len(particle_positions)-1):
            b = particle_positions[i+1]-particle_positions[i]
            bonds.append(b)
            l2 = t.dimensions[0]/2
            for i,b in enumerate(bonds):
                for j,x in enumerate(b):
                    if x>l2:
                        bonds[i][j] = x-l2*2
                    if x<-l2:
                        bonds[i][j] = x+l2*2
            a = b/np.linalg.norm(b)
            unit_bonds.append(a)
            length = np.linalg.norm(b)
            bond_lengths.append(length)
        bond_len.append(bond_lengths)


        for i in range(len(unit_bonds)-1):
            b1 = unit_bonds[0]
            b2 = unit_bonds[0+i]
            dot_product = np.dot(b1,b2)
            angles.append(dot_product)

        n=len(u.bonds)
        n_frames = u.trajectory.n_frames
        n_chains = 1
        norm = np.linspace(1,n- 1, n - 1)
        norm *= n_chains# * n_frames
        autocorrelation.append(angles)#/norm)

    '''average the data from trajectories together'''
    auto_average = []
    for j in range(len(autocorrelation[0])):
        k = []
        for i in range(len(autocorrelation)):
            k.append(autocorrelation[i][j])
        auto_average.append(np.mean(k))

    l_b = np.average(bond_len)
    x = [i for i in range(len(auto_average))]

    '''set negative results to 0'''
    for r in range(len(auto_average)):
        if auto_average[r] < 0:
            auto_average[r] = 0
    def expfunc(x, a):
        return np.exp(-x/a)

    exp_coeff = scipy.optimize.curve_fit(expfunc,x,auto_average)[0][0]

    l_p = exp_coeff * l_b

    fit = np.exp(-(x/exp_coeff))
    decorrelation = get_decorr(np.array(auto_average))

    return l_p, l_b, x, auto_average, fit, exp_coeff, autocorrelation, decorrelation, unit_bonds
