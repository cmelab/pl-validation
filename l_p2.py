import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm


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

def persistence_length(filepath, atom_index, monomer_count, start, stop, interval):
    """
    filepath needs to be a format in which you can
    create an mdanalysis universe from, we mostly use gsd files
    """
    u = mda.Universe(topology=filepath)

    n_monomers = monomer_count
    n_atoms_total = len(u.atoms)
    n_atoms_per_monomer = n_atoms_total // n_monomers
    monomer_atoms = u.atoms[atom_index::n_atoms_per_monomer]


    """create bonds list"""
    autocorrelation = []
    bond_len = []
    for t in u.trajectory[start:stop:interval]:
        particle_positions = []
        bonds = []
        unit_bonds = []
        bond_lengths = []
        angles = []

        for atom in monomer_atoms: # Looping through atoms in a monomer
            pos = atom.position
            particle_positions.append(pos)
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
