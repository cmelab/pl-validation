"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow.environment import DefaultSlurmEnvironment
import os
import logging
import argparse

logging.basicConfig(level=logging.DEBUG)



class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpu",
            help="Specify the partition to submit to."
        )


class R2(DefaultSlurmEnvironment):
    hostname_pattern = "r2"
    template = "r2.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpuq",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )


def sampled(job):
    return job.doc.get("done")

def initialized(job):
    pass

def initial_run_done(job):
    return job.isfile("trajectory.gsd")

def equilibrated(job):
    return job.doc.equilibrated


def sample(job):

    import flowermd
    from flowermd.utils import get_target_box_mass_density
    from flowermd.base import Simulation, Molecule
    from flowermd.library import FF_from_file
    from flowermd.base.system import Pack
    import unyt as u

    import mbuild as mb

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        if job.sp.remove_hydrogens:
            dt = 0.0003
        if not job.sp.remove_hydrogens:
            dt = 0.0001
        system_file = "/home/jacobbieri/repos/persistence-length/P3HT.mol2"
        mol_path = os.path.join(os.getcwd(), system_file)
        ff_file = "/home/jacobbieri/repos/persistence-length/P3HT.xml"
        ff_path = os.path.join(os.getcwd(), ff_file)

        def espaloma_mol(file_path):
            mol = mb.load(file_path)
            for p in mol.particles():
                p.name = f"_{p.name}"
            return mol

        esp_mol = espaloma_mol(mol_path)

        system = Pack(
                molecules=esp_mol,
                density=job.sp.density,
                packing_expand_factor=5
        )

        molff = FF_from_file(ff_path)
        system = Pack(molecules=esp_mol,density=job.sp.density * u.g/u.cm**3, packing_expand_factor=5)
        system.apply_forcefield(r_cut=2.5, force_field=molff, auto_scale=True,remove_charges=True, remove_hydrogens=True)
        
        job.doc.ref_distance = system.reference_distance
        job.doc.ref_mass = system.reference_mass
        job.doc.ref_energy = system.reference_energy

        gsd_path = os.path.join(job.ws, "trajectory.gsd")
        log_path = os.path.join(job.ws, "sim_data.txt")

        system_sim = Simulation(
            initial_state=system.hoomd_snapshot,
            gsd_write_freq=job.sp.n_steps/1000,
            gsd_file_name=gsd_path,
            log_file_name=log_path,
            log_write_freq=job.sp.n_steps/1000,
            dt=job.sp.dt
        )
        
        system_sim.run_NVT(
                kT=job.sp.kT,
                n_steps=job.sp.n_steps,
                tau_kt=job.sp.tau_kt
        )

if __name__ == '__main__':
    # Parse the command line arguments: python action.py --action <ACTION> [DIRECTORIES]
    parser = argparse.ArgumentParser()
    parser.add_argument('--action', required=True)
    parser.add_argument('directories', nargs='+')
    args = parser.parse_args()

    # Open the signac jobs
    project = signac.get_project()
    jobs = [project.open_job(id=directory) for directory in args.directories]

    # Call the action
    # globals()[args.action](*jobs)
