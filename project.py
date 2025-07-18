"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os
import logging

logging.basicConfig(level=logging.DEBUG)


class MyProject(FlowProject):
    pass


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

# Definition of project-related labels (classification)
@MyProject.label
def sampled(job):
    return job.doc.get("done")

@MyProject.label
def initialized(job):
    pass

@MyProject.label
def initial_run_done(job):
    return job.isfile("trajectory.gsd")

@MyProject.label
def equilibrated(job):
    return job.doc.equilibrated

@MyProject.post(initial_run_done)
@MyProject.operation(directives={"executable":"python -u","ngpu": "1"}, name="test2-sim1")

def sample(job):

    import flowermd
    from flowermd.utils import get_target_box_mass_density
    from flowermd.base import Simulation, Molecule
    from flowermd.library import FF_from_file
    from flowermd.base.system import Pack

    import mbuild as mb
    import foyer

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        if job.sp.remove_hydrogens:
            dt = 0.0003
        if not job.sp.remove_hydrogens:
            dt = 0.0001
        system_file = "/home/jacobbieri/projects/persistence-length/P3HT.mol2"
        mol_path = os.path.join(os.getcwd(), system_file)
        ff_file = "/home/jacobbieri/projects/persistence-length/P3HT.xml"
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
                mol_kwargs = {
                    "file_path": mol_path,
                },
                packing_expand_factor=5
        )

        system_ff = foyer.Forcefield(forcefield_files=ff_path)
        system.apply_forcefield(
                forcefield=system_ff,
                make_charge_neutral=True,
                remove_hydrogens=job.sp.remove_hydrogens,
                remove_charges=job.sp.remove_charges
        )

        job.doc.ref_distance = system.reference_distance
        job.doc.ref_mass = system.reference_mass
        job.doc.ref_energy = system.reference_energy

        gsd_path = os.path.join(job.ws, "trajectory.gsd")
        log_path = os.path.join(job.ws, "sim_data.txt")

        system_sim = Simulation(
            initial_state=system.hoomd_snapshot,
            forcefield=system.hoomd_forcefield,
            gsd_write_freq=job.sp.n_steps/1000,
            gsd_file_name=gsd_path,
            log_file_name=log_path,
            log_write_freq=100000,
            dt=dt
        )
        target_box = system.target_box*10/job.doc.ref_distance
        job.doc.target_box = target_box

        system_sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.shrink_steps,
                period=job.sp.shrink_period,
                tau_kt=job.sp.tau_kt,
                kT=job.sp.shrink_kT
        )
        system_sim.run_NVT(
                kT=job.sp.kT,
                n_steps=job.sp.n_steps,
                tau_kt=job.sp.tau_kt
        )

if __name__ == "__main__":
    MyProject(environment=Fry).main()
