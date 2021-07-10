#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.
"""
Refitting of model on the fly

This is a generic refitting function, specific ones and tweaks of
this one with the same interface are to be implemented here.
"""

import os.path
import shutil
import subprocess
from time import time

import ase.io
import numpy as np
from hybrid_md.state_objects import HybridMD


def refit(state: HybridMD):
    """Refit a GAP model, with in-place update

    This is a generic very simple solution, that should work as a
    first try starting from scratch. Change this function to your
    own system and fitting settings as needed.

    Parameters
    ----------
    state: HybridMD

    """
    # 2B + SOAP model
    gp_name = "GAP.xml"
    soap_n_sparse = 200
    frames_train = ase.io.read(state.xyz_filename, ":") + state.get_previous_data()

    # save the previous model
    if os.path.isfile(gp_name):
        shutil.move(gp_name, f"save__{time()}__{gp_name}")

    # training structures & delta
    ase.io.write("train.xyz", frames_train)
    delta = np.std([at.info["QM_energy"] / len(at) for at in frames_train]) / 4

    # descriptors
    desc_str_2b = (
        f"distance_Nb order=2 n_sparse=20 cutoff=4.5 cutoff_transition_width=1.0 "
        f"compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform "
        f"f0=0.0 add_species=T delta={delta}"
    )
    desc_str_soap = (
        f"soap n_sparse={soap_n_sparse} n_max=8 l_max=4 cutoff=4.0 cutoff_transition_width=1.0 "
        f"atom_sigma=0.5 add_species=True "
        f"delta={delta} covariance_type=dot_product zeta=4 sparse_method=cur_points"
    )

    # use lower kernel regularisation
    default_sigma = "0.001 0.020 0.1 1.0"

    # NEW descriptor str, with 2b and SOAP concatenated with the ":" character in the gap string.
    fit_str = (
        f"gap_fit at_file=train.xyz gp_file={gp_name} "
        f"energy_parameter_name=QM_energy force_parameter_name=QM_forces"
        f" virial_parameter_name=QM_virial "
        f"sparse_jitter=1.0e-8 do_copy_at_file=F sparse_separate_file=F "
        f"default_sigma={{ {default_sigma} }} e0_method=average "
        f"gap={{ {desc_str_soap} : {desc_str_2b} }}"
    )

    # fit the 2b+SOAP model
    proc = subprocess.run(
        fit_str, shell=True, capture_output=True, text=True, check=True
    )

    # print the outputs to file
    with open(f"stdout_{gp_name}_at_{time()}__.txt", "w") as file:
        file.write(proc.stdout)
    with open(f"stderr_{gp_name}_at_{time()}__.txt", "w") as file:
        file.write(proc.stderr)
