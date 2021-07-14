#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.
"""
Refitting of model on the fly

This is a generic refitting function, specific ones and tweaks of
this one with the same interface are to be implemented here.
"""

import importlib
import os.path
import shutil
import subprocess
from time import time

import ase.io
import numpy as np
from hybrid_md.state_objects import HybridMD


def refit(state: HybridMD):
    """Refit a GAP model, with in-place update

    This is a generic one, which can import the functio

    Parameters
    ----------
    state: HybridMD

    """
    if state.refit_function_name is None:
        return refit_generic(state, None, None)
    else:
        refit_function_import = state.refit_function_name

        # separate import path
        module_name = ".".join(refit_function_import.split(".")[:-1])
        function_name = refit_function_import.split(".")[-1]

        # import the module of the refit function
        try:
            module = importlib.import_module(module_name)
        except ModuleNotFoundError:
            raise RuntimeError(f"Refit function's module not found: {module_name}")

        # class of the calculator
        if hasattr(module, function_name):
            refit_function = getattr(module, function_name)
            assert callable(refit_function)
        else:
            raise RuntimeError(
                f"Refit function ({function_name}) not found in module {module_name}"
            )

        # YAY, all great now
        return refit_function(state)


def refit_turbo_si_c(state: HybridMD):
    # Refit function with TurboSOAP for the SiC example
    return refit_turbo_two_species(state, "6 14")


def refit_turbo_h_c(state: HybridMD):
    # Hydrogen and Carbon
    return refit_turbo_two_species(state, "1 6", 300)


def refit_fe_h(state: HybridMD):
    # Hydrogen and Iron

    frames_train = ase.io.read(state.xyz_filename, ":") + state.get_previous_data()
    delta = np.std([at.info["QM_energy"] / len(at) for at in frames_train]) / 4

    # soap
    soap_n_sparse = 300

    # there is no H-H interaction YET
    desc_str_2b = (
        "distance_Nb order=2 n_sparse=20 cutoff=4.5 cutoff_transition_width=1.0 "
        "compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform "
        f"f0=0.0 add_species=F delta={delta} "
        "Z={{1 26}} : "
        "distance_Nb order=2 n_sparse=20 cutoff=4.5 cutoff_transition_width=1.0 "
        "compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform "
        f"f0=0.0 add_species=F delta={delta} "
        "Z={{26 26}} "
    )

    # turbo SOAP
    soap_common = (
        f"soap n_sparse={soap_n_sparse} n_max=8 l_max=4 delta={delta} covariance_type=dot_product "
        f"zeta=4 sparse_method=cur_points n_species=2 "
    )
    desc_str_soap = (
        f"{soap_common} cutoff=3.0 cutoff_transition_width=0.6 atom_sigma=0.3 Z=1 "
        + "species_Z={{1 26}} : "
        f"{soap_common} cutoff=5.5 cutoff_transition_width=1.0 atom_sigma=0.5 Z=26 "
        + "species_Z={{1 26}}"
    )

    descriptor_strs = desc_str_2b + " : " + desc_str_soap

    # use lower kernel regularisation
    default_sigma = "0.002 0.050 1.0 1.0"

    return refit_generic(state, descriptor_strs, default_sigma)


def refit_turbo_two_species(state: HybridMD, species_str: str, soap_n_sparse=200):
    # refit with turbo-soap, given two species

    frames_train = ase.io.read(state.xyz_filename, ":") + state.get_previous_data()
    delta = np.std([at.info["QM_energy"] / len(at) for at in frames_train]) / 4

    # descriptors
    desc_str_2b = (
        f"distance_Nb order=2 n_sparse=20 cutoff=4.5 cutoff_transition_width=1.0 "
        f"compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform "
        f"f0=0.0 add_species=T delta={delta}"
    )

    # turbo SOAP
    desc_str_soap = (
        f"soap n_sparse=200 n_max=8 l_max=4 cutoff=4.5 cutoff_transition_width=1.0 "
        f"atom_sigma=0.5 add_species=True "
        f"delta={delta} covariance_type=dot_product zeta=4 sparse_method=cur_points"
    )

    descriptor_strs = desc_str_2b + " : " + desc_str_soap

    # use lower kernel regularisation
    default_sigma = "0.001 0.050 0.1 1.0"

    return refit_generic(state, descriptor_strs, default_sigma)


def refit_generic(
    state: HybridMD, descriptor_strs: str = None, default_sigma: str = None
):
    """Refit a GAP model, with in-place update

    This is a generic very simple solution, that should work as a
    first try starting from scratch. Change this function to your
    own system and fitting settings as needed.

    Parameters
    ----------
    state: HybridMD
    descriptor_strs : str
        descriptor strings, ':' separated, no brackets around them
    default_sigma : str
        default sigma, four numbers separated by ':'

    """
    if default_sigma is None:
        default_sigma = "0.005 0.050 0.1 1.0"

    # 2B + SOAP model
    gp_name = "GAP.xml"
    frames_train = ase.io.read(state.xyz_filename, ":") + state.get_previous_data()

    if descriptor_strs is None:
        # generic 2B+SOAP, need the frames for delta
        delta = np.std([at.info["QM_energy"] / len(at) for at in frames_train]) / 4
        desc_str_2b = (
            f"distance_Nb order=2 n_sparse=20 cutoff=4.5 cutoff_transition_width=1.0 "
            f"compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform "
            f"f0=0.0 add_species=T delta={delta}"
        )
        desc_str_soap = (
            f"soap n_sparse=200 n_max=8 l_max=4 cutoff=4.0 cutoff_transition_width=1.0 "
            f"atom_sigma=0.5 add_species=True "
            f"delta={delta} covariance_type=dot_product zeta=4 sparse_method=cur_points"
        )
        descriptor_strs = desc_str_2b + " : " + desc_str_soap

    # save the previous model
    if os.path.isfile(gp_name):
        shutil.move(gp_name, f"save__{time()}__{gp_name}")

    # training structures & delta
    ase.io.write("train.xyz", frames_train)
    if os.path.isfile("train.xyz.idx"):
        os.remove("train.xyz.idx")

    # assemble the fitting string
    fit_str = (
        f"gap_fit at_file=train.xyz gp_file={gp_name} "
        f"energy_parameter_name=QM_energy force_parameter_name=QM_forces"
        f" virial_parameter_name=QM_virial "
        f"sparse_jitter=1.0e-8 do_copy_at_file=F sparse_separate_file=T "
        f"default_sigma={{ {default_sigma} }} e0_method=average "
        f"gap={{ {descriptor_strs} }}"
    )

    with open("debug_output.txt", "w") as file:
        file.write(fit_str)

    # fit the 2b+SOAP model
    os.environ["OMP_NUM_THREADS"] = "32"
    proc = subprocess.run(
        fit_str, shell=True, capture_output=True, text=True, check=True
    )
    os.environ["OMP_NUM_THREADS"] = "1"

    # print the outputs to file
    with open(f"stdout_{gp_name}_at_{time()}__.txt", "w") as file:
        file.write(proc.stdout)
    with open(f"stderr_{gp_name}_at_{time()}__.txt", "w") as file:
        file.write(proc.stderr)
