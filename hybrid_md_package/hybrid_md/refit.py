#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021-2023.
"""
Refitting of model on the fly

This is a generic refitting function, specific ones and tweaks of
this one with the same interface are to be implemented here.
"""

import importlib
import os.path
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from time import time

import ase.io
import numpy as np

from hybrid_md.settings import SoapParamPreset
from hybrid_md.state_objects import HybridMD


def refit(state: HybridMD):
    """Refit a GAP model, with in-place update

    This is a generic one, which can import the function

    Parameters
    ----------
    state: HybridMD

    """
    refit_settings = state.settings.refit

    if refit_settings.function_name is None:
        return refit_generic(state)

    refit_function_import = state.settings.refit_function_name

    # separate import path
    module_name = ".".join(refit_function_import.split(".")[:-1])
    function_name = refit_function_import.split(".")[-1]

    # import the module of the refit function
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        raise RuntimeError(f"Refit function's module not found: {module_name}") from exc

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


def save_previous_model_files(state: HybridMD):
    """Handle Precious GAP-fit's files

    Saves them into a directory

    Parameters
    ----------
    state

    Returns
    -------

    """
    refit_settings = state.settings.refit

    main_xml = Path(refit_settings.gp_name)
    if not main_xml.is_file():
        # no file to deal with
        return

    # create new directory for saving the models to
    timestamp = datetime.utcnow().strftime("%Y-%m-%d_%H_%M_%S")
    save_to = Path(f"GAP_model-step{state.md_iteration:0>6}-at-{timestamp}")
    save_to.mkdir(exist_ok=True)

    for filename in sorted(Path("").glob(f"{refit_settings.gp_name}*")) + sorted(
        Path("").glob("train.xyz*")
    ):
        # all GAP files & training xyz
        if filename.is_file():
            shutil.move(filename, save_to / filename)


def refit_generic(state: HybridMD):
    """Refit a GAP model, with in-place update

    This is a generic very simple solution, that should work as a
    first try starting from scratch. Change this function to your
    own system and fitting settings as needed.

    Parameters
    ----------
    state: HybridMD
    """

    # extract main settings
    refit_settings = state.settings.refit

    # move any files leftover from the last fit
    save_previous_model_files(state)

    # deal with parameters

    # 2B + SOAP model
    frames_train = ase.io.read(state.xyz_filename, ":") + state.get_previous_data()

    if refit_settings.descriptor_str is not None:
        descriptor_str = refit_settings.descriptor_str
    else:
        # generic 2B+SOAP, need the frames for delta
        delta = np.std([at.info["QM_energy"] / len(at) for at in frames_train]) / 4
        desc_str_2b = (
            "distance_Nb order=2 n_sparse=20 cutoff=5.0 cutoff_transition_width=1.0 "
            "compact_clusters covariance_type=ard_se theta_uniform=1.0 "
            f"sparse_method=uniform f0=0.0 add_species=T delta={delta}"
        )

        # SOAP parameters
        if refit_settings.preset_soap_param is not None:
            desc_str_soap = refit_settings.preset_soap_param.soap_str_stub
        else:
            desc_str_soap = SoapParamPreset.fast.soap_str_stub

        desc_str_soap += f"delta={delta}"
        descriptor_str = desc_str_2b + " : " + desc_str_soap

    # save the previous model
    if os.path.isfile(refit_settings.gp_name):
        shutil.move(refit_settings.gp_name, f"save__{time()}__{refit_settings.gp_name}")

    # training structures & delta
    ase.io.write("train.xyz", frames_train)
    if os.path.isfile("train.xyz.idx"):
        os.remove("train.xyz.idx")

    # deal with e0 & e0_method together
    e0_args = ""
    if refit_settings.e0_method:
        e0_args += " e0_method=" + refit_settings.e0_method
    if refit_settings.e0:
        e0_args += " e0=" + refit_settings.e0

    # assemble the fitting string
    fit_str = (
        f"gap_fit at_file=train.xyz gp_file={refit_settings.gp_name}"
        f" energy_parameter_name=QM_energy"
        f" force_parameter_name=QM_forces"
        f" virial_parameter_name=QM_virial"
        f" do_copy_at_file=F sparse_separate_file=T"
        f" default_sigma={{ {refit_settings.default_sigma} }}"
        f" {e0_args}"
        f" gap={{ {descriptor_str} }} "
        f" {refit_settings.extra_gap_opts}"
    )

    # change to OMP-parallel mode for fitting
    omp_num_threads_before = os.environ.get("OMP_NUM_THREADS", "1")
    if refit_settings.use_omp:
        os.environ["OMP_NUM_THREADS"] = str(refit_settings.num_threads)

    # fit the model
    proc = subprocess.run(
        fit_str, shell=True, capture_output=True, text=True, check=True
    )

    # set OMP_NUM_THREADS back to normal
    if refit_settings.use_omp:
        os.environ["OMP_NUM_THREADS"] = omp_num_threads_before

    # print the outputs to file
    with open(f"{refit_settings.gp_name}.fit-stdout.{time()}.txt", "w") as file:
        file.write(proc.stdout)
    with open(f"{refit_settings.gp_name}.fit-stderr.{time()}.txt", "w") as file:
        file.write(proc.stderr)
