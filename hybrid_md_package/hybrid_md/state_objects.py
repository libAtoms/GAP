#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021-2023.
"""
Objects representing the state of a calculation
"""
import os
import shutil
from abc import ABC
from datetime import datetime
from enum import Enum, auto, unique
from pathlib import Path

import ase.io
import numpy as np
import yaml

from hybrid_md.settings import MainSettings


@unique
class StepKinds(Enum):
    """The kinds of steps we can have while the MD"""

    INITIAL = auto()
    LAST_INITIAL = auto()
    CHECK = auto()
    GENERIC = auto()


class HybridMD:
    """Main class for collecting state of the calculation and utilities.

    This can be refactored, this class may be doing too much at once.
    """

    def __init__(self, seed: str, md_iteration: int = None):
        self.seed = seed

        # files of interest
        self.log_filename = f"{self.seed}.hybrid-md-temporary-log"
        self.input_filename = f"{self.seed}.hybrid-md-input.yaml"
        self.xyz_filename = f"{self.seed}.hybrid-md.xyz"

        # associated sub-state objects
        self.xyz_this_run = SubStateXYZ(seed)
        self.carry = CarriedState(seed)
        self.settings = MainSettings.read_input(self.input_filename)

        # dummy arrays for results
        self.md_iteration = md_iteration

        # misc
        self._continuation_log_lines = None

    def handle_continuation(self):
        """Continuation of the calculation

        notes:
        - `self.md_iteration` is the step at the end of which the checkpoint was made
        - the last check step can be after the checkpoint was made, so we need to
        deal with that
        - in case of the adaptive method, we will try to do GAP up to the desired

        needs to do the following:
        - todo: back up the previous DFT data (xyz in new file)
        - figure out the checking interval & last check step
        """

        # make a backup of the .md file
        md_file = Path(f"{self.seed}.md")
        if md_file.is_file():
            timestamp = datetime.utcnow().strftime("%Y-%m-%d_%H_%M_%S")
            save_to = Path(f"{self.seed}.save-at-{timestamp}.md")
            shutil.copy(md_file, save_to)

        if os.path.isfile(self.carry.state_filename):
            # easy case: state has not disappeared
            self.carry.load()

            self._continuation_log_lines = [
                "CONTINUATION: read previous calculation's state",
                f"   MD iteration of checkpoint: {self.md_iteration}",
                "   MD iteration of last ab-inito calculation: "
                f"{self.carry.last_check_step}",
            ]

            # unpack the previous calculation's variables
            pre_last_check_step = self.carry.last_check_step
            pre_check_interval = self.carry.current_check_interval
            pre_restart_written = self.md_iteration
            pre_next_check = pre_last_check_step + pre_check_interval

            if pre_last_check_step <= pre_restart_written:
                # easy case: checkpoint was written after last ab-initio step
                new_check_interval = pre_check_interval
                new_last_check_step = pre_last_check_step - pre_restart_written
                self._continuation_log_lines.append(
                    "  - checkpoint was written after last ab-inito step, "
                    "no settings were changed"
                )
            else:
                # at least one ab-inito step was made after the checkpoint was written
                self._continuation_log_lines.append(
                    "  - checkpoint was written before last ab-inito step, "
                    "settings updated"
                )

                if self.settings.adaptive_method_parameters:
                    # adaptive method: cap the interval to the maximum
                    new_check_interval = pre_next_check - pre_restart_written
                    if (
                        self.settings.adaptive_method_parameters.n_max
                        < new_check_interval
                    ):
                        new_check_interval = (
                            self.settings.adaptive_method_parameters.n_max
                        )
                        self._continuation_log_lines.append(
                            "  - new check interval capped by n_max of adaptive method"
                        )
                else:
                    # otherwise unchanged check interval
                    new_check_interval = pre_check_interval

                # as if the last step was ab-inito, so we can utilise the whole
                # interval set right now
                new_last_check_step = 0

            # apply updates to the carried state
            self.carry.last_check_step = new_last_check_step
            self.carry.current_check_interval = new_check_interval

            self._continuation_log_lines.extend(
                [
                    f"   MD step of last check according to this run's counting: "
                    f"{new_last_check_step}",
                    f"   MD check interval: {new_check_interval}",
                ]
            )
        else:
            # just use the current settings & raise warning for the user
            self._continuation_log_lines = [
                "CONTINUATION: ",
                "WARNING Have not found previous state file "
                f"({self.carry.state_filename}), ",
                "so cannot perform clean restart. "
                "We are using the current input file's settings instead.",
            ]

        # --------------------------------------
        # deal with the number of initial steps:
        #  - this is not a common use-case, but needed for the logic to make sense
        #  - decrease the number of them in case the calculation stopped while in the
        #  initial stage
        #  - notice that one can restart a calculation in this state, and then restart
        #  again while still being in the initial stage (not expected to be a true
        #  use-case though)
        #
        if self.carry.continuation:
            previous_num_initial_steps = self.carry.continuation_initial_steps
        else:
            previous_num_initial_steps = self.settings.num_initial_steps

        if self.md_iteration < previous_num_initial_steps:
            self._continuation_log_lines.append(
                "INFO: Continuing calculation from before initial ab-initio "
                "steps were completed."
            )

            # decrease the number of initial steps
            self.carry.continuation_initial_steps = (
                previous_num_initial_steps - self.md_iteration
            )

        # --------------------------------------
        # keep track of this for the rest of the calculation
        self.carry.continuation = True

        self._continuation_log_lines.append(
            "! end of continuation info "
            "---------------------------------------------------"
        )

        # add line breaks to the ends of lines
        self._continuation_log_lines = ["\n"] + [
            x + "\n" for x in self._continuation_log_lines
        ]

    # --------------------------------------------------------------------------
    # step's IO

    def io_initial_step_banner(self):
        """Construct the banner shown at initialisation - written into the .castep"""
        lines = ["\n", "Using Hybrid MD\n", "! input file" + "-" * 66 + "\n"]

        with open(self.input_filename, "r") as file:
            lines.extend(file.readlines())

        lines.append("! end of input file" + "-" * 61)
        lines.append("")

        # add lines from continuation method
        if self._continuation_log_lines is not None:
            lines.extend(self._continuation_log_lines)

        self.write_to_tmp_log(lines, append=False)

    def write_to_tmp_log(self, lines, append):
        # write log to temporary file

        if append:
            mode = "a"
        else:
            mode = "w"

        with open(self.log_filename, mode) as file:
            file.writelines(lines)

    # --------------------------------------------------------------------------
    # Error table for IO

    def _tolerance_line(self, name: str, value: float, tolerance, unit: str):
        if tolerance is not None:
            # line to be printed
            yes_no = self._bool_to_str(tolerance > value)

            return (
                f" {name:>12} | {value:16.8f} | {tolerance:16.8f} |"
                f" {unit:10} | {yes_no:3} | <-- Hybrid-MD\n"
            )

        return (
            f" {name:>12} | {value:16.8f} |              Off | "
            f"{unit:10} |     | <-- Hybrid-MD\n"
        )

    @staticmethod
    def _tolerance_line_cumulative(name: str, value: float, unit: str):
        return (
            f" {name:>12}{value:30.8f}               | "
            f"{unit:10} | <-- Hybrid-MD-Cumul\n"
        )

    def check_tolerances(self):
        """Check if tolerances are met"""
        # todo: implement species specific checks
        # todo: implement passing X out of Y to pass rather than all()

        # pass if none of them are given
        checked_tolerances = [True]

        for val, tol in [
            (self.xyz_this_run.get_ediff(), self.get_tolerance("ediff")),
            (self.xyz_this_run.get_fmax(), self.get_tolerance("fmax")),
            (self.xyz_this_run.get_frmse(), self.get_tolerance("frmse")),
            (self.xyz_this_run.get_vmax(), self.get_tolerance("vmax")),
        ]:
            if tol is None:
                continue

            # perform the check
            checked_tolerances.append(val < tol)

        return all(checked_tolerances)

    def error_table(self):
        # line formatting
        separator = (
            " -------------+------------------+------------------"
            "+------------+-----+ <-- Hybrid-MD\n"
        )
        refit_str = "Refitting!" if self.carry.do_update_model else "No Refit"

        # container for table -> write to file later
        lines = [
            f"\n      Hybrid-MD: MD iteration{self.md_iteration:8}       "
            f"                             <-- Hybrid-MD\n"
            if self.md_iteration is not None
            else "\n",
            separator,
            "    Parameter |      value       |     tolerance    |"
            "    units   | OK? | <-- Hybrid-MD\n",
            separator,
            self._tolerance_line(
                "|Ediff|",
                self.xyz_this_run.get_ediff(),
                self.get_tolerance("ediff"),
                "eV/at",
            ),
            self._tolerance_line(
                "max |Fdiff|",
                self.xyz_this_run.get_fmax(),
                self.get_tolerance("fmax"),
                "eV/Å",
            ),
            self._tolerance_line(
                "force RMSE",
                self.xyz_this_run.get_frmse(),
                self.get_tolerance("frmse"),
                "eV/Å",
            ),
            self._tolerance_line(
                "max |Vdiff|",
                self.xyz_this_run.get_vmax(),
                self.get_tolerance("vmax"),
                "eV",
            ),
            separator,
            f"{refit_str:>72} <-- Hybrid-MD\n",
        ]

        self.write_to_tmp_log(lines, append=True)

    def cumulative_error_table(self):
        """Construct a table showing the errors, intended for the SEED.castep file"""
        postfix = " <-- Hybrid-MD-Cumul"

        separator = (
            f" -------------+------------------+------------------"
            f"+------------+-----+{postfix}\n"
        )

        lines = [
            separator,
            f"  Cumulative RMSE       value            "
            f"count: {self.xyz_this_run.get_count():>8}"
            f"  |    units   |{postfix}\n",
            separator,
            self._tolerance_line_cumulative(
                "Energy", self.xyz_this_run.get_cumulative_energy_rmse(), "eV/atom"
            ),
            self._tolerance_line_cumulative(
                "Forces", self.xyz_this_run.get_cumulative_force_rmse(), "eV/Å"
            ),
            self._tolerance_line_cumulative(
                "Virial", self.xyz_this_run.get_cumulative_virial_rmse(), "eV"
            ),
            separator,
        ]

        self.write_to_tmp_log(lines, append=True)

    def get_tolerance(self, key: str):
        """Get tolerance - wrapper for settings object to overwrite checking of
        virials if needed

        Parameters
        ----------
        key

        """
        if not self.xyz_this_run.use_virial and key in ["vmax"]:
            return None
        return self.settings.tolerances.get(key)

    # --------------------------------------------------------------------------
    # XYZ IO
    def read_xyz(self):
        """Wrapper for XYZ-file related part"""
        self.xyz_this_run.read_xyz()

    def get_previous_data(self):
        """Read the previous data from files given"""
        if self.settings.refit.previous_data is None:
            return []

        frames = []
        for fn in self.settings.refit.previous_data:
            frames.extend(ase.io.read(fn, ":"))
        return frames

    @staticmethod
    def _bool_to_str(value):
        if value:
            return "Yes"
        return " No"


class SeedAwareState(ABC):
    """State object with seed specified. Base-class for the other ones."""

    def __init__(self, seed: str):
        self.seed = seed


class CarriedState(SeedAwareState):
    """State of calculation carried between calls

    Notes
    -----

    * This is being written to file, so the next time this program is called we can
    continue from where we left off. We are using a YAML file named
    `<seed>.hybrid-md-state.yaml`.

    * In order to ensure the pre and post steps are called in the right order,
    we are saving if the next step if pre step or not.

    """

    do_comparison = False  # print the error table
    do_update_model = False  # REFIT
    next_ab_initio = False
    next_is_pre_step = True  # do them in order
    last_check_step = -1
    current_check_interval = -1
    # for continuation runs we need this, e.g. initial steps
    continuation = False
    continuation_initial_steps = 0

    def __init__(self, seed: str):
        super().__init__(seed)
        self.state_filename = f"{self.seed}.hybrid-md-state.yaml"

    def dump(self):
        """Save state to file"""
        with open(self.state_filename, "w") as file:
            yaml.dump(self.carry_dict(), file)

    def load(self):
        """Load state from file"""
        with open(self.state_filename, "r") as file:
            values = yaml.safe_load(file)
        self.unpack_dump(values)

    def carry_dict(self):
        """Dictionary representation of the carried data

        This could be improved with a dataclass instead
        """
        return {
            "do_comparison": self.do_comparison,
            "do_update_model": self.do_update_model,
            "next_is_pre_step": self.next_is_pre_step,
            "next_ab_initio": self.next_ab_initio,
            "last_check_step": self.last_check_step,
            "current_check_interval": self.current_check_interval,
            "continuation": self.continuation,
            "continuation_initial_steps": self.continuation_initial_steps,
        }

    def unpack_dump(self, values: dict):
        """Unpack values from dumped dictionary

        Parameters
        ----------
        values

        """
        self.do_comparison = values.get("do_comparison")
        self.do_update_model = values.get("do_update_model")
        self.next_is_pre_step = values.get("next_is_pre_step")
        self.next_ab_initio = values.get("next_ab_initio")
        self.last_check_step = values.get("last_check_step")
        self.current_check_interval = values.get("current_check_interval")
        self.continuation = values.get("continuation", False)
        self.continuation_initial_steps = values.get("continuation_initial_steps", 0)

    def reset(self):
        """reset the info"""
        self.do_comparison = False
        self.do_update_model = False
        self.next_is_pre_step = True


class SubStateXYZ(SeedAwareState):
    """Represents XYZ structural data gathered during the accelerated MD"""

    def __init__(self, seed: str):
        super().__init__(seed)
        self.xyz_filename = f"{self.seed}.hybrid-md.xyz"

        self.use_virial = False

        # arrays: filled with data from structures
        self.len_atoms = 1
        self.atomic_numbers = np.zeros(1)

        self.energy_ff = np.zeros(1)
        self.forces_pp = np.zeros(1)
        self.virial_pp = np.zeros(1)

        self.energy_qm = np.zeros(1)
        self.forces_pw = np.zeros(1)
        self.virial_pw = np.zeros(1)

    # --------------------------------------------------------------------------
    # XYZ IO
    def read_xyz(self):
        """Reads XYZ file with QM observations"""
        frames = ase.io.read(self.xyz_filename, ":")

        # global atoms information
        self.len_atoms = len(frames[0])
        self.atomic_numbers = frames[0].get_atomic_numbers()

        # energy : (n_frames)
        self.energy_qm = self._unpack_info(frames, "QM_energy")
        self.energy_ff = self._unpack_info(frames, "FF_energy")

        # forces : (n_frames, n_atoms, 3)
        self.forces_pw = self._unpack_arrays(frames, "QM_forces")
        self.forces_pp = self._unpack_arrays(frames, "FF_forces")

        # virial : (n_frames, 6)
        if "QM virial" in frames[0].info.keys():
            self.use_virial = True
            self.virial_pw = self._unpack_info(frames, "QM_virial")
            self.virial_pp = self._unpack_info(frames, "FF_virial")

    # --------------------------------------------------------------------------
    # last step's error measures
    def get_ediff(self):
        """energy difference per atom"""
        return np.abs(self.energy_ff[-1] - self.energy_qm[-1]) / self.len_atoms

    def get_fmax(self, atomic_number: int = None):
        """Maximum force component difference on the latest observation"""
        if atomic_number is None:
            # species agnostic maximum force component difference
            return np.max(np.abs(self.forces_pp[-1] - self.forces_pw[-1]))

        return self._max_abs(
            self.forces_pp[-1, self._z_mask(atomic_number)]
            - self.forces_pw[-1, self._z_mask(atomic_number)]
        )

    def get_frmse(self, atomic_number: int = None):
        """Force RMSE on the latest observation"""
        if atomic_number is None:
            # species agnostic force component RMSE
            return self._rmse(self.forces_pp[-1] - self.forces_pw[-1])

        return self._rmse(
            self.forces_pp[-1, self._z_mask(atomic_number)]
            - self.forces_pw[-1, self._z_mask(atomic_number)]
        )

    def get_vmax(self):
        """maximum virial component difference"""
        return np.max(np.abs(self.virial_pp[-1] - self.virial_pw[-1]))

    # -----------------------------------------------------------------------------------
    # cumulative error measures

    def get_count(self):
        """Count: number of QM observations"""
        return len(self.energy_qm)

    def get_cumulative_energy_rmse(self):
        """Calculate cumulative RMSE of energies"""
        # cumulative energy RMSE per atom
        return self._rmse((self.energy_ff - self.energy_qm) / self.len_atoms)

    def get_cumulative_force_rmse(self, atomic_number: int = None):
        """Calculate cumulative RMSE of forces"""
        if atomic_number is None:
            # species agnostic force component RMSE
            return self._rmse(self.forces_pp - self.forces_pw)

        return self._rmse(
            self.forces_pp[:, self._z_mask(atomic_number)]
            - self.forces_pw[:, self._z_mask(atomic_number)]
        )

    def get_cumulative_virial_rmse(self):
        """Calculate cumulative RMSE of virials"""
        return self._rmse(self.virial_pw - self.virial_pp)

    # -----------------------------------------------------------------------------------
    # helper functions

    @staticmethod
    def _rmse(array):
        return np.sqrt(np.mean(np.square(array)))

    @staticmethod
    def _max_abs(array):
        return np.max(np.abs(array))

    def _z_mask(self, atomic_number: int):
        return self.atomic_numbers == atomic_number

    @staticmethod
    def _unpack_info(frames, key: str):
        return np.array([at.info[key] for at in frames])

    @staticmethod
    def _unpack_arrays(frames, key: str):
        return np.array([at.arrays[key] for at in frames])
