import ase.io
import numpy as np
import yaml


class HybridMD:
    # info we want to carry
    do_comparison = False  # print the error table
    do_update_model = False  # REFIT
    next_ab_initio = False
    next_is_pre_step = True  # do them in order

    # tolerance
    tolerance_met = True
    tolerances = dict(
        ediff=None,  # in eV
        fmax=None,  # in eV/A
        frmse=None,  # in eV/A
        vmax=None,  # in eV -- not used this time
    )
    # if model can be updated, if False then we can only measure performance
    can_update = False

    # intervals
    check_interval = 1
    num_initial_steps = 0

    def __init__(self, seed: str, md_iteration: int = None):
        self.seed = seed

        self.state_filename = f"{self.seed}.hybrid-md-state.yaml"
        self.log_filename = f"{self.seed}.hybrid-md-log"
        self.input_filename = f"{self.seed}.hybrid-md-input.yaml"
        self.xyz_filename = f"{self.seed}.hybrid-md.xyz"

        self.use_virial = False  # if we are using virials in the calculations

        # read input -> tolerances, etc.
        self.read_input()

        # dummy arrays for results
        self.md_iteration = md_iteration
        self.len_atoms = 1
        self.atomic_numbers = np.zeros(1)

        self.energy_ff = np.zeros(1)
        self.forces_pp = np.zeros(1)
        self.virial_pp = np.zeros(1)

        self.energy_qm = np.zeros(1)
        self.forces_pw = np.zeros(1)
        self.virial_pw = np.zeros(1)

        # validation
        self.validate_settings()

    # -----------------------------------------------------------------------------------
    # IO for carried info
    def dump(self):
        # save state to file
        with open(self.state_filename, "w") as file:
            yaml.dump(self.carry_dict(), file)

    def load(self):
        # load state from file
        with open(self.state_filename, "r") as file:
            values = yaml.safe_load(file)
        self.unpack_dump(values)

    def carry_dict(self):
        return dict(
            do_comparison=self.do_comparison,
            do_update_model=self.do_update_model,
            next_is_pre_step=self.next_is_pre_step,
            next_ab_initio=self.next_ab_initio,
        )

    def unpack_dump(self, values: dict):
        self.do_comparison = values.get("do_comparison")
        self.do_update_model = values.get("do_update_model")
        self.next_is_pre_step = values.get("next_is_pre_step")
        self.next_ab_initio = values.get("next_ab_initio")

    def reset(self):
        # reset the info
        self.do_comparison = False
        self.do_update_model = False
        self.next_is_pre_step = True

    def read_input(self):
        # reads input settings of calculation
        with open(self.input_filename, "r") as file:
            data = yaml.safe_load(file)

        # unpack
        self.tolerances = data.get("tolerances", dict())
        self.can_update = data.get("can_update", False)
        self.check_interval = data.get("check_interval", 1)
        self.num_initial_steps = data.get("num_initial_steps", 0)

    def validate_settings(self):
        # any validation of the settings

        if self.num_initial_steps > 0 and not self.can_update:
            raise ValueError("Requesting initial DFT steps but cannot update model!")

    # -----------------------------------------------------------------------------------
    # Error table for IO

    def tolerance_line(self, name: str, value: float, tolerance, unit: str):
        if tolerance is not None:
            # checks tolerance as well
            self.tolerance_met = self.tolerances and tolerance > value

            # line to be printed
            ok = self.bool_to_str(tolerance > value)

            return f" {name:>12} | {value:16.8f} | {tolerance:16.8f} | {unit:10} | {ok:3} | <-- Hybrid-MD\n"
        else:
            return f" {name:>12} | {value:16.8f} |              Off | {unit:10} |     | <-- Hybrid-MD\n"

    def tolerance_line_cumulative(self, name: str, value: float, unit: str):

        return f" {name:>12}{value:30.8f}               | {unit:10} | <-- Hybrid-MD-Cumul\n"

    def check_tolerances(self):
        # todo: implement species specific check as well, and perhaps X out of Y to pass rather than all()

        # pass if none of them are given
        checked_tolerances = [True]

        for val, tol in [
            (self.get_ediff(), self.get_tolerance("ediff")),
            (self.get_fmax(), self.get_tolerance("fmax")),
            (self.get_frmse(), self.get_tolerance("frmse")),
            (self.get_vmax(), self.get_tolerance("vmax")),
        ]:
            if tol is None:
                continue

            # perform the check
            checked_tolerances.append(val < tol)

        return all(checked_tolerances)

    def error_table(self):
        # line formatting
        separator = " -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD\n"
        refit_str = "Refitting!" if self.do_update_model else "No Refit"

        # container for table -> write to file later
        lines = [
            f"\n      Hybrid-MD: MD iteration{self.md_iteration:8}                                    <-- Hybrid-MD\n"
            if self.md_iteration is not None
            else "\n",
            separator,
            "    Parameter |      value       |     tolerance    |    units   | OK? | <-- Hybrid-MD\n",
            separator,
            self.tolerance_line(
                "|Ediff|", self.get_ediff(), self.get_tolerance("ediff"), "eV/at"
            ),
            self.tolerance_line(
                "max |Fdiff|", self.get_fmax(), self.get_tolerance("fmax"), "eV/Å",
            ),
            self.tolerance_line(
                "force RMSE", self.get_frmse(), self.get_tolerance("frmse"), "eV/Å",
            ),
            self.tolerance_line(
                "max |Vdiff|", self.get_vmax(), self.get_tolerance("vmax"), "eV"
            ),
            separator,
            f"{refit_str:>72} <-- Hybrid-MD\n",
        ]

        with open(self.log_filename, "a") as file:
            file.writelines(lines)

    def cumulative_error_table(self):

        postfix = " <-- Hybrid-MD-Cumul"

        separator = f" -------------+------------------+------------------+------------+-----+{postfix}\n"

        lines = [
            separator,
            f"  Cumulative RMSE       value            count: {self.get_count():>8}  |    units   |{postfix}\n",
            separator,
            self.tolerance_line_cumulative(
                "Energy", self.get_cumulative_energy_rmse(), "eV/atom"
            ),
            self.tolerance_line_cumulative(
                "Forces", self.get_cumulative_force_rmse(), "eV/Å"
            ),
            self.tolerance_line_cumulative(
                "Virial", self.get_cumulative_virial_rmse(), "eV"
            ),
            separator,
        ]

        with open(self.log_filename, "a") as file:
            file.writelines(lines)

    def get_tolerance(self, key: str):
        if not self.use_virial and key in ["vmax"]:
            return None
        return self.tolerances.get(key, None)

    # -----------------------------------------------------------------------------------
    # XYZ IO
    def read_xyz(self):
        frames = ase.io.read(self.xyz_filename, ":")

        # global atoms information
        self.len_atoms = len(frames[0])
        self.atomic_numbers = frames[0].get_atomic_numbers()

        # energy : (n_frames)
        self.energy_qm = self.unpack_info(frames, "QM_energy")
        self.energy_ff = self.unpack_info(frames, "FF_energy")

        # forces : (n_frames, n_atoms, 3)
        self.forces_pw = self.unpack_arrays(frames, "QM_forces")
        self.forces_pp = self.unpack_arrays(frames, "FF_forces")

        # virial : (n_frames, 6)
        if "QM virial" in frames[0].info.keys():
            self.use_virial = True
            self.virial_pw = self.unpack_info(frames, "QM_virial")
            self.virial_pp = self.unpack_info(frames, "FF_virial")

    # -----------------------------------------------------------------------------------
    # last step's error measures
    def get_ediff(self):
        # energy difference per atom
        return np.abs(self.energy_ff[-1] - self.energy_qm[-1]) / self.len_atoms

    def get_fmax(self, z: int = None):
        if z is None:
            # species agnostic maximum force component difference
            return np.max(np.abs(self.forces_pp[-1] - self.forces_pw[-1]))
        else:
            return self.max_abs(
                self.forces_pp[-1, self.z_mask(z)] - self.forces_pw[-1, self.z_mask(z)]
            )

    def get_frmse(self, z: int = None):
        if z is None:
            # species agnostic force component RMSE
            return self.rmse(self.forces_pp[-1] - self.forces_pw[-1])
        else:
            return self.rmse(
                self.forces_pp[-1, self.z_mask(z)] - self.forces_pw[-1, self.z_mask(z)]
            )

    def get_vmax(self):
        # maximum virial component difference
        return np.max(np.abs(self.virial_pp[-1] - self.virial_pw[-1]))

    # -----------------------------------------------------------------------------------
    # cumulative error measures

    def get_count(self):
        return len(self.energy_qm)

    def get_cumulative_energy_rmse(self):
        # cumulative energy RMSE per atom
        return self.rmse((self.energy_ff - self.energy_qm) / self.len_atoms)

    def get_cumulative_force_rmse(self, z: int = None):
        if z is None:
            # species agnostic force component RMSE
            return self.rmse(self.forces_pp - self.forces_pw)
        else:
            return self.rmse(
                self.forces_pp[:, self.z_mask(z)] - self.forces_pw[:, self.z_mask(z)]
            )

    def get_cumulative_virial_rmse(self):
        return self.rmse(self.virial_pw - self.virial_pp)

    # -----------------------------------------------------------------------------------
    # helper functions

    @staticmethod
    def unpack_info(frames, key: str):
        return np.array([at.info[key] for at in frames])

    @staticmethod
    def unpack_arrays(frames, key: str):
        return np.array([at.arrays[key] for at in frames])

    @staticmethod
    def rmse(array):
        return np.sqrt(np.mean(np.square(array)))

    @staticmethod
    def max_abs(array):
        return np.max(np.abs(array))

    def z_mask(self, z: int):
        return self.atomic_numbers == z

    @staticmethod
    def bool_to_str(value):
        if value:
            return "Yes"
        else:
            return " No"
