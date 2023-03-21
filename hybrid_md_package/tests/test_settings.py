import os
import random
from abc import abstractmethod
from pathlib import Path
from string import ascii_letters

import pytest
from ase.io.extxyz import key_val_str_to_dict
from pytest import raises

from hybrid_md.settings import (
    Tolerances,
    SoapParamPreset,
    Refit,
    MainSettings,
    AdaptiveMethodSettings,
)


class TestTolerances:
    def test_init(self):
        Tolerances(ediff=0.1)

    def test_invalid_none_set(self):
        with raises(ValueError, match="At least one of the tolerances need to be set"):
            Tolerances()

    @pytest.mark.parametrize("key", ["ediff", "fmax", "frmse", "vmax"])
    @pytest.mark.parametrize("value", [-10000, -0.1, 0.0])
    def test_invalid_value(self, key, value):
        with raises(
            ValueError, match=f"Tolerance value `{key}` needs to be positive if set"
        ):
            Tolerances(**{key: value})

    def test_get_valid(self):
        # can get parameters which are valid
        tolerances = Tolerances(ediff=0.1)
        assert tolerances.get("ediff") == 0.1

    def test_get_none(self):
        tolerances = Tolerances(fmax=0.1)
        assert tolerances.get("ediff") is None
        assert tolerances.get("frmse") is None
        assert tolerances.get("vmax") is None

    def test_get_invalid_key(self):
        tolerances = Tolerances(fmax=0.1)
        with raises(KeyError):
            tolerances.get("dummy_key")

    def test_get_invalid_key_type(self):
        tolerances = Tolerances(fmax=0.1)
        with raises(KeyError):
            tolerances.get(1)  # type: ignore


class TestSoapParamPreset:
    def test_soap_str_fast(self):
        soap_params = SoapParamPreset.fast
        assert " n_sparse=500 " in soap_params.soap_str_stub
        assert " n_max=6 " in soap_params.soap_str_stub
        assert " l_max=3 " in soap_params.soap_str_stub

    def test_soap_str_medium(self):
        soap_params = SoapParamPreset.medium
        assert " n_sparse=1000 " in soap_params.soap_str_stub
        assert " n_max=8 " in soap_params.soap_str_stub
        assert " l_max=4 " in soap_params.soap_str_stub

    def test_soap_str_accurate(self):
        soap_params = SoapParamPreset.accurate
        assert " n_sparse=2000 " in soap_params.soap_str_stub
        assert " n_max=12 " in soap_params.soap_str_stub
        assert " l_max=6 " in soap_params.soap_str_stub

    @pytest.mark.parametrize(
        "soap_params",
        [SoapParamPreset.fast, SoapParamPreset.medium, SoapParamPreset.accurate],
    )
    def test_general_string(self, soap_params: SoapParamPreset):
        soap_str = soap_params.soap_str_stub
        assert soap_str.startswith("soap")

        soap_dict = key_val_str_to_dict(soap_params.soap_str_stub)

        # fixed keys
        for key, value in {
            "cutoff": 5.0,
            "cutoff_transition_width": 1.0,
            "atom_sigma": 0.5,
            "add_species": "True",
            "covariance_type": "dot_product",
            "zeta": 4,
            "sparse_method": "cur_points",
        }.items():
            assert key in soap_dict
            assert soap_dict[key] == value

        # non-fixed keys
        for key in {"n_max", "l_max", "n_sparse"}:
            assert key in soap_dict


class TestAdaptiveMethodSettings:
    def test_init(self):
        AdaptiveMethodSettings(n_min=10, n_max=100, factor=1.6)

    @pytest.mark.parametrize("n_min", [-1, 0])
    def test_invalid_n_min(self, n_min):
        with raises(
            ValueError, match="Adaptive method's minimum interval needs to be positive"
        ):
            AdaptiveMethodSettings(n_min=n_min, n_max=100, factor=1.6)

    @pytest.mark.parametrize("n_max", [-1, 0, 1])
    def test_invalid_n_max(self, n_max):
        with raises(
            ValueError,
            match="Adaptive method's maximum interval needs to be greater than 1",
        ):
            AdaptiveMethodSettings(n_min=1, n_max=n_max, factor=1.6)

    def test_invalid_min_max(self):
        with raises(
            ValueError,
            match="Adaptive method's maximum interval needs to be greater than it's "
            "minimum",
        ):
            AdaptiveMethodSettings(n_min=10, n_max=9, factor=1.6)

    @pytest.mark.parametrize("factor", [-10.0, 0.0, 1.0])
    def test_invalid_factor(self, factor):
        with raises(
            ValueError,
            match="Adaptive method's increment factor needs to be greater than 1",
        ):
            AdaptiveMethodSettings(n_min=1, n_max=10, factor=factor)


class TestRefit:
    def test_init(self):
        Refit()

    def test_defaults(self):
        refit = Refit()

        assert refit.previous_data == []
        assert refit.preset_soap_param == SoapParamPreset.fast
        assert refit.gp_name == "GAP.xml"
        assert refit.default_sigma == "0.005 0.050 0.1 1.0"
        assert refit.extra_gap_opts == "sparse_jitter=1.0e-8"
        assert refit.extra_gap_opts == "sparse_jitter=1.0e-8"
        assert refit.e0_method == "average"

        assert refit.descriptor_str is None
        assert refit.e0 is None
        assert refit.num_threads is None

    @pytest.mark.parametrize("num_threads", [2, 10, 999])
    def test_use_omp_true(self, num_threads):
        refit = Refit(num_threads=num_threads)
        assert refit.use_omp

    def test_use_omp_none(self):
        refit = Refit(num_threads=None)
        assert not refit.use_omp

    @pytest.mark.parametrize("num_threads", [-100, 0, 1])
    def test_use_omp_false(self, num_threads):
        refit = Refit(num_threads=num_threads)
        assert not refit.use_omp


# -----------------------------------------
# now a couple of tests for the main settings, re-using lots of boilerplate
# so arranged into a class hierarchy


class MainSettingsTestCase:
    @abstractmethod
    @pytest.fixture
    def input_file_content(self) -> str:
        """return input file (as YAML) and return as string"""

    @pytest.fixture(autouse=True)
    def tmp_input_file(self, tmp_path, input_file_content):
        """Sets up temporary directory & return written input file name"""

        # go to a temp dir
        os.chdir(tmp_path)

        # generate filename
        filename = Path(
            f"{''.join([random.choice(ascii_letters) for _ in range(10)])}.yaml"
        )

        # save content to file
        with filename.open(mode="w") as file:
            file.write(input_file_content)

        return filename

    @pytest.fixture
    def main_setting(self, tmp_input_file):
        # already read - use this for valid ones
        return MainSettings.read_input(tmp_input_file)

    @abstractmethod
    def test_inputs(self, *args):
        # implement this for checking what we read
        pass


class TestMainSettingsDefaults(MainSettingsTestCase):
    @pytest.fixture
    def input_file_content(self) -> str:
        return """
tolerances:
  ediff: 0.03  # eV
"""

    def test_inputs(self, main_setting):
        assert not main_setting.can_update
        assert main_setting.check_interval == 1
        assert main_setting.num_initial_steps == 0

        assert isinstance(main_setting.tolerances, Tolerances)
        assert main_setting.tolerances.ediff == 0.03
        assert main_setting.tolerances.fmax is None
        assert main_setting.tolerances.frmse is None
        assert main_setting.tolerances.vmax is None

        assert main_setting.adaptive_method_parameters is None
        assert isinstance(main_setting.refit, Refit)

        assert main_setting.refit.preset_soap_param == SoapParamPreset.fast


class TestMainSettingsDefaultsAllSpecified(TestMainSettingsDefaults):
    # running the same check as the one inherited from, though here all the
    # settings are set with their defaults
    @pytest.fixture
    def input_file_content(self) -> str:
        return """
    can_update: false
    check_interval: 1
    num_initial_steps: 0
    tolerances:
      ediff: 0.03
      fmax: null
      frmse: null
      vmax: null
    refit:
      function_name: null
      previous_data: []
      preset_soap_param: "fast"
      gp_name: "GAP.xml"
      default_sigma: "0.005 0.050 0.1 1.0"
      descriptor_str: null
      extra_gap_opts: "sparse_jitter=1.0e-8"
      e0: null
      e0_method: "average"
      num_threads: null
    """


class TestMainSettingsSimple(MainSettingsTestCase):
    @pytest.fixture
    def input_file_content(self) -> str:
        return """
can_update: true
check_interval: 5
num_initial_steps: 3
tolerances:
  ediff: 0.01  # eV
"""

    def test_inputs(self, main_setting):
        assert main_setting.can_update
        assert main_setting.check_interval == 5
        assert main_setting.num_initial_steps == 3

        assert isinstance(main_setting.tolerances, Tolerances)
        assert main_setting.tolerances.ediff == 0.01


class TestMainSettingsAdaptive(MainSettingsTestCase):
    @pytest.fixture
    def input_file_content(self) -> str:
        return """
can_update: true
check_interval: 10
tolerances:
  ediff: 0.01  # eV
adaptive_method_parameters:
  n_min: 10
  n_max: 5000
  factor: 1.3
"""

    def test_inputs(self, main_setting):
        assert main_setting.can_update
        assert main_setting.check_interval == 10
        assert main_setting.num_initial_steps == 0

        assert isinstance(main_setting.tolerances, Tolerances)
        assert main_setting.tolerances.ediff == 0.01

        assert isinstance(
            main_setting.adaptive_method_parameters, AdaptiveMethodSettings
        )

        assert main_setting.adaptive_method_parameters.n_min == 10
        assert main_setting.adaptive_method_parameters.n_max == 5000
        assert main_setting.adaptive_method_parameters.factor == 1.3


class TestMainSettingsRefit(MainSettingsTestCase):
    @pytest.fixture
    def input_file_content(self) -> str:
        return """
tolerances:
  ediff: 0.01  # eV
refit:
  function_name: null
  previous_data: [
    "filename1.xyz",
    "filename2.xyz",
]
  preset_soap_param: "fast"
  gp_name: "GAP888.xml"
  default_sigma: "0.01 0.09 1.5 1.09"
  descriptor_str: "Hello, this is a descriptor string"
  extra_gap_opts: "some extra GAP settings with numbers like 4.0"
  e0: "explicit E0 value"
  e0_method: "e0 method, any string OK here - GAP will need to validate"
  num_threads: 100
"""

    def test_inputs(self, main_setting):
        # same as otherwise
        assert not main_setting.can_update
        assert main_setting.check_interval == 1
        assert main_setting.num_initial_steps == 0
        assert isinstance(main_setting.tolerances, Tolerances)
        assert main_setting.tolerances.ediff == 0.01

        # refit settings
        refit = main_setting.refit
        assert refit.gp_name == "GAP888.xml"
        assert refit.previous_data == [
            "filename1.xyz",
            "filename2.xyz",
        ]
        assert refit.default_sigma == "0.01 0.09 1.5 1.09"
        assert refit.descriptor_str == "Hello, this is a descriptor string"
        assert refit.extra_gap_opts == "some extra GAP settings with numbers like 4.0"
        assert refit.e0 == "explicit E0 value"
        assert (
            refit.e0_method
            == "e0 method, any string OK here - GAP will need to validate"
        )
        assert refit.num_threads == 100
