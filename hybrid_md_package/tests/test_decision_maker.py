#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021-2023.
import os
from pathlib import Path

from hybrid_md.decision_making import SimpleDecisionMaker
from hybrid_md.state_objects import HybridMD, StepKinds

dummy_input_content = """
can_update: true
check_interval: 5
num_initial_steps: 3
tolerances:
  ediff: 0.01  # eV
refit:
  e0_method: "average"
"""


def test_simple_decision_maker(tmp_path):
    # go to a temp dir
    os.chdir(tmp_path)

    # create a dummy input file
    input_yaml = Path("dummy.hybrid-md-input.yaml")
    with input_yaml.open(mode="w") as file:
        file.write(dummy_input_content)

    # state object - reads the above
    state = HybridMD("dummy")

    for md_iteration, result in [
        # nb. 0 is the initial step and calls to this are from 1
        (1, StepKinds.INITIAL),
        (2, StepKinds.INITIAL),
        (3, StepKinds.LAST_INITIAL),
        (4, StepKinds.GENERIC),
        (8, StepKinds.CHECK),
        (14376, StepKinds.GENERIC),
        (6343, StepKinds.CHECK),
    ]:
        # print("num:", md_iteration, result)
        assert SimpleDecisionMaker(state).get_step_kind(md_iteration) == result
