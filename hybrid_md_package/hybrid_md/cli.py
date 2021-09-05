#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.
"""
CLI of hybrid MD implementation.

The QM calculators should call this with the subcommands at correct points in the calculation.
"""

import sys

import click
from hybrid_md.refit import refit
from hybrid_md.state_objects import HybridMD
from hybrid_md.decision_making import SimpleDecisionMaker, PreStepReturnNumber

VERBOSE = True


@click.group("hybrid-md")
def main():
    """
    Hybrid MD main CLI. See the subcommands for more details.
    """


@main.command("initialise")
@click.argument("seed", type=click.STRING)
def initialise(seed):
    """Initialisation of the Hybrid MD run.

    Answer: integer exit code (three bits encoded)
    - i_exit &  1: read log and put in code's log
    - i_exit &  2: start with ab-initio (on True), else GAP
    """

    # create the initial state object
    state = HybridMD(seed)
    state.next_is_pre_step = True

    # write state to disc, only `next_is_pre_step` relevant though
    state.dump()

    if VERBOSE:
        print(
            f"Hybrid-MD: INIT Step, exit: {0 if state.num_initial_steps == 0 else 1} "
            f" -- num_initial_steps {state.num_initial_steps}",
        )

    # write the log for the .castep file
    state.io_initial_step_banner()

    # exit status -- log reading always ON
    if state.num_initial_steps == 0:
        sys.exit(1)
    else:
        sys.exit(3)


@main.command("pre-step")
@click.argument("seed", type=click.STRING)
@click.argument("md-iteration", type=click.INT)
def pre_step(seed, md_iteration):
    """Start of the MD-Hybrid step, called by md step

    Answer: integer exit code (three bits encoded)
    - i_exit &  1: read log and put in code's log
    - i_exit &  2: do ab-initio now
    - i_exit &  4: do ab-initio in next iteration (rest of loop)
    - i_exit &  8: do update cell
    - i_exit & 16: use ab-initio forces in MD
    """

    # state of object
    state = HybridMD(seed, md_iteration)
    state.load()
    state.reset()
    if not state.next_is_pre_step:
        raise RuntimeError(
            "Hybrid MD steps called in the wrong order, expected post-step"
        )

    # decision making & update of state object
    decision = SimpleDecisionMaker(state).get_step_kind(md_iteration)
    converter = PreStepReturnNumber(state)
    return_value = converter.push_state(decision)

    # dump state
    state.dump()

    if VERBOSE:
        print(
            f"Hybrid-MD:  PRE Step, exit:{return_value:4}, md_iteration:{md_iteration:3}  -->",
            converter.do_ab_initio,
            converter.next_ab_initio,
            converter.do_update_cell,
            converter.use_qm_forces,
        )

    # pass the result back
    sys.exit(return_value)


@main.command("post-step")
@click.argument("seed", type=click.STRING)
@click.argument("md-iteration", type=click.INT)
def post_step(seed, md_iteration):
    """End of the MD-Hybrid step, called by md.f90

    Answer: integer exit code (three bits encoded)
    - i_exit &  1: read log and put in code's log
    - i_exit &  2: read the GAP model again
    """

    # state of object
    state = HybridMD(seed, md_iteration)
    state.load()
    if state.next_is_pre_step:
        raise RuntimeError(
            "Hybrid MD steps called in the wrong order, expected pre-step"
        )

    # main logic
    if state.do_comparison:
        # 1. read output
        # 2. calculate E, F, S errors for this frame and cumulatively, force by species
        state.read_xyz()

        # 3. decide if we are fitting or not
        tolerance_met = state.check_tolerances()
        if not tolerance_met and state.can_update:
            state.do_update_model = True

        # 4. IO: errors of this step and cumulative ones as well
        state.error_table()
        state.cumulative_error_table()

    # 5. update model -- triggered by either pre or post step
    if state.do_update_model:
        refit(state)

    # save state
    state.next_is_pre_step = True
    state.dump()

    if VERBOSE:
        print(
            f"Hybrid-MD: POST Step, exit:"
            f"{int(state.do_update_model):4}, md_iteration:{md_iteration:3}"
        )

    # exit status -- log reading always ON
    exit_code = int(state.do_comparison) + 2 * int(state.do_update_model)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
    # in case main() does not exit
    sys.exit(-999)
