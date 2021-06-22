import sys

import click

from hybrid_md.refit import refit
from hybrid_md.state_objects import HybridMD

VERBOSE = True


@click.group("hybrid-md")
def main():
    pass


@main.command("initialise")
@click.argument("seed", type=click.STRING)
def initialise(seed):
    """Initialisation of the Hybrid MD run.

    Answers:
    - 0: start with PP
    - 1: start with ab-initio
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

    # exit status
    if state.num_initial_steps == 0:
        sys.exit(0)
    else:
        sys.exit(1)


@main.command("pre-step")
@click.argument("seed", type=click.STRING)
@click.argument("md-iteration", type=click.INT)
def pre_step(seed, md_iteration):
    """Start of the MD-Hybrid step, called by md step

    Answer: integer exit code (three bits encoded)
    - i_exit & 1: do ab-initio now
    - i_exit & 2: do ab-initio in next iteration (rest of loop)
    - i_exit & 4: do update cell
    - i_exit & 8: use ab-initio forces in MD
    """

    # state of object
    state = HybridMD(seed, md_iteration)
    state.load()
    state.reset()
    if not state.next_is_pre_step:
        raise RuntimeError(
            "Hybrid MD steps called in the wrong order, expected post-step"
        )

    # mimicking the Fortran solution
    do_ab_initio = False  # 1
    next_ab_initio = False  # 2
    do_update_cell = False  # 4
    do_comparison = False  # print the error table
    do_update_model = False  # REFIT
    use_pw_forces = False  # 8

    if md_iteration < state.num_initial_steps:
        # initial steps with ab-initio, we have no PP model to work with
        do_ab_initio = True
        next_ab_initio = True
    elif md_iteration == state.num_initial_steps:
        # last initial ab-initio step, ask for model update now
        do_ab_initio = True
        do_update_model = True
    elif (md_iteration - state.num_initial_steps) % state.check_interval == 0:
        # check steps during the run
        do_ab_initio = True
        # this is the only case when the previous electronic state needs changing
        do_update_cell = True
        do_comparison = True
    else:
        # this is a generic step with PP
        pass

    # save values into state
    state.next_ab_initio = next_ab_initio
    state.next_is_pre_step = False
    state.do_comparison = do_comparison

    if do_update_model:
        if state.can_update:
            state.do_update_model = do_update_model
        else:
            raise RuntimeError(
                "Tried to update model, but input settings do not allow it!"
            )

    # dump state
    state.dump()

    return_value = 0
    for num, val in [
        (1, do_ab_initio),
        (2, next_ab_initio),
        (4, do_update_cell),
        (8, use_pw_forces),
    ]:
        return_value += num * int(val)

    if VERBOSE:
        print(
            f"Hybrid-MD:  PRE Step, exit:{return_value:4}, md_iteration:{md_iteration:3}  -->",
            do_ab_initio,
            next_ab_initio,
            do_update_cell,
            use_pw_forces,
        )

    # pass the result back
    sys.exit(return_value)


@main.command("post-step")
@click.argument("seed", type=click.STRING)
@click.argument("md-iteration", type=click.INT)
def post_step(seed, md_iteration):
    """End of the MD-Hybrid step, called by md.f90

    Answers:
    - 0: do nothing
    - 1: update PP model
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

    # exit status
    if state.do_update_model:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
    # in case main() does not exit
    sys.exit(-999)
