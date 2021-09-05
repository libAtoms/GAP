"""
Decision maker utilities for the Hybrid MD
"""

from abc import ABC, abstractmethod

from hybrid_md.state_objects import HybridMD, StepKinds


class DecisionMakerBase(ABC):
    """Base class for decision making

    used for pre-step decisions
    """

    def __init__(self, state: HybridMD):
        self.state = state

    @abstractmethod
    def get_step_kind(self, md_iteration: int) -> StepKinds:
        """ Perform the decision making in any way needed
        """
        ...


class SimpleDecisionMaker(DecisionMakerBase):
    """
    Simple uniform interval, like we had it implemented so far

    refactor this out of the cli module

    """

    def get_step_kind(self, md_iteration: int) -> StepKinds:
        """ Uniform checking, with optional initial DFT steps

        Parameters
        ----------
        md_iteration

        Returns
        -------
        step_kind
        """
        if md_iteration < self.state.num_initial_steps:
            return StepKinds.INITIAL

        if md_iteration == self.state.num_initial_steps:
            return StepKinds.LAST_INITIAL

        if (
            md_iteration - self.state.num_initial_steps
        ) % self.state.check_interval == 0:
            return StepKinds.CHECK

        return StepKinds.GENERIC


class AdaptiveDecisionMaker(DecisionMakerBase):
    """
    Adaptive decision making, with increasing and decreasing
    the interval based on accuracy target being met or not.
    """

    def get_step_kind(self, md_iteration: int):
        raise NotImplementedError


class PreStepReturnNumber:
    """
    Return number and changes to state object for Pre-step

    heavily mimicking the Fortran solution in terms of logic
    """

    do_ab_initio = False  # 2
    next_ab_initio = False  # 4
    do_update_cell = False  # 8
    do_comparison = False  # print the error table
    do_update_model = False  # REFIT
    use_qm_forces = False  # 16

    def __init__(self, state: HybridMD):
        self.state = state

    def push_state(self, value: StepKinds):
        """Actions to be done for the pre-step

        Parameters
        ----------
        value

        Returns
        -------

        """
        # set the internal booleans
        self._reset()
        self._set_internals(value)

        # save values into state
        self.state.next_is_pre_step = False
        self.state.next_ab_initio = self.next_ab_initio
        self.state.do_comparison = self.do_comparison

        if self.do_update_model:
            if self.state.can_update:
                self.state.do_update_model = self.do_update_model
            else:
                raise RuntimeError(
                    "Tried to update model, but input settings do not allow it!"
                )

        # return an integer
        return self._get_return_value()

    def _get_return_value(self):
        return_value = 0  # log reading off for this one
        for num, val in [
            (2, self.do_ab_initio),
            (4, self.next_ab_initio),
            (8, self.do_update_cell),
            (16, self.use_qm_forces),
        ]:
            return_value += num * int(val)

        return return_value

    def _reset(self):
        self.do_ab_initio = False  # 2
        self.next_ab_initio = False  # 4
        self.do_update_cell = False  # 8
        self.do_comparison = False  # print the error table
        self.do_update_model = False  # REFIT
        self.use_qm_forces = False  # 16

    def _set_internals(self, value: StepKinds):
        # sets internal values based on Step Kind

        if value == StepKinds.INITIAL:
            # initial steps with ab-initio, we have no PP model to work with
            self.do_ab_initio = True
            self.next_ab_initio = True
            self.use_qm_forces = True  # we don't have GAP forces yet
        elif value == StepKinds.LAST_INITIAL:
            # last initial ab-initio step, ask for model update now
            self.do_ab_initio = True
            self.do_update_model = True
            self.use_qm_forces = True  # we don't have GAP forces yet
        elif value == StepKinds.CHECK:
            # check steps during the run
            self.do_ab_initio = True
            # this is the only case when the previous electronic state needs changing
            self.do_update_cell = True
            self.do_comparison = True
        elif value == StepKinds.GENERIC:
            # this is a generic step with PP
            pass
        else:
            raise ValueError("Step kind given is not known")
