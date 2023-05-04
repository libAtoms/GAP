..
   Copyright (c) Tamas K. Stenczel, 2023.

.. _accelerated-aimd-decision-making:

Decision making
***************

Decision making methods are governing which forces are used at each step, and when to
refit a model. This is a key part of any acceleration and on-the-fly model improvement
cycle.

We currently supply two methods, which directly compare the GAP and ab-initio energy 
and forces. One is using a fixed interval for performing ab-initio calculations, while
the other one adapts to the accuracy of the model. Both make use of the tolerances set
by the user in the input file for deciding when to retrain the model.

For general calculations, particularly ones starting from scratch, **we recommend using 
the accuracy-adapted decision making method.**

Initial steps
=============

Regardless of the decision making method, the method allows for a number of initial 
ab-initio MD steps to be performed. This is useful for starting from scratch, where no
previous model is available.

In order to activate this set ``num_initial_steps > 0`` in the input file.

Fixed Interval Decision Making
==============================

.. _fixed interval:

Ab-initio calculations are performed at fixed intervals. Set the fixed interval with the
``check_interval`` keywords in the input file.

If any of the tolerance criteria are not met at this point, then the model is re-trained
with the available data.

Accuracy-Adapted Checking Interval
==================================

.. _adaptive interval:

This method adapts the interval between ab-initio steps according to a geometric factor
parameter, within user-specified upper and lower bounds. You need to specify the starting
interval, factor, and the bounds. 

Upon meeting the desired tolerances (i.e. being below the limits) at checking steps, the
interval until the next one is increased, otherwise the model is re-trained with the
newly available observations and the interval decreased by the factor specified.

This method is turned on by supplying the ``adaptive_method_parameters: ...`` input 
block in the decision making input file.

Further development
===================

Further developments are planned and are very welcome. There are lots of open avenues that can be explored with the framework provided here. Some of these the current developers have thought of, but that is surely a small subset of what people can and want to use this for. 

See a few examples below. If you are interested in collaborating on any of these, or have any other idea then please feel free to reach out, we are more than happy to work together or help others get started with the code.

* usage of bonding/structural information: trigger refitting when atomic connectivity is changed, or when the dynamics gets closer to potential bond breaking. Might be useful for modelling reactive systems, etc.

* use similarity to a database, say through SOAP descriptor distance, trigger data acquisition when this hits a threshold of moving further from the previous data

* integration of other ML and ab-initio codes (some in the making already)

