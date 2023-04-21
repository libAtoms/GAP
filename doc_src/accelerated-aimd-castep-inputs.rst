..
   Copyright (c) Tamas K. Stenczel, 2023.


********************
Castep Input File
********************

This section is about setting up castep's inputs to handle the acceleration.

Depending on which version you are using, the parameter names are a little different, see the two examples below.

Castep param file
=================

The settings governing the acceleration interface are located in the ``.param`` file's ``DEVEL_CODE`` block at
the moment, where the user needs to set a couple of things. However, this can be done once and copied over to
other calculations as-is as long as you keep your ``hybrid-md`` executable in your path rather than setting it
directly, and keeping the GAP model's name the same. Otherwise change these two things.

The rest of the parameters listed here are recommendations, which should make your life easy or were found to
be efficient.

DEVEL_CODE block
--------------------

Set the following:

``PP=T`` & ``MD: PP=T :ENDMD`` & ``PP_HYBRID=T`` which turn of the use of non ab-initio force calculators.

Turn on the usage of QUIP, initialisation arguments of a GAP model, and the path to it. Default is ``GAP.xml`` so
might just keep it as that.

.. code-block:: text

      pp:
          QUIP=T
          QUIP_PARAM_FILE=GAP.xml
          quip_init_args:IP GAP:endquip_init_args
      :endpp


Finally, set the decision maker program's path with the ``PP_HYBRID_EXEC`` block.
This is an "external executable" from Castep's point of view.

**Note on Castep v22:** The code's initial version was included in Castep's v22 academic release, where restarting
calculations was not possible yet, and parameter names had an ``MD_`` prefix. So if you are using v22, please
use ``MD_PP_HYBRID=T`` for turning on and ``MD_PP_HYBRID_EXEC`` for the executable.

Useful settings in .param
-------------------------

These are educated recommendations only, feel free to overwrite.

``NUM_BACKUP_ITER`` should be comparable to the number of expected MD steps between ab-initio calculations. This
controls the interval between checkpoints, from where restarting is possible. Anything after the checkpoint is
potentially lost if the job stops and needs restarting.

Use ``continuation: default`` for continuing a calculation where checkpoints were made and the seed name is the same.
Bear in mind, the MD steps is set to 0 after restart but time is incremented, so the directories with previous
GAP models will have indices accordingly.

``FINITE_BASIS_CORR: 0`` is recommended. This cannot be done at each ab-initio calculation, and the initial correction
value calculated at the beginning may not be applicable after a long MD run. If you know your system and your simulation
goals do disregard this.

Given that one is intending to use this method for longer MDs, ``md_sample_iter`` (interval of saving MD frames) should
be adjusted accordingly, since there is some IO cost associated with it as well as disc usage.

Use tight settings on energy and force convergence, since there provide the input data to your model. Worth starting
with looser settings for initial tests and then tightening for production runs.

In case of **variable cell calculations**, make sure to use very dense KPoint grid or tight spacing, otherwise the
training data given to GAP may not be consistent enough between frames and will yield sub-optimal model accuracy.

Example: Castep v23
======================

If you are using Castep v23 (academic release only), then then add the following to your `seed.param` file.
Note, the parameter prefixes are ``MD_PP_HYBRID`` unlike in v23.

.. code-block:: text

    %BLOCK DEVEL_CODE
      ! generally turns on PP, this is needed together with "PP_HYBRID=T"
      PP=T
      MD: PP=T :ENDMD

      ! settings of model called through QUIP
      pp:
          QUIP=T
          QUIP_PARAM_FILE=GAP.xml
          quip_init_args:IP GAP:endquip_init_args
      :endpp

      ! settings of PP Hybrid MD
      PP_HYBRID=T
      PP_HYBRID_EXEC:
        hybrid-md
      :ENDPP_HYBRID_EXEC
    %ENDBLOCK DEVEL_CODE



Example: Castep v22
======================

If you are using Castep v22 (academic release only), then then add the following to your `seed.param` file.
Note, the parameter prefixes are ``MD_PP_HYBRID`` unlike in v23.

.. code-block:: text

    %BLOCK DEVEL_CODE
      ! generally turns on PP, this is needed together with "MD_PP_HYBRID=T"
      PP=T
      MD: PP=T :ENDMD

      ! settings of model called through QUIP
      pp:
          QUIP=T
          QUIP_PARAM_FILE=GAP.xml
          quip_init_args:IP GAP:endquip_init_args
      :endpp

      ! settings of PP Hybrid MD
      MD_PP_HYBRID=T
      MD_PP_HYBRID_EXEC:
        hybrid-md
      :ENDMD_PP_HYBRID_EXEC
    %ENDBLOCK DEVEL_CODE



