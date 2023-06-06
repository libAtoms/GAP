..
   Copyright (c) Tamas K. Stenczel, 2023.

.. _accelerated-aimd-input-reference:

********************
Input file reference
********************

The acceleration code has a single input file, which defines the parameters of the 
decision making and refitting. Everything else is handled in the ab-initio code.

The file should be called as ``<SEED>.hybrid-md-input.yaml`` where the ``<SEED>`` is 
the calculation's filename prefix supplied by the ab-initio code as well.

See the possible contents and a couple of example below. For runnable example, see
the simulation examples.

Core & Mandatory parameters
===========================

A minimal example: 

.. code-block:: yaml
    
    # this is a complete seed.hybrid-md-input.yaml files
    can_update: true
    check_interval: 5
    num_initial_steps: 2
    tolerances:
      ediff: 0.005  # eV

This sets the following:

* ``can_update: true``: the model can be updated (default false)

* ``check_interval: 5``: interval between ab-initio calculations (see :ref:`Fixed interval <fixed interval>`)

* ``num_initial_steps: 2``: the first 2 steps are ab-initio and used to train the first model

* ``tolerances:`` block of tolerances for accepting the result
    - ``ediff: 0.005``: 5 meV/atom is the maximum allowed for continuations with the same model 


Tolerance section
=================

.. code-block:: yaml
    
    # nb. partial input only
    tolerances:
      ediff: 0.01   # eV
      fmax: null    # eV/A
      frmse: 0.100  # eV/A
      vmax: null    # eV (virial)


Tolerances can be specified for the following quantitities:

* ``ediff``: energy per atom difference

* ``fmax``: force-component maximum absolute difference (in eV/Å)

* ``frmse``: RMSE of force components (in eV/Å)

* ``vmax``: maximum virial stress difference, in eV

At least one of these needs to be set, specifying ``null`` and omitting the parameter
are equivalent. When checking them, the ones turned on are all used and the result needs
to be below all of them to pass the check.

Accuracy-Adapted method section
===============================

Section of parameters for the :ref:`Accuracy-adapted checking interval mode <adaptive interval>`. The existence
of the block turns this decision making method on.

.. code-block:: yaml
       
    # nb. partial input only
    check_interval: 10
    adaptive_method_parameters:
      n_min: 10
      n_max: 5000
      factor: 1.3


* ``check_interval: 10``: initial checking interval to start with

* ``n_min`` & ``n_max``: bounds for the number of steps between ab-initio steps

* ``factor: 1.3``: the geometric factor for increasing and decreasing the checking interval

Refitting-specific settings
===========================

Parameters for controlling the GAP model to be trained when the model is updated.

You can use the pre-set SOAP parameters, or supply your own entirely.

By default, the trained model includes a 2-body and a SOAP descriptor for all species, both with 5Å cutoff. *n.b. This may fail if you only have a single atom of any element, since there will be no data to train the self-interaction 2-body part for this element.* 

.. code-block:: yaml
    
    # nb. partial input only
    refit:
      e0_method: "isolated"
      num_threads: 128
    
* ``e0_method`` sets the method for GAP to choose the element-specific constant shift.
    - default is ``average``
    - see :ref:`gap_fit page <isolated atom>` for further information

* ``num_threads: 128`` tells the program to use 128 OMP threads for the fitting, this is set temporarily only and is reset to the previous value after fitting. Does not affect the parallelism of the ab-initio code. *n.b. this only single-node fitting is supported for GAP models at the moment for the acceleration program.*

One can additionally set ``e0`` in this section for passing explicit values, as in the ``gap_fit`` program. 

.. code-block:: yaml
    
    # nb. partial input only
    refit:
      e0_method: "isolated"
      num_threads: 128

Previous data
-------------

Previous data can be included in the model training, through specifying a list of xyz filenames under the ``previous_data`` parameter. These are expected to have their energy, force, and virial results saved with the same keywords as how the ab-initio code is saving them: ``QM_energy``, ``QM_forces``, ``QM_virial``, respectively.
 


Descriptor Parameters
---------------------

You can choose the descriptor parameters in multiple ways:

* using pre-set SOAP parameter values: ``preset_soap_param: "fast" | "medium" | "accurate"``, which then includes a 2-body as well 

* using an explicit descriptor string: ``descriptor_str: "..."`` In this case the 2-body descriptor is not added.

Other parameters
----------------

The rest of the parameters, which are mostly mimicking the keywords of ``gap_fit``. See the :ref:`Defaults <defaults>` section for the for the default values of these.

* ``gp_name`` name of model XML file

* ``default_sigma`` default kernel regularisation 

* ``extra_gap_opts`` anything else to pass to ``gap_fit``

* ``function_name`` allows the user to specify any python function installed in their environment by module import path for handling the refitting of the model. If this is set, then this function is used and everything else is ignored (though that function can see these parameters). This is an advanced developer feature.


Defaults
--------

.. _defaults:

The following input file contains the default values set explicitly. *n.b. this is not allowing refitting, due to* ``can_update: false``

.. code-block:: yaml
    
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


