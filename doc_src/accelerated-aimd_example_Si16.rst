******************************************************
Si16 - Walkthrough of a Simple Calculation
******************************************************

This is a simple walkthrough of an experiment that one can try out, even on a laptop.

We will be running MD of a small Si cell.

First calculation
******************************************************

We will start out by running a short MD, starting with 3 steps of ab-initio calculations, training a model
on that, and later only doing ab-initio calculations at every 10th step. We are running a total of 50 MD steps.

Write the following in the input file ``si16.hybrid-md-input.yaml``.

.. code-block:: yaml

    can_update: true
    check_interval: 10
    num_initial_steps: 3
    tolerances:
      ediff: 0.01  # eV

You need ``si16.cell`` and ``si16.param`` files as well. See examples here, of 8 atom si.

``si16.cell``:

.. code-block:: text

    %BLOCK LATTICE_CART
     7.679180  0.000000  0.000000
     0.000000  7.679180  0.000000
     0.000000  0.000000  5.430000
    %ENDBLOCK LATTICE_CART
    
    %BLOCK POSITIONS_ABS
    Si  0.000000  0.000000  0.000000
    Si  1.919795  0.000000  1.357500
    Si  1.919795  1.919795  2.715000
    Si  0.000000  1.919795  4.072500
    Si  0.000000  3.839590  0.000000
    Si  1.919795  3.839590  1.357500
    Si  1.919795  5.759385  2.715000
    Si  0.000000  5.759385  4.072500
    Si  3.839590  0.000000  0.000000
    Si  5.759385  0.000000  1.357500
    Si  5.759385  1.919795  2.715000
    Si  3.839590  1.919795  4.072500
    Si  3.839590  3.839590  0.000000
    Si  5.759385  3.839590  1.357500
    Si  5.759385  5.759385  2.715000
    Si  3.839590  5.759385  4.072500
    %ENDBLOCK POSITIONS_ABS



``si16.param``:

.. code-block:: text

    ###############################################################
    # functional & general settings
    ###############################################################
    xc_functional   LDA
    basis_precision PRECISE  ! basis set size from presets
    fix_occupancy   true
    opt_strategy    speed    ! faster runtime with more memory use

    # allow continuation
    NUM_BACKUP_ITER = 10     ! interval between checkpoints (.check file)
    continuation: default    ! continuation if .check file exists

    ###############################################################
    # MD settings
    ###############################################################
    task             = molecular dynamics

    md_ensemble      = NVT
    md_thermostat    = Langevin
    md_num_iter      = 50      ! number of MD steps we will take
    md_temperature   = 300 K
    md_sample_iter   = 10      ! interval between dumping MD frames
    md_delta_t       = 1 fs

    ###############################################################
    # Devel code: contains the `PP_HYBRID` method settings
    ###############################################################
    %BLOCK DEVEL_CODE
      ! turn on the PP & acceleration modules
      PP=T
      MD: PP=T :ENDMD
      PP_HYBRID=T

      ! settings of model called through QUIP
      pp:
          QUIP=T
          QUIP_PARAM_FILE=GAP.xml
          quip_init_args:IP GAP:endquip_init_args
      :endpp

      ! settings of PP Hybrid MD

      PP_HYBRID_EXEC:
        hybrid-md
      :endPP_HYBRID_EXEC
    %ENDBLOCK DEVEL_CODE

Let's run this!

.. code-block:: bash

    mpirun -n 4 castep.mpi si-0

The code should produce the following:

* standard Castep outputs:
    * ``si16.castep`` main output file (worth reading)
    * ``si16.md`` MD trajectory
    * a few more castep files, including a bibliography of the parts used
* ``GAP.xml`` model and associated ``GAP.xml.sparseX_<time>`` sparse point files. *This is the latest model trained*
* ``train.xyz``: the assembled training set of the last model
* ``GAP.xml.fit-stderr.<time>`` & ``GAP.xml.fit-stdout.<time>`` are the output streams of the ``gap_fit`` program. Worth taking a look at.
* ``GAP_model-step000..-at-<time>/`` directories with previous GAP models. The number in the name is when the model was replaced by the next one
* ``si16.hybrid-md.xyz``: all ab-initio observations from this calculation, where forces are gradients of the energy and don't include thermostat components.


Accuracy-Adapted Checking Interval
******************************************************

This was all well, but we can do much much better: once the model has some more training data then we likely don't need to
perform an ab-initio calculation every 10 steps. We can use our other decision making method, where the number of steps between
checks are adapting to the accuracy of the model. See more at :ref:`Accuracy-adapted checking interval mode <adaptive interval>`.

Let's move to a new directory, change the ``si16.hybrid-md-input.yaml`` to the following, increase the number of MD steps we are taking
(``md_num_iter = 200`` in ``si16.param``), and run the program again. Best done in a different directory.

.. code-block:: text

    can_update: true
    check_interval: 10
    num_initial_steps: 3
    tolerances:
      ediff: 0.01  # eV
    adaptive_method_parameters:   # <---
      n_min: 5                    # <---
      n_max: 100                  # <---
      factor: 1.5                 # <---

This should produce a more interesting output, where the number of steps between ab-initio calculations is not fixed, but lowe initially and
then increases towards the end of the calculation. We can look at the output, where an error table is shown:

.. code-block:: text

         Hybrid-MD: MD iteration      49                                    <-- Hybrid-MD
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD
       Parameter |      value       |     tolerance    |    units   | OK? | <-- Hybrid-MD
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD
         |Ediff| |       0.00228307 |       0.01000000 | eV/at      | Yes | <-- Hybrid-MD
     max |Fdiff| |       1.24702641 |              Off | eV/Å       |     | <-- Hybrid-MD
      force RMSE |       0.35238182 |              Off | eV/Å       |     | <-- Hybrid-MD
     max |Vdiff| |       0.00000000 |              Off | eV         |     | <-- Hybrid-MD
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD
                                                                   No Refit <-- Hybrid-MD
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD-Cumul
     Cumulative RMSE       value            count:       11  |    units   | <-- Hybrid-MD-Cumul
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD-Cumul
          Energy                    0.03409695               | eV/atom    | <-- Hybrid-MD-Cumul
          Forces                    0.82637217               | eV/Å       | <-- Hybrid-MD-Cumul
          Virial                    0.00000000               | eV         | <-- Hybrid-MD-Cumul
    -------------+------------------+------------------+------------+-----+ <-- Hybrid-MD-Cumul
                  Hybrid-MD: INCREASE interval to      7 at iter       49   <-- Hybrid-MD-Adapt

Here we can see a checking step, where the tolerance was met, so we are not refitting and are continuing with the previous model,
allowing for a longer interval until the next check. Otherwise we would refit and decrease the interval.

One can see all decisions made quickly by ``grep CREAS si16.castep``, which should produce something like:

.. code-block:: text

    Hybrid-MD: DECREASE interval to      6 at iter       13   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       19   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       24   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       29   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       34   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       39   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       44   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to      7 at iter       49   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     10 at iter       56   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     15 at iter       66   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     22 at iter       81   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to     14 at iter      103   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     21 at iter      117   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     31 at iter      138   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     46 at iter      169   <-- Hybrid-MD-Adapt

Where we can see how the interval quickly increased after some initial burn-in time. It is worth experimenting with the settings here,
seeing what is the maximum ``factor`` and ``n_max`` value which produces stable MD, so one can speed the calculation up as much as possible.

Continue a calculation
******************************************************

If we don't change anything in the directory where we ran our precious calculation, we can run the same again and
continue where we left off. Well, at the last checkpoint. You can increase the number of MD steps to say 1000 now and
see the results.

As you can see below, the interval has hit the ceiling, then the model was re-trained once at step 384 and slowed down
there, but hit the ceiling after again. Please note, that the number of steps is reset to 0 at the start of the MD run
by Castep and only time is incremented, while concatenating to the same output file.

.. code-block:: text

    Hybrid-MD: DECREASE interval to      6 at iter       13   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       19   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       24   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       29   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       34   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       39   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to      5 at iter       44   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to      7 at iter       49   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     10 at iter       56   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     15 at iter       66   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     22 at iter       81   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to     14 at iter      103   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     21 at iter      117   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     31 at iter      138   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     46 at iter      169   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     69 at iter       15   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter       84   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      184   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      284   <-- Hybrid-MD-Adapt
    Hybrid-MD: DECREASE interval to     66 at iter      384   <-- Hybrid-MD-Adapt   # n.b. refitting the model here
    Hybrid-MD: INCREASE interval to     99 at iter      450   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      549   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      649   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      749   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      849   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    100 at iter      949   <-- Hybrid-MD-Adapt

From here, you can restart the calculation again any time and change the settings of the adaptive interval. Try what
happens if you change the ceiling it to something much longer and let it run a longer calculation.

Re-use previous data & Set your own GAP parameters
******************************************************

We can see above how to go from a single structure and a few lines of simple settings to close to 100x accelerated ab-initio MD
plus a GAP model we can use for anything else as well.

Let's see how we can refine the model. We are running MD of Si here, we might want to tune the GAP model's settings to
the system. Along with that, we will allow the GAP model fitting to happen in parallel as well.

Finally, we will include the previous calculation's ab-initio observations, by copying the ``si16.hybrid-md.xyz`` file into
the new job directory as ``previous_si16.hybrid-md.xyz``.

Let's update the input file to the following:

.. code-block:: yaml

    can_update: true
    check_interval: 10
    num_initial_steps: 1
    tolerances:
      ediff: 0.01  # eV
    adaptive_method_parameters:
      n_min: 5
      n_max: 200
      factor: 1.5
    refit:
      e0_method: "average"
      num_threads: 4
      descriptor_str: "distance_Nb order=2 n_sparse=20 cutoff=5.5 cutoff_transition_width=1.0 compact_clusters covariance_type=ard_se theta_uniform=1.0 sparse_method=uniform f0=0.0 add_species=T delta=1.0 : soap n_sparse=1000 n_max=10 l_max=4 cutoff=5.0 cutoff_transition_width=1.0 atom_sigma=0.5 add_species=True covariance_type=dot_product zeta=4 sparse_method=cur_points delta=3.0 "
      previous_data: [
        "previous_si16.hybrid-md.xyz"
      ]
      default_sigma: "0.005 0.1 0.05 1.0"

Let's see what we have done in the refitting section

* increased the maximum number of steps between ab-initio evaluations to 200
* ``num_threads: 4`` allows the ``gap_fit`` program to be executed with 4 OMP threads, while the main Castep program's 4 MPI processes are waiting
* the ``descriptor_str`` line is settings the descriptors that GAP uses, these are cut down versions of the `GAP-18 Si model <http://doi.org/10.1103/PhysRevX.8.041048>`_
* the GP's kernel regularisation is updated with ``default_sigma``,  :ref:`see usage of gap_fit <gap_fit>`

Additionally, we can rattle the Si structure at the start of the calculation. Castep can do this by adding ``positions_noise 0.1 Ang``
to the ``si16.cell`` input file.

Having done all of these, run the calculation just as before and see the results.

You will notice in the ``gap_fit`` program outputs that 4 threads are used and that the previous data is included.
The calculation will start with a single ab-initio evaluation and model training, which is then followed by the
interval between checking steps quickly increasing.

See the relevant lines from ``si16.castep``:

.. code-block:: text

    Hybrid-MD: INCREASE interval to     15 at iter       11   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     22 at iter       26   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     33 at iter       48   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     49 at iter       81   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to     73 at iter      130   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    109 at iter      203   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    163 at iter      312   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    200 at iter      475   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    200 at iter      675   <-- Hybrid-MD-Adapt
    Hybrid-MD: INCREASE interval to    200 at iter      875   <-- Hybrid-MD-Adapt


What next
******************************************************

You have learnt how to use the method supplied and were able to perform some short MD with a model trained on the fly.
You can expand this, by running more meaningful longer MD, or trying different temperatures, or include defects in the
structure and investigate those.

The real test and valuable next step is trying this out for your own systems of interest. Think of what small simulation
you could do which you may already have the structures and settings for. Try plugging it into Castep and turn the
MD acceleration on, follow the ideas and steps learnt here, and see if anything useful comes out of it.

If you need any assistance, don't hesitate to contact the authors of this package.
