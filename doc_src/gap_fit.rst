*******************
The gap_fit program
*******************


In order to fit an interatomic potential to data with ``gap_fit``, you need

#. Data
#. A definition of descriptors and kernels that serve as basis functions
#. Some parameters that control the least squares fit

We will go over each of these in turn, and reference the relevant ``gap_fit`` options.


Data
****

The input file to ``gap_fit`` is a series of atomic structures (they
could be molecules, or periodic systems), with atomic numbers,
cartesian positions of the atoms, and some associated data. The
required file format is **extended XYZ**, which is similar to a
concatenation of ordinary XYZ, except that the second line (after the
first line that just contains the number of atoms) is not free format,
but has to be a series of ``key=value`` pairs, with some mandatory keys
and other restrictions. See the detailed definition [here].

The data associated with each structure is an energy, optionally
forces, and optionally virial stress. All of these are specified in
the extended XYZ file. Following the definition of the extended XYZ
file, the energy needs to be in units of eV, the forces in units of
eV/A and the virial stress in units of eV (this is the regular stress
multiplied by the volume). Watch out for the definition of the sign of
the stress! For example VASP uses the opposite sign convention to QUIP.
The atomic positions need to be in Angstroms, and the periodic unit cell is specified by three cartesian
lattice vectors.

An example structure in extended XYZ format::

  1
  config_type=isolated_atom gap_energy=-157.72725320 dft_virial="0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000" dft_energy=-158.54496821 nneightol=1.20000000 pbc="T T T" Lattice="20.00000000       0.00000000       0.00000000       0.00000000      20.00000000       0.00000000       0.00000000       0.00000000      20.00000000" Properties=species:S:1:pos:R:3:Z:I:1:map_shift:I:3:n_neighb:I:1:gap_force:R:3:dft_force:R:3
  Si             10.00000000     10.00000000     10.00000000      14      -1      -1      -1       0      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000```

The structure above has lots of data in it, not all of it is useful
for fitting (e.g. ``gap_energy`` happens to be a prediction by some
model, not relevant for making models, or ``nneightol`` which was
inserted by some earlier manupulation of the file). Note how the
value of the mandatory ``Properties`` key specifies the names and types
of the columns in the rest of the configuration. ``S`` stands for
'string', ``R`` stands for 'real', ``I`` is for 'integer' and ``L`` is for
logical. Columns ``species`` and ``pos`` (atom types and cartesian
positions, respectively) are mandatory.

**All structures in the extended XYZ format need a periodic lattice
unit cell (and defined by the ``Lattice`` key), even if a particular
direction (or any direction) is not supposed to be periodic. Just use
a lattice vector which is bigger than twice the cutoff of any
potential you plan to create or use.**

To get good fits, it is extremely helpful of the energy and the forces
are **consistent**, i.e. up to discretisation error (e.g. in k-space in
case of a DFT code), the forces are the negative derivatives of the
energy. E.g. if the data comes from density functional theory with a
finite electronic temperature (as is almost always used to help
convergence, even in the case of insulators, when there are defects),
it is the *free energy* (i.e. including the electronic entropy) that
is consistent with the Hellman-Feynman forces, rather than the
potential energy. In ASE, this is obtained by specifying the
``force_consistent=true`` option to ``Atoms.get_potential_energy()``.

Fits are tremendously improved by including forces. For 3N atoms,
there is only one energy per structure, but 3N force components
(obtained at almost the same cost as the single energy), so therefore
contributing a lot more data. Virial stresses are again really
important, especially for obtaining elastic constants for periodic
solids.  This is because the forces are always zero for an isotropic
change away from equilibrium, and give no information on the curvature
of the potential energy surface.

The type of data included can vary from configuration to
configuration, there is no restriction on having to include the same
type of data for all configurations.

Data field names
####################

Conceptually, the data the potential is fitted to are the energy and
its derivatives (forces and stresses), but the names with which these
data are referred to in the XYZ file are arbitrary (but have to be
consistent in the whole XYZ file). Unless otherwise indicated,
``gap_fit`` looks for data with keys ``energy``, ``force`` and ``virial``, and
this can be controlled using the command line options
``energy_parameter name``, ``force_parameter_name``, and
``virial_parameter_name``, respectively. This facilitates specifying
several types of energies and related quantities in the same XYZ, and
selected the one you want to fit to on the command line. For example,
energies computed with different DFT functionals can be stored in the
same XYZ file, or even with a completely different method such as
QMC.

The above example has the keys ``dft_energy`` for example, so when it is
used in fitting, the appropriate key names need to be specified.

Isolated atom
###################

Most of the time, it makes sense to fit an interatomic potential to
the **binding energy**, i.e. the total energy minus the energy of
isolated atoms.  There are several mechanisms in ``gap_fit`` to aid
this. The default behaviour is to look for a configuration among the
input structures in the input XYZ file that contains only one atom,
and has an energy specified for it. This energy will then be taken
and used to calculate the binding energies of all configurations
before the fit is made. We give it the name ``e0``, and it is also
stored in the XML file that specifies the potential. When predictions
are made, the ``e0`` value is added back to the prediction of the
binding energy. If the input XYZ contains multiple atom types, each
one needs to appear once by itself among the configurations.

Alternatively, if isolated atoms are not part of the input file, a
value for ``e0`` can be specified on the command line, in case of
multiple elements as a series of numbers:
``e0={H:123.456:O:789.10...}``. There is also an option to have an ``e0``
value be taken as the average per-atom energy of the input data. This
only makes sense if there is only one type of atom.

Note that if you create a dataset where the specified energies are
already binding energies, you should still specify 0.0 as the isolated atom
energy (either via a configuration, preferred, or on the command
line), otherwise the interatomic potential might (depending on the
descriptor) not formally be zero for isolated atoms.



Descriptors and Kernels
***********************


SOAP hyperparameters
########################

This section contains notes about choosing hyperparameters for the SOAP descriptor.

cutoff
************

Every finite range potential can be cast in the form of a sum over site
energies or atomic energies, and the cut-off radius defines the range of
this local term. The actual interaction range is of course twice the
cut-off radius, because atoms up to this distance can potentially
interact with one another via a many-body term centered on an atom in
between them.  When we approximate a quantum mechanical potential
energy (which is not formally local) using a local site energy with
cut-off radius, the error we necessarily incur can be
characterised in theform of a force variance.




Global fit parameters
*********************

blah blah

blah blah

Command line example
********************

Here is an annotated fitting example.

.. code-block:: bash

  gap_fit atoms_filename=database.xyz # input data in extended XYZ format
    gap={                              # start of descriptor and kernel spec
    distance_Nb                       # first descriptor is interatomic distance based
    order=2                           # descriptor is 2-body (i.e. a pair potential)
    cutoff=5.0                        # distance cutoff in the kernel, in Angstrom
    n_sparse=15                       # number of representative points, M in Sec. II
    covariance_type=ard_se            # form of kernel: squared exponential (Gaussian)
    delta=2.0                         # scaling of kernel, per descriptor, here it is per atom pair, in eV
    theta_uniform=2.5                 # length scale in Gaussian kernel in Angstrom
    sparse_method=uniform             # choice of representative points, here a uniform grid up to the cutoff
    compact_clusters=T                # how cutoff is applied, here a spherical manner around each atom
    :                                 # separator between descriptors
    soap                              # second descriptor is a SOAP
    l_max=6 n_max=12                  # number of angular and radial basis functions for SOAP
    atom_sigma=0.5                    # Gaussian smearing width of atom density for SOAP, in Angstrom
    cutoff=5.0                        # distance cutoff in the kernel, in Angstrom
    radial_scaling=-0.5               # exponent of atom density scaling, power of distance
    cutoff_transition_width=1.0       # distance across which kernel is smoothly taken to zero, in Angstrom
    central_weight=1.0                # relative weight of central atom in atom density for SOAP
    n_sparse=8000                     # number of representative points, M in Sec. II
    delta=0.2                         # scaling of kernel, per descriptor, here for SOAP it is per atom, in eV
    covariance_type=dot_product       # form of kernel
    zeta=4                            # power kernel is raised to - together with dot_product gives a polynomial kernel
    sparse_method=cur_points          # choice of representative points, here CUR decomposition of descriptor matrix
    }                                 # end of descriptor and kernel spec
   default_sigma={0.002 0.2 0.2 0.0}  # default regularisation corresponding to energy, force, virial, hessian
   config_type_sigma={                # start of per configuration-group regularisation spec, using groups defined in the input data file
    isolated_atom:0.0001:0.01:1.0:0.0:
    rss_rnd:0.03:0.4:0.5:0.0:
    rss_005:0.02:0.3:0.4:0.0:
    rss_200:0.01:0.2:0.2:0.0:
    rss_3c:0.005:0.1:0.1:0.00:
    cryst_dist:0.0003:0.03:0.05:0.00:
    cryst_dist_hp:0.005:0.1:0.1:0.0:
    liq_P4:0.003:0.3:0.5:0.0:
    liq_network:0.003:0.3:0.5:0.0:
    2D:0.001:0.03:0.05:0.0:
    ribbons:0.01:0.5:0.2:0.0
    }                                 # end of per configuration-group regularisation spec
   energy_parameter_name=energy       # name of the key in the input data file corresponding to the total energy
   force_parameter_name=forces        # name of the key in the input data file corresponding to the forces
   virial_parameter_name=virial       # name of the key in the input data file corresponding to the virial stress
   sparse_jitter=1.0e-8               # extra diagonal regulariser
   do_copy_at_file=F                  # copy input data into potential XML file?
   sparse_separate_file=T             # write representative point data into a separate file not in the main potential XML
   gp_file=gap.xml                    # name of output potential XML file
   core_param_file=P_r6_innercut.xml  # name of XML file containing the baseline potential (QUIP format)
   core_ip_args={IP Glue}             # initialisation string to call baseline potential


Command line options
********************

See `gap_fit --help` for description of options.

GAP options
***********

This does not exist any more: `quippy.gap_fit_parse_gap_str`.

``sparse_method`` options are:
 - RANDOM: default, chooses n_sparse random datapoints
 - PIVOT: based on the full covariance matrix finds the n_sparse "pivoting" points
 - CLUSTER: based on the full covariance matrix performs a k-medoid clustering into n_sparse clusters, returning the medoids
 - UNIFORM: makes a histogram of the data based on n_sparse and returns a data point from each bin
 - KMEANS: k-means clustering based on the data points
 - COVARIANCE: greedy data point selection based on the sparse covariance matrix, to minimise the GP variance of all datapoints
 - UNIQ: selects unique datapoints from the dataset
 - FUZZY: fuzzy k-means clustering
 - FILE: reads sparse points from a file
 - INDEX_FILE: reads indices of sparse points from a file
 - CUR_COVARIANCE: CUR, based on the full covariance matrix
 - CUR_POINTS: CUR, based on the datapoints
