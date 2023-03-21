.. _accelerated-aimd-installation:

Installation with CASTEP
************************

You will need the following program components:

#. QUIP with GAP
    - MPI-only library for linking to Castep
    - OpenMP ``gap_fit`` program
#. CASTEP: Academic release v22 and above
#. Python installation of ``hybrid-md`` - bundled with GAP

The two builds of QUIP are needed due to the different parallelism supported by the fitting and evaluation.

Linux - general setup
#####################

These are general instructions for a Linux machine, assuming the libraries needed
for the respective packages are installed. See documentation of QUIP & CASTEP for
these.

QUIP+GAP installations
----------------------

Clone the code, make sure to include the subpackages as well using the `--recursive` flag.

.. code-block:: bash

    # clone the QUIP repository
    git clone --recursive https://github.com/libAtoms/QUIP.git
    cd QUIP


Build the OpenMP version for fitting

.. code-block:: bash
    
    # arch: any compiler + OpenMP
    export QUIP_ARCH=linux_x86_64_gfortran_openmp

    # configure: build - Make sure to say yes to GAP
    make config

    # build: the gap_fit program
    make gap_programs

    # copy the executable to your desired installation directory
    cp ./build/$QUIP_ARCH/gap_fit ...your bin dir...

Build an MPI library version for linking to Castep. Make sure to use the same compiler
as for the Castep installation later on.

.. code-block:: bash

    # still in the QUIP directory created above 

    # remember this directory for configuring CASTEP
    export QUIP_ROOT=$(pwd)

    # arch: any compiler + MPI
    export QUIP_ARCH=linux_x86_64_gfortran_openmpi

    # configure build - Make sure to say yes to GAP
    make config

    # build: only the library
    make libquip.a

Castep
------

See the installation notes and documentation in your CASTEP source distribution,
the relevant part for this program is the inclusion of QUIP, which you can do by
specifying ``QUIP=system`` in the main ``Makefile`` plus providing the
``QUIP_ROOT`` (root directory of QUIP repo) and ``QUIP_ARCH`` (the MPI one).
These can be directly placed into the ``Makefile`` of CASTEP as well to be
exported for easier resume of build or rebuild later. 

.. code-block:: makefile
    
    # make sure these are specified
    COMMS_ARCH := mpi
    QUIP := system

    # example of explicitly setting the location of QUIP
    export QUIP_ROOT=...directory from above...
    export QUIP_ARCH=linux_x86_64_gfortran_openmpi

Having set these, just build and install CASTEP. 


``hybrid-md`` Python package
----------------------------

Simply install in the package from source, found in the ``GAP/hybrid_md_package/`` directory of this repo.

This is located at ``QUIP/src/GAP/hybrid_md_package/`` if you are looking at the QUIP source downloaded above.

.. code-block:: bash

    python -m pip install .

Archer2 cluster
###############

These are specific and tested instructions the UK's Archer2 https://www.archer2.ac.uk computer cluster.

An important gotcha on Archer2 is that the built-in maths libraries of the compiler are linking MPI by default, which 
breaks the setup for the ``gap_fit`` program, so we need to build that explicitly without MPI, see below.

Edit the CASTEP Makefile to include the following

.. code-block:: makefile
    
    COMMS_ARCH := mpi
    FFT := fftw3
    BUILD := fast
    MATHLIBS := mkl # optional

Full installation:

.. code-block:: bash

    # create bin directory for executables
    mkdir bin

    # Clone QUIP
    git clone --recursive https://github.com/libAtoms/QUIP.git --depth 1 --single-branch

    # load the correct modules
    module load cray-python
    module switch PrgEnv-cray PrgEnv-gnu/8.1.0
    module load cpe/22.04
    module load cray-fftw
    module load mkl/2023.0.0  # if using MKL for Castep

    # step 1: Python interpreter & installation of hybrid-md
    python -m virtualenv venv
    source venv/bin/activate
    
    python -m pip install ./QUIP/src/GAP/hybrid_md_package/

    # step 2: QUIP with MPI
    cd QUIP
    export QUIP_ROOT=$(pwd)
    export QUIP_ARCH=archer2_mpich
    make config   # configure: build - Make sure to say yes to GAP
    make libquip.a
    cd ../ # back to the starting dir

    # step 3. CASTEP with linking QUIP
    cd CASTEP/
    make -j8
    cp obj/linux_x86_64_gfortran10-XT--mpi/castep.mpi ../bin/
    cd ../ # back to the starting dir

    # step 4. Install gap_fit

    # IMPORTANT!! unload comms modules -> no MPI
    module load craype-network-none
    module remove cray-mpich
    
    cd QUIP
    export QUIP_ARCH=archer2_openmp
    make gap_programs
    cp build/archer2_openmp/gap_fit ../bin/
    cd ../ # back to the starting dir

