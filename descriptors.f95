! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   GAP (Gaussian Approximation Potental)
! HND X
! HND X
! HND X   Portions of GAP were written by Albert Bartok-Partay, Gabor Csanyi,
! HND X   Copyright 2006-2021.
! HND X
! HND X   Portions of GAP were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   GAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   GAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Gabor Csanyi or Albert Bartok-Partay. The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A. P. Bartok et al Physical Review Letters vol 104 p136403 (2010)
! HND X
! HND X   When using the SOAP kernel or its variants, please additionally cite:
! HND X
! HND X   A. P. Bartok et al Physical Review B vol 87 p184115 (2013)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module descriptors_module

   use error_module
   use system_module, only : dp, print, optional_default, system_timer, operator(//), split_string, string_to_int, split_string_simple, inoutput, OUTPUT, PRINT_VERBOSE, PRINT_NERD
   use linkedlist_module
   use units_module
   use periodictable_module
   use linearalgebra_module
   use dictionary_module
   use paramreader_module
   use atoms_module
   use atoms_types_module
   use topology_module
   use mpi_context_module
   use table_module
#ifdef DESCRIPTORS_NONCOMMERCIAL
   use permutation_maker_module
#endif
   use CInOutput_module
   use clusters_module
   use connection_module
   use angular_functions_module

   implicit none

   private
#ifdef GAP_VERSION
   integer, parameter :: gap_version = GAP_VERSION
#else
   integer, parameter :: gap_version = 0
#endif


   integer, parameter, public :: DT_NONE            =  0
   integer, parameter, public :: DT_BISPECTRUM_SO4  =  1
   integer, parameter, public :: DT_BISPECTRUM_SO3  =  2
   integer, parameter, public :: DT_BEHLER          =  3
   integer, parameter, public :: DT_DISTANCE_2B     =  4
   integer, parameter, public :: DT_COORDINATION    =  5
   integer, parameter, public :: DT_ANGLE_3B        =  6
   integer, parameter, public :: DT_CO_ANGLE_3B     =  7
   integer, parameter, public :: DT_CO_DISTANCE_2B  =  8
   integer, parameter, public :: DT_COSNX           =  9
   integer, parameter, public :: DT_TRIHIS          = 10
   integer, parameter, public :: DT_WATER_MONOMER   = 11
   integer, parameter, public :: DT_WATER_DIMER     = 12
   integer, parameter, public :: DT_A2_DIMER        = 13
   integer, parameter, public :: DT_AB_DIMER        = 14
   integer, parameter, public :: DT_BOND_REAL_SPACE = 15
   integer, parameter, public :: DT_ATOM_REAL_SPACE = 16
   integer, parameter, public :: DT_POWER_SO3       = 17
   integer, parameter, public :: DT_POWER_SO4       = 18
   integer, parameter, public :: DT_SOAP            = 19
   integer, parameter, public :: DT_AN_MONOMER      = 20
   integer, parameter, public :: DT_GENERAL_MONOMER = 21
   integer, parameter, public :: DT_GENERAL_DIMER   = 22
   integer, parameter, public :: DT_GENERAL_TRIMER  = 23
   integer, parameter, public :: DT_RDF             = 24
   integer, parameter, public :: DT_AS_DISTANCE_2B  = 25
   integer, parameter, public :: DT_MOLECULE_LO_D   = 26
   integer, parameter, public :: DT_alex            = 27
   integer, parameter, public :: DT_COM_DIMER       = 28
   integer, parameter, public :: DT_DISTANCE_NB     = 29
   integer, parameter, public :: DT_SOAP_EXPRESS    = 30
   integer, parameter, public :: DT_SOAP_TURBO      = 31

   integer, parameter :: NP_WATER_DIMER    = 8
   integer, parameter :: NP_A2_DIMER       = 8
   integer, parameter :: NP_AB_DIMER       = 2

   type transfer_parameters_type
      logical :: do_transfer
      real(dp) :: factor, r0, width
   endtype transfer_parameters_type

   type descriptor_data_mono
      real(dp), dimension(:), allocatable :: data
      real(dp), dimension(:,:,:), allocatable :: grad_data
      ! ci : atom indices amongst which to distribute energy of descriptor
      ! ii : all atoms involved in descriptor (for partial derivatives)
      integer, dimension(:), allocatable :: ci, ii
      real(dp), dimension(:,:), allocatable :: pos
      logical :: has_data
      logical, dimension(:), allocatable :: has_grad_data

      real(dp) :: covariance_cutoff = 1.0_dp
      real(dp), dimension(:,:), allocatable :: grad_covariance_cutoff
   endtype descriptor_data_mono

   type cplx_2d
      complex(dp), dimension(:,:), allocatable :: mm
   endtype cplx_2d

   type int_2d
      integer , dimension(:,:), allocatable :: mm
   endtype int_2d

   type real_2d
      real(dp), dimension(:,:), allocatable :: mm
   endtype real_2d

   type cplx_3d
      complex(dp), dimension(:,:,:), allocatable :: mm
   endtype cplx_3d

   !=======================================================================
   !==                begin descriptors
   !=======================================================================


   type RadialFunction_type
      integer :: n_max
      real(dp) :: cutoff, min_cutoff
      real(dp), dimension(:,:), allocatable :: RadialTransform
      real(dp), dimension(:), allocatable :: NormFunction

      logical :: initialised = .false.
   endtype RadialFunction_type

   type fourier_SO4_type
      real(dp) :: cutoff
      real(dp) :: z0_ratio
      real(dp) :: z0
      integer :: j_max, Z
      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype fourier_SO4_type

   type bispectrum_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w

      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype bispectrum_SO4

   type bispectrum_SO3

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype bispectrum_SO3

   type behler_g2
      integer :: Z_n = 0
      real(dp) :: eta
      real(dp) :: rs
      real(dp) :: rc
   endtype behler_g2

   type behler_g3
      integer,dimension(2) :: Z_n = 0
      real(dp) :: eta
      real(dp) :: lambda
      real(dp) :: zeta
      real(dp) :: rc
   endtype behler_g3

   type behler

      real(dp) :: cutoff = 0.0_dp
      logical :: initialised = .false.

      integer :: Z = 0
      integer :: n_g2, n_g3
      type(behler_g2), dimension(:), allocatable :: g2
      type(behler_g3), dimension(:), allocatable :: g3

   endtype behler

   type distance_2b
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: Z1, Z2
      character(STRING_LENGTH) :: resid_name
      logical :: only_intra, only_inter

      integer :: n_exponents, tail_exponent
      real(dp) :: tail_range
      integer, dimension(:), allocatable :: exponents

      logical :: has_tail
      logical :: initialised = .false.

   endtype distance_2b

   type coordination
      real(dp) :: cutoff
      real(dp) :: transition_width
      integer :: Z

      logical :: initialised = .false.

   endtype coordination

   type angle_3b
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype angle_3b

   type co_angle_3b
      real(dp) :: cutoff
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype co_angle_3b

   type co_distance_2b
      real(dp) :: cutoff
      real(dp) :: transition_width
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype co_distance_2b

   type cosnx

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype cosnx

   type trihis
      real(dp) :: cutoff
      integer :: n_gauss

      real(dp), dimension(:,:), allocatable :: gauss_centre
      real(dp), dimension(:,:), allocatable :: gauss_width

      logical :: initialised = .false.

   endtype trihis

   type water_monomer
      real(dp) :: cutoff

      logical :: initialised = .false.

   endtype water_monomer

   type water_dimer
      real(dp) :: cutoff, cutoff_transition_width
      real(dp) :: monomer_cutoff
      logical :: OHH_ordercheck
      real(dp) :: power,dist_shift

      logical :: initialised = .false.

   endtype water_dimer

   type A2_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number

      logical :: initialised = .false.

   endtype A2_dimer

   type AB_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number1, atomic_number2

      logical :: initialised = .false.

   endtype AB_dimer

   type atom_real_space
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: l_max
      real(dp) :: alpha
      real(dp) :: zeta

      logical :: initialised = .false.

   endtype atom_real_space

   type power_so3
      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype power_so3

   type power_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w

      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype power_SO4

   type soap
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      real(dp) :: alpha, atom_sigma, covariance_sigma0, central_weight

      integer :: cutoff_dexp
      real(dp) :: cutoff_scale
      real(dp) :: cutoff_rate
      integer :: l_max, n_max, n_Z, n_species

      integer ::  nu_R, nu_S
      integer, dimension(:), allocatable :: species_Z, Z
      real(dp), dimension(:), allocatable :: r_basis
      real(dp), dimension(:,:), allocatable :: transform_basis,cholesky_overlap_basis

      logical :: global = .false.
      logical :: central_reference_all_species = .false.
      logical :: diagonal_radial = .false.
      logical :: normalise = .true.
      logical :: initialised = .false.

      logical :: Z_mix = .false.
      logical :: R_mix = .false.
      logical :: sym_mix = .false.
      logical :: coupling = .false.
      integer :: K
   endtype soap


   type rdf
      real(dp) :: cutoff
      real(dp) :: transition_width, w_gauss
      integer :: Z, n_gauss
      real(dp), dimension(:), allocatable :: r_gauss

      logical :: initialised = .false.

   endtype rdf

   type as_distance_2b
      real(dp) :: min_cutoff, max_cutoff, as_cutoff, overlap_alpha
      real(dp) :: min_transition_width, max_transition_width, as_transition_width
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype as_distance_2b

   type alex

      integer :: Z, power_min, power_max
      real(dp) :: cutoff

      integer :: n_species
      integer, dimension(:), allocatable :: species_Z

      logical :: initialised = .false.
   endtype alex


   type distance_Nb
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: order
      integer, dimension(:), allocatable :: Z
      integer :: n_permutations
      integer, dimension(:,:), allocatable :: permutations
      logical, dimension(:,:,:), allocatable :: monomerConnectivities
      logical :: compact_clusters = .false.
      logical :: initialised = .false.
   endtype distance_Nb

   type soap_turbo
      ! User controllable parameters
      real(dp) :: rcut_hard, rcut_soft, nf
      integer  :: n_species, radial_enhancement, central_index, l_max, compress_P_nonzero
      character(len=STRING_LENGTH) :: basis, scaling_mode, compress_file, compress_mode

      real(dp), dimension(:), allocatable :: atom_sigma_r, atom_sigma_r_scaling, &
         atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, central_weight, compress_P_el
      integer, dimension(:), allocatable :: species_Z, alpha_max, compress_P_i, compress_P_j

      logical :: initialised = .false., compress = .false.
   endtype soap_turbo

#ifdef DESCRIPTORS_NONCOMMERCIAL
#include "descriptors_noncommercial_types.inc"
#endif

   !
   ! All the descriptors need to be public so that they are visible to the python wrapper
#ifdef DESCRIPTORS_NONCOMMERCIAL
   public :: soap, general_monomer, bispectrum_so4, bispectrum_so3, behler, distance_2b, &
        coordination, angle_3b, co_angle_3b, co_distance_2b, cosnx, trihis, water_monomer, &
        water_dimer, a2_dimer, bond_real_space, power_so3, power_so4, an_monomer, general_dimer, &
        general_trimer, rdf, as_distance_2b, molecule_lo_d, alex,  com_dimer,  distance_nb, &
        descriptor_data_mono, fourier_so4_type, radialfunction_type, transfer_parameters_type, &
        ab_dimer, atom_real_space, spherical_harmonics_type, behler_g2, behler_g3, soap_turbo, soap_express
#else
   public :: soap, bispectrum_so4, bispectrum_so3, behler, distance_2b, &
        coordination, angle_3b, co_angle_3b, co_distance_2b, cosnx, trihis, water_monomer, &
        water_dimer, a2_dimer, power_so3, power_so4, &
        rdf, as_distance_2b,  alex,  distance_nb, &
        descriptor_data_mono, fourier_so4_type, radialfunction_type, transfer_parameters_type, &
        ab_dimer, atom_real_space, spherical_harmonics_type, behler_g2, behler_g3, &
        soap_turbo
#endif

   !=======================================================================
   !==                end descriptors
   !=======================================================================

   type descriptor
      integer :: descriptor_type = DT_NONE

      type(bispectrum_SO4)  :: descriptor_bispectrum_SO4
      type(bispectrum_SO3)  :: descriptor_bispectrum_SO3
      type(behler)          :: descriptor_behler
      type(distance_2b)     :: descriptor_distance_2b
      type(coordination)    :: descriptor_coordination
      type(angle_3b)        :: descriptor_angle_3b
      type(co_angle_3b)     :: descriptor_co_angle_3b
      type(co_distance_2b)  :: descriptor_co_distance_2b
      type(cosnx)           :: descriptor_cosnx
      type(trihis)          :: descriptor_trihis
      type(water_monomer)   :: descriptor_water_monomer
      type(water_dimer)     :: descriptor_water_dimer
      type(A2_dimer)        :: descriptor_A2_dimer
      type(AB_dimer)        :: descriptor_AB_dimer
      type(atom_real_space) :: descriptor_atom_real_space
      type(power_so3)       :: descriptor_power_so3
      type(power_SO4)       :: descriptor_power_SO4
      type(soap)            :: descriptor_soap
      type(rdf)             :: descriptor_rdf
      type(as_distance_2b)  :: descriptor_as_distance_2b
      type(alex)            :: descriptor_alex
      type(distance_Nb)     :: descriptor_distance_Nb
      type(soap_turbo)      :: descriptor_soap_turbo
#ifdef DESCRIPTORS_NONCOMMERCIAL
      type(AN_monomer)      :: descriptor_AN_monomer
      type(general_monomer) :: descriptor_general_monomer
      type(general_dimer)   :: descriptor_general_dimer
      type(general_trimer)  :: descriptor_general_trimer
      type(molecule_lo_d)   :: descriptor_molecule_lo_d
      type(com_dimer)       :: descriptor_com_dimer
      type(soap_express)    :: descriptor_soap_express
      type(bond_real_space) :: descriptor_bond_real_space
#endif
   endtype

   type descriptor_data
      type(descriptor_data_mono), dimension(:), allocatable :: x
   endtype descriptor_data

   type cplx_1d
      complex(dp), dimension(:), allocatable :: m
   endtype cplx_1d

   type real_1d
      real(dp), dimension(:), allocatable :: m
   endtype real_1d

   type spherical_harmonics_type
      type(cplx_1d), dimension(:), allocatable :: spherical_harmonics
      type(cplx_2d), dimension(:), allocatable :: grad_spherical_harmonics
      real(dp) :: r
      real(dp), dimension(3) :: u
   endtype spherical_harmonics_type

   type neighbour_type
      type(spherical_harmonics_type), dimension(:), allocatable :: neighbour
   endtype neighbour_type

   type grad_spherical_harmonics_overlap_type
      type(cplx_3d), dimension(:), allocatable :: grad_integral
   endtype grad_spherical_harmonics_overlap_type

   public :: neighbour_type, real_space_fourier_coefficients, real_space_covariance_coefficient
   public :: SphericalYCartesian

   interface initialise
#ifdef DESCRIPTORS_NONCOMMERCIAL
      module procedure descriptor_initialise, RadialFunction_initialise, fourier_so4_initialise, &
      bispectrum_SO4_initialise, bispectrum_SO3_initialise, behler_initialise, distance_2b_initialise, &
      coordination_initialise, angle_3b_initialise, co_angle_3b_initialise, co_distance_2b_initialise, cosnx_initialise, trihis_initialise, &
      water_monomer_initialise, water_dimer_initialise, A2_dimer_initialise, AB_dimer_initialise, distance_Nb_initialise,  rdf_initialise, as_distance_2b_initialise, alex_initialise, &
      atom_real_space_initialise, power_so3_initialise, power_SO4_initialise, soap_initialise, soap_turbo_initialise, &
      general_monomer_initialise, general_dimer_initialise, general_trimer_initialise,  molecule_lo_d_initialise,  AN_monomer_initialise, &
      bond_real_space_initialise, transfer_initialise, com_dimer_initialise,  soap_express_initialise
#else
      module procedure descriptor_initialise, RadialFunction_initialise, fourier_so4_initialise, &
      bispectrum_SO4_initialise, bispectrum_SO3_initialise, behler_initialise, distance_2b_initialise, &
      coordination_initialise, angle_3b_initialise, co_angle_3b_initialise, co_distance_2b_initialise, cosnx_initialise, trihis_initialise, &
      water_monomer_initialise, water_dimer_initialise, A2_dimer_initialise, AB_dimer_initialise, distance_Nb_initialise,  rdf_initialise, as_distance_2b_initialise, alex_initialise, &
      atom_real_space_initialise, power_so3_initialise, power_SO4_initialise, soap_initialise, soap_turbo_initialise
#endif
   endinterface initialise
   public :: initialise

   interface finalise
#ifdef DESCRIPTORS_NONCOMMERCIAL
      module procedure descriptor_finalise, descriptor_data_finalise, RadialFunction_finalise, fourier_so4_finalise, cplx_2d_array1_finalise, cplx_3d_array2_finalise, &
      bispectrum_SO4_finalise, bispectrum_SO3_finalise, behler_finalise, distance_2b_finalise, coordination_finalise, angle_3b_finalise, co_angle_3b_finalise, &
      co_distance_2b_finalise, cosnx_finalise, trihis_finalise, water_monomer_finalise, water_dimer_finalise, rdf_finalise, as_distance_2b_finalise,  alex_finalise, &
      A2_dimer_finalise, AB_dimer_finalise, atom_real_space_finalise, power_so3_finalise, power_SO4_finalise, soap_finalise, distance_Nb_finalise,  soap_turbo_finalise, &
      AN_monomer_finalise, general_monomer_finalise, general_dimer_finalise, general_trimer_finalise,  molecule_lo_d_finalise, com_dimer_finalise, &
      bond_real_space_finalise, soap_express_finalise
#else
      module procedure descriptor_finalise, descriptor_data_finalise, RadialFunction_finalise, fourier_so4_finalise, cplx_2d_array1_finalise, cplx_3d_array2_finalise, &
      bispectrum_SO4_finalise, bispectrum_SO3_finalise, behler_finalise, distance_2b_finalise, coordination_finalise, angle_3b_finalise, co_angle_3b_finalise, &
      co_distance_2b_finalise, cosnx_finalise, trihis_finalise, water_monomer_finalise, water_dimer_finalise, rdf_finalise, as_distance_2b_finalise,  alex_finalise, &
      A2_dimer_finalise, AB_dimer_finalise, atom_real_space_finalise, power_so3_finalise, power_SO4_finalise, soap_finalise, distance_Nb_finalise, soap_turbo_finalise
#endif
   endinterface finalise
   public :: finalise

   interface calc
#ifdef DESCRIPTORS_NONCOMMERCIAL
      module procedure descriptor_calc, descriptor_calc_array, bispectrum_SO4_calc, bispectrum_SO3_calc, behler_calc, distance_2b_calc, coordination_calc, angle_3b_calc, co_angle_3b_calc, &
      co_distance_2b_calc, cosnx_calc, trihis_calc, water_monomer_calc, water_dimer_calc, A2_dimer_calc, AB_dimer_calc,  atom_real_space_calc, &
      power_so3_calc, power_SO4_calc, soap_calc, rdf_calc, as_distance_2b_calc, &
      distance_Nb_calc, alex_calc, soap_turbo_calc, &
      AN_monomer_calc,  soap_express_calc, general_monomer_calc, general_dimer_calc, general_trimer_calc,  molecule_lo_d_calc, com_dimer_calc, bond_real_space_calc
#else
      module procedure descriptor_calc, descriptor_calc_array, bispectrum_SO4_calc, bispectrum_SO3_calc, behler_calc, distance_2b_calc, coordination_calc, angle_3b_calc, co_angle_3b_calc, &
      co_distance_2b_calc, cosnx_calc, trihis_calc, water_monomer_calc, water_dimer_calc, A2_dimer_calc, AB_dimer_calc,  atom_real_space_calc, &
      power_so3_calc, power_SO4_calc, soap_calc, rdf_calc, as_distance_2b_calc, &
      distance_Nb_calc, alex_calc, soap_turbo_calc

#endif
   endinterface calc
   public :: calc

   interface cutoff
#ifdef DESCRIPTORS_NONCOMMERCIAL
      module procedure descriptor_cutoff, bispectrum_SO4_cutoff, bispectrum_SO3_cutoff, behler_cutoff, distance_2b_cutoff, coordination_cutoff, angle_3b_cutoff, co_angle_3b_cutoff, &
      co_distance_2b_cutoff, cosnx_cutoff, trihis_cutoff, water_monomer_cutoff, water_dimer_cutoff, A2_dimer_cutoff, AB_dimer_cutoff, atom_real_space_cutoff, &
      power_so3_cutoff, power_SO4_cutoff, soap_cutoff, alex_cutoff, distance_Nb_cutoff, rdf_cutoff, as_distance_2b_cutoff, soap_turbo_cutoff, &
      molecule_lo_d_cutoff, com_dimer_cutoff, soap_express_cutoff, AN_monomer_cutoff, general_monomer_cutoff, general_dimer_cutoff, general_trimer_cutoff,  bond_real_space_cutoff
#else
      module procedure descriptor_cutoff, bispectrum_SO4_cutoff, bispectrum_SO3_cutoff, behler_cutoff, distance_2b_cutoff, coordination_cutoff, angle_3b_cutoff, co_angle_3b_cutoff, &
      co_distance_2b_cutoff, cosnx_cutoff, trihis_cutoff, water_monomer_cutoff, water_dimer_cutoff, A2_dimer_cutoff, AB_dimer_cutoff, atom_real_space_cutoff, &
      power_so3_cutoff, power_SO4_cutoff, soap_cutoff, alex_cutoff, distance_Nb_cutoff, rdf_cutoff, as_distance_2b_cutoff, soap_turbo_cutoff
#endif
   endinterface cutoff
   public :: cutoff

   interface descriptor_sizes
#ifdef DESCRIPTORS_NONCOMMERCIAL
      module procedure descriptor_sizes, bispectrum_SO4_sizes, bispectrum_SO3_sizes, behler_sizes, distance_2b_sizes, coordination_sizes, angle_3b_sizes, co_angle_3b_sizes, &
      co_distance_2b_sizes, cosnx_sizes, trihis_sizes, water_monomer_sizes, water_dimer_sizes, A2_dimer_sizes, AB_dimer_sizes, atom_real_space_sizes, &
      power_so3_sizes, power_SO4_sizes, soap_sizes,  rdf_sizes, as_distance_2b_sizes, &
      alex_sizes, distance_Nb_sizes, soap_turbo_sizes, &
      molecule_lo_d_sizes, com_dimer_sizes,  soap_express_sizes, AN_monomer_sizes, general_monomer_sizes, general_dimer_sizes, general_trimer_sizes,  bond_real_space_sizes
#else
      module procedure descriptor_sizes, bispectrum_SO4_sizes, bispectrum_SO3_sizes, behler_sizes, distance_2b_sizes, coordination_sizes, angle_3b_sizes, co_angle_3b_sizes, &
      co_distance_2b_sizes, cosnx_sizes, trihis_sizes, water_monomer_sizes, water_dimer_sizes, A2_dimer_sizes, AB_dimer_sizes, atom_real_space_sizes, &
      power_so3_sizes, power_SO4_sizes, soap_sizes,  rdf_sizes, as_distance_2b_sizes, &
      alex_sizes, distance_Nb_sizes, soap_turbo_sizes
#endif
   endinterface descriptor_sizes
   public :: descriptor_sizes

   public :: descriptor_MPI_setup

   public :: descriptor, descriptor_data, descriptor_dimensions, descriptor_n_permutations, descriptor_permutations, descriptor_str_add_species
   public :: real_space_covariance
   public :: cplx_1d, cplx_2d

   contains


#ifdef DESCRIPTORS_NONCOMMERCIAL
#include "descriptors_noncommercial.inc"
#endif

   function get_descriptor_type(args_str,error)
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      integer :: get_descriptor_type

      type(Dictionary) :: params
      logical :: is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, &
         is_co_angle_3b, is_co_distance_2b, is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, &
         is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, is_soap, &
         is_AN_monomer, is_general_monomer, is_general_dimer, is_general_trimer, is_rdf, is_as_distance_2b, &
         is_molecule_lo_d, is_alex, is_com_dimer, is_distance_Nb, is_soap_express, is_soap_turbo

      INIT_ERROR(error)

      call initialise(params)
      call param_register(params, 'bispectrum_so4', 'false', is_bispectrum_so4, help_string="Type of descriptor is bispectrum_so4.")
      call param_register(params, 'bispectrum_so3', 'false', is_bispectrum_so3, help_string="Type of descriptor is bispectrum_so3.")
      call param_register(params, 'behler', 'false', is_behler, help_string="Type of descriptor is behler.")
      call param_register(params, 'distance_2b', 'false', is_distance_2b, help_string="Type of descriptor is distance_2b.")
      call param_register(params, 'coordination', 'false', is_coordination, help_string="Type of descriptor is coordination.")
      call param_register(params, 'angle_3b', 'false', is_angle_3b, help_string="Type of descriptor is angle_3b.")
      call param_register(params, 'co_angle_3b', 'false', is_co_angle_3b, help_string="Type of descriptor is co_angle_3b.")
      call param_register(params, 'co_distance_2b', 'false', is_co_distance_2b, help_string="Type of descriptor is co_distance_2b.")
      call param_register(params, 'cosnx', 'false', is_cosnx, help_string="Type of descriptor is cosnx.")
      call param_register(params, 'trihis', 'false', is_trihis, help_string="Type of descriptor is trihis.")
      call param_register(params, 'water_monomer', 'false', is_water_monomer, help_string="Type of descriptor is water_monomer.")
      call param_register(params, 'water_dimer', 'false', is_water_dimer, help_string="Type of descriptor is water_dimer.")
      call param_register(params, 'A2_dimer', 'false', is_A2_dimer, help_string="Type of descriptor is A2_dimer.")
      call param_register(params, 'AB_dimer', 'false', is_AB_dimer, help_string="Type of descriptor is AB_dimer.")
      call param_register(params, 'bond_real_space', 'false', is_bond_real_space, help_string="Type of descriptor is bond_real_space.")
      call param_register(params, 'atom_real_space', 'false', is_atom_real_space, help_string="Type of descriptor is atom_real_space.")
      call param_register(params, 'power_so3', 'false', is_power_so3, help_string="Type of descriptor is power_so3.")
      call param_register(params, 'power_so4', 'false', is_power_so4, help_string="Type of descriptor is power_so4.")
      call param_register(params, 'soap', 'false', is_soap, help_string="Type of descriptor is soap.")
      call param_register(params, 'AN_monomer', 'false', is_AN_monomer, help_string="Type of descriptor is AN_monomer.")
      call param_register(params, 'general_monomer', 'false', is_general_monomer, help_string="Type of descriptor is general_monomer.")
      call param_register(params, 'general_dimer', 'false', is_general_dimer, help_string="Type of descriptor is general_dimer.")
      call param_register(params, 'general_trimer', 'false', is_general_trimer, help_string="Type of descriptor is general_trimer.")
      call param_register(params, 'rdf', 'false', is_rdf, help_string="Type of descriptor is rdf.")
      call param_register(params, 'as_distance_2b', 'false', is_as_distance_2b, help_string="Type of descriptor is as_distance_2b.")
      call param_register(params, 'molecule_lo_d', 'false', is_molecule_lo_d, help_string="Type of descriptor is molecule_lo_d.")
      call param_register(params, 'alex', 'false', is_alex, help_string="Type of descriptor is alex.")
      call param_register(params, 'com_dimer', 'false', is_com_dimer, help_string="Type of descriptor is com_dimer.")
      call param_register(params, 'distance_Nb', 'false', is_distance_Nb, help_string="Type of descriptor is distance_Nb.")
      call param_register(params, 'soap_express', 'false', is_soap_express, help_string="Type of descriptor is soap_express.")
      call param_register(params, 'soap_turbo', 'false', is_soap_turbo, help_string="Type of descriptor is soap_turbo.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='descriptor_initialise args_str')) then
         RAISE_ERROR("descriptor_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (count( (/is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, is_co_angle_3b, is_co_distance_2b, &
      is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, &
      is_soap, is_AN_monomer, is_general_monomer, is_general_dimer, is_general_trimer, is_rdf, is_as_distance_2b, is_molecule_lo_d, is_alex, is_com_dimer, &
      is_distance_Nb, is_soap_express, is_soap_turbo /) ) /= 1) then
         RAISE_ERROR("descriptor_initialise found too few or too many IP Model types args_str='"//trim(args_str)//"'", error)
      endif

      get_descriptor_type = DT_NONE

      if( is_bispectrum_so4 ) then
         get_descriptor_type = DT_BISPECTRUM_SO4
      elseif( is_bispectrum_so3 ) then
         get_descriptor_type = DT_BISPECTRUM_SO3
      elseif( is_behler ) then
         get_descriptor_type = DT_BEHLER
      elseif( is_distance_2b ) then
         get_descriptor_type = DT_DISTANCE_2B
      elseif( is_coordination ) then
         get_descriptor_type = DT_COORDINATION
      elseif( is_angle_3b ) then
         get_descriptor_type = DT_ANGLE_3B
      elseif( is_co_angle_3b ) then
         get_descriptor_type = DT_CO_ANGLE_3B
      elseif( is_co_distance_2b ) then
         get_descriptor_type = DT_CO_DISTANCE_2B
      elseif( is_cosnx ) then
         get_descriptor_type = DT_COSNX
      elseif( is_trihis ) then
         get_descriptor_type = DT_TRIHIS
      elseif( is_water_monomer ) then
         get_descriptor_type = DT_WATER_MONOMER
      elseif( is_water_dimer ) then
         get_descriptor_type = DT_WATER_DIMER
      elseif( is_A2_dimer ) then
         get_descriptor_type = DT_A2_DIMER
      elseif( is_AB_dimer ) then
         get_descriptor_type = DT_AB_DIMER
      elseif( is_bond_real_space ) then
         get_descriptor_type = DT_BOND_REAL_SPACE
      elseif( is_atom_real_space ) then
         get_descriptor_type = DT_ATOM_REAL_SPACE
      elseif( is_power_so3 ) then
         get_descriptor_type = DT_POWER_SO3
      elseif( is_power_so4 ) then
         get_descriptor_type = DT_POWER_SO4
      elseif( is_soap ) then
         get_descriptor_type = DT_SOAP
      elseif( is_AN_monomer ) then
         get_descriptor_type = DT_AN_MONOMER
      elseif( is_general_monomer ) then
         get_descriptor_type = DT_GENERAL_MONOMER
      elseif( is_general_dimer ) then
         get_descriptor_type = DT_GENERAL_DIMER
      elseif( is_general_trimer ) then
         get_descriptor_type = DT_GENERAL_TRIMER
      elseif( is_rdf ) then
         get_descriptor_type = DT_RDF
      elseif( is_as_distance_2b ) then
         get_descriptor_type = DT_AS_DISTANCE_2B
      elseif( is_molecule_lo_d ) then
         get_descriptor_type = DT_MOLECULE_LO_D
      elseif( is_alex ) then
         get_descriptor_type = DT_ALEX
      elseif( is_com_dimer ) then
         get_descriptor_type = DT_COM_DIMER
      elseif( is_distance_Nb ) then
         get_descriptor_type = DT_DISTANCE_NB
      elseif( is_soap_express ) then
         get_descriptor_type = DT_SOAP_EXPRESS
      elseif( is_soap_turbo ) then
         get_descriptor_type = DT_SOAP_TURBO
      endif

   endfunction get_descriptor_type

   subroutine descriptor_initialise(this,args_str,error)
      type(descriptor), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      this%descriptor_type = get_descriptor_type(args_str,error)

      select case(this%descriptor_type)
      case(DT_BISPECTRUM_SO4)
         call initialise(this%descriptor_bispectrum_SO4,args_str,error)
      case(DT_BISPECTRUM_SO3)
         call initialise(this%descriptor_bispectrum_SO3,args_str,error)
      case(DT_BEHLER)
         call initialise(this%descriptor_behler,args_str,error)
      case(DT_DISTANCE_2B)
         call initialise(this%descriptor_distance_2b,args_str,error)
      case(DT_COORDINATION)
         call initialise(this%descriptor_coordination,args_str,error)
      case(DT_ANGLE_3B)
         call initialise(this%descriptor_angle_3b,args_str,error)
      case(DT_CO_ANGLE_3B)
         call initialise(this%descriptor_co_angle_3b,args_str,error)
      case(DT_CO_DISTANCE_2B)
         call initialise(this%descriptor_co_distance_2b,args_str,error)
      case(DT_COSNX)
         call initialise(this%descriptor_cosnx,args_str,error)
      case(DT_TRIHIS)
         call initialise(this%descriptor_trihis,args_str,error)
      case(DT_WATER_MONOMER)
         call initialise(this%descriptor_water_monomer,args_str,error)
      case(DT_WATER_DIMER)
         call initialise(this%descriptor_water_dimer,args_str,error)
      case(DT_A2_DIMER)
         call initialise(this%descriptor_A2_dimer,args_str,error)
      case(DT_AB_DIMER)
         call initialise(this%descriptor_AB_dimer,args_str,error)
      case(DT_ATOM_REAL_SPACE)
         call initialise(this%descriptor_atom_real_space,args_str,error)
      case(DT_POWER_SO3)
         call initialise(this%descriptor_power_so3,args_str,error)
      case(DT_POWER_SO4)
         call initialise(this%descriptor_power_so4,args_str,error)
      case(DT_SOAP)
         call initialise(this%descriptor_soap,args_str,error)
      case(DT_RDF)
         call initialise(this%descriptor_rdf,args_str,error)
      case(DT_AS_DISTANCE_2B)
         call initialise(this%descriptor_as_distance_2b,args_str,error)
      case(DT_ALEX)
         call initialise(this%descriptor_alex,args_str,error)
      case(DT_DISTANCE_NB)
         call initialise(this%descriptor_distance_Nb,args_str,error)
      case(DT_SOAP_TURBO)
         call initialise(this%descriptor_soap_turbo,args_str,error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
      case(DT_BOND_REAL_SPACE)
         call initialise(this%descriptor_bond_real_space,args_str,error)
      case(DT_AN_MONOMER)
         call initialise(this%descriptor_AN_monomer,args_str,error)
      case(DT_COM_DIMER)
         call initialise(this%descriptor_com_dimer,args_str,error)
      case(DT_MOLECULE_LO_D)
         call initialise(this%descriptor_molecule_lo_d,args_str,error)
      case(DT_GENERAL_MONOMER)
         call initialise(this%descriptor_general_monomer,args_str,error)
      case(DT_GENERAL_DIMER)
         call initialise(this%descriptor_general_dimer,args_str,error)
      case(DT_GENERAL_TRIMER)
         call initialise(this%descriptor_general_trimer,args_str,error)
      case(DT_SOAP_EXPRESS)
         call initialise(this%descriptor_soap_express,args_str,error)
#endif
      endselect

   endsubroutine descriptor_initialise

   subroutine descriptor_finalise(this,error)
      type(descriptor), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call finalise(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            call finalise(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            call finalise(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            call finalise(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            call finalise(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            call finalise(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            call finalise(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            call finalise(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            call finalise(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            call finalise(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            call finalise(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            call finalise(this%descriptor_water_dimer,error)
         case(DT_A2_dimer)
            call finalise(this%descriptor_A2_dimer,error)
         case(DT_AB_dimer)
            call finalise(this%descriptor_AB_dimer,error)
         case(DT_ATOM_REAL_SPACE)
            call finalise(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            call finalise(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            call finalise(this%descriptor_power_so4,error)
         case(DT_SOAP)
            call finalise(this%descriptor_soap,error)
         case(DT_RDF)
            call finalise(this%descriptor_rdf,error)
         case(DT_AS_DISTANCE_2b)
            call finalise(this%descriptor_as_distance_2b,error)
         case(DT_ALEX)
            call finalise(this%descriptor_alex,error)
         case(DT_DISTANCE_Nb)
            call finalise(this%descriptor_distance_Nb,error)
#ifdef DESCRIPTOR_NONCOMMERCIAL
         case(DT_COM_DIMER)
            call finalise(this%descriptor_com_dimer,error)
         case(DT_MOLECULE_LO_D)
            call finalise(this%descriptor_molecule_lo_d,error)
         case(DT_BOND_REAL_SPACE)
            call finalise(this%descriptor_bond_real_space,error)
         case(DT_GENERAL_MONOMER)
            call finalise(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            call finalise(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            call finalise(this%descriptor_general_trimer,error)
         case(DT_SOAP_EXPRESS)
            call finalise(this%descriptor_soap_express,error)
         case(DT_SOAP_TURBO)
            call finalise(this%descriptor_soap_turbo,error)
#endif
      endselect

      this%descriptor_type = DT_NONE

   endsubroutine descriptor_finalise

   subroutine descriptor_MPI_setup(this,at,mpi,mpi_mask,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(mpi%active) then
         select case(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_BISPECTRUM_SO3)
            RAISE_ERROR("descriptor_MPI_setup: bispectrum_so3 not MPI ready.", error)
         case(DT_BEHLER)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_DISTANCE_2B)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_COORDINATION)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_ANGLE_3B)
            RAISE_ERROR("descriptor_MPI_setup: angle_3b not MPI ready.", error)
         case(DT_CO_ANGLE_3B)
            RAISE_ERROR("descriptor_MPI_setup: co_angle_3b not MPI ready.", error)
         case(DT_CO_DISTANCE_2B)
            RAISE_ERROR("descriptor_MPI_setup: co_distance_2b not MPI ready.", error)
         case(DT_COSNX)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_TRIHIS)
            RAISE_ERROR("descriptor_MPI_setup: trihis not MPI ready.", error)
         case(DT_WATER_MONOMER)
            call descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_WATER_DIMER)
            call descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_A2_DIMER)
            RAISE_ERROR("descriptor_MPI_setup: A2_dimer not MPI ready.", error)
         case(DT_AB_DIMER)
            RAISE_ERROR("descriptor_MPI_setup: AB_dimer not MPI ready.", error)
         case(DT_ATOM_REAL_SPACE)
            RAISE_ERROR("descriptor_MPI_setup: atom_real_space not MPI ready.", error)
         case(DT_POWER_SO3)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_POWER_SO4)
            RAISE_ERROR("descriptor_MPI_setup: power_SO4 not MPI ready.", error)
         case(DT_SOAP)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_RDF)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_AS_DISTANCE_2B)
            RAISE_ERROR("descriptor_MPI_setup: as_distance_2b not MPI ready.", error)
         case(DT_ALEX)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_DISTANCE_NB)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
         case(DT_SOAP_TURBO)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_MOLECULE_LO_D)
            RAISE_ERROR("descriptor_MPI_setup: molecule_lo_d not MPI ready.", error)
         case(DT_BOND_REAL_SPACE)
            RAISE_ERROR("descriptor_MPI_setup: bond_real_space not MPI ready.", error)
         case(DT_AN_MONOMER)
            RAISE_ERROR("descriptor_MPI_setup: AN_monomer not MPI ready.", error)
         case(DT_GENERAL_MONOMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_GENERAL_DIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_GENERAL_TRIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_COM_DIMER)
            call descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
         case(DT_SOAP_EXPRESS)
            call descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
#endif
         case default
            RAISE_ERROR("descriptor_MPI_setup: descriptor type "//this%descriptor_type//" not recognised.",error)
         endselect
      else
         mpi_mask = .true.
      endif

   endsubroutine descriptor_MPI_setup

   subroutine descriptor_atomic_MPI_setup(at,mpi,mpi_mask,error)
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      mpi_mask = .false.
      do i = 1, at%N
         if( mod(i-1, mpi%n_procs) == mpi%my_proc ) mpi_mask(i) = .true.
      enddo

   endsubroutine descriptor_atomic_MPI_setup

   subroutine descriptor_water_monomer_dimer_MPI_setup(at,mpi,mpi_mask,error)
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      mpi_mask = .false.
      do i = 1, at%N
         if( at%Z(i) == 8 .and. mod(i-1, mpi%n_procs) == mpi%my_proc ) mpi_mask(i) = .true.
      enddo

   endsubroutine descriptor_water_monomer_dimer_MPI_setup


   subroutine descriptor_data_finalise(this,error)
      type(descriptor_data), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(allocated(this%x)) then
         do i = 1, size(this%x)
            if(allocated(this%x(i)%data)) deallocate(this%x(i)%data)
            if(allocated(this%x(i)%grad_data)) deallocate(this%x(i)%grad_data)
            if(allocated(this%x(i)%ci)) deallocate(this%x(i)%ci)
            if(allocated(this%x(i)%ii)) deallocate(this%x(i)%ii)
            if(allocated(this%x(i)%pos)) deallocate(this%x(i)%pos)
            if(allocated(this%x(i)%has_grad_data)) deallocate(this%x(i)%has_grad_data)
            if(allocated(this%x(i)%grad_covariance_cutoff)) deallocate(this%x(i)%grad_covariance_cutoff)
         enddo
         deallocate(this%x)
      endif

   endsubroutine descriptor_data_finalise

   subroutine RadialFunction_initialise(this,n_max,cutoff, min_cutoff,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, intent(in) :: n_max
      real(dp), intent(in) :: cutoff, min_cutoff
      integer, optional, intent(out) :: error

      real(dp), dimension(:,:), allocatable :: S, vS
      real(dp), dimension(:), allocatable :: eS
      integer :: i, j

      INIT_ERROR(error)

      call finalise(this)

      this%n_max = n_max
      this%cutoff = cutoff
      this%min_cutoff = min_cutoff

      allocate(this%RadialTransform(this%n_max,this%n_max),this%NormFunction(this%n_max))
      allocate(S(this%n_max,this%n_max), vS(this%n_max,this%n_max), eS(this%n_max))

      do i = 1, this%n_max
         this%NormFunction(i) = sqrt(this%cutoff**(2.0_dp*i+5.0_dp)/(2.0_dp*i+5.0_dp))
         do j = 1, this%n_max
            S(j,i) = sqrt((2.0_dp*i+5)*(2.0_dp*j+5))/(i+j+5.0_dp)
         enddo
      enddo

      call diagonalise(S,eS,vS)
      this%RadialTransform = matmul(matmul(vS,diag(1.0_dp/sqrt(eS))),transpose(vS))

      if(allocated(S)) deallocate(S)
      if(allocated(vS)) deallocate(vS)
      if(allocated(eS)) deallocate(eS)

      this%initialised = .true.

   endsubroutine RadialFunction_initialise

   subroutine RadialFunction_finalise(this,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%n_max = 0

      if(allocated(this%RadialTransform)) deallocate(this%RadialTransform)
      if(allocated(this%NormFunction)) deallocate(this%NormFunction)

      this%initialised = .false.

   endsubroutine RadialFunction_finalise

   subroutine cplx_2d_array1_finalise(this)
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: this
      integer :: j

      if(allocated(this)) then
         do j = lbound(this,1), ubound(this,1)
            if(allocated(this(j)%mm)) deallocate(this(j)%mm)
         enddo
         deallocate(this)
      endif
   endsubroutine cplx_2d_array1_finalise

   subroutine cplx_3d_array2_finalise(this)
      type(cplx_3d), dimension(:,:), allocatable, intent(inout) :: this
      integer :: i, j

      if(allocated(this)) then
         do j = lbound(this,2), ubound(this,2)
            do i = lbound(this,1), ubound(this,1)
               if(allocated(this(i,j)%mm)) deallocate(this(i,j)%mm)
            enddo
         enddo
         deallocate(this)
      endif

   endsubroutine cplx_3d_array2_finalise

   subroutine fourier_SO4_calc(this,at,i,U,dU,args_str,error)
      type(fourier_SO4_type), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(in) :: i
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: U
      type(cplx_3d), dimension(:,:), allocatable, intent(inout), optional :: dU
      integer, optional, intent(out) :: error
      character(len=*), intent(in), optional :: args_str

      complex(dp), dimension(:,:), allocatable :: Uc, Up
      complex(dp), dimension(:,:,:), allocatable :: dUc, dUp
      complex(dp) :: z0_pls_Iz, z0_min_Iz, x_pls_Iy, x_min_Iy
      complex(dp), dimension(3) :: dz0_pls_Iz, dz0_min_Iz, dx_pls_Iy, dx_min_Iy
      real(dp), dimension(3) :: diff, u_ij, dfcut, dz0, dr0
      real(dp) :: r0, r, fcut, z0, theta0
      integer :: n, n_i, ji, j, m1, m2
      integer, dimension(total_elements) :: species_map

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR('fourier_SO4_calc: object not initialised',error)
      endif

      species_map = 0
      do j = 1, size(this%species_Z)
         if(this%species_Z(j) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(j)) = j
         endif
      enddo


      if(allocated(U)) then
         if(lbound(U,1) /= 0 .or. ubound(U,1) /= this%j_max) call finalise(U)
      endif

      if(.not.allocated(U)) then
         allocate( U(0:this%j_max) )
         do j = 0, this%j_max
            allocate( U(j)%mm(-j:j,-j:j) )
            U(j)%mm = CPLX_ZERO
         enddo
      endif

      do j = 0, this%j_max
         U(j)%mm = CPLX_ZERO
         do m1 = -j, j, 2
            U(j)%mm(m1,m1) = CPLX_ONE
         enddo
      enddo

      allocate( Uc(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        Up(-this%j_max:this%j_max, -this%j_max:this%j_max) )

      Uc = CPLX_ZERO
      Up = CPLX_ZERO

      if(present(dU)) then
         if(allocated(dU)) call finalise(dU)

         ! dU is not allocated, allocate and zero it
         allocate( dU(0:this%j_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
         do j = 0, this%j_max
            allocate( dU(j,0)%mm(3,-j:j,-j:j) )
            dU(j,0)%mm = CPLX_ZERO
         enddo

         allocate( dUc(3,-this%j_max:this%j_max, -this%j_max:this%j_max), &
            dUp(3,-this%j_max:this%j_max, -this%j_max:this%j_max) )
         dUc = CPLX_ZERO
         dUp = CPLX_ZERO
      endif

      n_i = 0
      do n = 1, n_neighbours(at,i)
         ji = neighbour(at, i, n, distance=r, diff=diff, cosines=u_ij)
         if( r >= this%cutoff ) cycle

         n_i = n_i + 1

         theta0 = r / this%z0
         z0 = r / tan( theta0 )
         r0 = sin( theta0 ) / r

         z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) * r0
         z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) * r0
         x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) * r0
         x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) * r0

         fcut = cos_cutoff_function(r,this%cutoff) * this%w(species_map(at%Z(ji)))

         U(0)%mm(0,0) = U(0)%mm(0,0) + fcut
         Up(0:0,0:0) = CPLX_ONE

         if(present(dU)) then

            dfcut = -dcos_cutoff_function(r,this%cutoff)*u_ij * this%w(species_map(at%Z(ji)))
            dz0 = ( 1.0_dp / tan( theta0 ) - theta0 / sin(theta0)**2 ) * u_ij
            dr0 = ( cos( theta0 ) / (r*this%z0) - r0 / r ) * u_ij

            dz0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_pls_Iz(3) = dz0_pls_Iz(3) + CPLX_IMAG*r0

            dz0_min_Iz = ( z0 - CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_min_Iz(3) = dz0_min_Iz(3) - CPLX_IMAG*r0

            dx_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) )*dr0
            dx_pls_Iy(1) = dx_pls_Iy(1) + r0
            dx_pls_Iy(2) = dx_pls_Iy(2) + CPLX_IMAG*r0

            dx_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) )*dr0
            dx_min_Iy(1) = dx_min_Iy(1) + r0
            dx_min_Iy(2) = dx_min_Iy(2) - CPLX_IMAG*r0

            dUc = CPLX_ZERO
            dUp = CPLX_ZERO

            dU(0,0)%mm(:,0,0) = dU(0,0)%mm(:,0,0) + dfcut*CPLX_ONE

            allocate( dU(0,n_i)%mm(3,-0:0,-0:0) )

            dU(0,n_i)%mm(:,0,0) = - dfcut*CPLX_ONE
         endif

         do j = 1, this%j_max
            Uc(-j:j,-j:j) = CPLX_ZERO
            if(present(dU)) then
               dUc(:,-j:j,-j:j) = CPLX_ZERO
               allocate( dU(j,n_i)%mm(3,-j:j,-j:j) )
               dU(j,n_i)%mm = CPLX_ZERO
            endif

            do m1 = -j, j-2, 2
               do m2 = -j, j, 2
                  if( (j-m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) + &
                     sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                        sqrt( real(j-m2,dp)/real(j-m1,dp) ) * &
                        ( dz0_pls_Iz * Up(m2+1,m1+1) + z0_pls_Iz * dUp(:,m2+1,m1+1) )
                  endif

                  if( (j+m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) - &
                     CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * Up(m2-1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                        CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * &
                        ( dx_min_Iy * Up(m2-1,m1+1) + x_min_Iy * dUp(:,m2-1,m1+1) )

                  endif
               enddo
            enddo

            m1 = j
            do m2 = -j, j, 2
               if( (j+m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) + &
                  sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                     sqrt( real(j+m2,dp)/real(j+m1,dp) ) * &
                     ( dz0_min_Iz * Up(m2-1,m1-1) + z0_min_Iz * dUp(:,m2-1,m1-1) )
               endif

               if( (j-m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) - &
                  CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * Up(m2+1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                     CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * &
                     ( dx_pls_Iy * Up(m2+1,m1-1) + x_pls_Iy * dUp(:,m2+1,m1-1) )
               endif
            enddo

            U(j)%mm = U(j)%mm + Uc(-j:j,-j:j) * fcut
            Up(-j:j,-j:j) = Uc(-j:j,-j:j)
            if(present(dU)) then
               dUp(:,-j:j,-j:j) = dUc(:,-j:j,-j:j)
               dU(j,0)%mm = dU(j,0)%mm - dUc(:,-j:j,-j:j) * fcut
               dU(j,n_i)%mm = dU(j,n_i)%mm + dUc(:,-j:j,-j:j) * fcut
               do m1 = -j, j, 2
                  do m2 = -j, j, 2
                     dU(j,0)%mm(:,m2,m1) = dU(j,0)%mm(:,m2,m1) &
                        + Uc(m2,m1) * dfcut
                     dU(j,n_i)%mm(:,m2,m1) = dU(j,n_i)%mm(:,m2,m1) &
                        - Uc(m2,m1) * dfcut
                  enddo
               enddo
            endif

         enddo ! j
      enddo ! n

      if(allocated(Up)) deallocate(Up)
      if(allocated(Uc)) deallocate(Uc)
      if(allocated(dUp)) deallocate(dUp)
      if(allocated(dUc)) deallocate(dUc)

   endsubroutine fourier_SO4_calc

   subroutine fourier_so4_initialise(this,args_str,error)
      type(fourier_SO4_type), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '2.75', this%cutoff, help_string="Cutoff for SO4 bispectrum")
      call param_register(params, 'z0_ratio', '0.0', this%z0_ratio, help_string="Ratio of radius of 4D projection sphere times PI and the cutoff.")
      call param_register(params, 'j_max', '4', this%j_max, help_string="Max of expansion of bispectrum, i.e. resulution")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z")
      call param_register(params, 'n_Z_environment', '1', n_species, help_string="Number of species for the descriptor", altkey="n_species")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'Z_environment', '0', this%species_Z(1), help_string="Atomic number of species", altkey="species_Z")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'Z_environment', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species", altkey="species_Z")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%z0 = max(1.0_dp,this%z0_ratio) * this%cutoff/(PI-0.02_dp)

      this%initialised = .true.


   endsubroutine fourier_so4_initialise

   subroutine fourier_so4_finalise(this,error)
      type(fourier_so4_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%j_max = 0
      this%z0_ratio = 0.0_dp
      this%z0 = 0.0_dp
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      this%initialised = .false.

   endsubroutine fourier_so4_finalise

   subroutine bispectrum_so4_initialise(this,args_str,error)
      type(bispectrum_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w

      this%initialised = .true.

   endsubroutine bispectrum_so4_initialise

   subroutine bispectrum_so4_finalise(this,error)
      type(bispectrum_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine bispectrum_so4_finalise

   subroutine bispectrum_so3_initialise(this,args_str,error)
      type(bispectrum_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)

      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for bispectrum_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in bispectrum_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for bispectrum_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for bispectrum_so3-type descriptors")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z")
      call param_register(params, 'n_Z_environment', '1', n_species, help_string="Number of species for the descriptor", altkey="n_species")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'Z_environment', '0', this%species_Z(1), help_string="Atomic number of species", altkey="species_Z")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'Z_environment', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species", altkey="species_Z")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

      call print('Dimensions: '//bispectrum_so3_dimensions(this,error))

   endsubroutine bispectrum_so3_initialise

   subroutine bispectrum_so3_finalise(this,error)
      type(bispectrum_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine bispectrum_so3_finalise

   subroutine behler_initialise(this,args_str,error)
      type(behler), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(len=STRING_LENGTH) :: specification_str, specification_file
      logical :: has_specification,has_specification_file
      integer :: n_fields, i_field, i_g2, i_g3, n, sym_type
      character(len=128), dimension(:), allocatable :: specification
      character(len=16), dimension(7) :: sym_func

      type(inoutput) :: specification_inout

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params,"Z","0",this%Z, help_string="Central atom")
      call param_register(params,"specification","",specification_str, help_string="String to specify Parrinello-Behler descriptors", &
         has_value_target=has_specification)
      call param_register(params,"specification_file","",specification_file, help_string="File containing string to specify Parrinello-Behler descriptors", &
         has_value_target=has_specification_file)

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='behler_initialise args_str')) then
         RAISE_ERROR("behler_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      if( has_specification .or. has_specification_file ) then
         if( has_specification .and. has_specification_file ) then
             RAISE_ERROR("behler_initialise: both specification and specification_file specified",error)
         endif

         if(has_specification_file) then
            call initialise(specification_inout,trim(specification_file))
            read(specification_inout%unit,'(a)') specification_str
            call finalise(specification_inout)
         endif

         n_fields = num_fields_in_string_simple(specification_str,"|")
         allocate(specification(n_fields))
         call split_string_simple(specification_str,specification,n_fields,"|")

         this%n_g2 = 0
         this%n_g3 = 0

         do i_field = 1, n_fields
            call split_string_simple(specification(i_field),sym_func,n,":")
            sym_type = string_to_int(sym_func(1),error)
            select case(sym_type)
            case(2)
               this%n_g2 = this%n_g2 + 1
            case(3)
               this%n_g3 = this%n_g3 + 1
            case default
               RAISE_ERROR("behler_initialise: unknown symmetry function type "//sym_type,error)
            endselect
         enddo

         allocate(this%g2(this%n_g2))
         allocate(this%g3(this%n_g3))

         i_g2 = 0
         i_g3 = 0
         this%cutoff = 0.0_dp
         do i_field = 1, n_fields
            call split_string_simple(specification(i_field),sym_func,n,":")

            sym_type = string_to_int(sym_func(1),error)
            select case(sym_type)
            case(2)
               i_g2 = i_g2 + 1
               this%g2(i_g2)%Z_n = atomic_number_from_symbol(sym_func(2))
               this%g2(i_g2)%eta = string_to_real(sym_func(3)) / BOHR**2
               this%g2(i_g2)%rs = string_to_real(sym_func(4)) * BOHR
               this%g2(i_g2)%rc = string_to_real(sym_func(5)) * BOHR
               this%cutoff = max(this%cutoff,this%g2(i_g2)%rc)
            case(3)
               i_g3 = i_g3 + 1
               this%g3(i_g3)%Z_n(1) = atomic_number_from_symbol(sym_func(2))
               this%g3(i_g3)%Z_n(2) = atomic_number_from_symbol(sym_func(3))
               this%g3(i_g3)%eta = string_to_real(sym_func(4)) / BOHR**2
               this%g3(i_g3)%lambda = string_to_real(sym_func(5))
               this%g3(i_g3)%zeta = string_to_real(sym_func(6))
               this%g3(i_g3)%rc = string_to_real(sym_func(7)) * BOHR
               this%cutoff = max(this%cutoff,this%g3(i_g3)%rc)
            case default
               RAISE_ERROR("behler_initialise: unknown symmetry function type "//sym_type,error)
            endselect
         enddo
      else
         ! Default, for backwards compatibility
         this%n_g2 = 8
         this%n_g3 = 43

         allocate(this%g2(this%n_g2), this%g3(this%n_g3))
         do i_g2 = 1, this%n_g2
            this%g2(i_g2)%Z_n = 0
         enddo
         do i_g3 = 1, this%n_g3
            this%g3(i_g3)%Z_n = 0
         enddo

         this%g2(1)%eta = 0.001_dp / BOHR**2; this%g2(1)%rs = 0.000_dp * BOHR; this%g2(1)%rc = 11.338_dp * BOHR
         this%g2(2)%eta = 0.010_dp / BOHR**2; this%g2(2)%rs = 0.000_dp * BOHR; this%g2(2)%rc = 11.338_dp * BOHR
         this%g2(3)%eta = 0.020_dp / BOHR**2; this%g2(3)%rs = 0.000_dp * BOHR; this%g2(3)%rc = 11.338_dp * BOHR
         this%g2(4)%eta = 0.035_dp / BOHR**2; this%g2(4)%rs = 0.000_dp * BOHR; this%g2(4)%rc = 11.338_dp * BOHR
         this%g2(5)%eta = 0.060_dp / BOHR**2; this%g2(5)%rs = 0.000_dp * BOHR; this%g2(5)%rc = 11.338_dp * BOHR
         this%g2(6)%eta = 0.100_dp / BOHR**2; this%g2(6)%rs = 0.000_dp * BOHR; this%g2(6)%rc = 11.338_dp * BOHR
         this%g2(7)%eta = 0.200_dp / BOHR**2; this%g2(7)%rs = 0.000_dp * BOHR; this%g2(7)%rc = 11.338_dp * BOHR
         this%g2(8)%eta = 0.400_dp / BOHR**2; this%g2(8)%rs = 0.000_dp * BOHR; this%g2(8)%rc = 11.338_dp * BOHR

         this%g3( 1)%eta = 0.0001_dp / BOHR**2; this%g3( 1)%lambda = -1.000_dp; this%g3( 1)%zeta =  1.000_dp; this%g3( 1)%rc = 11.338_dp * BOHR
         this%g3( 2)%eta = 0.0001_dp / BOHR**2; this%g3( 2)%lambda =  1.000_dp; this%g3( 2)%zeta =  1.000_dp; this%g3( 2)%rc = 11.338_dp * BOHR
         this%g3( 3)%eta = 0.0001_dp / BOHR**2; this%g3( 3)%lambda = -1.000_dp; this%g3( 3)%zeta =  2.000_dp; this%g3( 3)%rc = 11.338_dp * BOHR
         this%g3( 4)%eta = 0.0001_dp / BOHR**2; this%g3( 4)%lambda =  1.000_dp; this%g3( 4)%zeta =  2.000_dp; this%g3( 4)%rc = 11.338_dp * BOHR
         this%g3( 5)%eta = 0.0030_dp / BOHR**2; this%g3( 5)%lambda = -1.000_dp; this%g3( 5)%zeta =  1.000_dp; this%g3( 5)%rc = 11.338_dp * BOHR
         this%g3( 6)%eta = 0.0030_dp / BOHR**2; this%g3( 6)%lambda =  1.000_dp; this%g3( 6)%zeta =  1.000_dp; this%g3( 6)%rc = 11.338_dp * BOHR
         this%g3( 7)%eta = 0.0030_dp / BOHR**2; this%g3( 7)%lambda = -1.000_dp; this%g3( 7)%zeta =  2.000_dp; this%g3( 7)%rc = 11.338_dp * BOHR
         this%g3( 8)%eta = 0.0030_dp / BOHR**2; this%g3( 8)%lambda =  1.000_dp; this%g3( 8)%zeta =  2.000_dp; this%g3( 8)%rc = 11.338_dp * BOHR
         this%g3( 9)%eta = 0.0080_dp / BOHR**2; this%g3( 9)%lambda = -1.000_dp; this%g3( 9)%zeta =  1.000_dp; this%g3( 9)%rc = 11.338_dp * BOHR
         this%g3(10)%eta = 0.0080_dp / BOHR**2; this%g3(10)%lambda =  1.000_dp; this%g3(10)%zeta =  1.000_dp; this%g3(10)%rc = 11.338_dp * BOHR
         this%g3(11)%eta = 0.0080_dp / BOHR**2; this%g3(11)%lambda = -1.000_dp; this%g3(11)%zeta =  2.000_dp; this%g3(11)%rc = 11.338_dp * BOHR
         this%g3(12)%eta = 0.0080_dp / BOHR**2; this%g3(12)%lambda =  1.000_dp; this%g3(12)%zeta =  2.000_dp; this%g3(12)%rc = 11.338_dp * BOHR
         this%g3(13)%eta = 0.0150_dp / BOHR**2; this%g3(13)%lambda = -1.000_dp; this%g3(13)%zeta =  1.000_dp; this%g3(13)%rc = 11.338_dp * BOHR
         this%g3(14)%eta = 0.0150_dp / BOHR**2; this%g3(14)%lambda =  1.000_dp; this%g3(14)%zeta =  1.000_dp; this%g3(14)%rc = 11.338_dp * BOHR
         this%g3(15)%eta = 0.0150_dp / BOHR**2; this%g3(15)%lambda = -1.000_dp; this%g3(15)%zeta =  2.000_dp; this%g3(15)%rc = 11.338_dp * BOHR
         this%g3(16)%eta = 0.0150_dp / BOHR**2; this%g3(16)%lambda =  1.000_dp; this%g3(16)%zeta =  2.000_dp; this%g3(16)%rc = 11.338_dp * BOHR
         this%g3(17)%eta = 0.0150_dp / BOHR**2; this%g3(17)%lambda = -1.000_dp; this%g3(17)%zeta =  4.000_dp; this%g3(17)%rc = 11.338_dp * BOHR
         this%g3(18)%eta = 0.0150_dp / BOHR**2; this%g3(18)%lambda =  1.000_dp; this%g3(18)%zeta =  4.000_dp; this%g3(18)%rc = 11.338_dp * BOHR
         this%g3(19)%eta = 0.0150_dp / BOHR**2; this%g3(19)%lambda = -1.000_dp; this%g3(19)%zeta = 16.000_dp; this%g3(19)%rc = 11.338_dp * BOHR
         this%g3(20)%eta = 0.0150_dp / BOHR**2; this%g3(20)%lambda =  1.000_dp; this%g3(20)%zeta = 16.000_dp; this%g3(20)%rc = 11.338_dp * BOHR
         this%g3(21)%eta = 0.0250_dp / BOHR**2; this%g3(21)%lambda = -1.000_dp; this%g3(21)%zeta =  1.000_dp; this%g3(21)%rc = 11.338_dp * BOHR
         this%g3(22)%eta = 0.0250_dp / BOHR**2; this%g3(22)%lambda =  1.000_dp; this%g3(22)%zeta =  1.000_dp; this%g3(22)%rc = 11.338_dp * BOHR
         this%g3(23)%eta = 0.0250_dp / BOHR**2; this%g3(23)%lambda = -1.000_dp; this%g3(23)%zeta =  2.000_dp; this%g3(23)%rc = 11.338_dp * BOHR
         this%g3(24)%eta = 0.0250_dp / BOHR**2; this%g3(24)%lambda =  1.000_dp; this%g3(24)%zeta =  2.000_dp; this%g3(24)%rc = 11.338_dp * BOHR
         this%g3(25)%eta = 0.0250_dp / BOHR**2; this%g3(25)%lambda = -1.000_dp; this%g3(25)%zeta =  4.000_dp; this%g3(25)%rc = 11.338_dp * BOHR
         this%g3(26)%eta = 0.0250_dp / BOHR**2; this%g3(26)%lambda =  1.000_dp; this%g3(26)%zeta =  4.000_dp; this%g3(26)%rc = 11.338_dp * BOHR
         this%g3(27)%eta = 0.0250_dp / BOHR**2; this%g3(27)%lambda = -1.000_dp; this%g3(27)%zeta = 16.000_dp; this%g3(27)%rc = 11.338_dp * BOHR
         this%g3(28)%eta = 0.0250_dp / BOHR**2; this%g3(28)%lambda =  1.000_dp; this%g3(28)%zeta = 16.000_dp; this%g3(28)%rc = 11.338_dp * BOHR
         this%g3(29)%eta = 0.0450_dp / BOHR**2; this%g3(29)%lambda = -1.000_dp; this%g3(29)%zeta =  1.000_dp; this%g3(29)%rc = 11.338_dp * BOHR
         this%g3(30)%eta = 0.0450_dp / BOHR**2; this%g3(30)%lambda =  1.000_dp; this%g3(30)%zeta =  1.000_dp; this%g3(30)%rc = 11.338_dp * BOHR
         this%g3(31)%eta = 0.0450_dp / BOHR**2; this%g3(31)%lambda = -1.000_dp; this%g3(31)%zeta =  2.000_dp; this%g3(31)%rc = 11.338_dp * BOHR
         this%g3(32)%eta = 0.0450_dp / BOHR**2; this%g3(32)%lambda =  1.000_dp; this%g3(32)%zeta =  2.000_dp; this%g3(32)%rc = 11.338_dp * BOHR
         this%g3(33)%eta = 0.0450_dp / BOHR**2; this%g3(33)%lambda = -1.000_dp; this%g3(33)%zeta =  4.000_dp; this%g3(33)%rc = 11.338_dp * BOHR
         this%g3(34)%eta = 0.0450_dp / BOHR**2; this%g3(34)%lambda =  1.000_dp; this%g3(34)%zeta =  4.000_dp; this%g3(34)%rc = 11.338_dp * BOHR
         this%g3(35)%eta = 0.0450_dp / BOHR**2; this%g3(35)%lambda = -1.000_dp; this%g3(35)%zeta = 16.000_dp; this%g3(35)%rc = 11.338_dp * BOHR
         this%g3(36)%eta = 0.0450_dp / BOHR**2; this%g3(36)%lambda =  1.000_dp; this%g3(36)%zeta = 16.000_dp; this%g3(36)%rc = 11.338_dp * BOHR
         this%g3(37)%eta = 0.0800_dp / BOHR**2; this%g3(37)%lambda = -1.000_dp; this%g3(37)%zeta =  1.000_dp; this%g3(37)%rc = 11.338_dp * BOHR
         this%g3(38)%eta = 0.0800_dp / BOHR**2; this%g3(38)%lambda =  1.000_dp; this%g3(38)%zeta =  1.000_dp; this%g3(38)%rc = 11.338_dp * BOHR
         this%g3(39)%eta = 0.0800_dp / BOHR**2; this%g3(39)%lambda = -1.000_dp; this%g3(39)%zeta =  2.000_dp; this%g3(39)%rc = 11.338_dp * BOHR
         this%g3(40)%eta = 0.0800_dp / BOHR**2; this%g3(40)%lambda =  1.000_dp; this%g3(40)%zeta =  2.000_dp; this%g3(40)%rc = 11.338_dp * BOHR
         this%g3(41)%eta = 0.0800_dp / BOHR**2; this%g3(41)%lambda = -1.000_dp; this%g3(41)%zeta =  4.000_dp; this%g3(41)%rc = 11.338_dp * BOHR
         this%g3(42)%eta = 0.0800_dp / BOHR**2; this%g3(42)%lambda =  1.000_dp; this%g3(42)%zeta =  4.000_dp; this%g3(42)%rc = 11.338_dp * BOHR
         this%g3(43)%eta = 0.0800_dp / BOHR**2; this%g3(43)%lambda =  1.000_dp; this%g3(43)%zeta = 16.000_dp; this%g3(43)%rc = 11.338_dp * BOHR

         this%cutoff = 11.338_dp * BOHR
      endif

      if( allocated(specification) ) deallocate(specification)

      this%initialised = .true.

   endsubroutine behler_initialise

   subroutine behler_finalise(this,error)
      type(behler), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%n_g2 = 0
      this%n_g3 = 0

      if(allocated(this%g2)) deallocate(this%g2)
      if(allocated(this%g3)) deallocate(this%g3)

      this%initialised = .false.

   endsubroutine behler_finalise

   subroutine distance_2b_initialise(this,args_str,error)
      type(distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      logical :: has_resid_name, has_exponents
      integer :: i

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for distance_2b-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.5', this%cutoff_transition_width, help_string="Transition width of cutoff for distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")
      call param_register(params, 'resid_name', '', this%resid_name, has_value_target=has_resid_name, help_string="Name of an integer property in the atoms object giving the residue id of the molecule to which the atom belongs.")
      call param_register(params, 'only_intra', 'F', this%only_intra, help_string="Only calculate INTRAmolecular pairs with equal residue ids (bonds)")
      call param_register(params, 'only_inter', 'F', this%only_inter, help_string="Only apply to INTERmolecular pairs with different residue ids (non-bonded)")

      call param_register(params, 'n_exponents', '1', this%n_exponents, help_string="Number of exponents")
      call param_register(params, 'tail_range', '1.0', this%tail_range, help_string="Tail order")
      call param_register(params, 'tail_exponent', '0', this%tail_exponent, &
         has_value_target = this%has_tail, help_string="Tail range")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_2b_initialise args_str')) then
         RAISE_ERROR("distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (this%only_intra .and. this%only_inter) then
         RAISE_ERROR("distance_2b_initialise: cannot specify both only_inter AND only_intra", error)
      end if
      if ((this%only_intra .or. this%only_inter) .and. (.not. has_resid_name)) then
         RAISE_ERROR("distance_2b_initialise: only_intra and only_inter require resid_name to be given as well", error)
      end if

      allocate(this%exponents(this%n_exponents))
      call initialise(params)
      if( this%n_exponents == 1 ) then
         call param_register(params, 'exponents',"1", this%exponents(1), &
            has_value_target=has_exponents,help_string="Exponents")
      else
         call param_register(params, 'exponents',repeat(" 1 ",this%n_exponents), this%exponents, &
            has_value_target=has_exponents,help_string="Exponents")
      endif
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_2b_initialise args_str')) then
         RAISE_ERROR("distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( .not. has_exponents .and. this%n_exponents > 1 ) then
         do i = 1, this%n_exponents
            this%exponents(i) = -i
         enddo
      endif

      this%initialised = .true.

   endsubroutine distance_2b_initialise

   subroutine distance_2b_finalise(this,error)
      type(distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.5_dp
      this%Z1 = 0
      this%Z2 = 0

      this%resid_name = ''
      this%only_intra = .false.
      this%only_inter = .false.

      this%tail_exponent = 0
      this%tail_range = 0.0_dp
      this%has_tail = .false.

      this%n_exponents = 0
      if(allocated(this%exponents)) deallocate(this%exponents)

      this%initialised = .false.

   endsubroutine distance_2b_finalise

   subroutine coordination_initialise(this,args_str,error)
      type(coordination), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)
      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for coordination-type descriptors")
      call param_register(params, 'transition_width', '0.20', this%transition_width, help_string="Width of transition region from 1 to 0")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z_center")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='coordination_initialise args_str')) then
         RAISE_ERROR("coordination_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine coordination_initialise

   subroutine coordination_finalise(this,error)
      type(coordination), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%Z = 0

      this%initialised = .false.

   endsubroutine coordination_finalise

   subroutine angle_3b_initialise(this,args_str,error)
      type(angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for angle_3b-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Cutoff transition width for angle_3b-type descriptors")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='angle_3b_initialise args_str')) then
         RAISE_ERROR("angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine angle_3b_initialise

   subroutine angle_3b_finalise(this,error)
      type(angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine angle_3b_finalise

   subroutine co_angle_3b_initialise(this,args_str,error)
      type(co_angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_angle_3b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_angle_3b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_angle_3b-type descriptors")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_angle_3b_initialise args_str')) then
         RAISE_ERROR("co_angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_angle_3b_initialise

   subroutine co_angle_3b_finalise(this,error)
      type(co_angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_angle_3b_finalise

   subroutine co_distance_2b_initialise(this,args_str,error)
      type(co_distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'transition_width', '0.50', this%transition_width, help_string="Transition width of cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_distance_2b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_distance_2b_initialise args_str')) then
         RAISE_ERROR("co_distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_distance_2b_initialise

   subroutine co_distance_2b_finalise(this,error)
      type(co_distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_distance_2b_finalise

   subroutine cosnx_initialise(this,args_str,error)
      type(cosnx), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for cosnx-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in cosnx-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for cosnx-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for cosnx-type descriptors")
      call param_register(params, 'Z_center', '0', this%Z, help_string="Atomic number of central atom", altkey="Z")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine cosnx_initialise

   subroutine cosnx_finalise(this,error)
      type(cosnx), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine cosnx_finalise

   subroutine trihis_initialise(this,args_str,error)
      type(trihis), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp), dimension(:), allocatable :: gauss_centre1D, gauss_width1D

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for trihis-type descriptors")
      call param_register(params, 'n_gauss', '0', this%n_gauss, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(gauss_centre1D(3*this%n_gauss),gauss_width1D(3*this%n_gauss))
      allocate(this%gauss_centre(3,this%n_gauss),this%gauss_width(3,this%n_gauss))

      call initialise(params)
      call param_register(params, 'trihis_gauss_centre', PARAM_MANDATORY, gauss_centre1D, help_string="Number of Gaussians for trihis-type descriptors")
      call param_register(params, 'trihis_gauss_width', PARAM_MANDATORY, gauss_width1D, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%gauss_centre = reshape(gauss_centre1D,(/3,this%n_gauss/))
      this%gauss_width = reshape(gauss_width1D,(/3,this%n_gauss/))

      deallocate(gauss_centre1D,gauss_width1D)

      this%initialised = .true.

   endsubroutine trihis_initialise

   subroutine trihis_finalise(this,error)
      type(trihis), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%n_gauss = 0

      if(allocated(this%gauss_centre)) deallocate(this%gauss_centre)
      if(allocated(this%gauss_width)) deallocate(this%gauss_width)

      this%initialised = .false.

   endsubroutine trihis_finalise

   subroutine water_monomer_initialise(this,args_str,error)
      type(water_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_monomer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_monomer_initialise args_str')) then
         RAISE_ERROR("water_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_monomer_initialise

   subroutine water_monomer_finalise(this,error)
      type(water_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp

      this%initialised = .false.

   endsubroutine water_monomer_finalise

   subroutine water_dimer_initialise(this,args_str,error)
      type(water_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for water_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for water_dimer-type descriptors")
      call param_register(params, 'OHH_ordercheck', 'T', this%OHH_ordercheck, help_string="T: find water molecules. F: use default order OHH")
      call param_register(params, 'power', '1.0', this%power, help_string="Power of distances to be used in the kernel")
      call param_register(params, 'dist_shift', '0.0', this%dist_shift, help_string="Distance shift for inverse distance descriptors.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_dimer_initialise args_str')) then
         RAISE_ERROR("water_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_dimer_initialise

   subroutine water_dimer_finalise(this,error)
      type(water_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%OHH_ordercheck = .true.
      this%power = 1.0_dp
      this%dist_shift = 0.0_dp

      this%initialised = .false.

   endsubroutine water_dimer_finalise

   subroutine A2_dimer_initialise(this,args_str,error)
      type(A2_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for A2_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for A2_dimer-type descriptors")
      call param_register(params, 'atomic_number', '1', this%atomic_number, help_string="Atomic number in A2_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='A2_dimer_initialise args_str')) then
         RAISE_ERROR("A2_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine A2_dimer_initialise

   subroutine A2_dimer_finalise(this,error)
      type(A2_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number = 0

      this%initialised = .false.

   endsubroutine A2_dimer_finalise

   subroutine AB_dimer_initialise(this,args_str,error)
      type(AB_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for AB_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for AB_dimer-type descriptors")
      call param_register(params, 'atomic_number1', '1', this%atomic_number1, help_string="Atomic number of atom 1 in AB_dimer-type descriptors")
      call param_register(params, 'atomic_number2', '9', this%atomic_number2, help_string="Atomic number of atom 2 in AB_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='AB_dimer_initialise args_str')) then
         RAISE_ERROR("AB_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( this%atomic_number1 == this%atomic_number2 ) then
         RAISE_ERROR("AB_dimer_initialise: AB_dimer_atomic_number1 = AB_dimer_atomic_number2 = "//this%atomic_number1//" which would require addtional permutational symmetries. Use A2_dimer descriptor instead.",error)
      endif

      this%initialised = .true.

   endsubroutine AB_dimer_initialise

   subroutine AB_dimer_finalise(this,error)
      type(AB_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number1 = 0
      this%atomic_number2 = 0

      this%initialised = .false.

   endsubroutine AB_dimer_finalise


   subroutine atom_real_space_initialise(this,args_str,error)
      type(atom_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for atom_real_space-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.00', this%cutoff_transition_width, help_string="Space transition width for atom_real_space-type descriptors")
      call param_register(params, 'l_max', '0', this%l_max, help_string="Cutoff for spherical harmonics expansion")
      call param_register(params, 'alpha', '1.0', this%alpha, help_string="Width of atomic Gaussians")
      call param_register(params, 'zeta', '1.0', this%zeta, help_string="Exponent of covariance function")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='atom_real_space_initialise args_str')) then
         RAISE_ERROR("atom_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine atom_real_space_initialise

   subroutine atom_real_space_finalise(this,error)
      type(atom_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp
      this%zeta = 0.0_dp

      this%initialised = .false.

   endsubroutine atom_real_space_finalise

   subroutine power_so3_initialise(this,args_str,error)
      type(power_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for power_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in power_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for power_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for power_so3-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine power_so3_initialise

   subroutine power_so3_finalise(this,error)
      type(power_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine power_so3_finalise

   subroutine power_so4_initialise(this,args_str,error)
      type(power_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w

      this%initialised = .true.

   endsubroutine power_so4_initialise

   subroutine power_so4_finalise(this,error)
      type(power_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine power_so4_finalise

   subroutine soap_initialise(this,args_str,error)
      type(soap), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp) :: alpha_basis, spacing_basis, cutoff_basis, basis_error_exponent
      real(dp), dimension(:,:), allocatable :: covariance_basis, overlap_basis, cholesky_overlap_basis
      integer :: i, j, xml_version

      type(LA_Matrix) :: LA_covariance_basis, LA_overlap_basis
      character(len=STRING_LENGTH) :: species_Z_str
      logical :: has_n_species, has_species_Z, has_central_reference_all_species


      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', PARAM_MANDATORY, this%cutoff, help_string="Cutoff for soap-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Cutoff transition width for soap-type descriptors")

      call param_register(params, 'cutoff_dexp', '0', this%cutoff_dexp, help_string="Cutoff decay exponent")
      call param_register(params, 'cutoff_scale', '1.0', this%cutoff_scale, help_string="Cutoff decay scale")
      call param_register(params, 'cutoff_rate', '1.0', this%cutoff_rate, help_string="Inverse cutoff decay rate")

      call param_register(params, 'l_max', PARAM_MANDATORY, this%l_max, help_string="L_max (spherical harmonics basis band limit) for soap-type descriptors")
      call param_register(params, 'n_max', PARAM_MANDATORY, this%n_max, help_string="N_max (number of radial basis functions) for soap-type descriptors")
      call param_register(params, 'atom_gaussian_width', PARAM_MANDATORY, this%atom_sigma, help_string="Width of atomic Gaussians for soap-type descriptors", altkey='atom_sigma')
      call param_register(params, 'central_weight', '1.0', this%central_weight, help_string="Weight of central atom in environment")
      call param_register(params, 'central_reference_all_species', 'F', this%central_reference_all_species, has_value_target=has_central_reference_all_species, &
           help_string="Place a Gaussian reference for all atom species densities."// &
           "By default (F) only consider when neighbour is the same species as centre")
      call param_register(params, 'average', 'F', this%global, help_string="Whether to calculate averaged SOAP - one descriptor per atoms object. If false (default) atomic SOAP is returned.")
      call param_register(params, 'diagonal_radial', 'F', this%diagonal_radial, help_string="Only return the n1=n2 elements of the power spectrum.")

      call param_register(params, 'covariance_sigma0', '0.0', this%covariance_sigma0, help_string="sigma_0 parameter in polynomial covariance function")
      call param_register(params, 'normalise', 'T', this%normalise, help_string="Normalise descriptor so magnitude is 1. In this case the kernel of two equivalent environments is 1.", altkey="normalize")
      call param_register(params, 'basis_error_exponent', '10.0', basis_error_exponent, help_string="10^(-basis_error_exponent) is the max difference between the target and the expanded function")

      call param_register(params, 'n_Z', '1', this%n_Z, help_string="How many different types of central atoms to consider")
      call param_register(params, 'n_species', '1', this%n_species, has_value_target=has_n_species, help_string="Number of species for the descriptor")
      call param_register(params, 'species_Z', '', species_Z_str, has_value_target=has_species_Z, help_string="Atomic number of species")
      call param_register(params, 'xml_version', '1426512068', xml_version, help_string="Version of GAP the XML potential file was created")

      call param_register(params, 'nu_R', '2', this%nu_R, help_string="radially sensitive correlation order")
      call param_register(params, 'nu_S', '2', this%nu_S, help_string="species sensitive correlation order")
      call param_register(params, 'Z_mix', 'F', this%Z_mix, help_string="mix Z channels together")
      call param_register(params, 'R_mix', 'F', this%R_mix, help_string="mix radial channels together")
      call param_register(params, 'sym_mix', 'F', this%sym_mix, help_string="symmetric mixing")
      call param_register(params, 'coupling', 'T', this%coupling, help_string="Full tensor product(=T) or Elementwise product(=F) between density channels")
      call param_register(params, 'K', '0', this%K, help_string="Number of mixing channels to create")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_initialise args_str')) then
         RAISE_ERROR("soap_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      ! backwards compatibility: the default used to be different before this version number
      if( xml_version < 1426512068 ) this%central_reference_all_species = .true.

      allocate(this%species_Z(0:this%n_species))
      allocate(this%Z(this%n_Z))
      this%species_Z(0)=0

      if( has_species_Z .and. .not. has_n_species ) then
         RAISE_ERROR("soap_initialise: is species_Z is present, n_species must be present, too.",error)
      endif

      call initialise(params)

      if( this%cutoff_dexp < 0 ) then
         RAISE_ERROR("soap_initialise: cutoff_dexp may not be less than 0",error)
      endif

      if( this%cutoff_scale <= 0.0_dp ) then
         RAISE_ERROR("soap_initialise: cutoff_scale must be greater than 0",error)
      endif

      if( this%cutoff_rate < 0.0_dp ) then
         RAISE_ERROR("soap_initialise: cutoff_rate may not be less than 0",error)
      endif

      if( has_n_species ) then
         if(this%n_species == 1) then
            call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         else
            call param_register(params, 'species_Z', '//MANDATORY//', this%species_Z(1:this%n_species), help_string="Atomic number of species")
         endif
      else
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
      endif

      if( .not. has_central_reference_all_species .and. this%n_species == 1 ) this%central_reference_all_species = .true.

      if( this%n_Z == 1 ) then
         call param_register(params, 'Z', '0', this%Z(1), help_string="Atomic number of central atom, 0 is the wild-card")
      else
         call param_register(params, 'Z', '//MANDATORY//', this%Z, help_string="Atomic numbers to be considered for central atom, must be a list")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_initialise args_str')) then
         RAISE_ERROR("soap_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)


      this%alpha = 0.5_dp / this%atom_sigma**2
      alpha_basis = this%alpha
      cutoff_basis = this%cutoff + this%atom_sigma * sqrt(2.0_dp * basis_error_exponent * log(10.0_dp))
      spacing_basis = cutoff_basis / this%n_max

      allocate(this%r_basis(this%n_max), this%transform_basis(this%n_max,this%n_max), &
         covariance_basis(this%n_max,this%n_max), overlap_basis(this%n_max,this%n_max), this%cholesky_overlap_basis(this%n_max,this%n_max))

      !this%r_basis(this%n_max) = cutoff_basis
      !do i = this%n_max-1, 1, -1
      !   this%r_basis(i)  = this%r_basis(i+1) - spacing_basis
      !enddo

      this%r_basis(1) = 0.0_dp
      do i = 2, this%n_max
         this%r_basis(i)  = this%r_basis(i-1) + spacing_basis
      enddo


      do i = 1, this%n_max
         do j = 1, this%n_max
            covariance_basis(j,i) = exp(-alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)
            !overlap_basis(j,i) = exp(-0.5_dp * alpha_basis* (this%r_basis(i) - this%r_basis(j))**2) * ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !print*, 'A', exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) )
            !print*, 'B', sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j))
            !print*, 'C', alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 )
            !print*, 'D', ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !overlap_basis(j,i) = exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
            !   ( sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
            !   alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
            !   ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )

            overlap_basis(j,i) = ( exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
               sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
               alpha_basis*exp(-0.5_dp * alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
               ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )
         enddo
      enddo

      !overlap_basis = overlap_basis * sqrt(pi / ( 8.0_dp * alpha_basis ) )
      overlap_basis = overlap_basis / sqrt(128.0_dp * alpha_basis**5)

      call initialise(LA_covariance_basis,covariance_basis)
      call initialise(LA_overlap_basis,overlap_basis)
      call LA_Matrix_Factorise(LA_overlap_basis, this%cholesky_overlap_basis)
      do i = 1, this%n_max
         do j = 1, i-1 !i + 1, this%n_max
            this%cholesky_overlap_basis(j,i) = 0.0_dp
         enddo
      enddo

      call Matrix_Solve(LA_covariance_basis,this%cholesky_overlap_basis,this%transform_basis)

      call finalise(LA_covariance_basis)
      call finalise(LA_overlap_basis)

      if(allocated(covariance_basis)) deallocate(covariance_basis)
      if(allocated(overlap_basis)) deallocate(overlap_basis)
      if(allocated(cholesky_overlap_basis)) deallocate(cholesky_overlap_basis)

      this%initialised = .true.

   endsubroutine soap_initialise

   subroutine soap_finalise(this,error)
      type(soap), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff_dexp = 0
      this%cutoff_scale = 1.0_dp
      this%cutoff_rate = 1.0_dp
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp
      this%central_weight = 0.0_dp
      this%central_reference_all_species = .false.
      this%global = .false.
      this%diagonal_radial = .false.
      this%covariance_sigma0 = 0.0_dp
      this%normalise = .true.

      this%n_max = 0
      this%n_Z = 0
      this%n_species = 0
      this%nu_R = 2
      this%nu_S = 2

      this%Z_mix = .false.
      this%R_mix = .false.
      this%sym_mix = .false.
      this%coupling = .true.
      this%K = 0

      if(allocated(this%r_basis)) deallocate(this%r_basis)
      if(allocated(this%transform_basis)) deallocate(this%transform_basis)
      if(allocated(this%cholesky_overlap_basis)) deallocate(this%cholesky_overlap_basis)
      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%Z)) deallocate(this%Z)

      this%initialised = .false.

   endsubroutine soap_finalise


   subroutine rdf_initialise(this,args_str,error)
      type(rdf), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: i
      real(dp) :: r_min, r_max
      logical :: has_r_max, has_w_gauss

      INIT_ERROR(error)

      call finalise(this)
      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for rdf-type descriptors")
      call param_register(params, 'transition_width', '0.20', this%transition_width, help_string="Width of transition region from 1 to 0")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'r_min', '0.0', r_min, help_string="Atomic number of central atom")
      call param_register(params, 'r_max', '0.0', r_max, has_value_target = has_r_max, help_string="Atomic number of central atom")
      call param_register(params, 'n_gauss', '10', this%n_gauss, help_string="Atomic number of central atom")
      call param_register(params, 'w_gauss', '0.0', this%w_gauss, has_value_target = has_w_gauss, help_string="Atomic number of central atom")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='rdf_initialise args_str')) then
         RAISE_ERROR("rdf_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%r_gauss(this%n_gauss))
      if(.not. has_w_gauss) this%w_gauss = this%cutoff / this%n_gauss * 2.0_dp
      if(.not. has_r_max) r_max = this%cutoff - this%w_gauss / 2.0_dp
      this%r_gauss = real( (/(i,i=1,this%n_gauss)/), kind=dp ) / real(this%n_gauss,kind=dp) * (r_max - r_min) + r_min

      this%initialised = .true.

   endsubroutine rdf_initialise

   subroutine rdf_finalise(this,error)
      type(rdf), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%Z = 0
      this%n_gauss = 0
      if( allocated(this%r_gauss) ) deallocate(this%r_gauss)

      this%initialised = .false.

   endsubroutine rdf_finalise

   subroutine as_distance_2b_initialise(this,args_str,error)
      type(as_distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Lower cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'max_cutoff', PARAM_MANDATORY, this%max_cutoff, help_string="Higher cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'as_cutoff', PARAM_MANDATORY, this%as_cutoff, help_string="Cutoff of asymmetricity")
      call param_register(params, 'overlap_alpha', '0.50', this%as_cutoff, help_string="Cutoff of asymmetricity")
      call param_register(params, 'min_transition_width', '0.50', this%min_transition_width, help_string="Transition width of lower cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'max_transition_width', '0.50', this%max_transition_width, help_string="Transition width of higher cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'as_transition_width', '0.10', this%as_transition_width, help_string="Transition width of asymmetricity cutoff for as_distance_2b-type descriptors")
      call param_register(params, 'coordination_cutoff', PARAM_MANDATORY, this%coordination_cutoff, help_string="Cutoff for coordination function in as_distance_2b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.50', this%coordination_transition_width, help_string="Transition width for as_distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='as_distance_2b_initialise args_str')) then
         RAISE_ERROR("as_distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine as_distance_2b_initialise

   subroutine as_distance_2b_finalise(this,error)
      type(as_distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%min_cutoff = 0.0_dp
      this%max_cutoff = 0.0_dp
      this%as_cutoff = 0.0_dp
      this%overlap_alpha = 0.0_dp
      this%min_transition_width = 0.0_dp
      this%max_transition_width = 0.0_dp
      this%as_transition_width = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine as_distance_2b_finalise


   subroutine alex_initialise(this,args_str,error)
      type(alex), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for alex-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'power_min', '5', this%power_min, help_string="Minimum power of radial basis for the descriptor")
      call param_register(params, 'power_max', '10', this%power_max, help_string="Maximum power of the radial basis for the descriptor")
      call param_register(params, 'n_species', '1', this%n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='alex_initialise args_str')) then
         RAISE_ERROR("alex_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(this%n_species))

      call initialise(params)
      if( this%n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='alex_initialise args_str')) then
         RAISE_ERROR("alex_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine alex_initialise

   subroutine alex_finalise(this,error)
      type(alex), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp

      if(allocated(this%species_Z)) deallocate(this%species_Z)

      this%initialised = .false.

   endsubroutine alex_finalise

   subroutine distance_Nb_initialise(this,args_str,error)
      type(distance_Nb), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(len=STRING_LENGTH) :: default_Z = ""
      integer :: i, j, k, i_p
      integer :: nEdges, nConnectivities, nMonomerConnectivities
      integer, dimension(:), allocatable :: n_permutations, connectivityList
      integer, dimension(:,:), allocatable :: atom_permutations, distance_matrix_index, edges
      integer :: xml_version
      logical :: has_compact_clusters

      logical, dimension(:,:,:), allocatable :: allConnectivities


      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', PARAM_MANDATORY, this%cutoff, help_string="Cutoff for distance_Nb-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.5', this%cutoff_transition_width, help_string="Transition width of cutoff for distance_Nb-type descriptors")
      call param_register(params, 'order', PARAM_MANDATORY, this%order, help_string="Many-body order, in terms of number of neighbours")
      call param_register(params, 'compact_clusters', "T", this%compact_clusters, help_string="If true, generate clusters where the atoms have at least one connection to the central atom. If false, only clusters where all atoms are connected are generated.", has_value_target=has_compact_clusters)
      call param_register(params, 'xml_version', '1596837814', xml_version, help_string="Version of GAP the XML potential file was created")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_Nb_initialise args_str')) then
         RAISE_ERROR("distance_Nb_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( this%order < 1 ) then
         RAISE_ERROR("distance_Nb_initialise: order must be greater than 0",error)
      endif

      if (.not. has_compact_clusters) then
        ! no compact_clusters specified explicitly, default depends on version
        if (xml_version < 1596837814) then
          ! before version where default was changed from false to true
          this%compact_clusters = .false.
        else
          ! after version where default was changed from false to true
          this%compact_clusters = .true.
        endif
      endif


      allocate(this%Z(this%order))
      default_Z = ""
      do i = 1, this%order
         default_Z = trim(default_Z) // " 0"
      enddo

      call initialise(params)
      if( this%order == 1 ) then
         call param_register(params, 'Z', trim(default_Z), this%Z(1), help_string="Atomic type of neighbours")
      else
         call param_register(params, 'Z', trim(default_Z), this%Z, help_string="Atomic type of neighbours")
      endif
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_Nb_initialise args_str')) then
         RAISE_ERROR("distance_Nb_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call sort_array(this%Z)
      call distance_Nb_n_permutations(this%Z, n_permutations)
      this%n_permutations = product(factorial_int(n_permutations))

      allocate(atom_permutations(this%order,this%n_permutations))
      call distance_Nb_permutations(n_permutations,atom_permutations)

      allocate(distance_matrix_index(this%order,this%order))
      allocate(this%permutations( max(1,(this%order - 1) * this%order / 2), this%n_permutations))

      if( this%order == 1 ) then
         this%permutations = 1
      else
         k = 0
         do i = 1, this%order
            do j = i+1, this%order
               k = k + 1
               distance_matrix_index(j,i) = k
               distance_matrix_index(i,j) = k
            enddo
         enddo

         do i_p = 1, this%n_permutations
            k = 0
            do i = 1, this%order
               do j = i+1, this%order
                  k = k + 1
                  this%permutations(k,i_p) = distance_matrix_index(atom_permutations(j,i_p), atom_permutations(i,i_p))
               enddo
            enddo
         enddo
      endif

      nEdges = this%order * (this%order - 1) / 2
      allocate( edges(2,nEdges))

      k = 0
      do i = 1, this%order
         do j = i+1, this%order
            k = k + 1
            edges(:,k) = (/i,j/)
         enddo
      enddo

      nConnectivities = 2**nEdges

      allocate(allConnectivities(this%order,this%order,nConnectivities))
      allocate(connectivityList(nEdges))

      nMonomerConnectivities = 0
      do i = 1, nConnectivities
         call integerDigits(i-1,2,connectivityList)
         allConnectivities(:,:,i) = .false.
         do j = 1, nEdges
            allConnectivities(edges(1,j),edges(2,j),i) = ( connectivityList(j) == 1 )
            allConnectivities(edges(2,j),edges(1,j),i) = ( connectivityList(j) == 1 )
         enddo

         if( graphIsConnected( allConnectivities(:,:,i) ) ) nMonomerConnectivities = nMonomerConnectivities + 1
      enddo

      allocate(this%monomerConnectivities(this%order,this%order,nMonomerConnectivities))
      j = 0
      do i = 1, nConnectivities
         if( graphIsConnected( allConnectivities(:,:,i) ) ) then
            j = j + 1
            this%monomerConnectivities(:,:,j) = allConnectivities(:,:,i)
         endif
      enddo

      if(allocated(n_permutations)) deallocate(n_permutations)
      if(allocated(atom_permutations)) deallocate(atom_permutations)
      if(allocated(distance_matrix_index)) deallocate(distance_matrix_index)
      if(allocated(edges)) deallocate(edges)
      if(allocated(allConnectivities)) deallocate(allConnectivities)
      if(allocated(connectivityList)) deallocate(connectivityList)
      this%initialised = .true.

   endsubroutine distance_Nb_initialise

   subroutine distance_Nb_finalise(this,error)
      type(distance_Nb), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.5_dp
      this%order = 0
      this%n_permutations = 0
      this%compact_clusters = .false.
      if(allocated(this%Z)) deallocate(this%Z)
      if(allocated(this%permutations)) deallocate(this%permutations)
      if(allocated(this%monomerConnectivities)) deallocate(this%monomerConnectivities)

      this%initialised = .false.

   endsubroutine distance_Nb_finalise


   subroutine distance_Nb_n_permutations(Z,n_permutations,error)
      integer, dimension(:), intent(in) :: Z
      integer, dimension(:), allocatable :: n_permutations
      integer, optional, intent(out) :: error

      integer :: i
      integer, dimension(:), allocatable :: uniq_Z

      INIT_ERROR(error)

      call uniq(Z,uniq_Z)
      call reallocate(n_permutations,size(uniq_Z))

      do i = 1, size(uniq_Z)
         n_permutations(i) = count( uniq_Z(i) == Z )
      enddo

      if(allocated(uniq_Z)) deallocate(uniq_Z)

   endsubroutine distance_Nb_n_permutations

   recursive subroutine distance_Nb_permutations(n_permutations,permutations)
      integer, dimension(:), intent(in) :: n_permutations
      integer, dimension(sum(n_permutations),product(factorial_int(n_permutations))), intent(inout) :: permutations

      integer, dimension(:), allocatable, save :: current_permutation
      integer :: i, j, n_lo, n_hi
      integer, save :: recursion_level = 0, i_current_permutation = 0

      recursion_level = recursion_level + 1


      if( recursion_level == 1 ) then
         i_current_permutation = 0
         allocate(current_permutation(sum(n_permutations)))
         current_permutation = 0
      endif


      do i = 1, size(n_permutations)
         if( i == 1 ) then
            n_lo = 1
         else
            n_lo = sum(n_permutations(1:i-1)) + 1
         endif
         n_hi = sum(n_permutations(1:i))
         do j = n_lo, n_hi
            if( i_current_permutation < size(permutations,2) ) then
               if( .not. any(j==current_permutation) .and. recursion_level >= n_lo .and. recursion_level <= n_hi ) then

                  current_permutation(recursion_level) = j
                  if( recursion_level == sum(n_permutations) ) then
                     i_current_permutation = i_current_permutation + 1
                     permutations(:,i_current_permutation) = current_permutation
                  else
                     call distance_Nb_permutations(n_permutations,permutations)
                  endif
               endif
            endif
         enddo
      enddo

      current_permutation(recursion_level) = 0

      recursion_level = recursion_level - 1

      if( recursion_level == 0 ) then
         deallocate(current_permutation)
      endif

   endsubroutine distance_Nb_permutations

   subroutine soap_turbo_initialise(this,args_str,error)
      use soap_turbo_compress_module

      type(soap_turbo), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      logical :: has_atom_sigma_angular

      integer :: l, k, i, j, m, n, n_nonzero
      real(dp) :: fact, fact1, fact2, ppi, atom_sigma_radial_normalised, cutoff_hard,&
         s2, I_n, N_n, N_np1, N_np2, I_np1, I_np2, C2
      character(len=64) :: compress_string

      type(LA_Matrix) :: LA_overlap
      real(dp), dimension(:), allocatable :: s
      real(dp), dimension(:,:), allocatable :: sqrt_overlap, u, v
      real(dp), parameter :: sqrt_two = sqrt(2.0_dp)

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'l_max', PARAM_MANDATORY, this%l_max, help_string="Angular basis resolution")
      call param_register(params, 'n_species', '1', this%n_species, help_string="Number of species for the descriptor")
      call param_register(params, 'rcut_hard', PARAM_MANDATORY, this%rcut_hard, help_string="Hard cutoff")
      call param_register(params, 'rcut_soft', PARAM_MANDATORY, this%rcut_soft, help_string="Soft cutoff")
      call param_register(params, 'nf', "4.0", this%nf, help_string="TODO")
      call param_register(params, 'radial_enhancement', "0", this%radial_enhancement, help_string="TODO")
      call param_register(params, 'basis', "poly3", this%basis, help_string="poly3 or poly3gauss")
      call param_register(params, 'scaling_mode', "polynomial", this%scaling_mode, help_string="TODO")
      call param_register(params, 'compress_file', "None", this%compress_file, help_string="TODO")
      call param_register(params, 'compress_mode', "None", this%compress_mode, help_string="TODO")
      call param_register(params, 'central_index', "1", this%central_index, help_string="Index of central atom species_Z in the >species< array")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_turbo_initialise args_str')) then
         RAISE_ERROR("soap_turbo_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      allocate(this%atom_sigma_r(this%n_species))
      allocate(this%atom_sigma_r_scaling(this%n_species))
      allocate(this%atom_sigma_t(this%n_species))
      allocate(this%atom_sigma_t_scaling(this%n_species))
      allocate(this%amplitude_scaling(this%n_species))
      allocate(this%central_weight(this%n_species))
      allocate(this%alpha_max(this%n_species))
      allocate(this%species_Z(this%n_species))

      call initialise(params)
      if(this%n_species == 1) then
         call param_register(params, 'alpha_max', PARAM_MANDATORY, this%alpha_max(1), &
            help_string="Radial basis resolution for each species")
         call param_register(params, 'atom_sigma_r', PARAM_MANDATORY, this%atom_sigma_r(1), &
            help_string="Width of atomic Gaussians for soap-type descriptors in the radial direction")
         call param_register(params, 'atom_sigma_r_scaling', PARAM_MANDATORY, this%atom_sigma_r_scaling(1), &
            help_string="Scaling rate of radial sigma: scaled as a function of neighbour distance")
         call param_register(params, 'atom_sigma_t', PARAM_MANDATORY, this%atom_sigma_t(1), &
            help_string="Width of atomic Gaussians for soap-type descriptors in the angular direction")
         call param_register(params, 'atom_sigma_t_scaling', PARAM_MANDATORY, this%atom_sigma_t_scaling(1), &
            help_string="Scaling rate of angular sigma: scaled as a function of neighbour distance")
         call param_register(params, 'amplitude_scaling', PARAM_MANDATORY, this%amplitude_scaling(1), &
            help_string="Scaling rate of amplitude: scaled as an inverse function of neighbour distance")
         call param_register(params, 'central_weight', PARAM_MANDATORY, this%central_weight(1), &
            help_string="Weight of central atom in environment")
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z(1), &
            help_string="Atomic number of species, including the central atom")
      else
         call param_register(params, 'alpha_max', '//MANDATORY//', this%alpha_max, &
            help_string="Radial basis resultion for each species")
         call param_register(params, 'atom_sigma_r', '//MANDATORY//', this%atom_sigma_r, &
            help_string="Width of atomic Gaussians for soap-type descriptors in the radial direction")
         call param_register(params, 'atom_sigma_r_scaling', '//MANDATORY//', this%atom_sigma_r_scaling, &
            help_string="Scaling rate of radial sigma: scaled as a function of neighbour distance")
         call param_register(params, 'atom_sigma_t', '//MANDATORY//', this%atom_sigma_t, &
            help_string="Width of atomic Gaussians for soap-type descriptors in the angular direction")
         call param_register(params, 'atom_sigma_t_scaling', '//MANDATORY//', this%atom_sigma_t_scaling, &
            help_string="Scaling rate of angular sigma: scaled as a function of neighbour distance")
         call param_register(params, 'amplitude_scaling', '//MANDATORY//', this%amplitude_scaling, &
                     help_string="Scaling rate of amplitude: scaled as an inverse function of neighbour distance")
         call param_register(params, 'central_weight', '//MANDATORY//', this%central_weight, &
            help_string="Weight of central atom in environment")
         call param_register(params, 'species_Z', '//MANDATORY//', this%species_Z, &
            help_string="Atomic number of species, including the central atom")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_turbo_initialise args_str')) then
         RAISE_ERROR("soap_turbo_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

!     Here we read in the compression information from a file (compress_file) or rely on a keyword provided
!     by the user (compress_mode) which leads to a predefined recipe to compress the soap_turbo descriptor
!     The file always takes precedence over the keyword.
      if( this%compress_file /= "None" )then
        this%compress = .true.
        open(unit=10, file=this%compress_file, status="old")
        read(10, *) (i, j=1,this%n_species), i, n
        read(10, '(A)') compress_string
        if( compress_string == "P_transformation" )then
          n_nonzero = -1
          do while( compress_string /= "end_transformation" )
            read(10, '(A)') compress_string
            n_nonzero = n_nonzero + 1
          end do
          this%compress_P_nonzero = n_nonzero
          allocate( this%compress_P_el(1:n_nonzero) )
          allocate( this%compress_P_i(1:n_nonzero) )
          allocate( this%compress_P_j(1:n_nonzero) )
          do i = 1, n_nonzero+1
            backspace(10)
          end do
          do i = 1, n_nonzero
            read(10,*) this%compress_P_i(i), this%compress_P_j(i), this%compress_P_el(i)
          end do
        else
!         Old way to handle compression for backcompatibility
          backspace(10)
          this%compress_P_nonzero = n
          allocate( this%compress_P_el(1:n) )
          allocate( this%compress_P_i(1:n) )
          allocate( this%compress_P_j(1:n) )
          do i = 1, n
            read(10, *) this%compress_P_j(i)
            this%compress_P_i(i) = i
            this%compress_P_el(i) = 1.0_dp
          end do
        end if
        close(10)
      else if( this%compress_mode /= "None" )then
        this%compress = .true.
        call get_compress_indices( this%compress_mode, this%alpha_max, this%l_max, n, this%compress_P_nonzero, &
                                   this%compress_P_i, this%compress_P_j, this%compress_P_el, "get_dim" )
        allocate( this%compress_P_i(1:this%compress_P_nonzero) )
        allocate( this%compress_P_j(1:this%compress_P_nonzero) )
        allocate( this%compress_P_el(1:this%compress_P_nonzero) )
        call get_compress_indices( this%compress_mode, this%alpha_max, this%l_max, n, this%compress_P_nonzero, &
                                   this%compress_P_i, this%compress_P_j, this%compress_P_el, "set_indices" )
      end if


      this%initialised = .true.

   endsubroutine soap_turbo_initialise

   subroutine soap_turbo_finalise(this,error)
      type(soap_turbo), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%rcut_hard = 0.0_dp
      this%rcut_soft = 0.0_dp
      this%nf = 0.0_dp
      this%n_species = 0
      this%radial_enhancement = 0
      this%central_index = 0
      this%l_max = 0

      if(allocated(this%alpha_max)) deallocate(this%alpha_max)
      if(allocated(this%atom_sigma_r)) deallocate(this%atom_sigma_r)
      if(allocated(this%atom_sigma_r_scaling)) deallocate(this%atom_sigma_r_scaling)
      if(allocated(this%atom_sigma_t)) deallocate(this%atom_sigma_t)
      if(allocated(this%atom_sigma_t_scaling)) deallocate(this%atom_sigma_t_scaling)
      if(allocated(this%amplitude_scaling)) deallocate(this%amplitude_scaling)
      if(allocated(this%central_weight)) deallocate(this%central_weight)
      if(allocated(this%species_Z)) deallocate(this%species_Z)

      this%initialised = .false.

   endsubroutine soap_turbo_finalise

   subroutine soap_turbo_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(soap_turbo), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_turbo_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%species_Z(this%central_index) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%rcut_hard) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine soap_turbo_sizes

   subroutine descriptor_str_add_species(this,species,descriptor_str,error)
      character(len=*), intent(in) :: this
      integer, dimension(:), intent(in) :: species
      character(len=STRING_LENGTH), dimension(:), allocatable, intent(out) :: descriptor_str
      integer, optional, intent(out) :: error

      integer :: my_descriptor_type, i, j, k, l, n_species, order, n
      integer, dimension(:,:), allocatable :: ZN
      real(dp), dimension(:), allocatable :: w
      type(Dictionary) :: params

      INIT_ERROR(error)

      if(allocated(descriptor_str)) deallocate(descriptor_str)

      my_descriptor_type = get_descriptor_type(this,error)
      n_species = size(species)

      select case(my_descriptor_type)
      case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_COSNX,DT_POWER_SO3,DT_POWER_SO4)
         allocate(w(n_species))
         allocate(descriptor_str(n_species))

         if( n_species == 1 ) then
            w = 1.0_dp
         else
            w = real( (/ (i, i=0, n_species-1) /), kind=dp ) / (n_species-1) * 0.5_dp + 0.5_dp
         endif

         do i = 1, n_species
            descriptor_str(i) = trim(this)//" n_species="//n_species//" Z="//species(i)//" species_Z={"//species//"} w={"//w//"}"
         enddo

         deallocate(w)
      case(DT_SOAP)
         allocate(descriptor_str(n_species))
         do i = 1, n_species
            descriptor_str(i) = trim(this)//" n_species="//n_species//" Z="//species(i)//" species_Z={"//species//"}"
         enddo
      case(DT_DISTANCE_2B,DT_CO_DISTANCE_2B,DT_AS_DISTANCE_2B)
         allocate(descriptor_str(n_species * (n_species+1) / 2))

         l = 0
         do i = 1, n_species
            do j = i, n_species
               l = l + 1
               descriptor_str(l) = trim(this)//" Z1="//species(i)//" Z2="//species(j)
            enddo
         enddo

      case(DT_COORDINATION,DT_RDF)
         allocate(descriptor_str(n_species))
         do i = 1, n_species
            descriptor_str(i) = trim(this)//" Z="//species(i)
         enddo
      case(DT_ANGLE_3B,DT_CO_ANGLE_3B)
         allocate(descriptor_str(n_species * n_species * (n_species+1) / 2))
         l = 0
         do i = 1, n_species
            do j = 1, n_species
               do k = j, n_species
                  l = l + 1
                  descriptor_str(l) = trim(this)//" Z="//species(i)//" Z1="//species(j)//" Z2="//species(k)
               enddo
            enddo
         enddo
      case(DT_GENERAL_MONOMER,DT_GENERAL_DIMER,DT_WATER_MONOMER,DT_WATER_DIMER,DT_A2_DIMER,DT_AB_DIMER,DT_TRIHIS,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_AN_MONOMER)
         allocate(descriptor_str(1))
         descriptor_str(1) = trim(this)
      case(DT_DISTANCE_NB)
         call initialise(params)
         call param_register(params, 'order', PARAM_MANDATORY, order, help_string="Many-body order, in terms of number of neighbours")
         if (.not. param_read_line(params, this, ignore_unknown=.true.,task='descriptor_str_add_species this')) then
            RAISE_ERROR("descriptor_str_add_species failed to parse descriptor string='"//trim(this)//"'", error)
         endif
         call finalise(params)

         n = 1
         do i = 1, order
            n = n * ( n_species + i - 1 ) / i ! avoids double counting
         enddo

         allocate(ZN(order,n),descriptor_str(n))

         call descriptor_str_add_species_distance_Nb(ZN,species,order)

         do i = 1, n
            descriptor_str(i) = trim(this)//" Z={"//ZN(:,i)//"}"
         enddo
         deallocate(ZN)

      case(DT_SOAP_EXPRESS)
         RAISE_ERROR("descriptor_str_add_species: no recipe for "//my_descriptor_type//" yet.",error)
      case(DT_SOAP_TURBO)
         RAISE_ERROR("descriptor_str_add_species: no recipe for "//my_descriptor_type//" yet.",error)
      case default
         RAISE_ERROR("descriptor_str_add_species: unknown descriptor type "//my_descriptor_type,error)
      endselect

   endsubroutine descriptor_str_add_species

   recursive subroutine descriptor_str_add_species_distance_Nb(ZN,species,order)
      integer, dimension(:,:), intent(inout) :: ZN
      integer, dimension(:), intent(in) :: species
      integer, intent(in) :: order

      integer :: i_species, n_species
      integer, save :: current_descriptor, current_order = 0
      integer, dimension(:), allocatable, save :: ZN_current

      n_species = size(species)

      if( current_order == 0 ) then                           ! first run, outermost order.
         current_descriptor = 0                               ! keeps track of descriptor
         current_order = 1                                    ! keeps track of order
         allocate(ZN_current(order))                          ! builds/updates atomic numbers gradually for each descriptor
      endif

      do i_species = 1, n_species
         if( current_order > 1 ) then                             ! no special atom, all atoms equivalent
            if( species(i_species) < ZN_current(current_order-1) ) cycle   ! avoids double-counting of neighbours
         endif

         ZN_current(current_order) = species(i_species)
         if( current_order < order ) then                                 ! calls recursively until we reach the last order
            current_order = current_order + 1
            call descriptor_str_add_species_distance_Nb(ZN,species,order)
         else                                                             ! when we reached the last order, fill the atomic numbers in the loop
            current_descriptor = current_descriptor + 1                   ! and add them to the output array
            ZN(:,current_descriptor) = ZN_current
         endif
      enddo

      current_order = current_order - 1                                   ! when the loop finished, step one level down

      if( current_order == 0 ) deallocate(ZN_current)                     ! when we reach zero, we finished.

   endsubroutine descriptor_str_add_species_distance_Nb

   subroutine descriptor_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(cutoff(this) > at%cutoff) then
         RAISE_ERROR("descriptor_calc: descriptor cutoff ("//cutoff(this)//") larger than atoms cutoff("//at%cutoff//")",error)
      endif

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call calc(this%descriptor_bispectrum_SO4,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_BISPECTRUM_SO3)
            call calc(this%descriptor_bispectrum_SO3,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error=error)
         case(DT_BEHLER)
            call calc(this%descriptor_behler,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error=error)
         case(DT_DISTANCE_2b)
            call calc(this%descriptor_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COORDINATION)
            call calc(this%descriptor_coordination,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ANGLE_3B)
            call calc(this%descriptor_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_CO_ANGLE_3B)
            call calc(this%descriptor_co_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_CO_DISTANCE_2b)
            call calc(this%descriptor_co_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COSNX)
            call calc(this%descriptor_cosnx,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_TRIHIS)
            call calc(this%descriptor_trihis,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_WATER_MONOMER)
            call calc(this%descriptor_water_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_WATER_DIMER)
            call calc(this%descriptor_water_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_A2_DIMER)
            call calc(this%descriptor_A2_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AB_DIMER)
            call calc(this%descriptor_AB_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ATOM_REAL_SPACE)
            call calc(this%descriptor_atom_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_POWER_SO3)
            call calc(this%descriptor_power_so3,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_POWER_SO4)
            call calc(this%descriptor_power_so4,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_SOAP)
            call calc(this%descriptor_soap,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_RDF)
            call calc(this%descriptor_rdf,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AS_DISTANCE_2b)
            call calc(this%descriptor_as_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_ALEX)
            call calc(this%descriptor_alex,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_DISTANCE_Nb)
            call calc(this%descriptor_distance_Nb,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_SOAP_TURBO)
            call calc(this%descriptor_soap_turbo,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_BOND_REAL_SPACE)
            call calc(this%descriptor_bond_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_AN_MONOMER)
            call calc(this%descriptor_AN_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_COM_DIMER)
            call calc(this%descriptor_com_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_MONOMER)
            call calc(this%descriptor_general_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_DIMER)
            call calc(this%descriptor_general_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_GENERAL_TRIMER)
            call calc(this%descriptor_general_trimer,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_MOLECULE_LO_D)
            call calc(this%descriptor_molecule_lo_d,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
         case(DT_SOAP_EXPRESS)
            call calc(this%descriptor_soap_express,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
#endif
         case default
            RAISE_ERROR("descriptor_calc: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_calc

   subroutine descriptor_calc_array(this,at,descriptor_out,covariance_cutoff,descriptor_index, &
         grad_descriptor_out,grad_descriptor_index,grad_descriptor_pos,grad_covariance_cutoff,args_str,error)

      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      real(dp), dimension(:,:), intent(out), optional :: descriptor_out
      real(dp), dimension(:), intent(out), optional :: covariance_cutoff
      integer, dimension(:,:), intent(out), optional :: descriptor_index
      real(dp), dimension(:,:,:), intent(out), optional :: grad_descriptor_out
      integer, dimension(:,:), intent(out), optional :: grad_descriptor_index
      real(dp), dimension(:,:), intent(out), optional :: grad_descriptor_pos
      real(dp), dimension(:,:), intent(out), optional :: grad_covariance_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(descriptor_data) :: my_descriptor_data
      type(Dictionary) :: params
      integer :: i, n, i_d, n_descriptors, n_cross, n_index
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer
      logical :: do_grad_descriptor, do_descriptor

      INIT_ERROR(error)

      do_descriptor = present(descriptor_out)
      do_grad_descriptor = present(grad_descriptor_out) .or. present(grad_descriptor_index) .or. present(grad_descriptor_pos)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='descriptor_calc_array args_str')) then
            RAISE_ERROR("descriptor_calc_array failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("descriptor_calc_array did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      if (associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      call calc(this,at,my_descriptor_data,do_descriptor=do_descriptor,do_grad_descriptor=do_grad_descriptor,args_str=args_str,error=error)

      if(present(descriptor_out)) &
         call check_size('descriptor_out',descriptor_out, (/descriptor_dimensions(this),n_descriptors/),'descriptor_calc_array',error)

      if(present(covariance_cutoff)) &
         call check_size('covariance_cutoff',covariance_cutoff,(/n_descriptors/),'descriptor_calc_array',error)

      if(present(descriptor_index)) &
         call check_size('descriptor_index',descriptor_index,(/n_index,n_descriptors/),'descriptor_calc_array',error)

      if(present(grad_descriptor_out)) &
         call check_size('grad_descriptor_out',grad_descriptor_out,(/descriptor_dimensions(this),3,n_cross/),'descriptor_calc_array',error)

      if(present(grad_descriptor_index)) &
         call check_size('grad_descriptor_index',grad_descriptor_index,(/2,n_cross/),'descriptor_calc_array',error)

      if(present(grad_descriptor_pos)) &
         call check_size('grad_descriptor_pos',grad_descriptor_pos,(/3,n_cross/),'descriptor_calc_array',error)

      if(present(grad_covariance_cutoff)) &
         call check_size('grad_covariance_cutoff',grad_covariance_cutoff,(/3,n_cross/),'descriptor_calc_array',error)

      if(do_descriptor) then
         do i = 1, n_descriptors
            descriptor_out(:,i) = my_descriptor_data%x(i)%data
            if(present(covariance_cutoff)) covariance_cutoff(i) = my_descriptor_data%x(i)%covariance_cutoff
            if(present(descriptor_index)) descriptor_index(:,i) = my_descriptor_data%x(i)%ci
         enddo
      endif

      if(do_grad_descriptor) then
         i_d = 0
         do i = 1, n_descriptors
            do n = lbound(my_descriptor_data%x(i)%ii,1),ubound(my_descriptor_data%x(i)%ii,1)
               i_d = i_d + 1
               if(present(grad_descriptor_index)) grad_descriptor_index(:,i_d) = (/i,my_descriptor_data%x(i)%ii(n)/)
               if(present(grad_descriptor_out)) grad_descriptor_out(:,:,i_d) = my_descriptor_data%x(i)%grad_data(:,:,n)
               if(present(grad_descriptor_pos)) grad_descriptor_pos(:,i_d) = my_descriptor_data%x(i)%pos(:,n)
               if(present(grad_covariance_cutoff)) grad_covariance_cutoff(:,i_d) = my_descriptor_data%x(i)%grad_covariance_cutoff(:,n)
            enddo
         enddo
      endif

      call finalise(my_descriptor_data,error=error)

   endsubroutine descriptor_calc_array

   subroutine bispectrum_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      complex(dp) :: sub
      complex(dp), dimension(3) :: dsub
      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r, tmp_cg
      integer :: i, n, n_i, ji, jn, j, m1, m2, j1, j2, m11, m12, m21, m22, &
         i_desc, i_bisp, d, n_descriptors, n_cross, l_n_neighbours, n_index
      integer, dimension(3) :: shift
      integer, dimension(total_elements) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('bispectrum_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bispectrum_SO4_calc args_str')) then
            RAISE_ERROR("bispectrum_SO4_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bispectrum_SO4_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call cg_initialise(this%j_max, 2)

      call finalise(descriptor_out)

      d = bispectrum_SO4_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            ! dU is not allocated, allocate and zero it
            allocate( dU(0:this%j_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do j = 0, this%j_max
               allocate( dU(j,0)%mm(3,-j:j,-j:j) )
               dU(j,0)%mm = CPLX_ZERO
            enddo

            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            ji = neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,args_str,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,args_str=args_str,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j1 = 0, this%j_max
               j2 = j1
               !do j2 = 0, this%j_max
                  do j = abs(j1-j2), min(this%j_max,j1+j2)
                     if( mod(j1+j2+j,2) == 1 ) cycle

                     i_bisp = i_bisp + 1

                     !do m1 = -j, j, 2
                     !   do m2 = -j, j, 2
                     !      sub = CPLX_ZERO
                     !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                     !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                     !            sub = sub + cg_array(j1,m11,j,m1,j1,m11+m1) &
                     !            * cg_array(j2,m21,j,m2,j2,m21+m2) &
                     !            * U(j1)%mm(m11,m11+m1) * U(j2)%mm(m21,m21+m2)
                     !         enddo
                     !      enddo
                     !      descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(-m2,m1))*(-1)**(m2/2)
                     !   enddo
                     !enddo

                     do m1 = -j, j, 2
                        do m2 = -j, j, 2
                           sub = CPLX_ZERO
                           do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                              do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                                 sub = sub + cg_array(j1,m11,j2,m1-m11,j,m1) &
                                 * cg_array(j1,m12,j2,m2-m12,j,m2) &
                                 * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                              enddo
                           enddo
                           descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(m1,m2))
                        enddo
                     enddo

                  enddo
               !enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 0, n_neighbours(at,i)
               if( n>0 ) then
                  ji = neighbour(at, i, n, distance=r)
                  if( r >= this%cutoff ) cycle
                  n_i = n_i + 1
               endif
               i_bisp = 0
               do j1 = 0, this%j_max
                  j2 = j1
                  !do j2 = 0, this%j_max
                     do j = abs(j1-j2), min(this%j_max,j1+j2)
                        if( mod(j1+j2+j,2) == 1 ) cycle

                        i_bisp = i_bisp + 1

                        !do m1 = -j, j, 2
                        !   do m2 = -j, j, 2
                        !      sub = CPLX_ZERO
                        !      dsub = CPLX_ZERO

                        !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                        !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                        !            tmp_cg =  cg_array(j1,m11,j,m1,j1,m11+m1) &
                        !              * cg_array(j2,m21,j,m2,j2,m21+m2)

                        !            sub = sub + tmp_cg &
                        !            * U(j1)%mm(m11,m1+m11) * U(j2)%mm(m21,m2+m21)
                        !            dsub = dsub + tmp_cg &
                        !            * ( dU(j1,n_i)%mm(:,m11,m1+m11) * U(j2)%mm(m21,m2+m21) + &
                        !            U(j1)%mm(m11,m1+m11) * dU(j2,n_i)%mm(:,m21,m2+m21) )
                        !         enddo
                        !      enddo
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                        !      ( dsub*conjg(U(j)%mm(-m2,m1)) + sub*conjg(dU(j,n_i)%mm(:,-m2,m1)) )*(-1)**(m2/2)
                        !   enddo
                        !enddo
                        do m1 = -j, j, 2
                           do m2 = -j, j, 2
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                                 do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2

                                    tmp_cg =  cg_array(j1,m11,j2,m1-m11,j,m1) &
                                    * cg_array(j1,m12,j2,m2-m12,j,m2)

                                    sub = sub + tmp_cg &
                                    * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                                    dsub = dsub + tmp_cg &
                                    * ( dU(j1,n_i)%mm(:,m11,m12) * U(j2)%mm(m1-m11,m2-m12) + &
                                    U(j1)%mm(m11,m12) * dU(j2,n_i)%mm(:,m1-m11,m2-m12) )
                                 enddo
                              enddo
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                              dsub*conjg(U(j)%mm(m1,m2)) + sub*conjg(dU(j,n_i)%mm(:,m1,m2))
                           enddo
                        enddo

                     enddo
                  !enddo
               enddo
            enddo
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('bispectrum_SO4_calc')

   endsubroutine bispectrum_SO4_calc

   subroutine bispectrum_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bispectrum_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, l1, l2, m1, i_desc, i_pow, l_n_neighbours, &
         n_i, n_descriptors, n_cross, n_index
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij

      complex(dp) :: sub, dsub(3)

      integer, dimension(total_elements) :: species_map

      INIT_ERROR(error)

      call system_timer('bispectrum_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call cg_initialise(this%l_max)

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bispectrum_SO3_calc args_str')) then
            RAISE_ERROR("bispectrum_SO3_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bispectrum_SO3_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("bispectrum_SO3_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = bispectrum_so3_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif


      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l1 = 0, this%l_max
                  l2 = l1
                  !do l2 = 0, this%l_max
                     do l = abs(l1-l2), min(this%l_max,l1+l2)
                        if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                        i_pow = i_pow + 1

                        do m = -l, l
                           sub = CPLX_ZERO
                           do m1 = max(-l1,m-l2),min(l1,m+l2)
                              sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                           enddo

                           descriptor_out%x(i_desc)%data(i_pow) = descriptor_out%x(i_desc)%data(i_pow) + fourier_so3(l,a)%m(m) * sub
                        enddo

                     enddo
                  !enddo
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l1 = 0, this%l_max
                     l2 = l1
                     !do l2 = 0, this%l_max
                        do l = abs(l1-l2), min(this%l_max,l1+l2)
                           if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                           i_pow = i_pow + 1

                           do m = -l, l
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m1 = max(-l1,m-l2),min(l1,m+l2)
                                 dsub = dsub + cg_array(l1,m1,l2,m-m1,l,m) * &
                                 ( conjg(dfourier_so3(l1,a,n)%mm(:,m1)) * conjg(fourier_so3(l2,a)%m(m-m1)) + &
                                   conjg(fourier_so3(l1,a)%m(m1)) * conjg(dfourier_so3(l2,a,n)%mm(:,m-m1)) )
                                 sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                              enddo

                              descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = descriptor_out%x(i_desc)%grad_data(i_pow,:,n) + &
                              fourier_so3(l,a)%m(m) * dsub + dfourier_so3(l,a,n)%mm(:,m) * sub
                           enddo
                        enddo
                     !enddo
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('bispectrum_so3_calc')

   endsubroutine bispectrum_so3_calc

   subroutine behler_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc_i, l_n_neighbours, &
         n_i, m_i, n_descriptors, n_cross, n_index
      integer, dimension(:), allocatable :: i_desc
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Ang, dAng, Rad, dRad_ij, dRad_ik, dRad_jk, f_cut_ij, f_cut_ik, f_cut_jk, df_cut_ij, df_cut_ik, df_cut_jk, g2, dg2
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('behler_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("behler_calc: descriptor object not initialised", error)
      endif

      if( at%cutoff < this%cutoff ) then
         RAISE_ERROR("behler_calc: cutoff of atoms object ("//at%cutoff//") less than cutoff of descriptor ("//this%cutoff//")", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='behler_calc args_str')) then
            RAISE_ERROR("behler_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("behler_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = behler_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      allocate(i_desc(at%N))

      i_desc = 0
      i_desc_i = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         if( this%Z /= 0 .and. this%Z /= at%Z(i) ) cycle
         i_desc_i = i_desc_i + 1
         i_desc(i) = i_desc_i

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc_i)%data(d))
            descriptor_out%x(i_desc_i)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc_i)%ci(n_index))
            descriptor_out%x(i_desc_i)%has_data = .false.
            descriptor_out%x(i_desc_i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc_i)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc_i)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc_i)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc_i)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc_i)%grad_data = 0.0_dp
            descriptor_out%x(i_desc_i)%ii = 0
            descriptor_out%x(i_desc_i)%pos = 0.0_dp
            descriptor_out%x(i_desc_i)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc_i)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc_i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

!$omp parallel do schedule(dynamic) default(none) shared(this,at,descriptor_out,my_do_descriptor, my_do_grad_descriptor, d, i_desc) &
!$omp private(i,j,k,i_desc_i,n_i,n,r_ij,u_ij,d_ij,shift_ij,f_cut_ij,df_cut_ij,g2,dg2,m_i,m,r_ik,u_ik,d_ik,d_jk,r_jk,u_jk,cos_ijk) &
!$omp private(dcosijk_ij,dcosijk_ik,a,b,f_cut_ik,f_cut_jk,df_cut_ik,df_cut_jk,Ang,Rad,dAng,dRad_ij,dRad_ik,dRad_jk)
      do i = 1, at%N
         if( this%Z /= 0 .and. this%Z /= at%Z(i) ) cycle

         if(i_desc(i) == 0) then
            cycle
         else
            i_desc_i = i_desc(i)
         endif

         if(my_do_descriptor) then
            descriptor_out%x(i_desc_i)%ci(1) = i
            descriptor_out%x(i_desc_i)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc_i)%ii(0) = i
            descriptor_out%x(i_desc_i)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc_i)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(n_i) = j
               descriptor_out%x(i_desc_i)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc_i)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_g2
               if ( r_ij >= this%g2(a)%rc .or. ( this%g2(a)%Z_n /=0 .and. this%g2(a)%Z_n /= at%Z(j) ) ) cycle

               f_cut_ij = cos_cutoff_function(r_ij,this%g2(a)%rc)
               if(my_do_grad_descriptor) df_cut_ij = dcos_cutoff_function(r_ij,this%g2(a)%rc)

               g2 = exp(-this%g2(a)%eta * (r_ij-this%g2(a)%rs)**2)
               if(my_do_descriptor) descriptor_out%x(i_desc_i)%data(a) = descriptor_out%x(i_desc_i)%data(a) + g2 * f_cut_ij
               if(my_do_grad_descriptor) then
                  dg2 = -2.0_dp * this%g2(a)%eta * (r_ij-this%g2(a)%rs) * g2
                  descriptor_out%x(i_desc_i)%grad_data(a,:,n_i) = ( dg2 * f_cut_ij + g2 * df_cut_ij ) * u_ij
                  descriptor_out%x(i_desc_i)%grad_data(a,:,0) = descriptor_out%x(i_desc_i)%grad_data(a,:,0) - descriptor_out%x(i_desc_i)%grad_data(a,:,n_i)
               endif
            enddo


            m_i = 0
            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               m_i = m_i + 1

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               u_jk = d_jk / r_jk

               cos_ijk = dot_product(u_ij,u_ik)

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do b = 1, this%n_g3
                  if( r_ik >= this%g3(b)%rc .or. r_jk >= this%g3(b)%rc) cycle
                  if( this%g3(b)%Z_n(1) /= 0 .and. this%g3(b)%Z_n(2) /= 0 ) then
                     if( .not. ( &
                        ( this%g3(b)%Z_n(1) == at%Z(j) .and. this%g3(b)%Z_n(2) == at%Z(k) ) .or. &
                        ( this%g3(b)%Z_n(1) == at%Z(k) .and. this%g3(b)%Z_n(2) == at%Z(j) ) ) ) cycle
                  endif

                  f_cut_ij = cos_cutoff_function(r_ij,this%g3(b)%rc)
                  f_cut_ik = cos_cutoff_function(r_ik,this%g3(b)%rc)
                  f_cut_jk = cos_cutoff_function(r_jk,this%g3(b)%rc)
                  if(my_do_grad_descriptor) then
                     df_cut_ij = dcos_cutoff_function(r_ij,this%g3(b)%rc)
                     df_cut_ik = dcos_cutoff_function(r_ik,this%g3(b)%rc)
                     df_cut_jk = dcos_cutoff_function(r_jk,this%g3(b)%rc)
                  endif

                  a = b + this%n_g2

                  Ang = (1.0_dp + this%g3(b)%lambda * cos_ijk)**this%g3(b)%zeta
                  Rad = exp( -this%g3(b)%eta * (r_ij**2 + r_ik**2 + r_jk**2) )
                  if(my_do_descriptor) descriptor_out%x(i_desc_i)%data(a) = descriptor_out%x(i_desc_i)%data(a) + 0.5_dp * Ang * Rad * f_cut_ij * f_cut_ik * f_cut_jk
                  if(my_do_grad_descriptor) then
                     dAng = this%g3(b)%zeta * (1.0_dp + this%g3(b)%lambda * cos_ijk)**(this%g3(b)%zeta -1.0_dp) * this%g3(b)%lambda
                     dRad_ij = -this%g3(b)%eta * 2.0_dp * r_ij * Rad
                     dRad_ik = -this%g3(b)%eta * 2.0_dp * r_ik * Rad
                     dRad_jk = -this%g3(b)%eta * 2.0_dp * r_jk * Rad

                     descriptor_out%x(i_desc_i)%grad_data(a,:,n_i) = descriptor_out%x(i_desc_i)%grad_data(a,:,n_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ij * Rad + Ang * ( dRad_ij * u_ij - dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ik * ( df_cut_ij * u_ij * f_cut_jk - f_cut_ij * df_cut_jk * u_jk ) )

                     descriptor_out%x(i_desc_i)%grad_data(a,:,m_i) = descriptor_out%x(i_desc_i)%grad_data(a,:,m_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ik * Rad + Ang * ( dRad_ik * u_ik + dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ij * ( df_cut_ik * u_ik * f_cut_jk + f_cut_ik * df_cut_jk * u_jk ) )

                     descriptor_out%x(i_desc_i)%grad_data(a,:,0) = descriptor_out%x(i_desc_i)%grad_data(a,:,0) - 0.5_dp * &
                     ( ( dAng * (dcosijk_ij+dcosijk_ik) * Rad + Ang * (dRad_ij * u_ij + dRad_ik * u_ik) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_jk * ( df_cut_ij * u_ij * f_cut_ik + f_cut_ij * df_cut_ik * u_ik ) )
                  endif


               enddo

            enddo
         enddo

         do b = 1, this%n_g3
            a = b + this%n_g2

            if(my_do_descriptor) descriptor_out%x(i_desc_i)%data(a) = descriptor_out%x(i_desc_i)%data(a) * 2.0_dp**(1.0_dp-this%g3(b)%zeta)
            if(my_do_grad_descriptor) descriptor_out%x(i_desc_i)%grad_data(a,:,:) = descriptor_out%x(i_desc_i)%grad_data(a,:,:) * 2.0_dp**(1.0_dp-this%g3(b)%zeta)
         enddo
      enddo
!$omp end parallel do

      if(allocated(i_desc)) deallocate(i_desc)

      call system_timer('behler_calc')

   endsubroutine behler_calc

   subroutine distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical :: needs_resid
      logical, dimension(:), pointer :: atom_mask_pointer
      integer, dimension(:), pointer :: resid_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n, n_index
      integer, dimension(3) :: shift
      real(dp) :: r_ij, covariance_cutoff, dcovariance_cutoff, tail, dtail
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)

      call system_timer('distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='distance_2b_calc args_str')) then
            RAISE_ERROR("distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      needs_resid = this%only_intra .or. this%only_inter
      if (needs_resid) then
         if (.not. assign_pointer(at, trim(this%resid_name), resid_pointer)) then
            RAISE_ERROR("distance_2b_calc did not find "//trim(this%resid_name)//" property (residue id) in the atoms object.", error)
         end if
      else
         resid_pointer => null()
      end if

      d = distance_2b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:1))
            allocate(descriptor_out%x(i)%ii(0:1))
            allocate(descriptor_out%x(i)%pos(3,0:1))
            allocate(descriptor_out%x(i)%has_grad_data(0:1))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:1))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if(associated(atom_mask_pointer)) then  ! skip if masked
            if(.not. atom_mask_pointer(i)) cycle ! skip if masked
         endif                                   ! skip if masked

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            if (needs_resid) then
               if (this%only_intra .and. resid_pointer(i) /= resid_pointer(j)) cycle
               if (this%only_inter .and. resid_pointer(i) == resid_pointer(j)) cycle
            end if

            i_desc = i_desc + 1

            covariance_cutoff = coordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
            if( this%has_tail .and. this%tail_exponent /= 0 ) then
               tail = ( erf(this%tail_range*r_ij) / r_ij )**this%tail_exponent
            else
               tail = 1.0_dp
            endif

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(:) = r_ij**this%exponents
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%covariance_cutoff = covariance_cutoff * tail
            endif
            if(my_do_grad_descriptor) then
               dcovariance_cutoff = dcoordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
               if( this%has_tail .and. this%tail_exponent /= 0 ) then
                  dtail = tail * this%tail_exponent * ( 2.0_dp*this%tail_range*exp(-this%tail_range**2*r_ij**2) / &
                     sqrt(pi) / erf(this%tail_range*r_ij) - 1.0_dp / r_ij )
               else
                  dtail = 0.0_dp
               endif

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(:,:,0) = -( this%exponents*r_ij**(this%exponents-1) ) .outer. u_ij
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -(dcovariance_cutoff*tail + covariance_cutoff*dtail)*u_ij

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(:,:,1) = - descriptor_out%x(i_desc)%grad_data(:,:,0)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

            endif
         enddo
      enddo

      call system_timer('distance_2b_calc')

   endsubroutine distance_2b_calc

   subroutine coordination_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, i_n, l_n_neighbours, i_desc, n_descriptors, n_cross, n_index
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, df_cut

      INIT_ERROR(error)

      call system_timer('coordination_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='coordination_calc args_str')) then
            RAISE_ERROR("coordination_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("coordination_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = coordination_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            if(my_do_descriptor) &
               descriptor_out%x(i_desc)%data(1) = descriptor_out%x(i_desc)%data(1) + coordination_function(r_ij,this%cutoff,this%transition_width)

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff,this%transition_width) * u_ij

               descriptor_out%x(i_desc)%grad_data(1,:,0) = descriptor_out%x(i_desc)%grad_data(1,:,0) - df_cut

               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,i_n) = df_cut
            endif
         enddo
      enddo

      call system_timer('coordination_calc')

   endsubroutine coordination_calc

   subroutine angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, n_index
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='angle_3b_calc args_str')) then
            RAISE_ERROR("angle_3b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("angle_3b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = angle_3b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
         endif

         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:2))
            allocate(descriptor_out%x(i)%ii(0:2))
            allocate(descriptor_out%x(i)%pos(3,0:2))
            allocate(descriptor_out%x(i)%has_grad_data(0:2))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:2))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         if(associated(atom_mask_pointer)) then   ! skip if masked
            if(.not. atom_mask_pointer(i)) cycle  ! skip if masked
         endif                                    ! skip if masked

         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
            dfc_j = dcoordination_function(r_ij,this%cutoff,this%cutoff_transition_width)

            do m = 1, n_neighbours(at,i)

               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,this%cutoff_transition_width)
               dfc_k = dcoordination_function(r_ik,this%cutoff,this%cutoff_transition_width)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%ci(1) = i
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik
               endif
            enddo
         enddo
      enddo

      call system_timer('angle_3b_calc')

   endsubroutine angle_3b_calc

   subroutine co_angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, &
         l_n_neighbours_coordination, n_index
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('co_angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='co_angle_3b_calc args_str')) then
            RAISE_ERROR("co_angle_3b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("co_angle_3b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("co_angle_3b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = co_angle_3b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         l_n_neighbours_coordination = n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)

               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               i_desc = i_desc + 1
               if(my_do_descriptor) then
                  allocate(descriptor_out%x(i_desc)%data(d))
                  descriptor_out%x(i_desc)%data = 0.0_dp
                  allocate(descriptor_out%x(i_desc)%ci(n_index))
                  descriptor_out%x(i_desc)%has_data = .false.
               endif

               if(my_do_grad_descriptor) then

                  allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%ii(0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%pos(3,0:2+l_n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%has_grad_data(0:2+l_n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_data = 0.0_dp
                  descriptor_out%x(i_desc)%ii = 0
                  descriptor_out%x(i_desc)%pos = 0.0_dp
                  descriptor_out%x(i_desc)%has_grad_data = .false.

                  allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:2+l_n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
               endif
            enddo
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' coordination_transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,do_descriptor,do_grad_descriptor,args_str,error)

      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,0.5_dp)
            dfc_j = dcoordination_function(r_ij,this%cutoff,0.5_dp)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,0.5_dp)
               dfc_k = dcoordination_function(r_ik,this%cutoff,0.5_dp)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%data(4) = descriptor_coordination%x(i)%data(1)
                  descriptor_out%x(i_desc)%ci(1) = i
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik
                  descriptor_out%x(i_desc)%grad_data(4,:,0) = descriptor_coordination%x(i)%grad_data(1,:,0)

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(3:) = descriptor_coordination%x(i)%ii(1:)
                  descriptor_out%x(i_desc)%pos(:,3:) = descriptor_coordination%x(i)%pos(:,1:)
                  descriptor_out%x(i_desc)%has_grad_data(3:) = descriptor_coordination%x(i)%has_grad_data(1:)
                  descriptor_out%x(i_desc)%grad_data(4,:,3:) = descriptor_coordination%x(i)%grad_data(1,:,1:)
               endif
            enddo
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_angle_3b_calc')

   endsubroutine co_angle_3b_calc

   subroutine co_distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n, &
         n_neighbours_coordination_i, n_neighbours_coordination_ij, n_index
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)
      call system_timer('co_distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='co_distance_2b_calc args_str')) then
            RAISE_ERROR("co_distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("co_distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("co_distance_2b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = co_distance_2b_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N

         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij)

            if(r_ij >= this%cutoff) cycle
!if(r_ij <3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(d))
               descriptor_out%x(i_desc)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc)%ci(n_index))
               descriptor_out%x(i_desc)%has_data = .false.
            endif

            if(my_do_grad_descriptor) then
               n_neighbours_coordination_ij = n_neighbours(at,i,max_dist=this%coordination_cutoff) + &
               n_neighbours(at,j,max_dist=this%coordination_cutoff) + 2

               allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%ii(0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%pos(3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%has_grad_data(0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_data = 0.0_dp
               descriptor_out%x(i_desc)%ii = 0
               descriptor_out%x(i_desc)%pos = 0.0_dp
               descriptor_out%x(i_desc)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
            endif
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,.true.,do_grad_descriptor,args_str,error)

      i_desc = 0
      do i = 1, at%N

         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij >= this%cutoff ) cycle
!if(r_ij <3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)

               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%data(2) = descriptor_coordination%x(i)%data(1) + descriptor_coordination%x(j)%data(1)
               descriptor_out%x(i_desc)%data(3) = (descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))**2

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij, this%cutoff,this%transition_width) !coordination_function(r_ij,3.5_dp, 0.5_dp, this%cutoff,this%transition_width)
            endif
            if(my_do_grad_descriptor) then
               n_neighbours_coordination_i = n_neighbours(at,i,max_dist=this%coordination_cutoff)

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%cutoff,this%transition_width)*u_ij !-dcoordination_function(r_ij,3.5_dp, 0.5_dp, this%cutoff,this%transition_width)*u_ij

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

               descriptor_out%x(i_desc)%ii(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%ii(:)
               descriptor_out%x(i_desc)%pos(:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,2:n_neighbours_coordination_i+2) = 2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(i)%grad_data(1,:,:)

               descriptor_out%x(i_desc)%ii(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%ii(:)
               descriptor_out%x(i_desc)%pos(:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,n_neighbours_coordination_i+3:) = -2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(j)%grad_data(1,:,:)

            endif
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_distance_2b_calc')

   endsubroutine co_distance_2b_calc

   subroutine cosnx_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc, i_cosnx, l_n_neighbours, n_i, &
         n_descriptors, n_cross, n_index
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, T_0_cos_ijk, T_1_cos_ijk, T_n_cos_ijk, U_0_cos_ijk, U_1_cos_ijk, U_n_cos_ijk, Ang
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, dAng_ij, dAng_ik
      real(dp), dimension(:), allocatable :: Rad_ij, Rad_ik, T_cos_ijk, U_cos_ijk
      real(dp), dimension(:,:), allocatable :: dRad_ij, dRad_ik
      integer, dimension(total_elements) :: species_map

      INIT_ERROR(error)

      call system_timer('cosnx_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='cosnx_calc args_str')) then
            RAISE_ERROR("cosnx_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("cosnx_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      print*, "jpd47", this%species_Z, size(this%species_Z)
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = cosnx_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(Rad_ij(this%n_max), Rad_ik(this%n_max))
      allocate(T_cos_ijk(0:this%l_max))
      if(my_do_grad_descriptor) then
         allocate(U_cos_ijk(-1:this%l_max))
         allocate(dRad_ij(3,this%n_max), dRad_ik(3,this%n_max))
      endif

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a) * this%w(species_map(at%Z(j)))
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij * this%w(species_map(at%Z(j)))
            enddo

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do a = 1, this%n_max
                  Rad_ik(a) = RadialFunction(this%Radial, r_ik, a) * this%w(species_map(at%Z(k)))
                  if(my_do_grad_descriptor) dRad_ik(:,a) = GradRadialFunction(this%Radial, r_ik, a) * u_ik * this%w(species_map(at%Z(k)))
               enddo

               if(this%l_max >= 0) then
                  T_cos_ijk(0) = 1.0_dp
                  T_0_cos_ijk = T_cos_ijk(0)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(-1) = 0.0_dp
                     U_cos_ijk(0) = 1.0_dp
                     U_0_cos_ijk = U_cos_ijk(0)
                  endif
               endif

               if(this%l_max >= 1) then
                  T_cos_ijk(1) = cos_ijk
                  T_1_cos_ijk = T_cos_ijk(1)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(1) = 2.0_dp*cos_ijk
                     U_1_cos_ijk = U_cos_ijk(1)
                  endif
               endif

               do b = 2, this%l_max
                  T_n_cos_ijk = 2*cos_ijk*T_1_cos_ijk - T_0_cos_ijk
                  T_0_cos_ijk = T_1_cos_ijk
                  T_1_cos_ijk = T_n_cos_ijk

                  T_cos_ijk(b) = T_n_cos_ijk

                  if(my_do_grad_descriptor) then
                     U_n_cos_ijk = 2*cos_ijk*U_1_cos_ijk - U_0_cos_ijk
                     U_0_cos_ijk = U_1_cos_ijk
                     U_1_cos_ijk = U_n_cos_ijk

                     U_cos_ijk(b) = U_n_cos_ijk
                  endif
               enddo

               i_cosnx = 0
               do a = 1, this%n_max
                  do b = 0, this%l_max
                     i_cosnx = i_cosnx + 1

                     Ang = T_cos_ijk(b)

                     if(my_do_descriptor) &
                        descriptor_out%x(i_desc)%data(i_cosnx) = descriptor_out%x(i_desc)%data(i_cosnx) + Rad_ij(a)*Rad_ik(a)*Ang*0.5_dp

                     if(my_do_grad_descriptor) then

                        dAng_ij = b*U_cos_ijk(b-1) * dcosijk_ij
                        dAng_ik = b*U_cos_ijk(b-1) * dcosijk_ik

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) - &
                        ( Rad_ij(a)*Rad_ik(a)*(dAng_ij+dAng_ik) + dRad_ij(:,a)*Rad_ik(a)*Ang + Rad_ij(a)*dRad_ik(:,a)*Ang ) * 0.5_dp

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) + &
                        (Rad_ij(a)*Rad_ik(a)*dAng_ij + dRad_ij(:,a)*Rad_ik(a)*Ang)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(Rad_ik)) deallocate(Rad_ik)
      if(allocated(T_cos_ijk)) deallocate(T_cos_ijk)
      if(allocated(U_cos_ijk)) deallocate(U_cos_ijk)
      if(allocated(dRad_ij)) deallocate(dRad_ij)
      if(allocated(dRad_ik)) deallocate(dRad_ik)

      call system_timer('cosnx_calc')

   endsubroutine cosnx_calc

   subroutine trihis_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, i_desc, n_index
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Sym_Cor_S, Sym_Cor_A, exp_desc
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, x, exp_arg, dexp_desc
      real(dp), dimension(3,3) :: dx_j, dx_k

      INIT_ERROR(error)

      call system_timer('trihis_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_calc: descriptor object not initialised", error)
      endif
      RAISE_ERROR("trihis_calc: ab686 noticed that this routine needs updating. Remove this line if you know what you are doing, then proceed.", error)

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='trihis_calc args_str')) then
            RAISE_ERROR("trihis_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("trihis_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("trihis_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = trihis_dimensions(this,error)

      allocate(descriptor_out%x(at%N))
      do i = 1, at%N
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%ii(0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%pos(3,0:n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%has_grad_data(0:n_neighbours(at,i)))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.
         endif
      enddo

      do i = 1, at%N

         if(my_do_descriptor) then
            descriptor_out%x(i)%ci(1) = i
            descriptor_out%x(i)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(0) = i
            descriptor_out%x(i)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i)%has_grad_data(0) = .true.
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            if(my_do_grad_descriptor) then
               descriptor_out%x(i)%ii(n) = j
               descriptor_out%x(i)%pos(:,n) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i)%has_grad_data(n) = .true.
            endif

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik >= this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               Sym_Cor_S = r_ij + r_ik
               Sym_Cor_A = (r_ij - r_ik)**2

               x = (/Sym_Cor_S, Sym_Cor_A, cos_ijk/)

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  dx_j(:,1) = u_ij
                  dx_j(:,2) = 2.0_dp*(r_ij - r_ik)*u_ij
                  dx_j(:,3) = dcosijk_ij

                  dx_k(:,1) = u_ik
                  dx_k(:,2) = -2.0_dp*(r_ij - r_ik)*u_ik
                  dx_k(:,3) = dcosijk_ik
               endif

               do i_desc = 1, this%n_gauss

                  exp_arg = (x - this%gauss_centre(:,i_desc))/this%gauss_width(:,i_desc)
                  exp_desc = exp(-0.5_dp*sum(exp_arg**2))

                  if(my_do_descriptor) &
                     descriptor_out%x(i)%data(i_desc) = descriptor_out%x(i)%data(i_desc) + exp_desc

                  if(my_do_grad_descriptor) then
                     dexp_desc = -exp_desc * exp_arg / this%gauss_width(:,i_desc)

                     descriptor_out%x(i)%grad_data(i_desc,:,0) = descriptor_out%x(i)%grad_data(i_desc,:,0) - &
                        matmul(dx_j+dx_k,dexp_desc)
                     descriptor_out%x(i)%grad_data(i_desc,:,n) = descriptor_out%x(i)%grad_data(i_desc,:,n) + &
                        2.0_dp*matmul(dx_j,dexp_desc)
                  endif
               enddo
            enddo
         enddo
      enddo

      call system_timer('trihis_calc')

   endsubroutine trihis_calc

   subroutine water_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i, iO, iH1, iH2, n_index
      integer :: i_desc, mpi_n_procs, mpi_my_proc
      integer, dimension(3) :: shift_1, shift_2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r1, r2
      real(dp), dimension(3) :: v1, v2, u1, u2

      INIT_ERROR(error)

      call system_timer('water_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='water_monomer_calc args_str')) then
            RAISE_ERROR("water_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("water_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = water_monomer_dimensions(this,error)
      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index, error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index, error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,3))
            allocate(descriptor_out%x(i)%ii(3))
            allocate(descriptor_out%x(i)%pos(3,3))
            allocate(descriptor_out%x(i)%has_grad_data(3))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,3))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(water_monomer_index(3,count(at%Z==8)))
      call find_water_monomer(at,water_monomer_index,error=error)

      i_desc = 0
      do i = 1, count(at%Z==8)

         iO = water_monomer_index(1,i)
         iH1 = water_monomer_index(2,i)
         iH2 = water_monomer_index(3,i)

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(iO)) cycle
         endif
         i_desc = i_desc + 1

         v1 = diff_min_image(at,iO,iH1,shift=shift_1)
         v2 = diff_min_image(at,iO,iH2,shift=shift_2)
         r1 = sqrt(dot_product(v1,v1))
         r2 = sqrt(dot_product(v2,v2))
         u1 = v1 / r1
         u2 = v2 / r2

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(:) = water_monomer_index(:,i)
            descriptor_out%x(i_desc)%has_data = .true.
            descriptor_out%x(i_desc)%data(1) = r1+r2
            descriptor_out%x(i_desc)%data(2) = (r1-r2)**2
            descriptor_out%x(i_desc)%data(3) = dot_product(v1,v2)
         endif

         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(:) = water_monomer_index(:,i)
            descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,iO)
            descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,iH1) + matmul(at%lattice,shift_1)
            descriptor_out%x(i_desc)%pos(:,3) = at%pos(:,iH2) + matmul(at%lattice,shift_2)
            descriptor_out%x(i_desc)%has_grad_data(:) = .true.

            descriptor_out%x(i_desc)%grad_data(1,:,1) = -u1-u2                  ! 1st descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(1,:,2) =  u1                     ! 1st descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(1,:,3) =  u2                     ! 1st descriptor wrt rH2
            descriptor_out%x(i_desc)%grad_data(2,:,1) =  2.0_dp*(r1-r2)*(u2-u1) ! 2nd descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(2,:,2) =  2.0_dp*(r1-r2)*u1      ! 2nd descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(2,:,3) = -2.0_dp*(r1-r2)*u2      ! 2nd descriptor wrt rH2
            descriptor_out%x(i_desc)%grad_data(3,:,1) =  -v1-v2                 ! 3rd descriptor wrt rO
            descriptor_out%x(i_desc)%grad_data(3,:,2) =  v2                     ! 3rd descriptor wrt rH1
            descriptor_out%x(i_desc)%grad_data(3,:,3) =  v1                     ! 3rd descriptor wrt rH2
         endif

      enddo

      deallocate(water_monomer_index)
      call system_timer('water_monomer_calc')

   endsubroutine water_monomer_calc

   subroutine water_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, n, &
         iAO, iAH1, iAH2, iBO, iBH1, iBH2, i_distance, n_index
      integer :: mpi_n_procs, mpi_my_proc
      integer, dimension(3) :: shift_AO_BO, shift_AO_AH1, shift_AO_AH2, shift_AO_BH1, shift_AO_BH2, &
         shift_BO_AH1, shift_BO_AH2, shift_BO_BH1, shift_BO_BH2, &
         shift_AH1_AH2, shift_AH1_BH1, shift_AH1_BH2, shift_AH2_BH1, shift_AH2_BH2, shift_BH1_BH2
      real(dp), dimension(3) :: diff_AO_BO, diff_AO_AH1, diff_AO_AH2, diff_AO_BH1, diff_AO_BH2, &
         diff_BO_AH1, diff_BO_AH2, diff_BO_BH1, diff_BO_BH2, &
         diff_AH1_AH2, diff_AH1_BH1, diff_AH1_BH2, diff_AH2_BH1, diff_AH2_BH2, diff_BH1_BH2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r_AO_BO, r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
         r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2
      integer, dimension(1) :: j_array
      real(dp), dimension(15) :: distances

      INIT_ERROR(error)

      call system_timer('water_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='water_dimer_calc args_str')) then
            RAISE_ERROR("water_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("water_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = water_dimer_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,6))
            allocate(descriptor_out%x(i)%ii(6))
            allocate(descriptor_out%x(i)%pos(3,6))
            allocate(descriptor_out%x(i)%has_grad_data(6))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,6))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp

         endif
      enddo

      n_monomers = 0
      do i = 1, at%N
         if(at%Z(i) == 8) n_monomers = n_monomers+1
      enddo

      allocate(water_monomer_index(3,n_monomers))
      call find_water_monomer(at,water_monomer_index,OHH_ordercheck=this%OHH_ordercheck,monomer_cutoff=this%monomer_cutoff,error=error)

      i_desc = 0
      do i = 1, n_monomers
         iAO = water_monomer_index(1,i)
         iAH1 = water_monomer_index(2,i)
         iAH2 = water_monomer_index(3,i)

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(iAO)) cycle
         endif

         diff_AO_AH1 = diff_min_image(at,iAO,iAH1,shift=shift_AO_AH1)
         diff_AO_AH2 = diff_min_image(at,iAO,iAH2,shift=shift_AO_AH2)
         diff_AH1_AH2 = diff_min_image(at,iAH1,iAH2,shift=shift_AH1_AH2)

         r_AO_AH1 = norm(diff_AO_AH1)
         r_AO_AH2 = norm(diff_AO_AH2)
         r_AH1_AH2 = norm(diff_AH1_AH2)

         do n = 1, n_neighbours(at,iAO)
            iBO = neighbour(at,iAO,n,distance=r_AO_BO, diff=diff_AO_BO, shift=shift_AO_BO )
            if(at%Z(iBO) /= 8) cycle
            if( r_AO_BO >= this%cutoff ) cycle
            i_desc = i_desc + 1
            j_array = find(water_monomer_index(1,:) == iBO)
            j = j_array(1)

            iBH1 = water_monomer_index(2,j)
            iBH2 = water_monomer_index(3,j)

            diff_BO_BH1 = diff_min_image(at,iBO,iBH1,shift=shift_BO_BH1)
            diff_BO_BH2 = diff_min_image(at,iBO,iBH2,shift=shift_BO_BH2)
            diff_BH1_BH2 = diff_min_image(at,iBH1,iBH2,shift=shift_BH1_BH2)

            r_BO_BH1 = norm(diff_BO_BH1)
            r_BO_BH2 = norm(diff_BO_BH2)
            r_BH1_BH2 = norm(diff_BH1_BH2)

            diff_AO_BH1 = diff_AO_BO + diff_BO_BH1
            diff_AO_BH2 = diff_AO_BO + diff_BO_BH2
            shift_AO_BH1 = shift_AO_BO + shift_BO_BH1
            shift_AO_BH2 = shift_AO_BO + shift_BO_BH2

            r_AO_BH1 = norm(diff_AO_BH1)
            r_AO_BH2 = norm(diff_AO_BH2)

            diff_BO_AH1 = -diff_AO_BO + diff_AO_AH1
            diff_BO_AH2 = -diff_AO_BO + diff_AO_AH2

            shift_BO_AH1 = -shift_AO_BO + shift_AO_AH1
            shift_BO_AH2 = -shift_AO_BO + shift_AO_AH2

            r_BO_AH1 = norm(diff_BO_AH1)
            r_BO_AH2 = norm(diff_BO_AH2)

            diff_AH1_BH1 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH1
            diff_AH1_BH2 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH2
            diff_AH2_BH1 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH1
            diff_AH2_BH2 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH2

            shift_AH1_BH1 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH1
            shift_AH1_BH2 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH2
            shift_AH2_BH1 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH1
            shift_AH2_BH2 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH2

            r_AH1_BH1 = norm(diff_AH1_BH1)
            r_AH1_BH2 = norm(diff_AH1_BH2)
            r_AH2_BH1 = norm(diff_AH2_BH1)
            r_AH2_BH2 = norm(diff_AH2_BH2)


            distances = (/r_AO_BO, &
                  r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
                  r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2/)

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ water_monomer_index(:,i),water_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (distances+this%dist_shift)**this%power

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_AO_BO, &
               this%cutoff,this%cutoff_transition_width)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ water_monomer_index(:,i),water_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,iAO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,iAH1) + matmul(at%lattice,shift_AO_AH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,3) = at%pos(:,iAH2) + matmul(at%lattice,shift_AO_AH2) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,4) = at%pos(:,iBO) + matmul(at%lattice,shift_AO_BO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,5) = at%pos(:,iBH1) + matmul(at%lattice,shift_AO_BH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,6) = at%pos(:,iBH2) + matmul(at%lattice,shift_AO_BH2) ! TODO: Have to figure out how to do this.

               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff_AO_BO / r_AO_BO     ! 1st descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(1,:,4) = -descriptor_out%x(i_desc)%grad_data(1,:,1)        ! 1st descriptor wrt OB

               descriptor_out%x(i_desc)%grad_data(2,:,1) = -diff_AO_AH1 / r_AO_AH1  ! 2nd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(2,:,2) = -descriptor_out%x(i_desc)%grad_data(2,:,1)        ! 2nd descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff_AO_AH2 / r_AO_AH2  ! 3rd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)        ! 3rd descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff_AO_BH1 / r_AO_BH1  ! 4th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(4,:,5) = -descriptor_out%x(i_desc)%grad_data(4,:,1)        ! 4th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff_AO_BH2 / r_AO_BH2  ! 5th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(5,:,6) = -descriptor_out%x(i_desc)%grad_data(5,:,1)        ! 5th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(6,:,4) = -diff_BO_AH1 / r_BO_AH1  ! 6th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,4)        ! 6th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(7,:,4) = -diff_BO_AH2 / r_BO_AH2  ! 7th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(7,:,3) = -descriptor_out%x(i_desc)%grad_data(7,:,4)        ! 7th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(8,:,4) = -diff_BO_BH1 / r_BO_BH1  ! 8th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(8,:,5) = -descriptor_out%x(i_desc)%grad_data(8,:,4)        ! 8th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(9,:,4) = -diff_BO_BH2 / r_BO_BH2  ! 9th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(9,:,6) = -descriptor_out%x(i_desc)%grad_data(9,:,4)        ! 9th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(10,:,2) = -diff_AH1_AH2 / r_AH1_AH2 ! 10th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(10,:,3) = -descriptor_out%x(i_desc)%grad_data(10,:,2)         ! 10th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(11,:,2) = -diff_AH1_BH1 / r_AH1_BH1 ! 11th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(11,:,5) = -descriptor_out%x(i_desc)%grad_data(11,:,2)         ! 11th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(12,:,2) = -diff_AH1_BH2 / r_AH1_BH2 ! 12th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(12,:,6) = -descriptor_out%x(i_desc)%grad_data(12,:,2)         ! 12th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(13,:,3) = -diff_AH2_BH1 / r_AH2_BH1 ! 13th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(13,:,5) = -descriptor_out%x(i_desc)%grad_data(13,:,3)         ! 13th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(14,:,3) = -diff_AH2_BH2 / r_AH2_BH2 ! 14th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(14,:,6) = -descriptor_out%x(i_desc)%grad_data(14,:,3)         ! 14th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(15,:,5) = -diff_BH1_BH2 / r_BH1_BH2 ! 15th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(15,:,6) = -descriptor_out%x(i_desc)%grad_data(15,:,5)         ! 15th descriptor wrt BH2

               do i_distance = 1, 15
                  descriptor_out%x(i_desc)%grad_data(i_distance,:,:) = descriptor_out%x(i_desc)%grad_data(i_distance,:,:) * &
                     (distances(i_distance)+this%dist_shift)**(this%power-1.0_dp) * this%power
               enddo

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -dcoordination_function(r_AO_BO,&
               this%cutoff,this%cutoff_transition_width) * diff_AO_BO / r_AO_BO
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,4) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)
            endif
         enddo
      enddo

      deallocate(water_monomer_index)
      call system_timer('water_dimer_calc')

   endsubroutine water_dimer_calc

   subroutine A2_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2, n_index
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('A2_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='A2_dimer_calc args_str')) then
            RAISE_ERROR("A2_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("A2_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("A2_dimer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = A2_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      n_monomers = count(at%Z == this%atomic_number) / 2

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index,error)

      i_desc = 0
      do i = 1, at%N
         iA1 = i
         iA2 = neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2,shift=shift_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2,shift=shift_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)

            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)

            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) >= this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ iA1, iA2, iB1, iB2 /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2/)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ iA1, iA2, iB1, iB2 /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 3rd descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 3rd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,1)         ! 4th descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,2) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(5,:,3) = -descriptor_out%x(i_desc)%grad_data(5,:,2)         ! 5th descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 6th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,4) = -descriptor_out%x(i_desc)%grad_data(6,:,2)         ! 6th descriptor wrt B2

            endif
         enddo
      enddo

      call system_timer('A2_dimer_calc')

   endsubroutine A2_dimer_calc

   subroutine AB_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2, n_index
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('AB_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='AB_dimer_calc args_str')) then
            RAISE_ERROR("AB_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("AB_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("AB_dimer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = AB_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_calc: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index,error)

      i_desc = 0
      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)


            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)

            r_A1_A2 = distance_min_image(at,iA1,iA2,shift=shift_A1_A2)
            r_B1_B2 = distance_min_image(at,iB1,iB2,shift=shift_B1_B2)

            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)
            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)

            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) >= this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(:) = (/ AB_monomer_index(:,i),AB_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_B1, r_A2_B2, r_A1_A2, r_B1_B2, r_A1_B2, r_A2_B1 /)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ AB_monomer_index(:,i),AB_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 2nd descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(4,:,2) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,2)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(5,:,4) = -descriptor_out%x(i_desc)%grad_data(5,:,1)         ! 4th descriptor wrt B2
               descriptor_out%x(i_desc)%grad_data(6,:,3) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,3)         ! 5th descriptor wrt B1

            endif
         enddo
      enddo

      deallocate(AB_monomer_index)
      call system_timer('AB_dimer_calc')

   endsubroutine AB_dimer_calc


   subroutine atom_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, grad_d, n_descriptors, n_cross, descriptor_mould_size, &
         i_desc, i_data, i, j, k, n, l, m, l_n_neighbours, i_n, n_index

      real(dp) :: r
      real(dp), dimension(3) :: diff
      real(dp), dimension(1) :: descriptor_mould
      integer, dimension(3) :: shift

      complex(dp), dimension(:), allocatable :: spherical_harmonics
      complex(dp), dimension(:,:), allocatable :: grad_spherical_harmonics

      INIT_ERROR(error)

      call system_timer('atom_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='atom_real_space_calc args_str')) then
            RAISE_ERROR("atom_real_space_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("atom_real_space_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("atom_real_space_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
         d = ( 2 * (this%l_max+1)**2 + 2 ) * l_n_neighbours

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            grad_d = 2 * (this%l_max+1)**2 + 2

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,1:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(1:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,1:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(spherical_harmonics(-this%l_max:this%l_max))
      if( my_do_grad_descriptor ) allocate(grad_spherical_harmonics(3,-this%l_max:this%l_max))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1
         i_data = 0
         i_n = 0

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            !descriptor_out%x(i_desc)%ii(0) = i
            !descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            !descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         do n = 1, n_neighbours(at,i)

            j = neighbour(at,i,n,distance = r, diff = diff, shift=shift)
            if(r >= this%cutoff) cycle
            i_n = i_n + 1

            i_data = i_data + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(i_data) = r
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = diff / r
            endif

            i_data = i_data + 1
            if(my_do_descriptor) descriptor_out%x(i_desc)%data(i_data) = real(i_n,dp)
            if(my_do_grad_descriptor) descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = real(i_n,dp)

            do l = 0, this%l_max
               descriptor_mould_size = size(transfer(spherical_harmonics(-l:l),descriptor_mould))

               do m = -l, l
                  if(my_do_descriptor) spherical_harmonics(m) = SphericalYCartesian(l,m,diff)
                  if(my_do_grad_descriptor) grad_spherical_harmonics(:,m) = GradSphericalYCartesian(l,m,diff)
               enddo

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(i_data+1:i_data+descriptor_mould_size) = transfer(spherical_harmonics(-l:l),descriptor_mould)
               endif

               if(my_do_grad_descriptor) then
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_data+1:i_data+descriptor_mould_size,k,i_n) = &
                     transfer(grad_spherical_harmonics(k,-l:l),descriptor_mould)
                  enddo
               endif

               i_data = i_data + descriptor_mould_size

            enddo
         enddo
      enddo

      if(allocated(spherical_harmonics)) deallocate(spherical_harmonics)
      if(allocated(grad_spherical_harmonics)) deallocate(grad_spherical_harmonics)

      call system_timer('atom_real_space_calc')

   endsubroutine atom_real_space_calc

   subroutine power_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, i_desc, i_pow, l_n_neighbours, n_i, &
         n_descriptors, n_cross, n_index
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij
      integer, dimension(total_elements) :: species_map

      INIT_ERROR(error)

      call system_timer('power_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='power_so3_calc args_str')) then
            RAISE_ERROR("power_so3_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("power_so3_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = power_so3_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross,&
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l = 0, this%l_max
                  i_pow = i_pow + 1

                  descriptor_out%x(i_desc)%data(i_pow) = dot_product(fourier_so3(l,a)%m,fourier_so3(l,a)%m)
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     i_pow = i_pow + 1

                     descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = 2.0_dp * matmul(conjg(dfourier_so3(l,a,n)%mm(:,:)),fourier_so3(l,a)%m(:))
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('power_so3_calc')

   endsubroutine power_so3_calc

   subroutine power_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r
      integer :: i, n, n_i, ji, jn, k, j, i_desc, i_bisp, d, &
         n_descriptors, n_cross, l_n_neighbours, n_index
      integer, dimension(3) :: shift
      integer, dimension(total_elements) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('power_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='power_SO4_calc args_str')) then
            RAISE_ERROR("power_SO4_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("power_SO4_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("power_SO4_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = power_SO4_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( associated(atom_mask_pointer) ) then
            if( .not. atom_mask_pointer(i) ) cycle
         endif

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            ji = neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r >= this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,args_str,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,args_str=args_str,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j = 0, this%j_max
               i_bisp = i_bisp + 1
               descriptor_out%x(i_desc)%data(i_bisp) =  sum( conjg(U(j)%mm)*U(j)%mm )
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 1, n_neighbours(at,i)
               ji = neighbour(at, i, n, distance=r)
               if( r >= this%cutoff ) cycle
               n_i = n_i + 1
               i_bisp = 0
               do j = 0, this%j_max
                  i_bisp = i_bisp + 1
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_bisp,k,n_i) = 2.0_dp * sum( conjg(U(j)%mm)*dU(j,n_i)%mm(k,:,:) )
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%grad_data(:,:,0) = -sum(descriptor_out%x(i_desc)%grad_data(:,:,:), dim=3)
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('power_SO4_calc')

   endsubroutine power_SO4_calc

   subroutine form_gs_index(this, gs_index, error)
      !replacement for the old rs_index
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: a, i, j, smin, smax, sdif, nmin, nmax, ndif, j_species, nu_R, nu_S
      type(int_2d), dimension(:), allocatable :: gs_index

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("form_gs_index: descriptor object not initialised", error)
      endif

      if (( this%nu_R > 2) .OR. (this%nu_R < 0)) then
         RAISE_ERROR("nu_R outside allowed range of 0-2", error)
      endif

      if (( this%nu_S > 2) .OR. (this%nu_S < 0)) then
         RAISE_ERROR("nu_S outside allowed range of 0-2", error)
      endif


      allocate(gs_index(2))
      nu_R = this%nu_R
      nu_S = this%nu_S

      !new loop
      do i = 1,2
         if (nu_R > 0) then
            nmin = 1
            nu_R  = nu_R - 1
         else
            nmin = 0
         endif
         nmax =this%n_max * nmin
         ndif = nmax-nmin

         if (nu_S > 0) then
            smin = 1
            nu_S  = nu_S - 1
         else
            smin = 0
         endif
         smax = this%n_species* smin
         sdif = smax-smin
         allocate(gs_index(i)%mm((ndif+1)*(sdif+1), 2))

         j = 0
         do j_species = smin, smax
            do a = nmin, nmax
               j = j +1
               gs_index(i)%mm(j,:) = (/ j_species, a /)
            enddo
         enddo
      enddo
   endsubroutine form_gs_index


   subroutine form_nu_W(this, W, sym_desc, error)
      !replacement for the old rs_index
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: K, nu_R, nu_S, i, dn, ds, ir, ic, n, s, s2, n2, n2_max, s2_max
      type(real_2d), dimension(:), allocatable :: W
      logical  :: sym_desc

      INIT_ERROR(error)

      if (( this%nu_R > 2) .OR. (this%nu_R < 0)) then
         RAISE_ERROR("nu_R outside allowed range of 0-2", error)
      endif

      if (( this%nu_S > 2) .OR. (this%nu_S < 0)) then
         RAISE_ERROR("nu_S outside allowed range of 0-2", error)
      endif

      ! decide if the l-slices are symmetric matricies
      if ((this%nu_R == 1) .OR. (this%nu_S == 1)) then
         sym_desc = .false.
      else
         sym_desc = .true.
      endif
      allocate(W(2))

      ! construct W(i) as required
      nu_R = this%nu_R
      nu_S = this%nu_S
      do i = 1,2
         ! determine size of W(i) and allocate
         K = 1
         ds = 0
         dn = 0
         n2_max = 1
         s2_max = 1

         if (nu_R > 0) then
            K = K * this%n_max
            nu_R = nu_R -1
            dn = 1
            n2_max = this%n_max
         endif
         if (nu_S > 0) then
            K = K * this%n_species
            nu_S = nu_S -1
            ds = 1
            s2_max = this%n_species
         endif
         allocate(W(i)%mm(this%n_max * this%n_species, K))
         W(i)%mm(:,:) = 0.0_dp
         !loop over S and N, populating W. 4 Loops but just looping over rows and columns of matrix
         ir = 0
         do s = 1, this%n_species
            do n = 1, this%n_max
               ir = ir + 1       ! row index in W
               do s2 = 1, s2_max
                  do n2 = 1, n2_max
                     if (ds*s == ds*s2 .and. dn*n == dn*n2) then
                        ic = 1 + (s2-1)*ds*n2_max + (n2-1)*dn
                        W(i)%mm(ir, ic) = 1.0
                     endif
                  enddo
               enddo
            enddo
         enddo

      enddo
   endsubroutine form_nu_W

   subroutine form_mix_W(this, W, sym_desc, error)
      !replacement for the old rs_index
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      type(real_2d), dimension(:), allocatable :: W
      logical  :: sym_desc
      integer  :: ik, in, is, ic, ir, j, n_seed
      real(dp), dimension(:,:), allocatable :: R
      integer, dimension(:), allocatable  :: seed

      sym_desc = this%sym_mix
      call random_seed(size=n_seed)
      allocate(seed(n_seed))
      seed = 0

      INIT_ERROR(error)
      print*, "In form_mix_W", this%Z_mix, this%R_mix, this%sym_mix
      allocate(W(2))
      if (this%R_mix .and. this%Z_mix) then
         do j = 1, 2
            allocate(W(j)%mm(this%n_species*this%n_max, this%K))
            if (this%sym_mix .and. j == 2) then
               W(2)%mm = W(1)%mm
            else
               do is = 1, this%n_species
                  seed = this%species_Z(is)
                  ir = (is-1)*this%n_max
                  call random_seed(put=seed)
                  call random_number(W(j)%mm(ir+1:ir+this%n_max, :))
               enddo
            endif
         enddo

      elseif (this%Z_mix) then
         do j = 1, 2
            allocate(W(j)%mm(this%n_species*this%n_max, this%K*this%n_max))
            if (this%sym_mix .and. j == 2) then
               W(2)%mm = W(1)%mm
            else
               allocate(R(this%n_species, this%K))
               do is = 1, this%n_species
                  seed = this%species_Z(is)
                  call random_seed(put=seed)
                  call random_number(R(is,:))
               enddo

               ir = 0
               do is = 1, this%n_species
                  do in = 1, this%n_max
                     ir = ir + 1
                     do ik = 1, this%k
                        ic = (ik-1)*this%n_max + in
                        W(j)%mm(ir, ic) = R(is, ik)
                     enddo
                  enddo
               enddo
               deallocate(R)
            endif
         enddo

      elseif (this%R_mix) then
         do j = 1, 2
            allocate(W(j)%mm(this%n_species*this%n_max, this%K*this%n_species))
            if (this%sym_mix .and. j == 2) then
               W(2)%mm = W(1)%mm
            else
               allocate(R(this%n_max, this%K))
               seed = this%n_max
               call random_seed(put=seed)
               call random_number(R)
               ir = 0
               do is = 1, this%n_species
                  do in = 1, this%n_max
                     ir = ir + 1
                     do ik = 1, this%k
                        ic = (is-1)*this%K + ik
                        W(j)%mm(ir, ic) = R(in, ik)
                     enddo
                  enddo
               enddo
               deallocate(R)
            endif
         enddo

      else
         RAISE_ERROR("form_mix_W: not mixing anything", error)
      endif

      !jpd47 center the uniform random numbers on zero
      do j = 1, 2
         W(j)%mm = W(j)%mm - 0.5
      enddo
   endsubroutine form_mix_W




   subroutine form_W(this, W, sym_desc, error)
      !replacement for the old rs_index
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      type(real_2d), dimension(:), allocatable :: W
      logical  :: sym_desc

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("form_W: descriptor object not initialised", error)
      endif

      if ((this%nu_R /= 2 .OR. this%nu_R /= 2) .and. (this%R_mix .or. this%Z_mix .or. this%sym_mix)) then
         RAISE_ERROR("(nu_R, nu_S) = (2,2) required to use channel mixing", error)
      endif

      if (this%R_mix .or. this%Z_mix .or. this%sym_mix) then
         call form_mix_W(this, W, sym_desc, error)
      else
         call form_nu_W(this, W, sym_desc, error)
      endif

   endsubroutine form_W




   ! main branch currently ~1000 lines long, would be nice not to blow this up
   subroutine soap_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

      type real_2d_array
         type(real_2d), dimension(:,:,:), allocatable :: x
      endtype real_2d_array

      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_1d), dimension(:), allocatable, save :: SphericalY_ij
      type(cplx_2d), dimension(:), allocatable, save :: grad_SphericalY_ij

      !SPEED type(cplx_1d), dimension(:,:,:), allocatable :: fourier_so3
      !SPEED type(cplx_2d), dimension(:,:,:), allocatable :: grad_fourier_so3
      type(real_1d), dimension(:,:,:), allocatable, save :: fourier_so3_r, fourier_so3_i, global_fourier_so3_r, global_fourier_so3_i
      type(real_2d), dimension(:,:,:), allocatable, save :: grad_fourier_so3_r, grad_fourier_so3_i
      real(dp), allocatable :: t_g_r(:,:), t_g_i(:,:), t_f_r(:,:), t_f_i(:,:), t_g_f_rr(:,:), t_g_f_ii(:,:)
      integer :: alpha

      logical :: my_do_descriptor, my_do_grad_descriptor, do_two_l_plus_one, sym_desc
      integer :: d, i, j, n, a, b, k, l, m, i_pow, i_coeff, l_n_neighbours, n_i, &
         n_descriptors, n_cross, i_species, j_species, ia, jb, i_desc_i, &
         xml_version, sum_l_n_neighbours, i_pair, i_pair_i, n_index, ub, ia_rs, jb_rs
      integer, dimension(3) :: shift_ij
      integer, dimension(:), allocatable :: i_desc
      integer, dimension(:,:), allocatable :: rs_index
      real(dp) :: r_ij, arg_bess, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, &
         exp_p, exp_m, f_cut, df_cut, norm_descriptor_i, radial_decay, dradial_decay, norm_radial_decay
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:,:), allocatable, save :: radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient, grad_descriptor_i
      real(dp), dimension(:), allocatable, save :: descriptor_i
      real(dp), dimension(:), allocatable :: global_fourier_so3_r_array, global_fourier_so3_i_array
      type(real_2d_array), dimension(:), allocatable :: global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array
      integer, dimension(total_elements) :: species_map
      type(int_2d), dimension(:), allocatable :: gs_index
      complex(dp), allocatable, save :: sphericalycartesian_all_t(:,:), gradsphericalycartesian_all_t(:,:,:)
      complex(dp) :: c_tmp(3)
      integer :: max_n_neigh

      ! jpd47 new variables
      type(real_2d), dimension(:), allocatable, save :: X_r, X_i, W
      type(real_2d), dimension(:, :), allocatable, save :: Y_r, Y_i, dT_i, dT_r
      type(real_2d), dimension(:, :, :), allocatable, save ::  dY_r, dY_i
      type(real_2d), dimension(:, :), allocatable, save :: dX_r, dX_i
      real(dp), dimension(:, :), allocatable, save :: Pl, Pl_g1, Pl_g2
      integer :: ic, K1, K2, ir, ig, ik
      !integer, dimension(:, :), allocatable, save :: neighbour_list
      real(dp) ::  tlpo
      real(dp) :: r_tmp(3)
      real, dimension(0:50) :: sc_times
      complex(dp), dimension(:), allocatable, save :: l_tmp
      logical :: original

!$omp threadprivate(radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient)
!$omp threadprivate(sphericalycartesian_all_t, gradsphericalycartesian_all_t)
!$omp threadprivate(fourier_so3_r, fourier_so3_i, X_i, X_r, Pl, Y_r, Y_i)
!$omp threadprivate(SphericalY_ij,grad_SphericalY_ij)
!$omp threadprivate(descriptor_i, grad_descriptor_i)
!$omp threadprivate(grad_fourier_so3_r, grad_fourier_so3_i, dY_r, dY_i, Pl_g1, Pl_g2, l_tmp, dX_r, dX_i)

      INIT_ERROR(error)
      sc_times = 0.0
      call cpu_time(sc_times(0))
      call system_timer('soap_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_calc: descriptor object not initialised", error)
      endif

      if (( this%nu_R > 2) .OR. (this%nu_R < 0)) then
         RAISE_ERROR("nu_R outside allowed range of 0-2", error)
      endif

      if (( this%nu_S > 2) .OR. (this%nu_S < 0)) then
         RAISE_ERROR("nu_S outside allowed range of 0-2", error)
      endif

      ! jpd47 logical for special routines to keep original power spectrum fast
      original = .false.
      if (this%coupling .and. this%nu_R == 2 .and. this%nu_S == 2 .and. .not. (this%R_mix .or. this%Z_mix)) original = .true.

      !jpd47 TODO replace this with form_WQ
      call cpu_time(sc_times(1))
      call form_W(this, W, sym_desc, error)
      call cpu_time(sc_times(2))
      !print*, "jpd47_timings: form_W took", sc_times(2)-sc_times(1)

      K1 = size(W(1)%mm(0,:))
      K2 = size(W(2)%mm(0,:))
      allocate(Pl(K1, K2))

      ! allocate(Wz(2, this%n_species))
      ! do ik = 1, 2
      !    do i_species = 1, this%n_species
      !       ic = (i_species-1)*this%n_max
      !       allocate(Wz(ik, i_species)%mm(this%n_max, size(W(ik)%mm(0,:))))
      !       Wz(ik, i_species)%mm = W(ik)%mm(ic+1:ic+this%n_max, :)
      !    enddo
      ! enddo


      ! print*, "K1 and K2 are", K1, K2
      ! print*, "K1 and K2 are", K1, K2
      ! print*, "shape of W(1) is", SHAPE(W(1)%mm), W(1)%mm(1,1)
      ! do ia = 1, SIZE(W(1)%mm(:,1))
      !    do jb = 1, SIZE(W(1)%mm(1,:))
      !       if (W(1)%mm(ia, jb) > 0.1 ) print*,"W is", ia, jb, W(1)%mm(ia, jb)
      !    enddo
      ! enddo
      ! print*, "shape of W(2) is", SHAPE(W(2)%mm), W(2)%mm(1,1)
      ! do ia = 1, SIZE(W(2)%mm(:,1))
      !    do jb = 1, SIZE(W(2)%mm(1,:))
      !       print*,"W is", ia, jb, W(2)%mm(ia, jb)
      !    enddo
      ! enddo


      call form_gs_index(this, gs_index, error)
      !jpd47 this would have been wrecking having with anything to do with form_mix_W!!
      !if ((this%nu_R == 1) .OR. (this%nu_S == 1)) then
      !   sym_desc = .false.
      !else
      !   sym_desc = .true.
      !endif

      do i_species = 1, this%n_species
         if(this%species_Z(i_species) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i_species)) = i_species
         endif
      enddo

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      has_atom_mask_name = .false. ! allow atom mask column in the atom table
      atom_mask_pointer => null()  ! allow atom mask column in the atom table
      xml_version = 1423143769     ! This is the version number where the 2l+1 normalisation of soap vectors was introduced
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
            help_string="Name of a logical property in the atoms object. For atoms where this property is " // &
            "true, descriptors are calculated.")

         call param_register(params, 'xml_version', '1423143769', xml_version, &
            help_string="Version of GAP the XML potential file was created")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='soap_calc args_str')) then
            RAISE_ERROR("soap_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("soap_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      if( this%cutoff_dexp > 0 ) then
         if( this%cutoff_rate == 0.0_dp ) then
            norm_radial_decay = 1.0_dp
         else
            norm_radial_decay = this%cutoff_rate / ( 1.0_dp + this%cutoff_rate )
         endif
      else
         norm_radial_decay = 1.0_dp
      endif

      do_two_l_plus_one = (xml_version >= 1423143769)

      allocate(rs_index(2,this%n_max*this%n_species))
      i = 0
      do i_species = 1, this%n_species
         do a = 1, this%n_max
            i = i + 1
            rs_index(:,i) = (/a,i_species/)
         enddo
      enddo

      call finalise(descriptor_out)

      d = soap_dimensions(this, error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      allocate(i_desc(at%N))

      max_n_neigh = 0
      do n_i = 1, at%N
         max_n_neigh = max(max_n_neigh, n_neighbours(at, n_i))
      end do

      call cpu_time(sc_times(2))
      !print*, "jpd47_timings: pre allocation has taken", sc_times(2)-sc_times(0)
      sc_times(1) = sc_times(2)

!$omp parallel default(none) shared(this,my_do_grad_descriptor,d,max_n_neigh, K1, K2) private(i_species, a, l, n_i, ub, ik)
      allocate(descriptor_i(d))
      if(my_do_grad_descriptor) allocate(grad_descriptor_i(d,3))

      allocate(radial_fun(0:this%l_max, this%n_max), radial_coefficient(0:this%l_max, this%n_max))
      !SPEED allocate(fourier_so3(0:this%l_max,this%n_max,this%n_species), SphericalY_ij(0:this%l_max))
      !allocate(fourier_so3_r(0:this%l_max,0:this%n_max,0:this%n_species), fourier_so3_i(0:this%l_max,0:this%n_max,0:this%n_species))
      allocate(SphericalY_ij(0:this%l_max))
      allocate(X_r(0:this%l_max), X_i(0:this%l_max))
      allocate(l_tmp(1:2*this%l_max + 1))
      do l = 0, this%l_max
         allocate(X_r(l)%mm(2*l+1, this%n_species*this%n_max))
         allocate(X_i(l)%mm(2*l+1, this%n_species*this%n_max))
      enddo

      allocate(Y_r(2, 0:this%l_max), Y_i(2, 0:this%l_max))
      do l = 0, this%l_max
         allocate(Y_r(1, l)%mm(2*l+1, K1))
         allocate(Y_i(1, l)%mm(2*l+1, K1))
         allocate(Y_r(2, l)%mm(2*l+1, K2))
         allocate(Y_i(2, l)%mm(2*l+1, K2))
      enddo


      if(my_do_grad_descriptor) then
         allocate(grad_radial_fun(0:this%l_max, this%n_max), grad_radial_coefficient(0:this%l_max, this%n_max))
         allocate(grad_SphericalY_ij(0:this%l_max))
      endif

      allocate(sphericalycartesian_all_t(0:this%l_max, -this%l_max:this%l_max))
      if(my_do_grad_descriptor) then
          allocate(gradsphericalycartesian_all_t(0:this%l_max, -this%l_max:this%l_max, 3))
      end if


      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
         if(my_do_grad_descriptor) allocate(grad_SphericalY_ij(l)%mm(3,-l:l))
      enddo

      if (my_do_grad_descriptor) then
         !jpd47 update general method to same form as original then this could be much neater
         if (original) then
            allocate(Pl_g1(K1, 3*this%n_max))
         else
            allocate(Pl_g1(K1, 3*K2), Pl_g2(3*K1, K2))
         endif

          ! jpd47 allocate new grad storage
         ! original only
         if (original) then
            allocate(dX_r(0:this%l_max, max_n_neigh), dX_i(0:this%l_max, max_n_neigh))
            do l = 0, this%l_max
               do n_i = 1, max_n_neigh
                  allocate(dX_r(l, n_i)%mm(2*l+1, 3*this%n_max))
                  allocate(dX_i(l, n_i)%mm(2*l+1, 3*this%n_max))
               enddo
            enddo
         ! general
         else
            allocate(dY_r(2, 0:this%l_max, max_n_neigh), dY_i(2, 0:this%l_max, max_n_neigh))
            do n_i = 1, max_n_neigh
               do ik = 1, SIZE(dY_r(:, 0, 1))
                  k = K1
                  if (ik == 2) k = K2
                  do l = 0, this%l_max
                     allocate(dY_r(ik, l, n_i)%mm(2*l+1, 3*k))
                     allocate(dY_i(ik, l, n_i)%mm(2*l+1, 3*k))
                  enddo
               enddo
            enddo
         endif

         !jpd47 temporary storage for the gradient cofficients before multiplication
         allocate(dT_r(0:2, 0:this%l_max), dT_i(0:2, 0:this%l_max))
         do l = 0, this%l_max
            allocate(dT_r(0, l)%mm(2*l+1, this%n_max), dT_i(0, l)%mm(2*l+1, this%n_max))
            allocate(dT_r(1, l)%mm(2*l+1, K1), dT_i(1, l)%mm(2*l+1, K1))
            allocate(dT_r(2, l)%mm(2*l+1, K2), dT_i(2, l)%mm(2*l+1, K2))
         enddo

      endif
!$omp end parallel
      call cpu_time(sc_times(2))
      !print*, "jpd47_timings: initial allocation took", sc_times(2)-sc_times(1)
      sc_times(1) = sc_times(2)

      i_desc = 0
      i_desc_i = 0
      do i = 1, at%N
         if( .not. any( at%Z(i) == this%Z ) .and. .not. any(this%Z == 0) ) cycle

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc_i = i_desc_i + 1
         i_desc(i) = i_desc_i

         if(.not. this%global) then ! atomic SOAP
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc_i)%data(d))
               !slow, no need
               !descriptor_out%x(i_desc_i)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc_i)%ci(n_index))
               descriptor_out%x(i_desc_i)%has_data = .false.
               descriptor_out%x(i_desc_i)%covariance_cutoff = 1.0_dp
            endif
            if(my_do_grad_descriptor) then
               l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

               allocate(descriptor_out%x(i_desc_i)%grad_data(d,3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%ii(0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%pos(3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%has_grad_data(0:l_n_neighbours))
               ! slow, no need
               ! descriptor_out%x(i_desc_i)%grad_data = 0.0_dp
               descriptor_out%x(i_desc_i)%grad_data(:,:,0) = 0.0_dp
               descriptor_out%x(i_desc_i)%ii = 0
               descriptor_out%x(i_desc_i)%pos = 0.0_dp
               descriptor_out%x(i_desc_i)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc_i)%grad_covariance_cutoff(3,0:l_n_neighbours))
               descriptor_out%x(i_desc_i)%grad_covariance_cutoff = 0.0_dp
            endif
         endif
      enddo

      allocate( &
         global_fourier_so3_r_array((this%l_max+1)**2 * (this%n_max+1) * (this%n_species+1)), &
         global_fourier_so3_i_array((this%l_max+1)**2 * (this%n_max+1) * (this%n_species+1)), &
         global_grad_fourier_so3_r_array( count(i_desc/=0) ), &
         global_grad_fourier_so3_i_array( count(i_desc/=0) ) )



      if(this%global) then
         if(my_do_descriptor) then
            allocate(descriptor_out%x(1)%data(d))
            allocate(descriptor_out%x(1)%ci(n_index))
            if( any(this%Z == 0) ) then
               descriptor_out%x(1)%ci(:) = (/ (i, i=1, at%N) /)
            else
               forall(i=1:at%N, any(at%Z(i) == this%Z)) descriptor_out%x(1)%ci(i_desc(i)) = i
            endif
            descriptor_out%x(1)%has_data = .true.
            descriptor_out%x(1)%covariance_cutoff = 1.0_dp
         endif ! my_do_descriptor
         if(my_do_grad_descriptor) then
            sum_l_n_neighbours = 0
            do i = 1, at%N

               if(i_desc(i) == 0) then
                  cycle
               else
                  i_desc_i = i_desc(i)
               endif

               l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
               sum_l_n_neighbours = sum_l_n_neighbours + l_n_neighbours + 1 ! include central atom as well!

               ! allocate( &
               ! global_grad_fourier_so3_r_array(i_desc_i)%x(0:this%l_max,0:this%n_max,l_n_neighbours), &
               !global_grad_fourier_so3_i_array(i_desc_i)%x(0:this%l_max,0:this%n_max,l_n_neighbours) )
               allocate( &
               global_grad_fourier_so3_r_array(i_desc_i)%x(0:this%l_max,0:this%n_max,max_n_neigh), &
               global_grad_fourier_so3_i_array(i_desc_i)%x(0:this%l_max,0:this%n_max,max_n_neigh) )

               do n_i = 1, l_n_neighbours
                  do a = 0, this%n_max
                     do l = 0, this%l_max
                        allocate( &
                           global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm(3,-l:l), &
                           global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm(3,-l:l) )
                        global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm = 0.0_dp
                        global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm = 0.0_dp
                     enddo ! l
                  enddo ! a
               enddo ! n_i
            enddo ! i

            allocate(descriptor_out%x(1)%grad_data(d,3,sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%ii(sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%pos(3,sum_l_n_neighbours))
            allocate(descriptor_out%x(1)%has_grad_data(sum_l_n_neighbours))

            allocate(descriptor_out%x(1)%grad_covariance_cutoff(3,sum_l_n_neighbours))
            descriptor_out%x(1)%grad_covariance_cutoff = 0.0_dp
         endif ! my_do_grad_descriptor

         global_fourier_so3_r_array = 0.0_dp
         global_fourier_so3_i_array = 0.0_dp
      endif ! this%global
      call cpu_time(sc_times(2))
      !print*, "jpd47_timings: second round of allocations took", sc_times(2)-sc_times(1)
      sc_times(1) = sc_times(2)

!$omp parallel do schedule(dynamic) default(none) shared(this, at, descriptor_out, my_do_descriptor, my_do_grad_descriptor, d, i_desc, species_map, rs_index, do_two_l_plus_one, gs_index, sym_desc, W, K1, K2, max_n_neigh) &
!$omp shared(global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array, norm_radial_decay) &
!$omp private(i, j, i_species, j_species, a, b, l, m, n, n_i, r_ij, u_ij, d_ij, shift_ij, i_pow, i_coeff, ia, jb, alpha, i_desc_i, ub, ia_rs, jb_rs, ic, ir, ig, ik) &
!$omp private(c_tmp, r_tmp) &
!$omp private(t_g_r, t_g_i, t_f_r, t_f_i, t_g_f_rr, t_g_f_ii) &
!$omp private(f_cut, df_cut, arg_bess, exp_p, exp_m, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lp, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, norm_descriptor_i) &
!$omp private(radial_decay, dradial_decay) &
!$omp reduction(+:global_fourier_so3_r_array,global_fourier_so3_i_array)


      do i = 1, at%N
         call cpu_time(sc_times(3))
         if(i_desc(i) == 0) then
            cycle
         else
            i_desc_i = i_desc(i)
         endif

         if(.not.this%global) then
            if(my_do_descriptor) then
               descriptor_out%x(i_desc_i)%ci(1) = i
               descriptor_out%x(i_desc_i)%has_data = .true.
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(0) = i
               descriptor_out%x(i_desc_i)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc_i)%has_grad_data(0) = .true.

               ! jpd47 zero the gradient contributions
               ! jpd47 original only
               if (original) then
                  do n_i= 1, max_n_neigh
                     do l = 0, this%l_max
                        dX_r(l, n_i)%mm = 0.0_dp
                        dX_i(l, n_i)%mm = 0.0_dp
                     enddo
                  enddo
               ! jpd47 general
               else
                  do n_i = 1, max_n_neigh
                     do l = 0, this%l_max
                        do k = 1, size(dY_r(:, 0, 1))
                           dY_r(k, l, n_i)%mm = 0.0_dp
                           dY_i(k, l, n_i)%mm = 0.0_dp
                        enddo
                     enddo
                  enddo
               endif
            endif
         endif

         radial_fun(0,:) = 0.0_dp
         radial_fun(0,1) = 1.0_dp
         radial_coefficient(0,:) = matmul( radial_fun(0,:), this%cholesky_overlap_basis)

         !jpd47 zero the coefficients and initialise counter
         do l = 0, this%l_max
            X_r(l)%mm = 0.0_dp
            X_i(l)%mm = 0.0_dp
         enddo

         do i_species = 0, this%n_species
            do a = 0, this%n_max

               if ((this%central_reference_all_species .or. this%species_Z(i_species) == at%Z(i) .or. this%species_Z(i_species) == 0) .and. i_species > 0 .and. a > 0) then
                  ic = (i_species-1) * this%n_max + a
                  X_r(0)%mm(1, ic) = this%central_weight * real(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)), dp)
                  X_i(0)%mm(1, ic) = this%central_weight * aimag(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)))
               endif
            enddo
         enddo



! soap_calc 20 takes 0.0052 s
         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij >= this%cutoff ) cycle

            n_i = n_i + 1

            i_species = species_map(at%Z(j))
            if( i_species == 0 ) cycle

            call cpu_time(sc_times(11))

            if(.not. this%global .and. my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(n_i) = j
               descriptor_out%x(i_desc_i)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc_i)%has_grad_data(n_i) = .true.
            endif

            f_cut = coordination_function(r_ij, this%cutoff, this%cutoff_transition_width)
            radial_decay = ( 1.0_dp + this%cutoff_rate ) / ( this%cutoff_rate + ( r_ij / this%cutoff_scale )**this%cutoff_dexp )
            radial_decay = norm_radial_decay * radial_decay

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
               dradial_decay = - this%cutoff_dexp * ( 1.0_dp + this%cutoff_rate ) * ( r_ij / this%cutoff_scale )**this%cutoff_dexp / &
                  ( r_ij * ( this%cutoff_rate + ( r_ij / this%cutoff_scale )**this%cutoff_dexp )**2  )
               dradial_decay = norm_radial_decay * dradial_decay

               df_cut = df_cut * radial_decay + f_cut * dradial_decay
            endif
            f_cut = f_cut * radial_decay

            do a = 1, this%n_max
               arg_bess = 2.0_dp * this%alpha * r_ij * this%r_basis(a)
               exp_p = exp( -this%alpha*( r_ij + this%r_basis(a) )**2 )
               exp_m = exp( -this%alpha*( r_ij - this%r_basis(a) )**2 )

               do l = 0, this%l_max
                  if( l == 0 ) then
                     if(arg_bess == 0.0_dp) then
                        !mo_spher_bess_fi_ki_l = 1.0_dp
                        mo_spher_bess_fi_ki_l = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) )
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        !mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                        !mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                        mo_spher_bess_fi_ki_lm = 0.5_dp * (exp_m + exp_p) / arg_bess
                        mo_spher_bess_fi_ki_l  = 0.5_dp * (exp_m - exp_p) / arg_bess
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     endif
                  else
                     if(arg_bess == 0.0_dp) then
                        mo_spher_bess_fi_ki_l = 0.0_dp
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                        mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                        if(my_do_grad_descriptor) then
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                           mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                        else
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                        endif
                     endif
                  endif

                  !radial_fun(l,a) = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) ) * mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  radial_fun(l,a) = mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  if(my_do_grad_descriptor) grad_radial_fun(l,a) = -2.0_dp * this%alpha * r_ij * mo_spher_bess_fi_ki_l + &
                     l*mo_spher_bess_fi_ki_l / r_ij + mo_spher_bess_fi_ki_lp * 2.0_dp * this%alpha * this%r_basis(a)

               enddo
            enddo

            radial_coefficient = matmul( radial_fun, this%transform_basis )
            if(my_do_grad_descriptor) grad_radial_coefficient = matmul( grad_radial_fun, this%transform_basis ) * f_cut + radial_coefficient * df_cut
            radial_coefficient = radial_coefficient * f_cut

            sphericalycartesian_all_t = SphericalYCartesian_all(this%l_max, d_ij)
            if(my_do_grad_descriptor) gradsphericalycartesian_all_t = GradSphericalYCartesian_all(this%l_max, d_ij)
            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian_all_t(l,m)
                  if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian_all_t(l,m,:)
               enddo
            enddo

            call cpu_time(sc_times(12))

            do l = 0, this%l_max
               do a = 1, this%n_max
                  ic = (i_species-1) * this%n_max + a
                  X_r(l)%mm(:, ic) = X_r(l)%mm(:, ic)  + radial_coefficient(l,a) * real(SphericalY_ij(l)%m(:))
                  X_i(l)%mm(:, ic) = X_i(l)%mm(:, ic)  + radial_coefficient(l,a) * aimag(SphericalY_ij(l)%m(:))
               enddo ! a
            enddo ! l

            if(my_do_grad_descriptor .and. original) then
               do k = 1, 3
                  do l = 0, this%l_max
                     !jpd47 special case for original power spectrum
                     do a = 1, this%n_max
                        ic = (a-1)*3 + k
                        l_tmp(1:2*l+1) =  grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(:) * u_ij(k) + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(k,:)
                        dX_r(l, n_i)%mm(:, ic) = real(l_tmp(1:2*l+1))
                        dX_i(l, n_i)%mm(:, ic) = aimag(l_tmp(1:2*l+1))
                     enddo ! a
                  enddo ! l
               enddo ! k
            endif ! my_do_grad_descriptor

            if(my_do_grad_descriptor .and. (.not. original)) then
               do k = 1, 3
                  do l = 0, this%l_max
                     do a = 1, this%n_max
                        l_tmp(1:2*l+1) =  grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(:) * u_ij(k) + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(k,:)
                        dT_r(0, l)%mm(:, a) =  real(l_tmp(1:2*l+1))
                        dT_i(0, l)%mm(:, a) = aimag(l_tmp(1:2*l+1))
                     enddo ! a

                     !jpd47 operate on the coefficients
                     ic = (i_species-1) * this%n_max
                     do ia = 1, 2
                        dT_r(ia, l)%mm = matmul(dT_r(0, l)%mm, W(ia)%mm(ic+1:ic+this%n_max, :))
                        dT_i(ia, l)%mm = matmul(dT_i(0, l)%mm, W(ia)%mm(ic+1:ic+this%n_max, :))
                     enddo

                     ! jpd47 package coefficients TODO combine these into a single loop
                     ! jpd47 put this into loop above to simplify
                     do ik = 1, K2
                        ir = (ik-1)*3 + k
                        dY_r(2, l, n_i)%mm(:, ir) = dT_r(2, l)%mm(:, ik)
                        dY_i(2, l, n_i)%mm(:, ir) = dT_i(2, l)%mm(:, ik)
                     enddo
                     !jpd47 TODO not doing this was cause of spurious factor of 2 because dY(1, .., ..) getting used when coupling=F and sym_desc = T
                     if (.not. sym_desc .or. .true.) then
                        do ik = 1, K1
                           ir = (ik-1)*3 + k
                           dY_r(1, l, n_i)%mm(:, ir) = dT_r(1, l)%mm(:, ik)
                           dY_i(1, l, n_i)%mm(:, ir) = dT_i(1, l)%mm(:, ik)
                        enddo
                     endif

                  enddo ! l
               enddo ! k
            endif ! my_do_grad_descriptor


            call cpu_time(sc_times(13))
            sc_times(14) = sc_times(14) + sc_times(12) - sc_times(11)
            sc_times(15) = sc_times(15) + sc_times(13) - sc_times(12)
         enddo ! n


         if(this%global .and. my_do_grad_descriptor) then
            global_grad_fourier_so3_r_array(i_desc_i)%x = grad_fourier_so3_r(:,:,1:n_neighbours(at,i,max_dist=this%cutoff))
            global_grad_fourier_so3_i_array(i_desc_i)%x = grad_fourier_so3_i(:,:,1:n_neighbours(at,i,max_dist=this%cutoff))
            !do n_i = lbound(grad_fourier_so3_r,3), ubound(grad_fourier_so3_r,3)
            !   do a = lbound(grad_fourier_so3_r,2), ubound(grad_fourier_so3_r,2)
            !      do l = lbound(grad_fourier_so3_r,1), ubound(grad_fourier_so3_r,1)
            !         global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_r(l,a,n_i)%mm
            !         global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_i(l,a,n_i)%mm
            !      enddo ! l
            !   enddo ! a
            !enddo ! n_i
         endif


         if(this%global) then
            i_coeff = 0
            do ia = 1, (this%n_species+1)*(this%n_max+1)
               a = rs_index(1,ia)
               i_species = rs_index(2,ia)
               do l = 0, this%l_max
                  global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_r(l,a,i_species)%m(:)
                  global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_i(l,a,i_species)%m(:)
                  i_coeff = i_coeff + 2*l+1
               enddo
            enddo
         endif

         call cpu_time(sc_times(4))
         sc_times(5) = sc_times(5) + sc_times(4) - sc_times(3)

         !jpd47 new power spectrum calculation
         call cpu_time(sc_times(6))

         if (this%coupling) then
            !jpd47 standard full tensor product coupling between density channels.
            i_pow = 0
            do l = 0, this%l_max
               tlpo = 1.0_dp
               if (do_two_l_plus_one) tlpo = 1.0_dp / sqrt(2.0_dp * l + 1.0_dp)
               ! special case for regular power spectrum
               if (original) then
                  ! Pl = matmul(tranpose(X_r(l)%mm), X_r(l)%mm)
                  call dgemm('T', 'N', K1, K1, 2*l+1, tlpo, X_r(l)%mm, 2*l+1, X_r(l)%mm, 2*l+1, 0.0_dp, Pl, K1)
                  call dgemm('T', 'N', K1, K1, 2*l+1, tlpo, X_i(l)%mm, 2*l+1, X_i(l)%mm, 2*l+1, 1.0_dp, Pl, K1)
               ! everything else
               else
                  !Y_r(1, l)%mm = matmul(X_r(l)%mm, W(1)%mm)
                  call dgemm('N', 'N', 2*l+1, K1, this%n_max*this%n_species, 1.0_dp, X_r(l)%mm, 2*l+1, W(1)%mm, this%n_max*this%n_species, 0.0_dp, Y_r(1, l)%mm, 2*l+1)
                  call dgemm('N', 'N', 2*l+1, K1, this%n_max*this%n_species, 1.0_dp, X_i(l)%mm, 2*l+1, W(1)%mm, this%n_max*this%n_species, 0.0_dp, Y_i(1, l)%mm, 2*l+1)

                  !skipping this for regular power spec saves 1e-3
                  if (sym_desc) then
                     Y_r(2, l)%mm = Y_r(1, l)%mm
                     Y_i(2, l)%mm = Y_i(1, l)%mm
                  else
                     Y_r(2, l)%mm = matmul(X_r(l)%mm, W(2)%mm)
                     Y_i(2, l)%mm = matmul(X_i(l)%mm, W(2)%mm)
                  endif

                  !Pl = matmul(transpose(Y_r(1, l)%mm), Y_r(2, l)%mm) + matmul(transpose(Y_i(1, l)%mm), Y_i(2, l)%mm)
                  call dgemm('T', 'N', K1, K2, 2*l+1, tlpo, Y_r(1, l)%mm, 2*l+1, Y_r(2, l)%mm, 2*l+1, 0.0_dp, Pl, K1)
                  call dgemm('T', 'N', K1, K2, 2*l+1, tlpo, Y_i(1, l)%mm, 2*l+1, Y_i(2, l)%mm, 2*l+1, 1.0_dp, Pl, K1)
               endif

               ! jpd47 unpack l-slice
               i_pow = l + 1
               do ia = 1, K1
                  ub = K2
                  if (sym_desc) then
                     ub = ia
                  endif
                  a = rs_index(1,ia)
                  do jb = 1, ub
                     b = rs_index(1,jb)
                     if (this%diagonal_radial .and. a /= b) cycle
                     descriptor_i(i_pow) = Pl(ia, jb)
                     if( ia /= jb .and. sym_desc) descriptor_i(i_pow) = descriptor_i(i_pow) * SQRT_TWO
                     i_pow = i_pow + this%l_max+1
                  enddo
               enddo

            enddo
            ! deallocate(Pl)
         else
            !jpd74 elementwise coupling between density channels. For use with tensor-reduced compression.
            do l = 0, this%l_max
               tlpo = 1.0_dp
               if (do_two_l_plus_one) tlpo = 1.0_dp / sqrt(2.0_dp * l + 1.0_dp)
               Y_r(1, l)%mm = matmul(X_r(l)%mm, W(1)%mm)
               Y_i(1, l)%mm = matmul(X_i(l)%mm, W(1)%mm)
               if (this%sym_mix) then
                  Y_r(2, l)%mm = Y_r(1, l)%mm
                  Y_i(2, l)%mm = Y_i(1, l)%mm
               else
                  Y_r(2, l)%mm = matmul(X_r(l)%mm, W(2)%mm)
                  Y_i(2, l)%mm = matmul(X_i(l)%mm, W(2)%mm)
               endif

               i_pow = l + 1
               do ik = 1, K1
                  descriptor_i(i_pow) =  dot_product(Y_r(1, l)%mm(:, ik), Y_r(2, l)%mm(:, ik)) + dot_product(Y_i(1, l)%mm(:, ik), Y_i(2, l)%mm(:, ik))
                  if (do_two_l_plus_one) descriptor_i(i_pow)  = descriptor_i(i_pow) * tlpo
                  i_pow = i_pow + this%l_max + 1
               enddo
            enddo
         endif

         !jpd47 normalise the descriptor
         descriptor_i(d) = 0.0_dp
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))
         if(.not. this%global .and. my_do_descriptor) then
            if(this%normalise) then
               descriptor_out%x(i_desc_i)%data = descriptor_i / norm_descriptor_i
            else
               descriptor_out%x(i_desc_i)%data = descriptor_i
            endif
            descriptor_out%x(i_desc_i)%data(d) = this%covariance_sigma0
         endif

         call cpu_time(sc_times(7))
         !jpd47 new gradients calcuation
         if (my_do_grad_descriptor) then
            n_i = 0
            do n = 1, n_neighbours(at,i)
               j = neighbour(at, i, n, distance = r_ij)
               if( r_ij >= this%cutoff ) cycle
               n_i = n_i + 1
               if( species_map(at%Z(j)) == 0 ) cycle
               grad_descriptor_i = 0.0_dp

               if (this%coupling .and. (.not. original)) then
                  do l = 0, this%l_max
                     ! call dgemm(transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
                     tlpo = 1.0_dp
                     if (do_two_l_plus_one) tlpo = 1.0_dp / sqrt(2.0_dp * l + 1.0_dp)
                     call dgemm('T','N', K1, 3 * K2, 2*l+1, tlpo, Y_r(1, l)%mm, 2*l+1, dY_r(2, l, n_i)%mm, 2*l+1, 0.0_dp, Pl_g1, K1)
                     call dgemm('T','N', K1, 3 * K2, 2*l+1, tlpo, Y_i(1, l)%mm, 2*l+1, dY_i(2, l, n_i)%mm, 2*l+1, 1.0_dp, Pl_g1, K1)
                     ! K1 x 3 K2

                     if (.not. sym_desc) then
                        Pl_g2 = matmul(transpose(dY_r(1, l, n_i)%mm), Y_r(2,l)%mm) + matmul(transpose(dY_i(1, l, n_i)%mm), Y_i(2,l)%mm)
                        ! 3 K1 x K2
                        if(do_two_l_plus_one) Pl_g2 = Pl_g2 / sqrt(2.0_dp * l + 1.0_dp)
                     endif

                     call cpu_time(sc_times(18))
                     ! jpd47 loop over neighbour atoms "unravelling" matrix form of gradients

                     i_pow = l + 1
                     do ia = 1, K1
                        ub = K2
                        if (sym_desc) ub = ia
                        a = modulo(ia, this%n_max)    !jpd47 never checked, probably wrong
                        do jb = 1, ub
                           b = modulo(jb, this%n_max)
                           if (this%diagonal_radial .and. a /= b) cycle
                           ic = (jb-1) * 3
                           ir = (ia-1) * 3
                           if (sym_desc) then
                              r_tmp = Pl_g1(ia, ic+1:ic+3) + Pl_g1(jb, ir+1:ir+3)
                           else
                              r_tmp = Pl_g1(ia, ic+1:ic+3) + Pl_g2(ir+1:ir+3, jb)
                           endif

                           if(ia /= jb .and. sym_desc ) r_tmp = r_tmp * SQRT_TWO
                           !descriptor_out%x(i_desc_i)%grad_data(i_pow,:,n_i) = r_tmp
                           grad_descriptor_i(i_pow, :) = r_tmp
                           i_pow = i_pow + this%l_max+1
                        enddo
                     enddo
                  enddo !l
                  call cpu_time(sc_times(19))
                  sc_times(20) = sc_times(20) + sc_times(19) - sc_times(18)

            !coupling=F and sym_mix=T leads to spurious factor of 2 where gradients are 2x actual
            elseif(.not. this%coupling ) then
               do l = 0, this%l_max
                  tlpo = 1.0_dp
                  if (do_two_l_plus_one) tlpo = 1.0_dp / sqrt(2.0_dp * l + 1.0_dp)
                  i_pow = l + 1
                  do ik = 1, K1
                     ir = (ik-1) * 3
                     !jpd47 doing 3 gradient directions in one shot via matmul
                     r_tmp = matmul(transpose(dY_r(1, l, n_i)%mm(:, ir+1:ir+3)), Y_r(2, l)%mm(:, ik)) + matmul(transpose(dY_i(1, l, n_i)%mm(:, ir+1:ir+3)), Y_i(2, l)%mm(:, ik) )
                     r_tmp = r_tmp + matmul(transpose(dY_r(2, l, n_i)%mm(:, ir+1:ir+3)), Y_r(1, l)%mm(:, ik) ) + matmul(transpose(dY_i(2, l, n_i)%mm(:, ir+1:ir+3)), Y_i(1, l)%mm(:, ik) )
                     grad_descriptor_i(i_pow, :) = r_tmp * tlpo
                     i_pow = i_pow + this%l_max + 1
                  enddo
               enddo

            !jpd47 original power spectrum gradients as special case to exploit sparsity of dX_r w.r.t the neighbour species
            else
               do l = 0, this%l_max
                  tlpo = 1.0_dp
                  if (do_two_l_plus_one) tlpo = 1.0_dp / sqrt(2.0_dp * l + 1.0_dp)
                  !jpd47 TODO try swapping order here... might matter which one is transposed given big size difference
                  call dgemm('T','N', K1, 3 * this%n_max, 2*l+1, tlpo, X_r(l)%mm, 2*l+1, dX_r(l, n_i)%mm, 2*l+1, 0.0_dp, Pl_g1, K1)
                  call dgemm('T','N', K1, 3 * this%n_max, 2*l+1, tlpo, X_i(l)%mm, 2*l+1, dX_i(l, n_i)%mm, 2*l+1, 1.0_dp, Pl_g1, K1)

                  call cpu_time(sc_times(18))
                  ! jpd47 loop over neighbour atoms "unravelling" matrix form of gradients
                  i_pow = l + 1
                  do ia = 1, K1
                     a = rs_index(1,ia)
                     i_species = rs_index(2,ia)
                     do jb = 1, ia
                        b = rs_index(1,jb)
                        j_species = rs_index(2,jb)
                        if (this%diagonal_radial .and. a /= b) cycle
                        if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) then
                           ic = (a-1) * 3
                           grad_descriptor_i(i_pow, :) =  grad_descriptor_i(i_pow, :) + Pl_g1(jb, ic+1:ic+3)
                        endif
                        if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) then
                           ic = (b-1) * 3
                           grad_descriptor_i(i_pow, :) =  grad_descriptor_i(i_pow, :) + Pl_g1(ia, ic+1:ic+3)
                        endif

                        if(ia /= jb) grad_descriptor_i(i_pow, :) = grad_descriptor_i(i_pow, :) * SQRT_TWO
                        i_pow = i_pow + this%l_max+1
                     enddo !jb
                  enddo !ia
               enddo !l
            endif
               !jpd47 now normalise the gradients - checked this bit
               grad_descriptor_i(d, 1:3) = 0.0_dp
               if(.not. this%global) then
                  if( this%normalise ) then
                     grad_descriptor_i = grad_descriptor_i / norm_descriptor_i
                     c_tmp = matmul(descriptor_i,grad_descriptor_i) / norm_descriptor_i**2
                     do k = 1, 3
                        descriptor_out%x(i_desc_i)%grad_data(:,k,n_i) = grad_descriptor_i(:,k) - descriptor_i * c_tmp(k)
                     enddo
                  else
                     descriptor_out%x(i_desc_i)%grad_data(:,:,n_i) = grad_descriptor_i
                  endif
                  descriptor_out%x(i_desc_i)%grad_data(:,:,0) = descriptor_out%x(i_desc_i)%grad_data(:,:,0) - descriptor_out%x(i_desc_i)%grad_data(:,:,n_i)
               endif
            enddo !n_i
         endif


         call cpu_time(sc_times(8))
         sc_times(9) = sc_times(9) + sc_times(7)-sc_times(6)
         sc_times(10) = sc_times(10) + sc_times(8) - sc_times(7)


      enddo ! i
!$omp end parallel do
      print*, "jpd47_timings"
      print*, "jpd47_timings: (5) total coefficients and gradients took", sc_times(5)
      print*, "jpd47_timings: (14) geometry for coefs and gradients took", sc_times(14)
      print*, "jpd47_timings: (15) packing and operating on coefs and gradients took", sc_times(15)
      print*, "jpd47_timings: (23) packing gradients took", sc_times(23)
      print*, "jpd47_timings: (9) assembling power spectrum took", sc_times(9)
      print*, "jpd47_timings: (10) assembling gradients took", sc_times(10)
      print*, "jpd47_timings: (20) unpacking  gradients took", sc_times(20)
      print*, "jpd47_timings: (26) normalising gradients took", sc_times(26)

      call cpu_time(sc_times(2))
      print*, "jpd47_timings: (2-1) looping over atoms took", sc_times(2)-sc_times(1)
      sc_times(1) = sc_times(2)
      !SPEED if(allocated(fourier_so3)) then
      !SPEED    do i_species = 1, this%n_species
      !SPEED       do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
      !SPEED          do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
      !SPEED             deallocate(fourier_so3(l,a,i_species)%m)
      !SPEED          enddo
      !SPEED       enddo
      !SPEED    enddo
      !SPEED    deallocate(fourier_so3)
      !SPEED endif

!$omp parallel default(none) shared(this, max_n_neigh) private(i_species, a, l, n_i, ub)
      if(allocated(fourier_so3_r)) then
         do i_species = lbound(fourier_so3_r,3), ubound(fourier_so3_r,3)
            do a = lbound(fourier_so3_r,2), ubound(fourier_so3_r,2)
               do l = lbound(fourier_so3_r,1), ubound(fourier_so3_r,1)
                  deallocate(fourier_so3_r(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(fourier_so3_r)
      endif
      if(allocated(fourier_so3_i)) then
         do i_species = lbound(fourier_so3_i,3), ubound(fourier_so3_i,3)
            do a = lbound(fourier_so3_i,2), ubound(fourier_so3_i,2)
               do l = lbound(fourier_so3_i,1), ubound(fourier_so3_i,1)
                  deallocate(fourier_so3_i(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(fourier_so3_i)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(grad_SphericalY_ij)) then
         do l = lbound(grad_SphericalY_ij,1), ubound(grad_SphericalY_ij,1)
            deallocate(grad_SphericalY_ij(l)%mm)
         enddo
         deallocate(grad_SphericalY_ij)
      endif

      if (allocated(sphericalycartesian_all_t)) deallocate(sphericalycartesian_all_t)
      if (allocated(gradsphericalycartesian_all_t)) deallocate(gradsphericalycartesian_all_t)

      if(allocated(radial_fun)) deallocate(radial_fun)
      if(allocated(radial_coefficient)) deallocate(radial_coefficient)
      if(allocated(grad_radial_fun)) deallocate(grad_radial_fun)
      if(allocated(grad_radial_coefficient)) deallocate(grad_radial_coefficient)
      if(allocated(descriptor_i)) deallocate(descriptor_i)

      !print *, "about to deallocate grad_descriptor_i"
      if(allocated(grad_descriptor_i)) deallocate(grad_descriptor_i)

        if (allocated(grad_fourier_so3_r)) then ! should really check for grad_fourier_so3_i also
            do n_i = 1, max_n_neigh
               do a = 0, this%n_max
                  do l = 0, this%l_max
                     !SPEED deallocate(grad_fourier_so3(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_r(l,a,n_i)%mm)) deallocate(grad_fourier_so3_r(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_i(l,a,n_i)%mm)) deallocate(grad_fourier_so3_i(l,a,n_i)%mm)
                  enddo
               enddo
            enddo
        endif
        !SPEED deallocate(grad_fourier_so3)
        if (allocated(grad_fourier_so3_r)) deallocate(grad_fourier_so3_r)
        if (allocated(grad_fourier_so3_i)) deallocate(grad_fourier_so3_i)

         !jpd47 trsoap make this parallel in the future
      if (allocated(X_r)) then
         do l = 0, this%l_max
            if (allocated(X_r(l)%mm)) deallocate(X_r(l)%mm)
            if (allocated(X_i(l)%mm)) deallocate(X_i(l)%mm)
         enddo
         if (allocated(X_r)) deallocate(X_r)
         if (allocated(X_i)) deallocate(X_i)
      endif

      if (allocated(Y_r)) then
         do k = 1, 2
            do l = 0, this%l_max
               if (allocated(Y_r(k, l)%mm)) deallocate(Y_r(k, l)%mm)
               if (allocated(Y_i(k, l)%mm)) deallocate(Y_i(k, l)%mm)
            enddo
         enddo
         if (allocated(Y_r)) deallocate(Y_r)
         if (allocated(Y_i)) deallocate(Y_i)
      endif

      if (allocated(dX_r)) then
         do l = 0, this%l_max
            do n_i = 1, max_n_neigh
               if (allocated(dX_r(l, n_i)%mm)) deallocate(dX_r(l, n_i)%mm)
               if (allocated(dX_i(l, n_i)%mm)) deallocate(dX_i(l, n_i)%mm)
            enddo
         enddo
         if (allocated(dX_r)) deallocate(dX_r)
         if (allocated(dX_i)) deallocate(dX_i)
      endif

      if (allocated(dY_R)) then
         do n_i = 1, size(dY_R(1, 0, :))
            do ik = 1, size(dY_R(:, 0, 1))
               do l = 0, this%l_max
                  if (allocated(dY_i(ik, l, n_i)%mm)) deallocate(dY_i(ik, l, n_i)%mm)
                  if (allocated(dY_r(ik, l, n_i)%mm)) deallocate(dY_r(ik, l, n_i)%mm)
               enddo
            enddo
         enddo
         if (allocated(dY_i)) deallocate(dY_i)
         if (allocated(dY_r)) deallocate(dY_r)
      endif

      if (allocated(dT_r)) then
         do k = 0, 2
            do l = 0, this%l_max
               deallocate(dT_r(k, l)%mm, dT_i(k, l)%mm)
            enddo
         enddo
         deallocate(dT_r, dT_i)
      endif


      if (allocated(Pl_g1)) deallocate(Pl_g1)
      if (allocated(Pl_g2)) deallocate(Pl_g2)

      if (allocated(l_tmp)) deallocate(l_tmp)

!$omp end parallel
      if (allocated(W)) deallocate(W)
      ! if (allocated(Wz)) then
      !    do k = 1,2
      !       do i_species = 1, this%n_species
      !          if (allocated(Wz(k, i_species)%mm)) deallocate(Wz(k, i_species)%mm)
      !       enddo
      !    enddo
      !    deallocate(Wz)
      ! endif
      if (allocated(Pl)) deallocate(Pl)




      if(this%global) then
         allocate(global_fourier_so3_r(0:this%l_max,0:this%n_max,0:this%n_species), global_fourier_so3_i(0:this%l_max,0:this%n_max,0:this%n_species), &
            descriptor_i(d) )

         i_coeff = 0
         do ia = 1, (this%n_species+1)*(this%n_max+1)
            a = rs_index(1,ia)
            i_species = rs_index(2,ia)
            do l = 0, this%l_max
               allocate(global_fourier_so3_r(l,a,i_species)%m(-l:l))
               allocate(global_fourier_so3_i(l,a,i_species)%m(-l:l))
               global_fourier_so3_r(l,a,i_species)%m(:) = global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1)
               global_fourier_so3_i(l,a,i_species)%m(:) = global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1)
               i_coeff = i_coeff + 2*l+1
            enddo
         enddo

         i_pow = 0
         do ia = 1, size(gs_index(1)%mm(:,0))
            a = gs_index(1)%mm(ia,2)
            i_species = gs_index(1)%mm(ia,1)

            !set upper bound for the second loop
            ub = size(gs_index(2)%mm(:,0))
            if (sym_desc) then
               ub = ia
            endif
            do jb = 1, ub
               b = gs_index(2)%mm(jb, 2)
               j_species = gs_index(2)%mm(jb, 1)

               if(this%diagonal_radial .and. a /= b) cycle

               do l = 0, this%l_max
                  i_pow = i_pow + 1
                  !SPEED descriptor_i(i_pow) = real( dot_product(fourier_so3(l,a,i_species)%m, fourier_so3(l,b,j_species)%m) )
                  descriptor_i(i_pow) = dot_product(global_fourier_so3_r(l,a,i_species)%m, global_fourier_so3_r(l,b,j_species)%m) + dot_product(global_fourier_so3_i(l,a,i_species)%m, global_fourier_so3_i(l,b,j_species)%m)
                  if(do_two_l_plus_one) descriptor_i(i_pow) = descriptor_i(i_pow) / sqrt(2.0_dp * l + 1.0_dp)
                  if( ia /= jb .and. sym_desc) then
                     descriptor_i(i_pow) = descriptor_i(i_pow) * SQRT_TWO
                  endif
               enddo !l
            enddo !jb
         enddo !ia


         descriptor_i(d) = 0.0_dp
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))
         if( norm_descriptor_i .feq. 0.0_dp ) norm_descriptor_i = tiny(1.0_dp)
         if(my_do_descriptor) then
            if(this%normalise) then
               descriptor_out%x(1)%data = descriptor_i / norm_descriptor_i
            else
               descriptor_out%x(1)%data = descriptor_i
            endif
            descriptor_out%x(1)%data(d) = this%covariance_sigma0
         endif

         if(my_do_grad_descriptor .and. .false.) then
       allocate(t_g_r((this%n_max+1)*3, 2*this%l_max+1), t_g_i((this%n_max+1)*3, 2*this%l_max+1))
	    allocate(t_f_r((this%n_max+1)*(this%n_species+1), 2*this%l_max+1), t_f_i((this%n_max+1)*(this%n_species+1), 2*this%l_max+1))
	    allocate(t_g_f_rr((this%n_max+1)*3, (this%n_max+1)*(this%n_species+1)), t_g_f_ii((this%n_max+1)*3, (this%n_max+1)*(this%n_species+1)))
            allocate(grad_descriptor_i(d,3))

            i_pair = 0
            do i = 1, at%N

               if(i_desc(i) == 0) then
                  cycle
               else
                  i_desc_i = i_desc(i)
               endif

               i_pair = i_pair + 1
               i_pair_i = i_pair ! accumulates \frac{ \partial p^{(j)} }{ \partial r_{ji\alpha} }

               descriptor_out%x(1)%ii(i_pair_i) = i
               descriptor_out%x(1)%pos(:,i_pair_i) = 0.0_dp
               descriptor_out%x(1)%has_grad_data(i_pair_i) = .true.
               descriptor_out%x(1)%grad_data(:,:,i_pair_i) = 0.0_dp

               n_i = 0
               do n = 1, n_neighbours(at,i)
                  j = neighbour(at, i, n, distance = r_ij, diff = d_ij)
                  if( r_ij >= this%cutoff ) cycle

                  n_i = n_i + 1
                  i_pair = i_pair + 1 ! \frac{ \partial p^{(i)} }{ \partial r_{ij\alpha} }

                  descriptor_out%x(1)%ii(i_pair) = j
                  descriptor_out%x(1)%pos(:,i_pair) = d_ij
                  descriptor_out%x(1)%has_grad_data(i_pair) = .true.

                  i_pow = 0
                  grad_descriptor_i = 0.0_dp

                  !global gradient loop
                  do l=0, this%l_max
                     do a = 0, this%n_max
                        do alpha=1, 3
                           t_g_r(3*a+alpha, 1:2*l+1) = global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm(alpha,-l:l)
                           t_g_i(3*a+alpha, 1:2*l+1) = global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm(alpha,-l:l)
                        enddo
                     enddo

                     do ia = 1, (this%n_species+1)*(this%n_max+1)
                        a = rs_index(1,ia)
                        i_species = rs_index(2,ia)

                        t_f_r(ia, 1:2*l+1) = global_fourier_so3_r(l,a,i_species)%m(-l:l)
                        t_f_i(ia, 1:2*l+1) = global_fourier_so3_i(l,a,i_species)%m(-l:l)
                     enddo

                     call dgemm('N','T',(this%n_max+1)*3, (this%n_max+1)*(this%n_species+1), 2*l+1, 1.0_dp, &
                     t_g_r(1,1), size(t_g_r,1), t_f_r(1,1), size(t_f_r,1), 0.0_dp, t_g_f_rr(1,1), size(t_g_f_rr, 1))
                     call dgemm('N','T',(this%n_max+1)*3, (this%n_max+1)*(this%n_species+1), 2*l+1, 1.0_dp, &
                     t_g_i(1,1), size(t_g_i,1), t_f_i(1,1), size(t_f_i,1), 0.0_dp, t_g_f_ii(1,1), size(t_g_f_ii, 1))
                     !t_g_f_rr = matmul(t_g_r,transpose(t_f_r))
                     !t_g_f_ii = matmul(t_g_i,transpose(t_f_i))

                     i_pow = l+1
                     do ia = 1, size(gs_index(1)%mm(:,0))
                        a = gs_index(1)%mm(ia,2)
                        i_species = gs_index(1)%mm(ia,1)
                        ia_rs = i_species*(this%n_max+1)+a+1

                        !set upper bound for the second loop
                        ub = size(gs_index(2)%mm(:,0))
                        if (sym_desc) then
                           ub = ia
                        endif
                        do jb = 1, ub
                           b = gs_index(2)%mm(jb, 2)
                           j_species = gs_index(2)%mm(jb, 1)
                           jb_rs = j_species*(this%n_max+1)+b+1

                           if(this%diagonal_radial .and. a /= b) cycle

                           if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*a+1:3*(a+1),jb_rs) + t_g_f_ii(3*a+1:3*(a+1),jb_rs)
                           if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*b+1:3*(b+1),ia_rs) + t_g_f_ii(3*b+1:3*(b+1),ia_rs)

                           if(do_two_l_plus_one) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) / sqrt(2.0_dp * l + 1.0_dp)
                           if(ia /= jb .and. sym_desc ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                           i_pow = i_pow + this%l_max+1
                        enddo
                     enddo

                  end do !l



                  grad_descriptor_i(d, 1:3) = 0.0_dp
                  if( this%normalise ) then
                     descriptor_out%x(1)%grad_data(:,:,i_pair) = grad_descriptor_i / norm_descriptor_i
                     do k = 1, 3
                        descriptor_out%x(1)%grad_data(:,k,i_pair) = descriptor_out%x(1)%grad_data(:,k,i_pair) - descriptor_i * dot_product(descriptor_i,grad_descriptor_i(:,k)) / norm_descriptor_i**3
                     enddo
                  else
                     descriptor_out%x(1)%grad_data(:,:,i_pair) = grad_descriptor_i
                  endif

                  descriptor_out%x(1)%grad_data(:,:,i_pair_i) = descriptor_out%x(1)%grad_data(:,:,i_pair_i) - descriptor_out%x(1)%grad_data(:,:,i_pair)
               enddo ! n/n_i

            enddo ! i

            deallocate(grad_descriptor_i)
            deallocate(t_f_r, t_f_i)
            deallocate(t_g_r, t_g_i)
            deallocate(t_g_f_rr, t_g_f_ii)
         endif ! my_do_grad_descriptor

         do i_species = lbound(global_fourier_so3_r,3), ubound(global_fourier_so3_r,3)
            do a = lbound(global_fourier_so3_r,2), ubound(global_fourier_so3_r,2)
               do l = lbound(global_fourier_so3_r,1), ubound(global_fourier_so3_r,1)
                  deallocate(global_fourier_so3_r(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(global_fourier_so3_r)
         do i_species = lbound(global_fourier_so3_i,3), ubound(global_fourier_so3_i,3)
            do a = lbound(global_fourier_so3_i,2), ubound(global_fourier_so3_i,2)
               do l = lbound(global_fourier_so3_i,1), ubound(global_fourier_so3_i,1)
                  deallocate(global_fourier_so3_i(l,a,i_species)%m)
               enddo
            enddo
         enddo
         deallocate(global_fourier_so3_i)

         if(allocated(descriptor_i)) deallocate(descriptor_i)
      endif ! this%global

      if(allocated(global_fourier_so3_r_array)) deallocate(global_fourier_so3_r_array)
      if(allocated(global_fourier_so3_i_array)) deallocate(global_fourier_so3_i_array)

      if(allocated(global_grad_fourier_so3_r_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_r_array,1), ubound(global_grad_fourier_so3_r_array,1)
            if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_r_array)
      endif
      if(allocated(global_grad_fourier_so3_i_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_i_array,1), ubound(global_grad_fourier_so3_i_array,1)
            if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_i_array)
      endif

      if(allocated(rs_index)) deallocate(rs_index)
      if(allocated(i_desc)) deallocate(i_desc)
      if (allocated(gs_index)) deallocate(gs_index)



      call system_timer('soap_calc')

      call cpu_time(sc_times(2))
      print*, "jpd47_timings: in total soap_calc took", sc_times(2)-sc_times(0)
   endsubroutine soap_calc





   subroutine rdf_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(rdf), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, i_n, l_n_neighbours, i_desc, n_descriptors, n_cross, n_index
      integer, dimension(3) :: shift
      real(dp) :: r_ij, f_cut, df_cut
      real(dp), dimension(3) :: u_ij
      real(dp), dimension(:), allocatable :: rdf_ij

      INIT_ERROR(error)

      call system_timer('rdf_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='rdf_calc args_str')) then
            RAISE_ERROR("rdf_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("rdf_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = rdf_dimensions(this,error)
      allocate(rdf_ij(d))

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            rdf_ij = exp( -0.5_dp * (r_ij - this%r_gauss)**2 / this%w_gauss**2 )
            f_cut = coordination_function(r_ij,this%cutoff,this%transition_width)

            if(my_do_descriptor) &
               descriptor_out%x(i_desc)%data = descriptor_out%x(i_desc)%data + rdf_ij * f_cut

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff,this%transition_width)

               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.

               descriptor_out%x(i_desc)%grad_data(:,:,i_n) = ( - ( rdf_ij * (r_ij - this%r_gauss) / this%w_gauss**2 ) * f_cut + rdf_ij * df_cut ) .outer. u_ij
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,i_n)
            endif
         enddo
      enddo

      if(allocated(rdf_ij)) deallocate(rdf_ij)

      call system_timer('rdf_calc')

   endsubroutine rdf_calc

   subroutine as_distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(as_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, &
         n_neighbours_coordination_i, n_neighbours_coordination_ij, n_index
      integer, dimension(3) :: shift
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, cos_jik, f_cut_i, f_cut_j, f_cut_ij, f_cut_ik, f_cut_jk, f_cut_as_i, f_cut_as_j, rho_i, rho_j
      real(dp), dimension(3) :: u_ij, u_ik, u_jk

      INIT_ERROR(error)
      call system_timer('as_distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='as_distance_2b_calc args_str')) then
            RAISE_ERROR("as_distance_2b_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("as_distance_2b_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("as_distance_2b_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = as_distance_2b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)


      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij, cosines=u_ij)

            if(r_ij > this%max_cutoff .or. r_ij < this%min_cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            rho_i = 0.0_dp
            f_cut_i = 0.0_dp

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance=r_ik, cosines=u_ik)

               if(r_ik > this%coordination_cutoff) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               f_cut_ik = coordination_function(r_ik,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_i = f_cut_i + f_cut_ik
               rho_i = rho_i + 0.5_dp * ( erf(cos_ijk/this%overlap_alpha) + 1.0_dp ) * f_cut_ik**2
            enddo

            rho_i = rho_i / f_cut_i

            if(rho_i > this%as_cutoff) cycle

            rho_j = 0.0_dp
            f_cut_j = 0.0_dp

            do m = 1, n_neighbours(at,j)
               k = neighbour(at, j, m, distance=r_jk, cosines=u_jk)

               if(r_jk > this%coordination_cutoff) cycle

               cos_jik = dot_product(-u_ij,u_jk)
               f_cut_jk = coordination_function(r_jk,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_j = f_cut_j + f_cut_jk
               rho_j = rho_j + 0.5_dp * ( erf(cos_jik/this%overlap_alpha) + 1.0_dp ) * f_cut_jk**2
            enddo

            if(rho_j > this%as_cutoff) cycle
            ! all three conditions fulfilled: pair within lower and upper cutoff, asymmetricity lower than threshold

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(d))
               descriptor_out%x(i_desc)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc)%ci(n_index))
               descriptor_out%x(i_desc)%has_data = .false.
            endif

            if(my_do_grad_descriptor) then
               n_neighbours_coordination_ij = n_neighbours(at,i,max_dist=this%coordination_cutoff) + &
               n_neighbours(at,j,max_dist=this%coordination_cutoff) + 2

               allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%ii(0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%pos(3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%has_grad_data(0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_data = 0.0_dp
               descriptor_out%x(i_desc)%ii = 0
               descriptor_out%x(i_desc)%pos = 0.0_dp
               descriptor_out%x(i_desc)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
            endif
         enddo
      enddo

      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if(r_ij > this%max_cutoff .or. r_ij < this%min_cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            rho_i = 0.0_dp
            f_cut_i = 0.0_dp

            do m = 1, n_neighbours(at,i)
               k = neighbour(at, i, m, distance=r_ik, cosines=u_ik)

               if(r_ik > this%coordination_cutoff) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               f_cut_ik = coordination_function(r_ik,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_i = f_cut_i + f_cut_ik
               rho_i = rho_i + 0.5_dp * ( erf(cos_ijk/this%overlap_alpha) + 1.0_dp ) * f_cut_ik**2
            enddo

            rho_i = rho_i / f_cut_i

            if(rho_i > this%as_cutoff) cycle

            rho_j = 0.0_dp
            f_cut_j = 0.0_dp

            do m = 1, n_neighbours(at,j)
               k = neighbour(at, j, m, distance=r_jk, cosines=u_jk)

               if(r_jk > this%coordination_cutoff) cycle

               cos_jik = dot_product(-u_ij,u_jk)
               f_cut_jk = coordination_function(r_jk,this%coordination_cutoff,this%coordination_transition_width)

               f_cut_j = f_cut_j + f_cut_jk
               rho_j = rho_j + 0.5_dp * ( erf(cos_jik/this%overlap_alpha) + 1.0_dp ) * f_cut_jk**2
            enddo

            if(rho_j > this%as_cutoff) cycle
            ! all three conditions fulfilled: pair within lower and upper cutoff, asymmetricity lower than threshold

            i_desc = i_desc + 1

            f_cut_ij = coordination_function(r_ij,this%max_cutoff,this%max_transition_width,this%min_cutoff,this%min_transition_width)
            f_cut_as_i = coordination_function(rho_i,this%as_cutoff, this%as_transition_width)
            f_cut_as_j = coordination_function(rho_j,this%as_cutoff, this%as_transition_width)

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%ci(1:2) = (/i,j/)

               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%data(2) = f_cut_i + f_cut_j
               descriptor_out%x(i_desc)%data(3) = (f_cut_i - f_cut_j)**2

               descriptor_out%x(i_desc)%covariance_cutoff = f_cut_ij * f_cut_as_i * f_cut_as_j
            endif
            if(my_do_grad_descriptor) then
               n_neighbours_coordination_i = n_neighbours(at,i,max_dist=this%coordination_cutoff)

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%coordination_cutoff,this%coordination_transition_width)*u_ij

               !descriptor_out%x(i_desc)%ii(1) = j
               !descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               !descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               !descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               !descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

               !descriptor_out%x(i_desc)%ii(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%ii(:)
               !descriptor_out%x(i_desc)%pos(:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%pos(:,:)
               !descriptor_out%x(i_desc)%has_grad_data(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%has_grad_data(:)
               !descriptor_out%x(i_desc)%grad_data(2,:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%grad_data(1,:,:)
               !descriptor_out%x(i_desc)%grad_data(3,:,2:n_neighbours_coordination_i+2) = 2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
               !   descriptor_coordination%x(i)%grad_data(1,:,:)

               !descriptor_out%x(i_desc)%ii(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%ii(:)
               !descriptor_out%x(i_desc)%pos(:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%pos(:,:)
               !descriptor_out%x(i_desc)%has_grad_data(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%has_grad_data(:)
               !descriptor_out%x(i_desc)%grad_data(2,:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%grad_data(1,:,:)
               !descriptor_out%x(i_desc)%grad_data(3,:,n_neighbours_coordination_i+3:) = -2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
               !   descriptor_coordination%x(j)%grad_data(1,:,:)

            endif
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('as_distance_2b_calc')

   endsubroutine as_distance_2b_calc


   subroutine alex_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(alex), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor

      integer :: i, j, n, d, p, q, r, a, b, c, n_i, n_radial, pp, i_desc, &
         l_n_neighbours, desc_index, n_cross, n_descriptors, n_index
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: d_ij
      real(dp), dimension(:), allocatable :: neighbour_dists
      real(dp), dimension(:,:), allocatable :: neighbour_vecs
      integer, dimension(total_elements) :: species_map
      real(dp), allocatable :: S0(:), S1(:,:), S2(:,:,:), S0der(:,:,:), S1der(:,:,:,:), S2der(:,:,:,:,:)

      INIT_ERROR(error)

      call system_timer('alex_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("alex_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='alex_calc args_str')) then
            RAISE_ERROR("alex_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("alex_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = alex_dimensions(this,error)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))

      n_radial = this%power_max-this%power_min+1

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         ! number of neighbours for the current atom within the descriptor cutoff
         l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)
         allocate(neighbour_vecs(3,l_n_neighbours), neighbour_dists(l_n_neighbours))
         allocate(S0(n_radial), S1(3,n_radial), S2(3,3,n_radial))
         if(my_do_grad_descriptor) then
            allocate( &
               S0der(n_radial,l_n_neighbours,3), &
               S1der(3,n_radial,l_n_neighbours,3), &
               S2der(3,3,n_radial,l_n_neighbours,3) )
         endif

         n_i = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, diff=d_ij)
            if( r_ij >= this%cutoff ) cycle
            n_i = n_i + 1
            neighbour_vecs(:,n_i) = d_ij
            neighbour_dists(n_i) = r_ij
         end do

         do p = 1,n_radial
            pp = -(p+this%power_min-1)
            S0(p) = sum(neighbour_dists**pp)
            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  S0der(p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(:,n_i)
               enddo
            endif

            S1(:,p) = matmul(neighbour_vecs, neighbour_dists**pp)
            !do a = 1,3
               !S1(a, p) = sum(neighbour_vecs(a,:)*neighbour_dists**pp)
            !end do
            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  do a = 1,3
                     S1der(a,p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(a,n_i) * neighbour_vecs(:,n_i)
                     S1der(a,p,n_i,a) = S1der(a,p,n_i,a) + neighbour_dists(n_i)**pp
                  end do
               enddo
            endif

            !do a=1,3
            do b=1,3
               S2(:,b,p) = matmul(neighbour_vecs, neighbour_vecs(b,:)*neighbour_dists**pp)
               !S2(a,b,p) = sum(neighbour_vecs(a,:)*neighbour_vecs(b,:)*neighbour_dists**pp)
            end do
            !end do

            if(my_do_grad_descriptor) then
               do n_i = 1, l_n_neighbours
                  do a = 1,3
                     do b = 1,3
                        S2der(a,b,p,n_i,:) = pp * neighbour_dists(n_i)**(pp-2) * neighbour_vecs(a,n_i) * neighbour_vecs(b,n_i) * neighbour_vecs(:,n_i)
                     end do
                  end do

                  do a = 1,3
                     do b = 1,3
                        S2der(a,b,p,n_i,b) = S2der(a,b,p,n_i,b) + neighbour_dists(n_i)**pp * neighbour_vecs(a,n_i)
                        S2der(a,b,p,n_i,a) = S2der(a,b,p,n_i,a) + neighbour_dists(n_i)**pp * neighbour_vecs(b,n_i)
                     end do
                  end do
               enddo
            endif
         end do

         descriptor_out%x(i_desc)%data(1:n_radial) = S0
         descriptor_out%x(i_desc)%data(n_radial+1:n_radial+n_radial**2) = reshape(matmul(transpose(S1), S1), (/n_radial**2/))
         desc_index = n_radial+n_radial**2+1
         do p = 1,n_radial
            do q = 1,n_radial
               descriptor_out%x(i_desc)%data(desc_index) = sum(S2(:,:,p) * S2(:,:,q))
               desc_index = desc_index + 1
            end do
         end do

         do p = 1,n_radial
            do q = 1,n_radial
               do r = 1,n_radial
                  descriptor_out%x(i_desc)%data(desc_index) =  dot_product(S1(:,p), matmul(S2(:,:,q), S1(:,r)))
                  desc_index = desc_index + 1
               end do
            end do
         end do

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 1, n_neighbours(at,i)
               j = neighbour(at, i, n, distance = r_ij, shift=shift_ij)
               if( r_ij >= this%cutoff ) cycle

               n_i = n_i + 1

               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.

               descriptor_out%x(i_desc)%grad_data(1:n_radial,:,n_i) = S0der(:,n_i,:)

               desc_index = n_radial + 1
               do p = 1,n_radial
                  do q = 1,n_radial
                     do a = 1, 3
                        do c = 1, 3
                           descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                              S1der(a,p,n_i,c)*S1(a,q) + S1(a,p)*S1der(a,q,n_i,c)
                        enddo
                     enddo
                     desc_index = desc_index + 1
                  enddo
               enddo

               do p = 1, n_radial
                  do q = 1, n_radial
                     do a = 1, 3
                        do b = 1, 3
                           do c = 1, 3
                              descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                                 S2der(a,b,p,n_i,c)*S2(a,b,q) + S2(a,b,p)*S2der(a,b,q,n_i,c)
                           enddo
                        enddo
                     enddo
                     desc_index = desc_index + 1
                  enddo
               enddo

               do p = 1, n_radial
                  do q = 1, n_radial
                     do r = 1, n_radial
                        do a = 1, 3
                           do b = 1, 3
                              do c = 1, 3
                                 descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) = descriptor_out%x(i_desc)%grad_data(desc_index,c,n_i) + &
                                    S1der(a,p,n_i,c) * S2(a,b,q)        * S1(b,r) + &
                                    S1(a,p)        * S2der(a,b,q,n_i,c) * S1(b,r) + &
                                    S1(a,p)        * S2(a,b,q)        * S1der(b,r,n_i,c)
                              enddo
                           enddo
                        enddo
                        desc_index = desc_index + 1
                     enddo
                  enddo
               enddo
            enddo

            descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n_i)
            deallocate(S0der, S1der, S2der)
         endif

         deallocate(neighbour_vecs, neighbour_dists, S0, S1, S2)
      enddo


      call system_timer('alex_calc')

   endsubroutine alex_calc

   subroutine distance_Nb_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i_desc, i_data, i, j, ii, jj, kk, ll, &
         iConnectivity, n_index
      integer, dimension(3) :: s_i, s_j
      real(dp) :: r_ij, fcut_connectivity
      real(dp), dimension(3) :: dfcut_connectivity
      real(dp), dimension(3) :: d_ij
      integer, dimension(:,:,:), allocatable :: atoms_in_descriptors
      real(dp), dimension(:,:), allocatable :: fcut_pair, dfcut_pair
      real(dp), dimension(:,:,:), allocatable :: directions

      logical, dimension(:), pointer :: atom_mask_pointer => null()

      INIT_ERROR(error)

      call system_timer('distance_Nb_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='distance_Nb_calc args_str')) then
            RAISE_ERROR("distance_Nb_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         atom_mask_pointer => null()

         if( has_atom_mask_name ) then
            if( .not. this%compact_clusters ) then
               RAISE_ERROR("distance_Nb_calc: MPI/LAMMPS ready only for compact_clusters=T type of distance_Nb.", error)
            endif

            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("distance_Nb_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      d = distance_Nb_dimensions(this,error)
      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask = atom_mask_pointer,n_index=n_index, error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,this%order))
            allocate(descriptor_out%x(i)%ii(this%order))
            allocate(descriptor_out%x(i)%pos(3,this%order))
            allocate(descriptor_out%x(i)%has_grad_data(this%order))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,this%order))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if(associated(atom_mask_pointer)) then
         call distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors=atoms_in_descriptors,mask=atom_mask_pointer,error=error)
      else
         call distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors=atoms_in_descriptors,error=error)
      endif

      allocate(fcut_pair(this%order,this%order))
      if( my_do_grad_descriptor ) then
         allocate(dfcut_pair(this%order,this%order), directions(3,this%order,this%order))
      endif

      do i_desc = 1, n_descriptors
         if( this%order == 1 ) then
            descriptor_out%x(i_desc)%data = 0.0_dp
            if( my_do_grad_descriptor ) descriptor_out%x(i_desc)%grad_data = 0.0_dp
         else
            i_data = 0
            do ii = 1, this%order
               i = atoms_in_descriptors(1,ii,i_desc)
               s_i = atoms_in_descriptors(2:4,ii,i_desc)
               do jj = ii+1, this%order
                  i_data = i_data + 1
                  j = atoms_in_descriptors(1,jj,i_desc)
                  s_j = atoms_in_descriptors(2:4,jj,i_desc)
                  d_ij = at%pos(:,j) - at%pos(:,i) + matmul(at%lattice,s_j-s_i)
                  r_ij = sqrt(sum(d_ij**2))

                  fcut_pair(jj,ii) = coordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
                  fcut_pair(ii,jj) = fcut_pair(jj,ii)

                  descriptor_out%x(i_desc)%data(i_data) = r_ij
                  if( my_do_grad_descriptor ) then
                     dfcut_pair(ii,jj) = dcoordination_function(r_ij,this%cutoff,this%cutoff_transition_width)
                     dfcut_pair(jj,ii) = dfcut_pair(ii,jj)

                     directions(:,ii,jj) = d_ij / r_ij
                     directions(:,jj,ii) =  - directions(:,ii,jj)
                     descriptor_out%x(i_desc)%grad_data(i_data,:,jj) = directions(:,ii,jj)
                     descriptor_out%x(i_desc)%grad_data(i_data,:,ii) = &
                       - descriptor_out%x(i_desc)%grad_data(i_data,:,jj)
                  endif
               enddo
            enddo

            descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp
            if ( this%compact_clusters ) then

               descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp

               do jj = 2, this%order
                  descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff * fcut_pair(jj,1)
               enddo

               if( my_do_grad_descriptor ) then
                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = 0.0_dp
                  do kk = 2, this%order
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = 1.0_dp
                     do jj = 2, this%order
                        if( jj == kk ) then
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = &
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) * dfcut_pair(jj,1) * (-directions(:,jj,1))
                        else
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = &
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) * fcut_pair(jj,1)
                        endif
                     enddo
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = &
                     descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk)
                  enddo
               endif


            else
               do iConnectivity = 1, size(this%monomerConnectivities,3)

                  fcut_connectivity = 1.0_dp

                  do ii = 1, this%order
                     do jj = ii+1, this%order
                        if( this%monomerConnectivities(jj,ii,iConnectivity) ) then
                           fcut_connectivity = fcut_connectivity * fcut_pair(jj,ii)
                        else
                           fcut_connectivity = fcut_connectivity * ( 1.0_dp - fcut_pair(jj,ii) )
                        endif
                     enddo
                  enddo
                  descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff + fcut_connectivity

                  if( my_do_grad_descriptor ) then
                     do kk = 1, this%order
                        do ll = kk+1, this%order
                           dfcut_connectivity = 1.0_dp
                           do ii = 1, this%order
                              do jj = ii+1, this%order
                                 if( this%monomerConnectivities(jj,ii,iConnectivity) ) then
                                    if( kk == ii .and. ll == jj ) then
                                       dfcut_connectivity = dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    elseif( kk == jj .and. ll == ii ) then
                                       dfcut_connectivity = dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    else
                                       dfcut_connectivity = dfcut_connectivity * fcut_pair(jj,ii)
                                    endif
                                 else
                                    if( kk == ii .and. ll == jj ) then
                                       dfcut_connectivity = - dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    elseif( kk == jj .and. ll == ii) then
                                       dfcut_connectivity = - dfcut_connectivity * dfcut_pair(jj,ii) * directions(:,ll,kk)
                                    else
                                       dfcut_connectivity = dfcut_connectivity * ( 1.0_dp - fcut_pair(jj,ii) )
                                    endif
                                 endif
                              enddo !jj
                           enddo !ii
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,kk) + &
                              dfcut_connectivity
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,ll) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,ll) - &
                              dfcut_connectivity
                        enddo !ll
                     enddo !kk
                  endif

               enddo
            endif

         endif

         descriptor_out%x(i_desc)%ci = atoms_in_descriptors(1,:,i_desc)
         descriptor_out%x(i_desc)%has_data = .true.
         if( my_do_grad_descriptor ) then
            descriptor_out%x(i_desc)%ii = descriptor_out%x(i_desc)%ci
            descriptor_out%x(i_desc)%pos = at%pos(:,descriptor_out%x(i_desc)%ii) + &
               matmul(at%lattice,atoms_in_descriptors(2:4,:,i_desc))
            descriptor_out%x(i_desc)%has_grad_data = .true.
         endif

      enddo

      if(allocated(atoms_in_descriptors)) deallocate(atoms_in_descriptors)
      if(allocated(fcut_pair)) deallocate(fcut_pair)
      if(allocated(dfcut_pair)) deallocate(dfcut_pair)
      if(allocated(directions)) deallocate(directions)

      call system_timer('distance_Nb_calc')

   endsubroutine distance_Nb_calc

   subroutine soap_turbo_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      use soap_turbo_desc

      type(soap_turbo), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, do_timing
      integer :: d, i, j, k, n, i_n, l_n_neighbours, &
         i_desc, n_descriptors, n_cross, n_index, n_atom_pairs
      real(dp) :: r_ij
      real(dp), dimension(3) :: d_ij, u_ij
      real(dp), dimension(:), allocatable :: rjs, thetas, phis, rcut_hard, rcut_soft, nf, global_scaling
      real(dp), dimension(:,:), allocatable :: descriptor_i
      real(dp), dimension(:,:,:), allocatable :: grad_descriptor_i
      integer, dimension(:), allocatable :: species_map
      integer, dimension(3) :: shift_ij
      logical, dimension(:,:), allocatable :: mask

      INIT_ERROR(error)

      call system_timer('soap_turbo_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_turbo_calc: descriptor object not initialised", error)
      endif

!     This is to make the code compatible with the newest TurboGAP (which as multisoap support)
      allocate( rcut_hard(this%n_species) )
      allocate( rcut_soft(this%n_species) )
      allocate( nf(this%n_species) )
      allocate( global_scaling(this%n_species) )
      rcut_hard = this%rcut_hard
      rcut_soft = this%rcut_soft
      nf = this%nf
      global_scaling = 1.0_dp

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

!      allocate(species_map(maxval(this%species_Z)))
      allocate(species_map(1:118))
      species_map = 0
      species_map(this%species_Z) = (/(i, i = 1, this%n_species)/)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)

         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         call param_register(params, 'do_timing', 'F', do_timing, help_string="Do timing or not")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='coordination_calc args_str')) then
            RAISE_ERROR("soap_turbo_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif

         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("soap_turbo_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = soap_turbo_dimensions(this,error)

      allocate(descriptor_i(d,1))

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
        if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         if( at%Z(i) /= this%species_Z(this%central_index) ) cycle

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            allocate(descriptor_out%x(i_desc)%ci(n_index))
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
!           l_n_neighbours = n_neighbours(at,i,max_dist=this%rcut_hard)
            l_n_neighbours = 0
            do n = 1, n_neighbours(at, i)
               j = neighbour(at, i, n, distance = r_ij)
!              The neighbors list past to the soap_turbo library must only contained the "seen" species
               if( r_ij < this%rcut_hard .and. species_map(at%Z(j)) > 0)then
                  l_n_neighbours = l_n_neighbours + 1
               endif
            enddo
            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:l_n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:l_n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif
         if( at%Z(i) /= this%species_Z(this%central_index) ) cycle

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(1) = i
            descriptor_out%x(i_desc)%has_data = .true.
         endif
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

!         n_atom_pairs = n_neighbours(at,i, max_dist = this%rcut_hard) + 1 !Including the central atom
         n_atom_pairs = 1 !Including the central atom
         do n = 1, n_neighbours(at, i)
            j = neighbour(at, i, n, distance = r_ij)
!           The neighbors list past to the soap_turbo library must only contained the "seen" species
            if( r_ij < this%rcut_hard .and. species_map(at%Z(j)) > 0)then
               n_atom_pairs = n_atom_pairs + 1
            endif
         enddo
         allocate( rjs(n_atom_pairs) )
         allocate( thetas(n_atom_pairs) )
         allocate( phis(n_atom_pairs) )
         allocate( mask(n_atom_pairs,this%n_species) )
         mask = .false.

         i_n = 1 ! Start with central atom
         rjs(i_n) = 0.0_dp
         thetas(i_n) = 0.0_dp
         phis(i_n) = 0.0_dp
         mask(i_n,species_map(at%Z(i))) = .true.
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, diff = d_ij, cosines = u_ij)

            if( r_ij >= this%rcut_hard .or. species_map(at%Z(j)) == 0 ) cycle
            i_n = i_n + 1

            rjs(i_n) = r_ij

            thetas(i_n) = dacos( u_ij(3) )
            phis(i_n) = datan2( d_ij(2), d_ij(1) )
            mask(i_n,species_map(at%Z(j))) = .true.
         enddo

         if( my_do_grad_descriptor ) then
            i_n = 0
            do n = 1, n_neighbours(at,i)
               j = neighbour(at, i, n, distance = r_ij, shift = shift_ij)
               if( r_ij >= this%rcut_hard .or. species_map(at%Z(j)) == 0 ) cycle
               i_n = i_n + 1
               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
            enddo
         endif

         descriptor_i = 0.0_dp
         if( my_do_grad_descriptor ) then
            allocate(grad_descriptor_i(3,d,n_atom_pairs))
            grad_descriptor_i = 0.0_dp
         endif
         call get_soap(1, (/n_atom_pairs/), this%n_species, reshape( (/species_map(at%Z(i))/), (/1,1/)), (/1/), &
            n_atom_pairs, mask, rjs, thetas, phis, this%alpha_max, this%l_max, rcut_hard, rcut_soft, nf, &
            global_scaling, this%atom_sigma_r, this%atom_sigma_r_scaling, &
            this%atom_sigma_t, this%atom_sigma_t_scaling, this%amplitude_scaling, this%radial_enhancement, this%central_weight, &
            this%basis, this%scaling_mode, .false., my_do_grad_descriptor, this%compress, this%compress_P_nonzero, this%compress_P_i, &
            this%compress_P_j, this%compress_P_el, descriptor_i, grad_descriptor_i)

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%data = descriptor_i(:,1)
         endif

         if(my_do_grad_descriptor) then
            do k = 1, 3
               descriptor_out%x(i_desc)%grad_data(:,k,0:n_atom_pairs-1) = grad_descriptor_i(k,:,1:n_atom_pairs)
            enddo
         endif

         deallocate(rjs)
         deallocate(thetas)
         deallocate(phis)
         deallocate(mask)
         if(allocated(grad_descriptor_i)) deallocate(grad_descriptor_i)
      enddo

      deallocate(descriptor_i)
      deallocate(species_map)
      deallocate(rcut_hard)
      deallocate(rcut_soft)
      deallocate(nf)
      deallocate(global_scaling)

      call system_timer('soap_turbo_calc')

   endsubroutine soap_turbo_calc

   subroutine distance_Nb_calc_get_clusters(this,at,atoms_in_descriptors,n_descriptors,mask,error)

      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, dimension(:,:,:), intent(out), allocatable, optional :: atoms_in_descriptors
      integer, intent(out), optional :: n_descriptors
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: error

      integer, dimension(:,:,:), allocatable :: my_atoms_in_descriptors

      if( present(atoms_in_descriptors) ) then
         call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors=n_descriptors,mask=mask,error=error)
      else
         call distance_Nb_calc_neighbour_loop(this,at,my_atoms_in_descriptors,n_descriptors=n_descriptors,mask=mask,error=error)
         if(allocated(my_atoms_in_descriptors)) deallocate(my_atoms_in_descriptors)
      endif

   endsubroutine distance_Nb_calc_get_clusters

!   recursive subroutine distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors,error)
!
!      type(distance_Nb), intent(in) :: this
!      type(atoms), intent(in) :: at
!      integer, dimension(:,:,:), intent(inout), allocatable :: atoms_in_descriptors
!      integer, intent(out), optional :: n_descriptors
!      integer, intent(out), optional :: error
!
!      integer, save :: current_order = 0
!      integer :: i, j, n, order, i_desc, d
!      real(dp) :: r_ij
!      integer, dimension(3) :: shift_i, shift_j, shift_ij
!      integer, dimension(:,:), allocatable :: current_descriptor
!
!      type(LinkedList_i2d), pointer :: LL_atoms_in_descriptors => null()
!
!      INIT_ERROR(error)
!
!      current_order = current_order + 1
!
!      if( current_order == 1 ) then
!         allocate(current_descriptor(4,1))
!
!         do i = 1, at%N
!            if( any( at%Z(i) == this%Z ) .or. any( 0 == this%Z ) ) then
!               current_descriptor(:,1) = (/i,0,0,0/)
!               call append(LL_atoms_in_descriptors,current_descriptor,error)
!            endif
!         enddo
!
!         deallocate(current_descriptor)
!         call retrieve(LL_atoms_in_descriptors,atoms_in_descriptors)
!         call finalise(LL_atoms_in_descriptors)
!         if( this%order > 1 ) &
!            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
!
!         if( present(n_descriptors) ) n_descriptors = size(atoms_in_descriptors,3)
!      else
!         if( .not. allocated(atoms_in_descriptors) ) then
!            RAISE_ERROR("distance_Nb_calc_neighbour_loop: atoms_in_descriptors must be allocated",error)
!         endif
!
!         allocate(current_descriptor(4,current_order))
!         do i_desc = 1, size(atoms_in_descriptors,3)
!            do order = 1, size(atoms_in_descriptors,2)
!               i = atoms_in_descriptors(1,order,i_desc)
!               shift_i = atoms_in_descriptors(2:4,order,i_desc)
!               loop_n: do n = 1, n_neighbours(at,i)
!                  j = neighbour(at,i,n,distance = r_ij, shift = shift_ij)
!
!                  if( r_ij > this%cutoff ) cycle
!                  if( .not. is_subset(this%Z, at%Z( (/j,atoms_in_descriptors(1,:,i_desc)/) ), error) .and. all(this%Z /= 0) ) cycle
!
!                  shift_j = shift_ij + shift_i
!
!                  current_descriptor(:,1:current_order-1) = atoms_in_descriptors(:,:,i_desc)
!                  current_descriptor(:,current_order) = (/j, shift_j/)
!                  if( order_and_check_for_duplicates(current_descriptor,at) ) then
!                     do d = current_order, 1, -1
!                        current_descriptor(2:4,d) = current_descriptor(2:4,d) - current_descriptor(2:4,1)
!                     enddo
!                     if( .not. is_in_LinkedList(LL_atoms_in_descriptors,current_descriptor,error) ) &
!                       call append(LL_atoms_in_descriptors,current_descriptor,error)
!                  endif
!               enddo loop_n
!            enddo
!         enddo
!
!         deallocate(current_descriptor)
!         call retrieve(LL_atoms_in_descriptors,atoms_in_descriptors)
!         call finalise(LL_atoms_in_descriptors)
!         if( current_order < this%order ) &
!            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
!      endif
!
!      current_order = current_order - 1
!
!   endsubroutine distance_Nb_calc_neighbour_loop

   recursive subroutine distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors,n_descriptors,mask,error)

      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, dimension(:,:,:), intent(inout), allocatable :: atoms_in_descriptors
      integer, intent(out), optional :: n_descriptors
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: error

      integer, save :: current_order = 0
      integer :: i, j, n, order, i_desc, d
      real(dp) :: r_ij
      integer, dimension(3) :: shift_i, shift_j, shift_ij
      integer, dimension(:,:), allocatable :: current_descriptor

      type(Table)  :: Table_atoms_in_descriptors, Table_atoms_in_descriptors_uniq

      INIT_ERROR(error)

      current_order = current_order + 1

      if( current_order == 1 ) then
         call initialise(Table_atoms_in_descriptors, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)
         allocate(current_descriptor(4,1))

         do i = 1, at%N
            if( any( at%Z(i) == this%Z ) .or. any( 0 == this%Z ) ) then
               if( present(mask) ) then
                  if( .not. mask(i) ) cycle
               endif

               current_descriptor(:,1) = (/i,0,0,0/)
               call append(Table_atoms_in_descriptors,current_descriptor(:,1))
            endif
         enddo

         deallocate(current_descriptor)

         allocate(atoms_in_descriptors(4,1,Table_atoms_in_descriptors%N))
         atoms_in_descriptors = reshape(Table_atoms_in_descriptors%int(:,1:Table_atoms_in_descriptors%N),(/4,1,Table_atoms_in_descriptors%N/))

         call finalise(Table_atoms_in_descriptors)

         if( this%order > 1 ) &
            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)

         if( present(n_descriptors) ) n_descriptors = size(atoms_in_descriptors,3)
      else
         if( .not. allocated(atoms_in_descriptors) ) then
            RAISE_ERROR("distance_Nb_calc_neighbour_loop: atoms_in_descriptors must be allocated",error)
         endif

         call initialise(Table_atoms_in_descriptors, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)
         allocate(current_descriptor(4,current_order))

         do i_desc = 1, size(atoms_in_descriptors,3)
            do order = 1, merge(1,size(atoms_in_descriptors,2),this%compact_clusters) !size(atoms_in_descriptors,2)
               ! if compact_clusters == T, only neighbours of the first (central) atom is considered
               i = atoms_in_descriptors(1,order,i_desc)
               shift_i = atoms_in_descriptors(2:4,order,i_desc)
               loop_n: do n = 1, n_neighbours(at,i)
                  j = neighbour(at,i,n,distance = r_ij, shift = shift_ij)

                  if( r_ij > this%cutoff ) cycle
                  if( .not. is_subset(this%Z, at%Z( (/j,atoms_in_descriptors(1,:,i_desc)/) ), error) .and. all(this%Z /= 0) ) cycle

                  shift_j = shift_ij + shift_i

                  current_descriptor(:,1:current_order-1) = atoms_in_descriptors(:,:,i_desc)
                  current_descriptor(:,current_order) = (/j, shift_j/)
                  if( order_and_check_for_duplicates(current_descriptor(:,merge(2,1,this%compact_clusters):),at) ) then
                     ! if compact_clusters == T, leave first atom alone
                     do d = current_order, 1, -1
                        current_descriptor(2:4,d) = current_descriptor(2:4,d) - current_descriptor(2:4,1)
                     enddo
                     call append(Table_atoms_in_descriptors,reshape(current_descriptor,(/4*current_order/)))

                     !if( .not. is_in_LinkedList(LL_atoms_in_descriptors,current_descriptor,error) ) &
                     !  call append(LL_atoms_in_descriptors,current_descriptor,error)
                  endif
               enddo loop_n
            enddo
         enddo

         deallocate(current_descriptor,atoms_in_descriptors)
         call initialise(Table_atoms_in_descriptors_uniq, Nint = 4*current_order, Nreal = 0, Nstr = 0, Nlogical = 0, error=error)

         if( Table_atoms_in_descriptors%N > 0 ) then
            call heap_sort(Table_atoms_in_descriptors%int(:,1:Table_atoms_in_descriptors%N))
            call append(Table_atoms_in_descriptors_uniq,Table_atoms_in_descriptors%int(:,1))
            do i_desc = 2, Table_atoms_in_descriptors%N
               if( .not. all( Table_atoms_in_descriptors%int(:,i_desc) == Table_atoms_in_descriptors%int(:,i_desc-1) ) ) &
                   call append(Table_atoms_in_descriptors_uniq,Table_atoms_in_descriptors%int(:,i_desc))
            enddo
         endif

         allocate(atoms_in_descriptors(4,current_order,Table_atoms_in_descriptors_uniq%N))
         atoms_in_descriptors = reshape(Table_atoms_in_descriptors_uniq%int(:,1:Table_atoms_in_descriptors_uniq%N),(/4,current_order,Table_atoms_in_descriptors_uniq%N/))

         call finalise(Table_atoms_in_descriptors)
         call finalise(Table_atoms_in_descriptors_uniq)

         if( current_order < this%order ) &
            call distance_Nb_calc_neighbour_loop(this,at,atoms_in_descriptors = atoms_in_descriptors,n_descriptors=n_descriptors,error=error)
      endif

      current_order = current_order - 1

   endsubroutine distance_Nb_calc_neighbour_loop

   function order_and_check_for_duplicates(array,at)
      integer, dimension(:,:), intent(inout) :: array
      type(atoms), intent(in) :: at
      logical :: order_and_check_for_duplicates

      integer :: ii, jj, n
      integer, dimension(size(array,1)) :: tmp
      logical :: do_swap

      integer, dimension(size(array,1)+1,size(array,2)) :: Z_array

      Z_array(1,:) = at%Z(array(1,:))
      Z_array(2:,:) = array(:,:)

      call heap_sort(Z_array)

      do ii = 2, size(Z_array,2)
         if( all( Z_array(:,ii-1) == Z_array(:,ii) ) ) then
            order_and_check_for_duplicates = .false.
            return
         endif
      enddo

      array(:,:) = Z_array(2:,:)

      order_and_check_for_duplicates = .true.

   endfunction order_and_check_for_duplicates

   function is_subset(set,subset,error)
      logical :: is_subset
      integer, dimension(:), intent(in) :: set, subset
      integer, optional, intent(out) :: error

      logical, dimension(size(set)) :: found
      integer :: i, j

      INIT_ERROR(error)
      if( size(set) < size(subset) ) then
         RAISE_ERROR("is_subset: size of set must be greater than or equal to the size of subset",error)
      endif

      found = .false.
      loop_i: do i = 1, size(subset)
         do j = 1, size(set)
            if(set(j) == subset(i) .and. .not. found(j)) then
               found(j) = .true.
               cycle loop_i
            endif
         enddo
      enddo loop_i

      is_subset = ( count(found) == size(subset) )

   endfunction is_subset


   function descriptor_dimensions(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: descriptor_dimensions

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_dimensions = bispectrum_SO4_dimensions(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_dimensions = bispectrum_SO3_dimensions(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_dimensions = behler_dimensions(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_dimensions = distance_2b_dimensions(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_dimensions = coordination_dimensions(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_dimensions = angle_3b_dimensions(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_dimensions = co_angle_3b_dimensions(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_dimensions = co_distance_2b_dimensions(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_dimensions = cosnx_dimensions(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_dimensions = trihis_dimensions(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_dimensions = water_monomer_dimensions(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_dimensions = water_dimer_dimensions(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_dimensions = A2_dimer_dimensions(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_dimensions = AB_dimer_dimensions(this%descriptor_AB_dimer,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_dimensions = atom_real_space_dimensions(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_dimensions = power_so3_dimensions(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_dimensions = power_so4_dimensions(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_dimensions = soap_dimensions(this%descriptor_soap,error)
         case(DT_RDF)
            descriptor_dimensions = rdf_dimensions(this%descriptor_rdf,error)
         case(DT_AS_DISTANCE_2b)
            descriptor_dimensions = as_distance_2b_dimensions(this%descriptor_as_distance_2b,error)
         case(DT_ALEX)
            descriptor_dimensions = alex_dimensions(this%descriptor_alex,error)
         case(DT_DISTANCE_Nb)
            descriptor_dimensions = distance_Nb_dimensions(this%descriptor_distance_Nb,error)
         case(DT_SOAP_TURBO)
            descriptor_dimensions = soap_turbo_dimensions(this%descriptor_soap_turbo,error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_BOND_REAL_SPACE)
            descriptor_dimensions = bond_real_space_dimensions(this%descriptor_bond_real_space,error)
         case(DT_AN_MONOMER)
            descriptor_dimensions = AN_monomer_dimensions(this%descriptor_AN_monomer,error)
         case(DT_COM_DIMER)
            descriptor_dimensions = com_dimer_dimensions(this%descriptor_com_dimer,error)
         case(DT_GENERAL_MONOMER)
            descriptor_dimensions = general_monomer_dimensions(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            descriptor_dimensions = general_dimer_dimensions(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            descriptor_dimensions = general_trimer_dimensions(this%descriptor_general_trimer,error)
         case(DT_MOLECULE_LO_D)
            descriptor_dimensions = molecule_lo_d_dimensions(this%descriptor_molecule_lo_d,error)
         case(DT_SOAP_EXPRESS)
            descriptor_dimensions = soap_express_dimensions(this%descriptor_soap_express,error)
#endif
         case default
            RAISE_ERROR("descriptor_dimensions: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_dimensions

   function bispectrum_SO4_dimensions(this,error) result(i)
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: j, j1, j2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do j1 = 0, this%j_max
         j2 = j1
         !do j2 = 0, this%j_max
            do j = abs(j1-j2), min(this%j_max,j1+j2)
               if( mod(j1+j2+j,2) == 1 ) cycle
               i = i + 1
            enddo
         !enddo
      enddo

   endfunction bispectrum_SO4_dimensions

   function bispectrum_SO3_dimensions(this,error) result(i)
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: a, l1, l2, l

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do a = 1, this%n_max
         do l1 = 0, this%l_max
            l2 = l1
            !do l2 = 0, this%l_max
               do l = abs(l1-l2), min(this%l_max,l1+l2)
                  if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                  i = i + 1
               enddo
            !enddo
         enddo
      enddo

   endfunction bispectrum_SO3_dimensions

   function behler_dimensions(this,error) result(i)
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_g2 + this%n_g3

   endfunction behler_dimensions

   function distance_2b_dimensions(this,error) result(i)
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_exponents

   endfunction distance_2b_dimensions

   function coordination_dimensions(this,error) result(i)
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction coordination_dimensions

   function angle_3b_dimensions(this,error) result(i)
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction angle_3b_dimensions

   function co_angle_3b_dimensions(this,error) result(i)
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 4

   endfunction co_angle_3b_dimensions

   function co_distance_2b_dimensions(this,error) result(i)
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction co_distance_2b_dimensions

   function cosnx_dimensions(this,error) result(i)
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction cosnx_dimensions

   function trihis_dimensions(this,error) result(i)
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_gauss

   endfunction trihis_dimensions

   function water_monomer_dimensions(this,error) result(i)
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction water_monomer_dimensions

   function water_dimer_dimensions(this,error) result(i)
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 15

   endfunction water_dimer_dimensions

   function A2_dimer_dimensions(this,error) result(i)
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction A2_dimer_dimensions

   function AB_dimer_dimensions(this,error) result(i)
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction AB_dimer_dimensions


   function atom_real_space_dimensions(this,error) result(i)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 2 * (this%l_max+1)**2 + 2

   endfunction atom_real_space_dimensions

   function power_so3_dimensions(this,error) result(i)
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction power_so3_dimensions

   function power_SO4_dimensions(this,error) result(i)
      type(power_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = this%j_max + 1

   endfunction power_SO4_dimensions

   function soap_dimensions(this,error) result(i)
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i, K1, K2
      logical :: sym_desc
      type(real_2d), dimension(:), allocatable :: W

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_dimensions: descriptor object not initialised", error)
      endif

      call form_W(this, W, sym_desc, error)
      K1 = size(W(1)%mm(0,:))
      K2 = size(W(2)%mm(0,:))

      if (this%diagonal_radial) then
         if (this%Z_mix .or. this%R_mix .or. this%nu_R /= 2 .or. this%nu_S /= 2 .or. (.not. this%coupling)) then
            RAISE_ERROR("soap_dimensions: can't combine diagonal radial with any other compression strategies", error)
         endif
         i = (this%l_max+1) * this%n_max * this%n_species * (this%n_species+1) / 2 + 1
      elseif (this%coupling) then
         if (sym_desc) then
            i = (this%l_max+1) * (K1 * (K1+1)) /2 + 1
         else
            i = (this%l_max+1) * K1 * K2 + 1
         endif
      else
         if (K1 /= K2) then
            RAISE_ERROR("require K1=K2 to use elementwise coupling", error)
         endif
         i = K1 * (this%l_max + 1) + 1
      endif

      if (allocated(W)) deallocate(W)
   endfunction soap_dimensions



   function rdf_dimensions(this,error) result(i)
      type(rdf), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_gauss

   endfunction rdf_dimensions

   function as_distance_2b_dimensions(this,error) result(i)
      type(as_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction as_distance_2b_dimensions


   function alex_dimensions(this,error) result(i)
      type(alex), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i, nradial

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_dimensions: descriptor object not initialised", error)
      endif

      nradial = this%power_max-this%power_min + 1
      i = nradial+2*nradial**2+nradial**3

   endfunction alex_dimensions

   function distance_Nb_dimensions(this,error) result(i)
      type(distance_Nb), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_dimensions: descriptor object not initialised", error)
      endif

      i = max(1,this%order * ( this%order - 1 ) / 2)

   endfunction distance_Nb_dimensions

   function descriptor_cutoff(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: descriptor_cutoff

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_cutoff = cutoff(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_cutoff = cutoff(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_cutoff = cutoff(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_cutoff = cutoff(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_AB_dimer,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_cutoff = cutoff(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_cutoff = cutoff(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_cutoff = cutoff(this%descriptor_soap,error)
         case(DT_RDF)
            descriptor_cutoff = cutoff(this%descriptor_rdf,error)
         case(DT_ALEX)
            descriptor_cutoff = cutoff(this%descriptor_alex,error)
         case(DT_DISTANCE_Nb)
            descriptor_cutoff = cutoff(this%descriptor_distance_Nb,error)
         case(DT_SOAP_TURBO)
            descriptor_cutoff = cutoff(this%descriptor_soap_turbo,error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_BOND_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_bond_real_space,error)
         case(DT_MOLECULE_LO_D)
            descriptor_cutoff = cutoff(this%descriptor_molecule_lo_d,error)
         case(DT_AN_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_AN_monomer,error)
         case(DT_GENERAL_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_general_monomer,error)
         case(DT_GENERAL_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_general_dimer,error)
         case(DT_GENERAL_TRIMER)
            descriptor_cutoff = cutoff(this%descriptor_general_trimer,error)
         case(DT_COM_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_com_dimer,error)
         case(DT_SOAP_EXPRESS)
            descriptor_cutoff = cutoff(this%descriptor_soap_express,error)
#endif
         case default
            RAISE_ERROR("descriptor_cutoff: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_cutoff

   function bispectrum_SO4_cutoff(this,error)
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO4_cutoff = this%cutoff

   endfunction bispectrum_SO4_cutoff

   function bispectrum_SO3_cutoff(this,error)
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO3_cutoff = this%cutoff

   endfunction bispectrum_SO3_cutoff

   function behler_cutoff(this,error)
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: behler_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_cutoff: descriptor object not initialised", error)
      endif

      behler_cutoff = this%cutoff

   endfunction behler_cutoff

   function distance_2b_cutoff(this,error)
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_cutoff: descriptor object not initialised", error)
      endif

      distance_2b_cutoff = this%cutoff

   endfunction distance_2b_cutoff

   function co_distance_2b_cutoff(this,error)
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_cutoff: descriptor object not initialised", error)
      endif

      co_distance_2b_cutoff = this%cutoff

   endfunction co_distance_2b_cutoff

   function coordination_cutoff(this,error)
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: coordination_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_cutoff: descriptor object not initialised", error)
      endif

      coordination_cutoff = this%cutoff

   endfunction coordination_cutoff

   function angle_3b_cutoff(this,error)
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_cutoff: descriptor object not initialised", error)
      endif

      angle_3b_cutoff = this%cutoff

   endfunction angle_3b_cutoff

   function co_angle_3b_cutoff(this,error)
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_cutoff: descriptor object not initialised", error)
      endif

      co_angle_3b_cutoff = this%cutoff

   endfunction co_angle_3b_cutoff

   function cosnx_cutoff(this,error)
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: cosnx_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_cutoff: descriptor object not initialised", error)
      endif

      cosnx_cutoff = this%cutoff

   endfunction cosnx_cutoff

   function trihis_cutoff(this,error)
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: trihis_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_cutoff: descriptor object not initialised", error)
      endif

      trihis_cutoff = this%cutoff

   endfunction trihis_cutoff

   function water_monomer_cutoff(this,error)
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_cutoff: descriptor object not initialised", error)
      endif

      water_monomer_cutoff = this%cutoff

   endfunction water_monomer_cutoff

   function water_dimer_cutoff(this,error)
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_cutoff: descriptor object not initialised", error)
      endif

      water_dimer_cutoff = this%cutoff

   endfunction water_dimer_cutoff

   function A2_dimer_cutoff(this,error)
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: A2_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_cutoff: descriptor object not initialised", error)
      endif

      A2_dimer_cutoff = this%cutoff

   endfunction A2_dimer_cutoff

   function AB_dimer_cutoff(this,error)
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: AB_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_cutoff: descriptor object not initialised", error)
      endif

      AB_dimer_cutoff = this%cutoff

   endfunction AB_dimer_cutoff


   function atom_real_space_cutoff(this,error)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: atom_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_cutoff: descriptor object not initialised", error)
      endif

      atom_real_space_cutoff = this%cutoff

   endfunction atom_real_space_cutoff

   function power_so3_cutoff(this,error)
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_cutoff: descriptor object not initialised", error)
      endif

      power_so3_cutoff = this%cutoff

   endfunction power_so3_cutoff

   function power_so4_cutoff(this,error)
      type(power_so4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so4_cutoff: descriptor object not initialised", error)
      endif

      power_so4_cutoff = this%cutoff

   endfunction power_so4_cutoff

   function soap_cutoff(this,error)
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: soap_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_cutoff: descriptor object not initialised", error)
      endif

      soap_cutoff = this%cutoff

   endfunction soap_cutoff


   function rdf_cutoff(this,error)
      type(rdf), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: rdf_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_cutoff: descriptor object not initialised", error)
      endif

      rdf_cutoff = this%cutoff

   endfunction rdf_cutoff

   function as_distance_2b_cutoff(this,error)
      type(as_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: as_distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_cutoff: descriptor object not initialised", error)
      endif

      as_distance_2b_cutoff = this%max_cutoff

   endfunction as_distance_2b_cutoff


   function alex_cutoff(this,error)
      type(alex), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: alex_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_cutoff: descriptor object not initialised", error)
      endif

      alex_cutoff = this%cutoff

   endfunction alex_cutoff

   function distance_Nb_cutoff(this,error)
      type(distance_Nb), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: distance_Nb_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_cutoff: descriptor object not initialised", error)
      endif

      distance_Nb_cutoff = this%cutoff

   endfunction distance_Nb_cutoff

   function soap_turbo_cutoff(this,error)
      type(soap_turbo), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: soap_turbo_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_turbo_cutoff: descriptor object not initialised", error)
      endif

      soap_turbo_cutoff = this%rcut_hard

   endfunction soap_turbo_cutoff


   subroutine descriptor_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call bispectrum_SO4_sizes(this%descriptor_bispectrum_SO4,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_BISPECTRUM_SO3)
            call bispectrum_SO3_sizes(this%descriptor_bispectrum_SO3,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_BEHLER)
            call behler_sizes(this%descriptor_behler,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_DISTANCE_2b)
            call distance_2b_sizes(this%descriptor_distance_2b,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_COORDINATION)
            call coordination_sizes(this%descriptor_coordination,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_ANGLE_3B)
            call angle_3b_sizes(this%descriptor_angle_3b,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_CO_ANGLE_3B)
            call co_angle_3b_sizes(this%descriptor_co_angle_3b,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_CO_DISTANCE_2b)
            call co_distance_2b_sizes(this%descriptor_co_distance_2b,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_COSNX)
            call cosnx_sizes(this%descriptor_cosnx,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_TRIHIS)
            call trihis_sizes(this%descriptor_trihis,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_WATER_MONOMER)
            call water_monomer_sizes(this%descriptor_water_monomer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_WATER_DIMER)
            call water_dimer_sizes(this%descriptor_water_dimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_A2_DIMER)
            call A2_dimer_sizes(this%descriptor_A2_dimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_AB_DIMER)
            call AB_dimer_sizes(this%descriptor_AB_dimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_ATOM_REAL_SPACE)
            call atom_real_space_sizes(this%descriptor_atom_real_space,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_POWER_SO3)
            call power_so3_sizes(this%descriptor_power_so3,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_POWER_SO4)
            call power_so4_sizes(this%descriptor_power_so4,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_SOAP)
            call soap_sizes(this%descriptor_soap,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_RDF)
            call rdf_sizes(this%descriptor_rdf,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_AS_DISTANCE_2b)
            call as_distance_2b_sizes(this%descriptor_as_distance_2b,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_ALEX)
            call alex_sizes(this%descriptor_alex,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_DISTANCE_Nb)
            call distance_Nb_sizes(this%descriptor_distance_Nb,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_SOAP_TURBO)
            call soap_turbo_sizes(this%descriptor_soap_turbo,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_BOND_REAL_SPACE)
            call bond_real_space_sizes(this%descriptor_bond_real_space,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_MOLECULE_LO_D)
            call molecule_lo_d_sizes(this%descriptor_molecule_lo_d,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_AN_MONOMER)
            call AN_monomer_sizes(this%descriptor_AN_monomer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_GENERAL_MONOMER)
            call general_monomer_sizes(this%descriptor_general_monomer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_GENERAL_DIMER)
            call general_dimer_sizes(this%descriptor_general_dimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_GENERAL_TRIMER)
            call general_trimer_sizes(this%descriptor_general_trimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_COM_DIMER)
            call com_dimer_sizes(this%descriptor_com_dimer,at, &
               n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
         case(DT_SOAP_EXPRESS)
            call soap_express_sizes(this%descriptor_soap_express,at, &
                 n_descriptors,n_cross,mask=mask,n_index=n_index,error=error)
#endif
         case default
            RAISE_ERROR("descriptor_sizes: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_sizes

   subroutine bispectrum_SO4_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine bispectrum_SO4_sizes

   subroutine bispectrum_SO3_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(bispectrum_SO3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine bispectrum_SO3_sizes

   subroutine behler_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         if( this%Z /= 0 .and. this%Z /= at%Z(i) ) cycle

         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine behler_sizes

   subroutine distance_2b_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2
      real(dp) :: r_ij

      logical :: needs_resid
      integer, dimension(:), pointer :: resid_pointer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_sizes: descriptor object not initialised", error)
      endif

      needs_resid = this%only_intra .or. this%only_inter
      if (needs_resid) then
         if (.not. assign_pointer(at, trim(this%resid_name), resid_pointer)) then
            RAISE_ERROR("distance_2b_sizes did not find "//trim(this%resid_name)//" property (residue id) in the atoms object.", error)
         end if
      else
         resid_pointer => null()
      end if

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij)
            if(r_ij >= this%cutoff) cycle
!if(r_ij < 3.5_dp) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            if (needs_resid) then
               if (this%only_intra .and. resid_pointer(i) /= resid_pointer(j)) cycle
               if (this%only_inter .and. resid_pointer(i) == resid_pointer(j)) cycle
            end if

            n_descriptors = n_descriptors + 1
         enddo
      enddo

      n_cross = n_descriptors*2

      if( present(n_index) ) n_index = 2

   endsubroutine distance_2b_sizes

   subroutine coordination_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine coordination_sizes

   subroutine angle_3b_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle

               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
            enddo
         enddo
      enddo
      n_cross = n_descriptors * 3

      if( present(n_index) ) n_index = 1

   endsubroutine angle_3b_sizes

   subroutine co_angle_3b_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m, n_neighbours_coordination
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif

         n_neighbours_coordination = n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)
            if( r_ij >= this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, n_neighbours(at,i)
               if( n == m ) cycle
               k = neighbour(at, i, m, distance = r_ik)
               if( r_ik >= this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 3 + n_neighbours_coordination
            enddo
         enddo
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine co_angle_3b_sizes

   subroutine co_distance_2b_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      real(dp) :: r_ij
      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij)
            if( r_ij >= this%cutoff ) cycle
!if( r_ij < 3.5_dp ) cycle


            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 4 + n_neighbours(at,i,max_dist=this%coordination_cutoff) + n_neighbours(at,j,max_dist=this%coordination_cutoff)
         enddo
      enddo

      if( present(n_index) ) n_index = 2

   endsubroutine co_distance_2b_sizes

   subroutine cosnx_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine cosnx_sizes

   subroutine trihis_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N

      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_cross = n_cross + n_neighbours(at,i) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine trihis_sizes

   subroutine water_monomer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(at%Z(i) == 8) then
            if(present(mask)) then
               if(.not. mask(i)) cycle
            endif
            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 3
         endif
      enddo

      if( present(n_index) ) n_index = 3

   endsubroutine water_monomer_sizes

   subroutine water_dimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, n
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
call print("mask present ? "//present(mask))
      do i = 1, at%N
         if(at%Z(i) == 8) then
            if(present(mask)) then
               if(.not. mask(i)) cycle
            endif
            do n = 1, n_neighbours(at,i)
               j = neighbour(at,i,n,distance=r_ij)
               if(at%Z(j) == 8 .and. r_ij < this%cutoff) then
                  n_descriptors = n_descriptors + 1
                  n_cross = n_cross + 6
               endif
            enddo
         endif
      enddo

      if( present(n_index) ) n_index = 6

   endsubroutine water_dimer_sizes

   subroutine A2_dimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, iA1, iA2, iB1, iB2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         iA1 = i
         iA2 = neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A1_B2 = distance_min_image(at,iA1,iB2)

            r_A2_B1 = distance_min_image(at,iA2,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)

            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

      if( present(n_index) ) n_index = 4

   endsubroutine A2_dimer_sizes

   subroutine AB_dimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, n_monomers, iA1, iA2, iB1, iB2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_sizes: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)

            r_A1_A2 = distance_min_image(at,iA1,iA2)
            r_B1_B2 = distance_min_image(at,iB1,iB2)

            r_A1_B2 = distance_min_image(at,iA1,iB2)
            r_A2_B1 = distance_min_image(at,iA2,iB1)

            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

      deallocate(AB_monomer_index)

      if( present(n_index) ) n_index = 4

   endsubroutine AB_dimer_sizes


   subroutine atom_real_space_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff)*2
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine atom_real_space_sizes

   subroutine power_so3_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine power_so3_sizes

   subroutine power_SO4_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine power_SO4_sizes

   subroutine soap_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( .not. any( at%Z(i) == this%Z ) .and. .not. any(this%Z == 0) ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if(this%global) then
         n_descriptors = 1
         if( present(n_index) ) then
            if( any(this%Z == 0) ) then
               n_index = at%N
            else
               n_index = count( (/(any(at%Z(i)==this%Z),i=1,at%N)/) )
            endif
         endif
      else
         if( present(n_index) ) n_index = 1
      endif

   endsubroutine soap_sizes

   subroutine rdf_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(rdf), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("rdf_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine rdf_sizes

   subroutine as_distance_2b_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(as_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      real(dp) :: r_ij
      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("as_distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij)
            if( r_ij > this%max_cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 4 + n_neighbours(at,i,max_dist=this%coordination_cutoff) + n_neighbours(at,j,max_dist=this%coordination_cutoff)
         enddo
      enddo

      if( present(n_index) ) n_index = 2

   endsubroutine as_distance_2b_sizes


   subroutine alex_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(alex), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("alex_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine alex_sizes

   subroutine distance_Nb_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(distance_Nb), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_Nb_sizes: descriptor object not initialised", error)
      endif

      call distance_Nb_calc_get_clusters(this,at,n_descriptors=n_descriptors,mask=mask,error=error)
      n_cross = n_descriptors * this%order

      if( present(n_index) ) n_index = this%order

   endsubroutine distance_Nb_sizes

   function soap_turbo_dimensions(this,error) result(i)
      type(soap_turbo), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: n_max

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_turbo_dimensions: descriptor object not initialised", error)
      endif

      if( this%compress )then
         i = maxval(this%compress_P_i)
      else
         n_max = sum(this%alpha_max)
         i = ( this%l_max+1 ) * ( n_max*(n_max+1) ) / 2
      endif

   endfunction soap_turbo_dimensions

   function descriptor_n_permutations(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error

      integer :: descriptor_n_permutations, i

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,&
            DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP,DT_RDF, DT_ALEX, DT_COM_DIMER, &
            DT_SOAP_EXPRESS,DT_SOAP_TURBO)

            descriptor_n_permutations = 1

         case(DT_WATER_DIMER)
            descriptor_n_permutations = NP_WATER_DIMER
         case(DT_A2_DIMER)
            descriptor_n_permutations = NP_A2_DIMER
         case(DT_AB_DIMER)
            descriptor_n_permutations = NP_AB_DIMER
#ifdef DESCRIPTORS_NONCOMMERCIAL
         case(DT_AN_MONOMER)
            if(this%descriptor_AN_monomer%do_atomic) then
               descriptor_n_permutations = factorial(this%descriptor_AN_monomer%N-1)
            else
               descriptor_n_permutations = factorial(this%descriptor_AN_monomer%N)
            endif
         case(DT_GENERAL_MONOMER)
            if (.not. this%descriptor_general_monomer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_monomer%permutation_data%n_perms
         case(DT_GENERAL_DIMER)
            if (.not. this%descriptor_general_dimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_dimer%permutation_data%n_perms
         case(DT_GENERAL_TRIMER)
            if (.not. this%descriptor_general_trimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_general_trimer%permutation_data%n_perms
         case(DT_MOLECULE_LO_D)
            if (.not. this%descriptor_molecule_lo_d%permutation_data%initialised)then
              RAISE_ERROR("descriptor_n_permutations: permutation_data not initialised "//this%descriptor_type,error)
            end if
            descriptor_n_permutations = this%descriptor_molecule_lo_d%permutation_data%n_perms
#endif
         case(DT_DISTANCE_NB)
            descriptor_n_permutations = this%descriptor_distance_Nb%n_permutations
         case default
            RAISE_ERROR("descriptor_n_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_n_permutations

   subroutine descriptor_permutations(this,permutations,error)
      type(descriptor), intent(in) :: this
#ifdef DESCRIPTORS_NONCOMMERCIAL
      type(permutation_data_type) :: my_permutation_data
#endif
      integer, dimension(:,:), intent(out) :: permutations
      integer, optional, intent(out) :: error

      integer :: i, d, np, n, m, ip, j
      integer,dimension(1) :: unit_vec
      integer, dimension(:), allocatable :: this_perm
      integer, dimension(:,:), allocatable :: distance_matrix, atom_permutations, sliced_permutations

      INIT_ERROR(error)

      d = descriptor_dimensions(this,error)
      np = descriptor_n_permutations(this,error)
      call check_size('permutations',permutations, (/d,np/),'descriptor_permutations',error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,&
            DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP,DT_RDF, DT_ALEX, DT_COM_DIMER,&
            DT_SOAP_EXPRESS,DT_SOAP_TURBO)

            permutations(:,1) = (/ (i, i = 1, size(permutations,1)) /)
         case(DT_WATER_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/) ! original order
            permutations(:,2) = (/1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 13, 14, 11, 12, 15/) ! swap Hs on monomer A
            permutations(:,3) = (/1, 2, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 14, 13, 15/) ! swap Hs on monomer B
            permutations(:,4) = (/1, 3, 2, 5, 4, 7, 6, 9, 8, 10, 14, 13, 12, 11, 15/) ! swap Hs on both monomers
            permutations(:,5) = (/1, 8, 9, 6, 7, 4, 5, 2, 3, 15, 11, 13, 12, 14, 10/) ! swap monomers A and B
            permutations(:,6) = (/1, 9, 8, 6, 7, 5, 4, 2, 3, 15, 12, 14, 11, 13, 10/) ! swap monomers and Hs on monomer A
            permutations(:,7) = (/1, 8, 9, 7, 6, 4, 5, 3, 2, 15, 13, 11, 14, 12, 10/) ! swap monomers and Hs on monomer B
            permutations(:,8) = (/1, 9, 8, 7, 6, 5, 4, 3, 2, 15, 14, 12, 13, 11, 10/) ! swap monomers and Hs on both monomers

         case(DT_A2_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/1, 2, 5, 6, 3, 4/) ! swap atoms on monomer A
            permutations(:,3) = (/1, 2, 4, 3, 6, 5/) ! swap atoms on monomer B
            permutations(:,4) = (/1, 2, 6, 5, 4, 3/) ! swap atoms on both monomers
            permutations(:,5) = (/2, 1, 3, 5, 4, 6/) ! swap monomers A and B
            permutations(:,6) = (/2, 1, 5, 3, 6, 4/) ! swap monomers and atoms on monomer A
            permutations(:,7) = (/2, 1, 4, 6, 3, 5/) ! swap monomers and atoms on monomer B
            permutations(:,8) = (/2, 1, 6, 4, 5, 3/) ! swap monomers and atoms on both monomers

         case(DT_AB_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/2, 1, 3, 4, 6, 5/) ! swap monomers
#ifdef DESCRIPTORS_NONCOMMERCIAL
#include "descriptors_noncommercial_permutations.inc"
#endif
         case(DT_DISTANCE_NB)
            permutations = this%descriptor_distance_Nb%permutations
         case default
            RAISE_ERROR("descriptor_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_permutations


   subroutine real_space_fourier_coefficients(at,l_max,atom_coefficient)
      type(atoms), intent(in) :: at
      integer, intent(in) :: l_max
      type(neighbour_type), dimension(:), allocatable :: atom_coefficient

      integer :: i, j, n, l, m
      real(dp) :: r
      real(dp), dimension(3) :: d

      if(.not.allocated(atom_coefficient)) allocate(atom_coefficient(at%N))

      do i = 1, at%N
         if(.not. allocated(atom_coefficient(i)%neighbour)) allocate(atom_coefficient(i)%neighbour(n_neighbours(at,i)))
         do n = 1, n_neighbours(at,i)

            j = neighbour(at,i,n,distance = r, diff = d)
            atom_coefficient(i)%neighbour(n)%r = r
            atom_coefficient(i)%neighbour(n)%u = d / r

            if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics)) allocate( atom_coefficient(i)%neighbour(n)%spherical_harmonics(0:l_max), &
            atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(0:l_max) )
            do l = 0, l_max
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m)) &
               allocate(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(-l:l))
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm)) &
               allocate(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(3,-l:l))

               atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m = CPLX_ZERO
               atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm = CPLX_ZERO

               do m = -l, l
                  atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(m) = SphericalYCartesian(l,m,d)
                  atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(:,m) = GradSphericalYCartesian(l,m,d)
               enddo
            enddo
         enddo
      enddo

   endsubroutine real_space_fourier_coefficients

   function real_space_covariance_coefficient(anc1,anc2,i1,i2,alpha,l_max,f1,f2)
      type(neighbour_type), dimension(:), intent(in) :: anc1, anc2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(out), optional :: f1, f2

      real(dp) :: real_space_covariance_coefficient

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: n1, n2, l, m1, m2, k
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, grad_mo_spher_bess_fi_ki_l
      real(dp), dimension(3) :: u1, u2, grad_arg_bess1, grad_fac_exp1, grad_arg_bess2, grad_fac_exp2
      type(cplx_2d), dimension(:), allocatable :: integral_r
      type(grad_spherical_harmonics_overlap_type), dimension(:), allocatable :: grad_integral_r1, grad_integral_r2

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2))

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      if(present(f1)) then
         allocate(grad_integral_r1(0:size(anc1(i1)%neighbour)))
         do n1 = 0, size(anc1(i1)%neighbour)
            allocate(grad_integral_r1(n1)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r1(n1)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r1(n1)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif

      if(present(f2)) then
         allocate(grad_integral_r2(0:size(anc2(i2)%neighbour)))
         do n2 = 0, size(anc2(i2)%neighbour)
            allocate(grad_integral_r2(n2)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r2(n2)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r2(n2)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif
      do n1 = 1, size(anc1(i1)%neighbour)
         r1 = anc1(i1)%neighbour(n1)%r
         u1 = anc1(i1)%neighbour(n1)%u
         do n2 = 1, size(anc2(i2)%neighbour)
            r2 = anc2(i2)%neighbour(n2)%r

            u2 = anc2(i2)%neighbour(n2)%u

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            if(present(f1)) then
               grad_arg_bess1 = alpha*r2*u1
               grad_fac_exp1 = -fac_exp*alpha*r1*u1
            endif

            if(present(f2)) then
               grad_arg_bess2 = alpha*r1*u2
               grad_fac_exp2 = -fac_exp*alpha*r2*u2
            endif

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                  if(do_derivative) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
               else
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  if(do_derivative) then
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                  else
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                  endif

               endif


               if(do_derivative) grad_mo_spher_bess_fi_ki_l = 0.5_dp * (mo_spher_bess_fi_ki_lp - mo_spher_bess_fi_ki_l / arg_bess + mo_spher_bess_fi_ki_lm)

               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                     if(present(f1)) then
                        grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * &
                        ( conjg(anc1(i1)%neighbour(n1)%grad_spherical_harmonics(l)%mm(:,m1)) * mo_spher_bess_fi_ki_l*fac_exp + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp1 ) )
                     endif

                     if(present(f2)) then
                        grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * &
                        ( anc2(i2)%neighbour(n2)%grad_spherical_harmonics(l)%mm(:,m2) * mo_spher_bess_fi_ki_l*fac_exp + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp2 ) )
                     endif

                  enddo
               enddo
            enddo
         enddo
      enddo

      if(present(f1)) then
         f1 = 0.0_dp
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f1(k,n1+1) = f1(k,n1+1) + real(sum(conjg(grad_integral_r1(n1)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f1 = 2.0_dp * f1
      endif

      if(present(f2)) then
         f2 = 0.0_dp
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f2(k,n2+1) = f2(k,n2+1) + real(sum(conjg(grad_integral_r2(n2)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f2 = 2.0_dp * f2
      endif

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance_coefficient = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

      if(present(f1)) then
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r1(n1)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r1(n1)%grad_integral)
         enddo
         deallocate(grad_integral_r1)
      endif

      if(present(f2)) then
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r2(n2)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r2(n2)%grad_integral)
         enddo
         deallocate(grad_integral_r2)
      endif

   endfunction real_space_covariance_coefficient

   function real_space_covariance(at1,at2,i1,i2,alpha,l_max,f1,f2)
      type(atoms), intent(in) :: at1, at2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(inout), optional :: f1, f2

      real(dp) :: real_space_covariance

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: j1, j2, n1, n2, l, m1, m2
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm
      real(dp), dimension(3) :: d1, d2
      type(cplx_2d), dimension(:), allocatable :: integral_r

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2))

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      do n1 = 1, n_neighbours(at1,i1)
         j1 = neighbour(at1,i1,n1,distance = r1, diff = d1)
         do n2 = 1, n_neighbours(at2,i2)
            j2 = neighbour(at2,i2,n2,distance = r2, diff = d2)

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lmm = sinh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm
               elseif( l == 1 ) then
                  mo_spher_bess_fi_ki_lm = ( arg_bess*cosh(arg_bess) - sinh(arg_bess) ) / arg_bess**2
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lm
               else
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l+1)*mo_spher_bess_fi_ki_lm / arg_bess
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
               endif

               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(SphericalYCartesian(l,m1,d1)) * SphericalYCartesian(l,m2,d2)*mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                  enddo
               enddo
            enddo
         enddo
      enddo

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

   endfunction real_space_covariance

   function RadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: RadialFunction

      real(dp), dimension(this%n_max) :: h
      integer :: j

      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = (this%cutoff-r)**(j+2) / this%NormFunction(j)
         enddo
         RadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         RadialFunction = 0.0_dp
      endif

   endfunction RadialFunction

   function GradRadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: GradRadialFunction

      real(dp), dimension(this%n_max) :: h
      integer :: j

      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = - (j+2) * (this%cutoff-r)**(j+1) / this%NormFunction(j)
         enddo
         GradRadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         GradRadialFunction = 0.0_dp
      endif

   endfunction GradRadialFunction



   function graphIsConnected(connectivityMatrix,error)

      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(out), optional :: error
      logical :: graphIsConnected

      logical, dimension(:), allocatable :: visitedVertices

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphIsConnected: not square matrix",error)
      endif

      allocate(visitedVertices(size(connectivityMatrix,1)))

      call graphBFS(connectivityMatrix,1,visitedVertices=visitedVertices,error=error)
      graphIsConnected = all(visitedVertices)

      deallocate(visitedVertices)

   endfunction graphIsConnected

   subroutine graphBFS(connectivityMatrix,startVertex,visitedVertices,tree,error)

      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(in) :: startVertex
      logical, dimension(:), target, intent(out), optional :: visitedVertices
      integer, dimension(:,:), allocatable, intent(out), optional :: tree
      integer, intent(out), optional :: error

      type(LinkedList_i1d), pointer :: LL_edges => null(), LL_remove => null(), LL_tree => null()

      logical, dimension(:), pointer :: my_visitedVertices
      integer, dimension(:), pointer :: edge
      integer, dimension(2) :: vw

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphBFS: not square matrix",error)
      endif

      if( present( visitedVertices ) ) then
         my_visitedVertices => visitedVertices
      else
         allocate(my_visitedVertices(size(connectivityMatrix,1)))
      endif

      my_visitedVertices = .false.
      call graphSearch(connectivityMatrix,startVertex,LL_edges,my_visitedVertices,error)
      do while( associated(LL_edges) )
         LL_remove => LL_edges
         edge => retrieve_node(LL_remove)
         vw = edge
         call delete_node(LL_edges,LL_remove)
         if( .not. my_visitedVertices(vw(2)) ) then

            if(present(tree)) call append(LL_tree,vw)
            call graphSearch(connectivityMatrix, vw(2), LL_edges, my_visitedVertices,error)
         endif
      enddo

      if( .not. present( visitedVertices ) ) deallocate(my_visitedVertices)

      if (present(tree)) then
         call retrieve(LL_tree,tree)
         call finalise(LL_tree)
      endif

   endsubroutine graphBFS

   subroutine graphSearch(connectivityMatrix, vertex, LL_edges, visitedVertices,error)
      logical, dimension(:,:), intent(in) :: connectivityMatrix
      integer, intent(in) :: vertex
      type(LinkedList_i1d), pointer, intent(inout) :: LL_edges
      logical, dimension(:), intent(inout) :: visitedVertices
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( .not. is_square(connectivityMatrix) ) then
         RAISE_ERROR("graphSearch: not square matrix",error)
      endif

      visitedVertices(vertex) = .true.

      do i = 1, size(connectivityMatrix,1)
         if( connectivityMatrix(i,vertex) ) call append(LL_edges,(/vertex,i/))
      enddo

   endsubroutine graphSearch


endmodule descriptors_module
