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

!
!
! This file contains descriptors written by others, not A. P. Bartok and Gabor Csanyi
! and the code here is owned by their respective authors, as indicated below. 
!

!!!! Author: Wojciech Slachta
   type bond_real_space
      real(dp) :: bond_cutoff
      real(dp) :: bond_transition_width
      real(dp) :: cutoff
      real(dp) :: transition_width
      real(dp) :: atom_sigma
      integer :: max_neighbours

      logical :: initialised = .false.

   endtype bond_real_space

!!!! Authors: Alan Nichol and S. T. John

   type AN_monomer
      real(dp) :: cutoff
      integer :: atomic_number
      integer :: N

      logical :: initialised = .false.
      logical :: do_atomic = .false.

   endtype AN_monomer

   type general_monomer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature
      real(dp) :: cutoff, cutoff_transition_width
      logical  :: atom_ordercheck, internal_swaps_only
      logical  :: strict
      real(dp) :: power
      logical  :: initialised = .false.
   endtype general_monomer

   type general_dimer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature_one, signature_two
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff
      logical  :: atom_ordercheck, internal_swaps_only, use_smooth_cutoff, monomers_identical,double_count
      logical  :: strict, use_com, mpifind, strict_mask
      type(transfer_parameters_type) :: transfer_parameters
      logical  :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
      real(dp) :: power,dist_shift
   endtype general_dimer

   type general_trimer
      type(permutation_data_type) :: permutation_data
      integer, dimension(:), allocatable :: signature_one, signature_two, signature_three
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff, monomer_three_cutoff
      logical :: atom_ordercheck, internal_swaps_only, use_smooth_cutoff, one_two_identical, one_three_identical, two_three_identical
      logical :: strict, use_com, mpifind
      logical :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
      real(dp) :: power,dist_shift
   endtype general_trimer


   type molecule_lo_d
      type(permutation_data_type) :: permutation_data
      type(Atoms) :: template_atoms
      integer :: n_atoms, max_dimension ! max_dimension is descriptor dimension if include all interatomic distances
      integer, dimension(:), allocatable :: signature, included_components
      integer, dimension(:,:), allocatable :: component_atoms
      real(dp) :: cutoff, cutoff_transition_width
      integer :: neighbour_graph_depth
      logical :: atom_ordercheck, use_smooth_cutoff
      logical :: initialised = .false.
      type(Table) :: bonds, atom_pairs
      integer :: distance_transform
   endtype molecule_lo_d


   type com_dimer
      integer, dimension(:), allocatable :: signature_one, signature_two
      real(dp) :: cutoff, cutoff_transition_width, monomer_one_cutoff, monomer_two_cutoff
      logical :: atom_ordercheck, use_smooth_cutoff, monomers_identical
      logical :: strict, mpifind
      type(transfer_parameters_type) :: transfer_parameters
      logical :: initialised = .false.
      logical, dimension(:), allocatable :: is_intermolecular, cutoff_contributor
   endtype com_dimer

!!!!! Authors: A. P. Bartok and Miguel Caro
   type soap_express
      ! User controllable parameters
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      real(dp) :: cutoff_decay_rate

      real(dp) :: atom_sigma_radial, atom_sigma_angular, central_weight, covariance_sigma0, &
         atom_sigma_scaling_radial, atom_sigma_scaling_angular, amplitude_scaling

      integer :: l_max, n_max

      ! Internal pre-initialised variables
      real(dp), dimension(:,:), allocatable :: basis_transformation_coefficients, &
         overlap
      real(dp), dimension(:), allocatable :: semifactorial_table, Y_lm_prefactor

      integer :: angular_array_size

      logical :: initialised = .false.
   endtype soap_express


   type soap_turbo
      ! User controllable parameters
      real(dp) :: rcut_hard, rcut_soft, nf
      integer  :: n_species, radial_enhancement, central_index, l_max
      character(len=STRING_LENGTH) :: basis, scaling_mode, compress_file

      real(dp), dimension(:), allocatable :: atom_sigma_r, atom_sigma_r_scaling, &
         atom_sigma_t, atom_sigma_t_scaling, amplitude_scaling, central_weight
      integer, dimension(:), allocatable :: species_Z, alpha_max, compress_indices

      logical :: initialised = .false., compress = .false.
   endtype soap_turbo
