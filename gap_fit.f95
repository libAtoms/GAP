! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   GAP (Gaussian Approximation Potental)
! HND X   
! HND X
! HND X   Portions of GAP were written by Albert Bartok-Partay, Gabor Csanyi, 
! HND X   and Sascha Klawohn. Copyright 2006-2021.
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

program gap_fit_program

  use libatoms_module
  use gp_predict_module
  use gap_fit_module
  use task_manager_module

  implicit none

  type(gap_fit) :: main_gap_fit

  call system_initialise(verbosity=PRINT_NORMAL, enable_timing=.false.)
  call gap_fit_init_mpi_scalapack(main_gap_fit)

  call gap_fit_parse_command_line(main_gap_fit)
  call gap_fit_parse_gap_str(main_gap_fit)

  call gap_fit_read_core_param_file(main_gap_fit)

  call add_template_string(main_gap_fit) ! if descriptor requires a template xyz file and this is provided, write to a string and add to descriptor_str

  call read_descriptors(main_gap_fit) ! initialises descriptors from the descriptor_str and sets max_cutoff according to that.
  call read_fit_xyz(main_gap_fit)   ! reads in xyz into an array of atoms objects. sets cutoff and does calc_connect on each frame
  call print('XYZ file read')

  call gap_fit_init_task_manager(main_gap_fit)

  call get_species_xyz(main_gap_fit) ! counts the number of species present in the xyz file.
  call add_multispecies_gaps(main_gap_fit)

  call get_n_sparseX_for_files(main_gap_fit)
  call parse_config_type_sigma(main_gap_fit)
  call parse_config_type_n_sparseX(main_gap_fit)

  if(any(main_gap_fit%add_species)) then ! descriptor_str might have changed. reinitialises descriptors from the descriptor_str and sets max_cutoff according to that.
     call read_descriptors(main_gap_fit)
  endif
  call print('Multispecies support added where requested')

  call fit_n_from_xyz(main_gap_fit) ! counts number of energies, forces, virials. computes number of descriptors and gradients.
  call gap_fit_distribute_tasks(main_gap_fit)
  if (main_gap_fit%task_manager%n_workers > 1) call fit_n_from_xyz(main_gap_fit)
  call gap_fit_set_mpi_blocksizes(main_gap_fit)
  call gap_fit_estimate_memory(main_gap_fit)

  if (main_gap_fit%dryrun) then
     call print('Exit before major allocations because dryrun is true.')
     call system_finalise()
     stop
  end if

  call set_baselines(main_gap_fit) ! sets e0 etc.

  call fit_data_from_xyz(main_gap_fit) ! converts atomic neighbourhoods (bond neighbourhoods etc.) do descriptors, and feeds those to the GP
  call print('Cartesian coordinates transformed to descriptors')

  if(main_gap_fit%sparsify_only_no_fit) then
     if (gap_fit_is_root(main_gap_fit)) then
        call initialise(main_gap_fit%gp_sp, main_gap_fit%my_gp)
        call gap_fit_print_xml(main_gap_fit, main_gap_fit%gp_file, main_gap_fit%sparseX_separate_file)
     end if
     call system_finalise()
     stop
  end if

  call enable_timing()
  call system_timer('GP sparsify')

  call gp_covariance_sparse(main_gap_fit%my_gp)
  call gap_fit_print_linear_system_dump_file(main_gap_fit)
  call gpSparse_fit(main_gap_fit%gp_sp, main_gap_fit%my_gp, main_gap_fit%task_manager, main_gap_fit%condition_number_norm)

  if (gap_fit_is_root(main_gap_fit)) call gap_fit_print_xml(main_gap_fit, main_gap_fit%gp_file, main_gap_fit%sparseX_separate_file)

  call system_timer('GP sparsify')
  call system_finalise()

end program gap_fit_program
