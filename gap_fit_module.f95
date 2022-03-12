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

#include "error.inc"

module gap_fit_module

  use error_module
  use libatoms_module
  use descriptors_module
  use gp_predict_module
  use gp_fit_module
  use fox_wxml
  use potential_module
  use ScaLAPACK_module
  use task_manager_module

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: THETA_LENGTH = 10000

  integer, parameter :: E0_ISOLATED = 1
  integer, parameter :: E0_AVERAGE = 2
  integer, parameter :: EXCLUDE_LOC = -1

#ifdef GAP_VERSION
  integer, parameter, private :: gap_version = GAP_VERSION
#else
  integer, parameter, private :: gap_version = huge(1)
#endif

  type gap_fit
  !% everything from the command line
     type(Atoms), dimension(:), allocatable :: at
     
     character(len=STRING_LENGTH) :: at_file='', core_ip_args = '', e0_str, local_property0_str, &
     energy_parameter_name, local_property_parameter_name, force_parameter_name, virial_parameter_name, &
     stress_parameter_name, hessian_parameter_name, config_type_parameter_name, sigma_parameter_name, &
     config_type_sigma_string, core_param_file, gp_file, template_file, force_mask_parameter_name, &
     condition_number_norm, linear_system_dump_file

     character(len=10240) :: command_line = ''
     real(dp), dimension(total_elements) :: e0, local_property0
     real(dp) :: max_cutoff
     real(dp), dimension(4) :: default_sigma
     real(dp) :: default_local_property_sigma
     real(dp) :: sparse_jitter, e0_offset, hessian_delta
     integer :: e0_method = E0_ISOLATED
     logical :: do_core = .false., do_copy_at_file, has_config_type_sigma, sigma_per_atom = .true.
     logical :: sparsify_only_no_fit = .false.
     integer :: n_frame = 0
     integer :: n_coordinate = 0
     integer :: n_ener = 0
     integer :: n_force = 0
     integer :: n_virial = 0
     integer :: n_hessian = 0
     integer :: n_local_property = 0
     integer :: n_species = 0
     integer :: min_save
     integer :: mpi_blocksize
     type(extendable_str) :: quip_string
     type(Potential) :: core_pot

     type(gpFull) :: my_gp
     type(gpSparse) :: gp_sp

     type(MPI_Context) :: mpi_obj
     type(ScaLAPACK) :: ScaLAPACK_obj
     type(task_manager_type) :: task_manager

     type(descriptor), dimension(:), allocatable :: my_descriptor
     character(len=STRING_LENGTH), dimension(200) :: gap_str

     real(dp), dimension(:), allocatable :: delta, f0, theta_uniform, zeta, unique_hash_tolerance, unique_descriptor_tolerance !, theta
     real(dp), dimension(:,:), allocatable :: sigma
     integer, dimension(:), allocatable :: n_sparseX, sparse_method, target_type, n_cross, n_descriptors, species_Z, covariance_type
     integer, dimension(:,:), allocatable :: config_type_n_sparseX
     character(len=STRING_LENGTH), dimension(:), allocatable :: theta_file, sparse_file, theta_fac_string, config_type, config_type_n_sparseX_string, print_sparse_index
     logical, dimension(:), allocatable :: mark_sparse_atoms, add_species, has_theta_fac, has_theta_uniform, has_theta_file, has_zeta

     logical :: sparseX_separate_file
     logical :: sparse_use_actual_gpcov
     logical :: has_template_file, has_e0, has_local_property0, has_e0_offset, has_linear_system_dump_file

  endtype gap_fit
     
  private

  public :: fit_n_from_xyz
  public :: fit_data_from_xyz
  public :: e0_from_xyz
  public :: w_Z_from_xyz
  public :: gap_fit
  public :: gap_fit_print_xml
  public :: file_print_xml
!  public :: print_sparse
  public :: set_baselines
  public :: parse_config_type_sigma
  public :: parse_config_type_n_sparseX
  public :: read_fit_xyz
  public :: read_descriptors
  public :: get_species_xyz
  public :: add_multispecies_gaps
  public :: add_template_string
  public :: gap_fit_parse_command_line
  public :: gap_fit_parse_gap_str
  public :: gap_fit_read_core_param_file

  public :: gap_fit_init_mpi_scalapack
  public :: gap_fit_init_task_manager
  public :: gap_fit_distribute_tasks

  public :: gap_fit_is_root

  public :: gap_fit_print_linear_system_dump_file
  public :: gap_fit_estimate_memory

contains

  subroutine gap_fit_parse_command_line(this)
  !% This subroutine parses the main command line options.
     type(gap_fit), intent(inout), target :: this
     type(Dictionary) :: params

     character(len=STRING_LENGTH), pointer :: at_file, e0_str, local_property0_str, &
          core_param_file, core_ip_args, &
          energy_parameter_name, local_property_parameter_name, force_parameter_name, &
          virial_parameter_name, stress_parameter_name, hessian_parameter_name, &
          config_type_parameter_name, sigma_parameter_name, config_type_sigma_string, &
          gp_file, template_file, force_mask_parameter_name, condition_number_norm, &
          linear_system_dump_file

     character(len=STRING_LENGTH) ::  gap_str, verbosity, sparse_method_str, covariance_type_str, e0_method, &
        parameter_name_prefix

     logical, pointer :: sigma_per_atom, do_copy_at_file, sparseX_separate_file, sparse_use_actual_gpcov
     logical :: do_ip_timing, has_sparse_file, has_theta_uniform, has_at_file, has_gap, has_default_sigma
     logical, pointer :: sparsify_only_no_fit
     
     real(dp), pointer :: e0_offset, sparse_jitter, hessian_delta
     real(dp), dimension(:), pointer :: default_sigma
     real(dp), pointer :: default_local_property_sigma

     integer :: rnd_seed
     integer, pointer :: mpi_blocksize

     at_file => this%at_file
     e0_str => this%e0_str
     local_property0_str => this%local_property0_str
     e0_offset => this%e0_offset
     default_sigma => this%default_sigma
     default_local_property_sigma => this%default_local_property_sigma
     sparse_jitter => this%sparse_jitter
     hessian_delta => this%hessian_delta
     core_param_file => this%core_param_file
     core_ip_args => this%core_ip_args
     energy_parameter_name => this%energy_parameter_name
     local_property_parameter_name => this%local_property_parameter_name
     force_parameter_name => this%force_parameter_name
     virial_parameter_name => this%virial_parameter_name
     stress_parameter_name => this%stress_parameter_name
     hessian_parameter_name => this%hessian_parameter_name
     config_type_parameter_name => this%config_type_parameter_name
     sigma_parameter_name => this%sigma_parameter_name
     force_mask_parameter_name => this%force_mask_parameter_name
     config_type_sigma_string => this%config_type_sigma_string
     sigma_per_atom => this%sigma_per_atom
     do_copy_at_file => this%do_copy_at_file
     sparseX_separate_file => this%sparseX_separate_file
     sparse_use_actual_gpcov => this%sparse_use_actual_gpcov
     gp_file => this%gp_file
     template_file => this%template_file
     sparsify_only_no_fit => this%sparsify_only_no_fit
     condition_number_norm => this%condition_number_norm
     linear_system_dump_file => this%linear_system_dump_file
     mpi_blocksize => this%mpi_blocksize
     
     call initialise(params)
     
     call param_register(params, 'atoms_filename', '//MANDATORY//', at_file, has_value_target = has_at_file, help_string="XYZ file with fitting configurations", altkey="at_file")
     call param_register(params, 'gap', '//MANDATORY//', gap_str, has_value_target = has_gap, help_string="Initialisation string for GAPs")
     call param_register(params, 'e0', '0.0', e0_str, has_value_target = this%has_e0, &
          help_string="Atomic energy value to be subtracted from energies before fitting (and added back on after prediction). &
          & Specifiy a single number (used for all species) or by species: {Ti:-150.0:O:-320...}. energy = baseline + GAP + e0")
     
     call param_register(params, 'local_property0', '0.0', local_property0_str, has_value_target = this%has_local_property0, &
          help_string="Local property value to be subtracted from the local property before fitting (and added back on after prediction). &
          & Specifiy a single number (used for all species) or by species: {H:20.0:Cl:35.0...}.")
     
     call param_register(params, 'e0_offset', '0.0', e0_offset, has_value_target = this%has_e0_offset, &
          help_string="Offset of baseline. If zero, the offset is the average atomic energy of the input data or the e0 specified manually.")
   
     call param_register(params, 'e0_method','isolated',e0_method, &
        help_string="Method to determine e0, if not explicitly specified. Possible options: isolated (default, each atom &
        present in the XYZ needs to have an isolated representative, with a valid energy), average (e0 is the average of &
        all total energies across the XYZ)")

     call param_register(params, 'default_kernel_regularisation', '//MANDATORY//', default_sigma, has_value_target = has_default_sigma, &
         help_string="error in [energies forces virials hessians]", altkey="default_sigma")
   
     call param_register(params, 'default_kernel_regularisation_local_property', '0.001', default_local_property_sigma, &
         help_string="error in local_property", altkey="default_local_property_sigma")

     call param_register(params, 'sparse_jitter', "1.0e-10", sparse_jitter, &
         help_string="Extra regulariser used to regularise the sparse covariance matrix before it is passed to the linear solver. Use something small, it really shouldn't affect your results, if it does, your sparse basis is still very ill-conditioned.")
     
     call param_register(params, 'hessian_displacement', "1.0e-2", hessian_delta, &
         help_string="Finite displacement to use in numerical differentiation when obtaining second derivative for the Hessian covariance", altkey="hessian_delta")
     
     call param_register(params, 'baseline_param_filename', 'quip_params.xml', core_param_file, &
         help_string="QUIP XML file which contains a potential to subtract from data (and added back after prediction)", altkey="core_param_file")
     
     call param_register(params, 'baseline_ip_args', '', core_ip_args, has_value_target = this%do_core, &
          help_string=" QUIP init string for a potential to subtract from data (and added back after prediction)", altkey="core_ip_args")
     
     call param_register(params, 'energy_parameter_name', 'energy', energy_parameter_name, &
          help_string="Name of energy property in the input XYZ file that describes the data")
     
     call param_register(params, 'local_property_parameter_name', 'local_property', local_property_parameter_name, &
          help_string="Name of local_property (column) in the input XYZ file that describes the data")
     
     call param_register(params, 'force_parameter_name', 'force', force_parameter_name, &
          help_string="Name of force property (columns) in the input XYZ file that describes the data")
     
     call param_register(params, 'virial_parameter_name', 'virial', virial_parameter_name, &
          help_string="Name of virial property in the input XYZ file that describes the data")

     call param_register(params, 'stress_parameter_name', 'stress', stress_parameter_name, &
          help_string="Name of stress property (6-vector or 9-vector) in the input XYZ file that describes the data - stress values only used if virials are not available (opposite sign, standard Voigt order)")
     
     call param_register(params, 'hessian_parameter_name', 'hessian', hessian_parameter_name, &
          help_string="Name of hessian property (column) in the input XYZ file that describes the data")
     
     call param_register(params, 'config_type_parameter_name', 'config_type', config_type_parameter_name, &
          help_string="Allows grouping on configurations into. This option is the name of the key that indicates the configuration type in the input XYZ file. With the default, the key-value pair config_type=blah would place that configuration into the group blah.")
     
     call param_register(params, 'kernel_regularisation_parameter_name', 'sigma', sigma_parameter_name, &
          help_string="kernel regularisation parameters for a given configuration in the database. &
          Overrides the command line values (both defaults and config-type-specific values). In the input XYZ file, it must be prepended by energy_, force_, virial_ or hessian_", altkey="sigma_parameter_name")
     
     call param_register(params, 'force_mask_parameter_name', 'force_mask', force_mask_parameter_name, &
          help_string="To exclude forces on specific atoms from the fit. In the XYZ, it must be a logical column.")
     
     call param_register(params, 'parameter_name_prefix', '', parameter_name_prefix, &
          help_string="Prefix that gets uniformly appended in front of {energy,local_property,force,virial,...}_parameter_name")
     
     call param_register(params, 'config_type_kernel_regularisation', '', config_type_sigma_string, has_value_target = this%has_config_type_sigma, &
          help_string="What kernel regularisation values to choose for each type of data, when the configurations are grouped into config_types. Format: {configtype1:energy:force:virial:hessian:config_type2:energy:force:virial:hessian...}", altkey="config_type_sigma")

     call param_register(params, 'kernel_regularisation_is_per_atom', 'T', sigma_per_atom, &
          help_string="Interpretation of the energy and virial sigmas specified in >>default_kernel_regularisation<< and >>config_type_kernel_regularisation<<. &
          If >>T<<, they are interpreted as per-atom errors, and the variance will be scaled according to the number of atoms in the configuration. &
          If >>F<< they are treated as absolute errors and no scaling is performed. &
          NOTE: values specified on a per-configuration basis (see >>kernel_regularisation_parameter_name<<) are always absolute, not per-atom.", altkey="sigma_per_atom")
   
     call param_register(params, 'do_copy_atoms_file', 'T', do_copy_at_file, &
          help_string="Copy the input XYZ file into the GAP XML file (should be set to False for NetCDF input).", altkey="do_copy_at_file")
   
     call param_register(params, 'sparse_separate_file', 'T', sparseX_separate_file, &
          help_string="Save sparse point data in separate file in binary (use it for large datasets)")
   
     call param_register(params, 'sparse_use_actual_gpcov', 'F', sparse_use_actual_gpcov, &
          help_string="Use actual GP covariance for sparsification methods")
   
     call param_register(params, 'gap_file', 'gap_new.xml', gp_file, &
          help_string="Name of output XML file that will contain the fitted potential", altkey="gp_file")
   
     call param_register(params, 'verbosity', 'NORMAL', verbosity, &
          help_string="Verbosity control. Options: NORMAL, VERBOSE, NERD, ANALYSIS.") ! changed name to ANALYSIS now that we are grown up
   
     call param_register(params, "rnd_seed", "-1", rnd_seed, &
          help_string="Random seed.")
   
     call param_register(params, "openmp_chunk_size", "1", openmp_chunk_size, &
          help_string="Chunk size in OpenMP scheduling")
   
     call param_register(params, 'do_ip_timing', 'F', do_ip_timing, &
          help_string="To enable or not timing of the interatomic potential.")
   
     call param_register(params, 'template_file', 'template.xyz', template_file, has_value_target=this%has_template_file, &
          help_string="Template XYZ file for initialising object")

     call param_register(params, 'sparsify_only_no_fit', 'F', sparsify_only_no_fit, &
          help_string="If true, sparsification is done, but no fitting. print the sparse index by adding print_sparse_index=file.dat to the descriptor string.")
     
     call param_register(params, 'condition_number_norm', ' ', condition_number_norm, &
          help_string="Norm for condition number of matrix A; O: 1-norm, I: inf-norm, <space>: skip calculation (default)")

     call param_register(params, 'linear_system_dump_file', '', linear_system_dump_file, has_value_target=this%has_linear_system_dump_file, &
          help_string="Basename prefix of linear system dump files. Skipped if empty (default).")

     call param_register(params, 'mpi_blocksize', '0', mpi_blocksize, &
          help_string="Blocksize of MPI distributed matrices. Affects efficiency and memory usage. Max if 0 (default).")

     if (.not. param_read_args(params, command_line=this%command_line)) then
        call print("gap_fit")
        call system_abort('Exit: Mandatory argument(s) missing...')
     endif
     call print_title("Input parameters")
     call param_print(params)
     call print_title("")
     call finalise(params)
     

     if (len_trim(parameter_name_prefix) > 0) then
        energy_parameter_name = trim(parameter_name_prefix) // trim(energy_parameter_name)
        local_property_parameter_name = trim(parameter_name_prefix) // trim(local_property_parameter_name)
        force_parameter_name = trim(parameter_name_prefix) // trim(force_parameter_name)
        virial_parameter_name = trim(parameter_name_prefix) // trim(virial_parameter_name)
        hessian_parameter_name = trim(parameter_name_prefix) // trim(hessian_parameter_name)
        stress_parameter_name = trim(parameter_name_prefix) // trim(stress_parameter_name)
        config_type_parameter_name = trim(parameter_name_prefix) // trim(config_type_parameter_name)
        sigma_parameter_name = trim(parameter_name_prefix) // trim(sigma_parameter_name)
        force_mask_parameter_name = trim(parameter_name_prefix) // trim(force_mask_parameter_name)
     endif

     if (sparsify_only_no_fit) then
        force_parameter_name = '//IGNORE//'
        virial_parameter_name = '//IGNORE//'
        hessian_parameter_name = '//IGNORE//'
        stress_parameter_name = '//IGNORE//'
        call print_warning("sparsify_only_no_fit == T: force, virial, hessian, stress parameters are ignored.")
     end if
   
     if( len_trim(this%gp_file) > 216 ) then    ! The filename's length is limited to 255 char.s in some filesystem. 
                                        ! Without this check, the fit would run but produce a core file and only a temporary xml file. 
                                        ! The limit is set to 216 as the sparse file can be 39 characters longer.
       call system_abort("gap_file's name "//this%gp_file//" is too long. Please start the fit again with a shorter name.")
     endif

     if(do_ip_timing) call enable_timing()

     select case(verbosity)
       case ("NORMAL")
         call verbosity_push(PRINT_NORMAL)
       case ("VERBOSE")
         call verbosity_push(PRINT_VERBOSE)
       case ("NERD")
         call verbosity_push(PRINT_NERD)
       case ("ANALYSIS")                     ! changed name now that we are grown up
         call verbosity_push(PRINT_ANALYSIS) ! changed name now that we are grown up
       case default
         call system_abort("confused by verbosity " // trim(verbosity))
     end select

     select case(lower_case(e0_method))
       case ("isolated")
         this%e0_method = E0_ISOLATED
       case ("average")
         this%e0_method = E0_AVERAGE
       case default
         call system_abort("confused by e0_method " // trim(e0_method))
     end select

     if (rnd_seed >= 0) call system_set_random_seeds(rnd_seed)

     call print_title('Gaussian Approximation Potentials - Database fitting')
     call print('')
     call print('Initial parsing of command line arguments finished.')

     call split_string(gap_str,':;','{}',this%gap_str(:),this%n_coordinate,matching=.true.)

     call print('Found '//this%n_coordinate//' GAPs.')

  endsubroutine gap_fit_parse_command_line

  subroutine set_baselines(this)
     type(gap_fit), intent(inout) :: this

     integer :: i

     this%e0 = 0.0_dp

     if( count( (/this%has_e0, this%has_e0_offset/) ) > 1 ) then
        call print_warning('Both e0 and e0_offset has been specified. That means your atomic energy is e0 + e0_offset')
     endif

     if( this%has_e0 ) then
        call parse_atomtype_value_str(this%e0_str,this%e0)
     else
        call e0_from_xyz(this) ! calculates the average atomic energy so it can be subtracted later.
     endif

     if( this%has_e0_offset ) this%e0 = this%e0 + this%e0_offset

     if( .not. this%has_e0 ) then
        do i = 1, size(this%e0)
           if( all(i/=this%species_Z) ) this%e0(i) = 0.0_dp
        enddo
        call print('E0/atom = '//this%e0)
     endif

     if( this%has_local_property0 ) then
        call parse_atomtype_value_str(this%local_property0_str,this%local_property0)
        this%e0 = 0.0_dp
     else
        this%local_property0 = 0.0_dp
     endif

  endsubroutine set_baselines

  subroutine parse_atomtype_value_str(this,values,error)

     character(len=STRING_LENGTH), intent(in) :: this
     real(dp), dimension(total_elements), intent(out) :: values
     integer, intent(out), optional :: error

     integer :: n_string_array, i, z
     character(len=STRING_LENGTH), dimension(2*total_elements) :: string_array

     INIT_ERROR(error)

     call split_string(this,':','{}',string_array(:),n_string_array,matching=.true.)
     if(n_string_array == 1) then
        values = string_to_real(trim(string_array(1)))
     elseif(mod(n_string_array,2) == 0) then
        values = 0.0_dp
        do i = 1, n_string_array / 2
           z = atomic_number(trim( string_array((i-1)*2+1) ))
           if( z==0 )  then
              RAISE_ERROR("parse_atomtype_value_str: invalid atomic symbol "//trim(string_array((i-1)*2+1)),error)
           endif
           values(z) = string_to_real(trim( string_array(2*i) ))
        enddo
     else
        RAISE_ERROR("parse_atomtype_value_str: number of fields is an odd number. It must be a list of pairs of values, such as {Ti:-150.4:O:-345.1}",error)
     endif

  endsubroutine parse_atomtype_value_str

  subroutine gap_fit_parse_gap_str(this)
  !% This subroutine parses the options given in the gap string, for each GAP.
     type(gap_fit), intent(inout), target :: this
     type(Dictionary) :: params

     integer :: i_coordinate

     real(dp) :: delta, f0, theta_uniform, zeta, unique_hash_tolerance, unique_descriptor_tolerance
     integer :: n_sparseX, sparse_method, covariance_type
     character(len=STRING_LENGTH) :: config_type_n_sparseX_string, theta_fac_string, theta_file, sparse_file, print_sparse_index, &
                                     covariance_type_str, sparse_method_str
     logical :: mark_sparse_atoms, add_species, has_sparse_file

     allocate(this%delta(this%n_coordinate))
     allocate(this%f0(this%n_coordinate))
     allocate(this%n_sparseX(this%n_coordinate))
     allocate(this%config_type_n_sparseX_string(this%n_coordinate))
     allocate(this%theta_fac_string(this%n_coordinate))
     allocate(this%theta_uniform(this%n_coordinate))
     allocate(this%theta_file(this%n_coordinate))
     allocate(this%has_theta_fac(this%n_coordinate))
     allocate(this%has_theta_uniform(this%n_coordinate))
     allocate(this%has_theta_file(this%n_coordinate))
     allocate(this%sparse_file(this%n_coordinate))
     allocate(this%mark_sparse_atoms(this%n_coordinate))
     allocate(this%sparse_method(this%n_coordinate))
     allocate(this%add_species(this%n_coordinate))
     allocate(this%covariance_type(this%n_coordinate))
     allocate(this%zeta(this%n_coordinate))
     allocate(this%has_zeta(this%n_coordinate))
     allocate(this%print_sparse_index(this%n_coordinate))
     allocate(this%unique_hash_tolerance(this%n_coordinate))
     allocate(this%unique_descriptor_tolerance(this%n_coordinate))

     do i_coordinate = 1, this%n_coordinate
        call initialise(params)

        call param_register(params, 'energy_scale', "//MANDATORY//", delta, &
             help_string="Set the typical scale of the function you are fitting (or the specific term if you use multiple descriptors). It is equivalent to the standard deviation of the Gaussian process in the probabilistic view, and typically this would be &
             set to the standard deviation (i.e. root mean square) of the function &
             that is approximated with the Gaussian process. ", altkey="delta")

        call param_register(params, 'f0', '0.0', f0, &
             help_string="Set the mean of the Gaussian process. Defaults to 0.")

        call param_register(params, 'n_sparse', "0", n_sparseX, &
             help_string="Number of sparse points to use in the sparsification of the Gaussian process")

        call param_register(params, 'config_type_n_sparse', '', config_type_n_sparseX_string, &
             help_string="Number of sparse points in each config type. Format: {type1:50:type2:100}")

        call param_register(params, 'sparse_method', 'RANDOM', sparse_method_str, &
             help_string="Sparsification method. RANDOM(default), PIVOT, CLUSTER, UNIFORM, KMEANS, COVARIANCE, NONE, FUZZY, FILE, &
             INDEX_FILE, CUR_COVARIANCE, CUR_POINTS")

        call param_register(params, 'lengthscale_factor', '1.0', theta_fac_string, has_value_target = this%has_theta_fac(i_coordinate), &
             help_string="Set the width of Gaussians for the Gaussian and PP kernel by multiplying the range of each descriptor by lengthscale_factor. &
             Can be a single number or different for each dimension. For multiple theta_fac separate each value by whitespaces.", altkey="theta_fac")

        call param_register(params, 'lengthscale_uniform', '0.0', theta_uniform, has_value_target = this%has_theta_uniform(i_coordinate), &
             help_string="Set the width of Gaussians for the Gaussian and PP kernel, same in each dimension.", altkey="theta_uniform")

        call param_register(params, 'lengthscale_file', '', theta_file, has_value_target = this%has_theta_file(i_coordinate), &
             help_string="Set the width of Gaussians for the Gaussian kernel from a file. &
             There should be as many real numbers as the number of dimensions, in a single line", altkey="theta_file")

        call param_register(params, 'sparse_file', '', sparse_file, has_value_target = has_sparse_file, &
             help_string="Sparse points from a file. If sparse_method=FILE, descriptor values (real) listed in a text file, one &
             & >>element<< per line. If sparse_method=INDEX_FILE, 1-based index of sparse points, one per line.")

        call param_register(params, 'mark_sparse_atoms', 'F', mark_sparse_atoms, &
             help_string="Reprints the original xyz file after sparsification process. &
             sparse propery added, true for atoms associated with a sparse point.")

        call param_register(params, 'add_species', 'T', add_species, &
             help_string="Create species-specific descriptor, using the descriptor string as a template.")

        call param_register(params, 'covariance_type', "//MANDATORY//", covariance_type_str, &
             help_string="Type of covariance function to use. Available: Gaussian, DOT_PRODUCT, BOND_REAL_SPACE, PP (piecewise polynomial)")

        !call param_register(params, 'theta', '1.0', main_gap_fit%theta(i_coordinate), &
        !help_string="Width of Gaussians for use with bond real space covariance.")

        call param_register(params, 'soap_exponent', '1.0', zeta, has_value_target = this%has_zeta(i_coordinate), &
             help_string="Exponent of soap type dot product covariance kernel", altkey="zeta")

        call param_register(params, 'print_sparse_index', '', print_sparse_index, &
             help_string="If given, after determinining the sparse points, their 1-based indices are appended to this file")

        call param_register(params, 'unique_hash_tolerance', '1.0e-10', unique_hash_tolerance, &
             help_string="Hash tolerance when filtering out duplicate data points")

        call param_register(params, 'unique_descriptor_tolerance', '1.0e-10', unique_descriptor_tolerance, &
             help_string="Descriptor tolerance when filtering out duplicate data points")

        if (.not. param_read_line(params, this%gap_str(i_coordinate), ignore_unknown=.true., task='main program gap_str('//i_coordinate//')')) then
           call system_abort("main program failed to parse gap string ("//i_coordinate//")='"//trim(this%gap_str(i_coordinate))//"'")
        endif
        call finalise(params)

        this%delta(i_coordinate) = delta
        this%f0(i_coordinate) = f0
        this%n_sparseX(i_coordinate) = n_sparseX
        this%config_type_n_sparseX_string(i_coordinate) = config_type_n_sparseX_string
        this%theta_fac_string(i_coordinate) = theta_fac_string
        this%theta_uniform(i_coordinate) = theta_uniform
        this%theta_file(i_coordinate) = theta_file
        this%sparse_file(i_coordinate) = sparse_file
        this%mark_sparse_atoms(i_coordinate) = mark_sparse_atoms
        this%add_species(i_coordinate) = add_species
        this%zeta(i_coordinate) = zeta
        this%print_sparse_index(i_coordinate) = print_sparse_index
        this%unique_hash_tolerance(i_coordinate) = unique_hash_tolerance
        this%unique_descriptor_tolerance(i_coordinate) = unique_descriptor_tolerance

        select case(lower_case(trim(sparse_method_str)))
        case('random')
           this%sparse_method(i_coordinate) = GP_SPARSE_RANDOM
        case('pivot')
           this%sparse_method(i_coordinate) = GP_SPARSE_PIVOT
        case('cluster')
           this%sparse_method(i_coordinate) = GP_SPARSE_CLUSTER
        case('uniform')
           this%sparse_method(i_coordinate) = GP_SPARSE_UNIFORM
        case('kmeans')
           this%sparse_method(i_coordinate) = GP_SPARSE_KMEANS
        case('covariance')
           this%sparse_method(i_coordinate) = GP_SPARSE_COVARIANCE
        case('uniq')
           call system_abort("sparse method UNIQ is no longer in use. Use NONE instead." )
        case('fuzzy')
           this%sparse_method(i_coordinate) = GP_SPARSE_FUZZY
        case('file')
           this%sparse_method(i_coordinate) = GP_SPARSE_FILE
        case('index_file')
           this%sparse_method(i_coordinate) = GP_SPARSE_INDEX_FILE
        case('cur_covariance')
           this%sparse_method(i_coordinate) = GP_SPARSE_CUR_COVARIANCE
        case('cur_points')
           this%sparse_method(i_coordinate) = GP_SPARSE_CUR_POINTS
        case('none')
           this%sparse_method(i_coordinate) = GP_SPARSE_NONE
        case default
           call system_abort("unknown sparse method "//trim(sparse_method_str))
        endselect

        if( has_sparse_file ) then
           if( this%sparse_method(i_coordinate) /= GP_SPARSE_FILE .and. &
               this%sparse_method(i_coordinate) /= GP_SPARSE_INDEX_FILE ) then
              call system_abort('"sparse_file" specified in command line, but sparse method not "file" or "index_file"')
           endif
        endif

        select case(lower_case(trim(covariance_type_str)))
        case('none')
           call system_abort("covariance type cannot be"//trim(covariance_type_str))
           this%covariance_type(i_coordinate) = COVARIANCE_NONE
        case('gaussian')
           this%covariance_type(i_coordinate) = COVARIANCE_ARD_SE
        case('ard_se') ! backwards compatibility
           this%covariance_type(i_coordinate) = COVARIANCE_ARD_SE
        case('dot_product')
           this%covariance_type(i_coordinate) = COVARIANCE_DOT_PRODUCT
        case('bond_real_space')     
           this%covariance_type(i_coordinate) = COVARIANCE_BOND_REAL_SPACE
        case('pp')
           this%covariance_type(i_coordinate) = COVARIANCE_PP
        case default
           call system_abort("unknown covariance type"//trim(covariance_type_str)//". Available: Gaussian, DOT_PRODUCT, BOND_REAL_SPACE, PP (piecewise polynomial)")
        endselect

     enddo

     call print('Descriptors have been parsed')

  endsubroutine gap_fit_parse_gap_str
  
  subroutine read_fit_xyz(this)

    type(gap_fit), intent(inout) :: this

    type(cinoutput) :: xyzfile
    integer :: n_con
    logical :: file_exists

    if( allocated(this%at) ) then
       do n_con = 1, this%n_frame
          call finalise(this%at(n_con))
       enddo
       deallocate(this%at)
       this%n_frame = 0
    endif

    inquire(file=this%at_file, exist=file_exists)
    if( .not. file_exists ) then
       call system_abort("read_fit_xyz: at_file "//this%at_file//" could not be found")
    endif

    call initialise(xyzfile,this%at_file,mpi=this%mpi_obj)
    this%n_frame = xyzfile%n_frame

    allocate(this%at(this%n_frame))

    do n_con = 1, this%n_frame
       call read(xyzfile,this%at(n_con),frame=n_con-1)
       call set_cutoff(this%at(n_con), this%max_cutoff)
       call calc_connect(this%at(n_con))
    enddo

    call finalise(xyzfile)

    if(this%n_frame <= 0) then
      call system_abort("read_fit_xyz: "//this%n_frame//" frames read from "//this%at_file//".")
   endif

  endsubroutine read_fit_xyz

  subroutine read_descriptors(this)

    type(gap_fit), intent(inout) :: this

    integer :: i

    this%max_cutoff = 0.0_dp

    if(allocated(this%my_descriptor)) then
       do i = 1, size(this%my_descriptor)
          call finalise(this%my_descriptor(i))
       enddo
       deallocate(this%my_descriptor)
    endif

    allocate(this%my_descriptor(this%n_coordinate))
    do i = 1, this%n_coordinate
       call initialise(this%my_descriptor(i),this%gap_str(i))
       if( this%max_cutoff < cutoff(this%my_descriptor(i)) ) this%max_cutoff = cutoff(this%my_descriptor(i))
    enddo

  endsubroutine read_descriptors

  subroutine fit_n_from_xyz(this)

    type(gap_fit), intent(inout) :: this

    logical :: do_collect_tasks, do_filter_tasks

    type(Atoms) :: at

    integer :: n_con
    logical :: has_ener, has_force, has_virial, has_stress_3_3, has_stress_voigt, has_hessian, has_local_property, has_force_mask, exclude_atom
    real(dp) :: ener, virial(3,3), stress_3_3(3,3)
    real(dp) :: stress_voigt(6)
    real(dp), pointer, dimension(:,:) :: f, hessian_eigenvector_j
    real(dp), pointer, dimension(:) :: local_property
    logical, pointer, dimension(:) :: force_mask
    integer :: i, j, k
    integer :: n_descriptors, n_cross, n_hessian
    integer :: n_current, n_last

    do_collect_tasks = (this%task_manager%active .and. .not. this%task_manager%distributed)
    do_filter_tasks = (this%task_manager%active .and. this%task_manager%distributed)

    if (allocated(this%n_cross)) deallocate(this%n_cross)
    if (allocated(this%n_descriptors)) deallocate(this%n_descriptors)
    allocate(this%n_cross(this%n_coordinate))
    allocate(this%n_descriptors(this%n_coordinate))

    this%n_cross = 0
    this%n_descriptors = 0
    this%n_ener = 0
    this%n_force = 0
    this%n_virial = 0
    this%n_hessian = 0
    this%n_local_property = 0
    n_last = 0

    do n_con = 1, this%n_frame
       if (do_filter_tasks) then
          if (this%task_manager%tasks(n_con)%worker_id /= this%task_manager%my_worker_id) cycle
       end if

       has_ener = get_value(this%at(n_con)%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(this%at(n_con),this%force_parameter_name, f)
       has_virial = get_value(this%at(n_con)%params,this%virial_parameter_name,virial)
       has_stress_voigt = get_value(this%at(n_con)%params,this%stress_parameter_name,stress_voigt)
       has_stress_3_3 = get_value(this%at(n_con)%params,this%stress_parameter_name,stress_3_3)
       has_hessian = get_value(this%at(n_con)%params,"n_"//this%hessian_parameter_name,n_hessian)
       has_local_property = assign_pointer(this%at(n_con),this%local_property_parameter_name, local_property)
       has_force_mask = assign_pointer(this%at(n_con),trim(this%force_mask_parameter_name),force_mask)

       if( has_ener ) then
          this%n_ener = this%n_ener + 1
       endif

       if( has_force ) then
          do i = 1, this%at(n_con)%N
             exclude_atom = .false.
             if(has_force_mask) exclude_atom = force_mask(i)

             if( .not. exclude_atom ) this%n_force = this%n_force + 3
          enddo
       endif

       if( has_stress_voigt .or. has_stress_3_3 ) then
          if( has_stress_voigt .and. has_stress_3_3 ) then
              call system_abort("fit_n_from_xyz: conflict in stress between 6-vector and 9-vector (really 3x3 matrix)")
          endif
          ! if has_stress is true, virial is available whether or not virial
          ! field has been detected
          has_virial = .true.
       endif

       if( has_virial ) then
          this%n_virial = this%n_virial + 6
       endif

       if( has_hessian ) then
          this%n_hessian = this%n_hessian + n_hessian
          at = this%at(n_con)
       endif

       if( has_local_property ) then
          this%n_local_property = this%n_local_property + this%at(n_con)%N
       endif

       if( has_local_property .and. ( has_ener .or. has_force .or. has_virial .or. has_hessian ) ) then
          call system_abort("fit_n_from_xyz: local_property and (energy or force or virial or hessian) present in configuration, currently not allowed.")
       endif

       do i = 1, this%n_coordinate
          call descriptor_sizes(this%my_descriptor(i),this%at(n_con),n_descriptors,n_cross)

          if( has_force ) then
             this%n_cross(i) = this%n_cross(i) + n_cross*3
          endif

          if( has_virial ) then
             this%n_cross(i) = this%n_cross(i) + n_cross*6
          endif

          this%n_descriptors(i) = this%n_descriptors(i) + n_descriptors

          if( has_hessian ) then
             do j = 1, n_hessian
                if( .not. assign_pointer(this%at(n_con),trim(this%hessian_parameter_name)//j, hessian_eigenvector_j) ) &
                   call system_abort("fit_n_from_xyz: could not find the "//j//"th of "//n_hessian//" hessian eigenvector")

                hessian_eigenvector_j = hessian_eigenvector_j / sqrt( sum(hessian_eigenvector_j**2) )

                do k = -1, 1, 2
                   at%pos = this%at(n_con)%pos + k * this%hessian_delta * hessian_eigenvector_j
                   call set_cutoff(at,this%max_cutoff)
                   call calc_connect(at)
                   call descriptor_sizes(this%my_descriptor(i),at,n_descriptors,n_cross)

                   this%n_descriptors(i) = this%n_descriptors(i) + n_descriptors
                   this%n_cross(i) = this%n_cross(i) + n_descriptors
                enddo

             enddo
          endif
       enddo

       if (do_collect_tasks) then
         n_current = this%n_ener + this%n_local_property + this%n_force + this%n_virial + this%n_hessian
         call task_manager_add_task(this%task_manager, n_current - n_last)
         n_last = n_current
       end if

       call finalise(at)
    enddo

    if (.not. do_filter_tasks) then
      call print_title("Report on number of descriptors found")
      do i = 1, this%n_coordinate
         call print("---------------------------------------------------------------------")
         call print("Descriptor "//i//": "//this%gap_str(i))
         call print("Number of descriptors:                        "//this%n_descriptors(i))
         call print("Number of partial derivatives of descriptors: "//this%n_cross(i))
      enddo
      call print_title("")
    end if

  end subroutine fit_n_from_xyz

  subroutine fit_data_from_xyz(this,error)

    type(gap_fit), intent(inout) :: this
    integer, optional, intent(out) :: error

    logical :: do_filter_tasks

    type(inoutput) :: theta_inout
    type(descriptor_data) :: my_descriptor_data

    type(Atoms) :: at
    integer :: d
    integer :: n_con
    logical :: has_ener, has_force, has_virial, has_stress_voigt, has_stress_3_3, has_hessian, has_local_property, &
       has_config_type, has_energy_sigma, has_force_sigma, has_virial_sigma, has_hessian_sigma, &
       has_force_atom_sigma, has_force_component_sigma, has_local_property_sigma, has_force_mask, exclude_atom
    real(dp) :: ener, ener_core, my_cutoff, energy_sigma, force_sigma, virial_sigma, hessian_sigma, local_property_sigma, &
       grad_covariance_cutoff, use_force_sigma
    real(dp), dimension(3) :: pos
    real(dp), dimension(3,3) :: virial, virial_core, stress_3_3
    real(dp), dimension(6) :: stress_voigt
    real(dp), dimension(:), allocatable :: theta, theta_fac, hessian, hessian_core, grad_data
    real(dp), dimension(:), pointer :: force_atom_sigma
    real(dp), dimension(:,:), pointer :: f, hessian_eigenvector_i, f_hessian, force_component_sigma
    real(dp), dimension(:), pointer :: local_property
    logical, dimension(:), pointer :: force_mask
    real(dp), dimension(:,:), allocatable :: f_core
    integer, dimension(:,:), allocatable :: force_loc, permutations
    integer :: ie, i, j, n, k, l, i_coordinate, n_hessian, n_energy_sigma, n_force_sigma, n_force_atom_sigma, &
    n_force_component_sigma, n_hessian_sigma, n_virial_sigma, n_local_property_sigma, n_descriptors
    integer, dimension(:), allocatable :: xloc, hessian_loc, local_property_loc
    integer, dimension(3,3) :: virial_loc

    integer :: i_config_type, n_config_type, n_theta_fac
    character(len=STRING_LENGTH) :: config_type
    character(len=THETA_LENGTH) :: theta_string
    character(len=STRING_LENGTH), dimension(:), allocatable :: theta_string_array

    INIT_ERROR(error)

    do_filter_tasks = (this%task_manager%active .and. this%task_manager%distributed)

    my_cutoff = 0.0_dp
    call gp_setParameters(this%my_gp,this%n_coordinate,this%n_ener+this%n_local_property,this%n_force+this%n_virial+this%n_hessian,this%sparse_jitter)

    do i_coordinate = 1, this%n_coordinate
       d = descriptor_dimensions(this%my_descriptor(i_coordinate))

       call gp_setParameters(this%my_gp,i_coordinate, d, this%n_descriptors(i_coordinate), this%n_cross(i_coordinate), this%delta(i_coordinate), this%f0(i_coordinate), &
                      covariance_type=this%covariance_type(i_coordinate) )
       call gp_addDescriptor(this%my_gp,i_coordinate,trim(this%gap_str(i_coordinate)))

       allocate(permutations(d,descriptor_n_permutations(this%my_descriptor(i_coordinate))))
       call descriptor_permutations(this%my_descriptor(i_coordinate),permutations)
       call gp_setPermutations(this%my_gp,i_coordinate,permutations)
       deallocate(permutations)

       my_cutoff = max(my_cutoff,cutoff(this%my_descriptor(i_coordinate)))
    enddo

    call print_title("Report on number of target properties found in training XYZ:")
    call print("Number of target energies (property name: "//trim(this%energy_parameter_name)//") found: "//sum(this%task_manager%MPI_obj, this%n_ener))
    call print("Number of target local_properties (property name: "//trim(this%local_property_parameter_name)//") found: "//sum(this%task_manager%MPI_obj, this%n_local_property))
    call print("Number of target forces (property name: "//trim(this%force_parameter_name)//") found: "//sum(this%task_manager%MPI_obj, this%n_force))
    call print("Number of target virials (property name: "//trim(this%virial_parameter_name)//") found: "//sum(this%task_manager%MPI_obj, this%n_virial))
    call print("Number of target Hessian eigenvalues (property name: "//trim(this%hessian_parameter_name)//") found: "//sum(this%task_manager%MPI_obj, this%n_hessian))
    call print_title("End of report")

    if( this%do_core ) call Initialise(this%core_pot, args_str=this%core_ip_args, param_str=string(this%quip_string))

    n_energy_sigma = 0
    n_force_sigma = 0
    n_force_atom_sigma = 0
    n_force_component_sigma = 0
    n_hessian_sigma = 0
    n_virial_sigma = 0
    n_local_property_sigma = 0

    do n_con = 1, this%n_frame
       if (do_filter_tasks) then
          if (this%task_manager%tasks(n_con)%worker_id /= this%task_manager%my_worker_id) cycle
       end if

       has_ener = get_value(this%at(n_con)%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(this%at(n_con),this%force_parameter_name, f)
       has_virial = get_value(this%at(n_con)%params,this%virial_parameter_name,virial)
       has_stress_voigt = get_value(this%at(n_con)%params,this%stress_parameter_name,stress_voigt)
       has_stress_3_3 = get_value(this%at(n_con)%params,this%stress_parameter_name,stress_3_3)
       has_hessian = get_value(this%at(n_con)%params,"n_"//this%hessian_parameter_name,n_hessian)
       has_config_type = get_value(this%at(n_con)%params,this%config_type_parameter_name,config_type)
       has_local_property = assign_pointer(this%at(n_con),this%local_property_parameter_name,local_property)

       has_energy_sigma = get_value(this%at(n_con)%params,'energy_'//trim(this%sigma_parameter_name),energy_sigma)
       has_force_sigma = get_value(this%at(n_con)%params,'force_'//trim(this%sigma_parameter_name),force_sigma)
       has_virial_sigma = get_value(this%at(n_con)%params,'virial_'//trim(this%sigma_parameter_name),virial_sigma)
       has_hessian_sigma = get_value(this%at(n_con)%params,'hessian_'//trim(this%sigma_parameter_name),hessian_sigma)
       has_force_atom_sigma = assign_pointer(this%at(n_con),'force_atom_'//trim(this%sigma_parameter_name),force_atom_sigma)
       has_force_component_sigma = assign_pointer(this%at(n_con),'force_component_'//trim(this%sigma_parameter_name),force_component_sigma)
       has_local_property_sigma = get_value(this%at(n_con)%params,'local_property_'//trim(this%sigma_parameter_name),local_property_sigma)
       has_force_mask = assign_pointer(this%at(n_con),trim(this%force_mask_parameter_name),force_mask)

       if ((.not. has_virial) .and. (has_stress_3_3 .or. has_stress_voigt)) then
          if (has_stress_voigt) then
             virial(1,1) = stress_voigt(1)
             virial(2,2) = stress_voigt(2)
             virial(3,3) = stress_voigt(3)
             virial(2,3) = stress_voigt(4)
             virial(3,1) = stress_voigt(5)
             virial(1,2) = stress_voigt(6)
             virial(3,2) = virial(2,3)
             virial(1,3) = virial(3,1)
             virial(2,1) = virial(1,2)
          else if (has_stress_3_3) then
             virial = stress_3_3
          else
             call system_abort("Frame "//n_con//" has no virial and stress that is neither a 9-vector (3x3)"// &
                " nor 6-vector (Voigt)")
          endif
          virial = -virial * cell_volume(this%at(n_con))
          has_virial = .true.
       endif

       if( has_force_atom_sigma .and. has_force_component_sigma ) then
          call print_warning("Frame "//n_con//" contains both force_atom_"//trim(this%sigma_parameter_name)// &
             " and force_component_"//trim(this%sigma_parameter_name)//" parameters. Per-component values will be used.")
       endif

       if( has_hessian ) then
          allocate(hessian(n_hessian))
          do i = 1, n_hessian
             if( .not. get_value(this%at(n_con)%params,trim(this%hessian_parameter_name)//i,hessian(i)) ) &
             call system_abort("fit_data_from_xyz: did not find "//i//"th of "//n_hessian//" hessian element" )
          enddo
       endif

       if( has_config_type ) then
          config_type = trim(config_type)
       else
          config_type = "default"
       endif

       if( .not. allocated(this%config_type) ) call system_abort('config_type not allocated')
       n_config_type = 0
       do i_config_type = 1, size(this%config_type)
          if( trim(this%config_type(i_config_type)) == trim(config_type) ) n_config_type = i_config_type
       enddo

       if( n_config_type == 0 ) then ! get the number of the "default" type as default
          do i_config_type = 1, size(this%config_type)
             if( trim(this%config_type(i_config_type)) == "default" ) n_config_type = i_config_type
          enddo
       endif

       if( this%do_core ) then
          allocate( f_core(3,this%at(n_con)%N) )
          ener_core = 0.0_dp
          f_core = 0.0_dp
          virial_core = 0.0_dp

          if( this%at(n_con)%cutoff < max(cutoff(this%core_pot),my_cutoff) ) then
             call set_cutoff(this%at(n_con), max(cutoff(this%core_pot),my_cutoff))
             call calc_connect(this%at(n_con))
          endif

          if(has_virial .and. has_force) then
             call calc(this%core_pot,this%at(n_con),energy=ener_core,force=f_core,virial=virial_core)
          elseif(has_force) then
             call calc(this%core_pot,this%at(n_con),energy=ener_core,force=f_core)
          elseif(has_virial) then
             call calc(this%core_pot,this%at(n_con),energy=ener_core,virial=virial_core)
          else
             call calc(this%core_pot,this%at(n_con),energy=ener_core)
          end if

          if(has_hessian) then
             allocate( hessian_core(n_hessian), f_hessian(3,this%at(n_con)%N) )
             hessian_core = 0.0_dp
             at = this%at(n_con)
             call set_cutoff(at, cutoff(this%core_pot))
             do i = 1, n_hessian
                if( .not. assign_pointer(this%at(n_con),trim(this%hessian_parameter_name)//i, hessian_eigenvector_i) ) &
                call system_abort("fit_data_from_xyz: could not find "//i//"th of "//n_hessian//" hessian eigenvector.")

                hessian_eigenvector_i = hessian_eigenvector_i / sqrt( sum(hessian_eigenvector_i**2) )

                do j = -1, 1, 2
                   at%pos = this%at(n_con)%pos + j * this%hessian_delta * hessian_eigenvector_i
                   call calc_connect(at)
                   call calc(this%core_pot,at,force = f_hessian)
                   hessian_core(i) = hessian_core(i) + j * sum(f_hessian*hessian_eigenvector_i) / 2.0_dp / this%hessian_delta
                enddo
             enddo
             call finalise(at)

             hessian = hessian - hessian_core
             deallocate(hessian_core, f_hessian)
          endif

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core

          deallocate(f_core)
       endif

       if(has_ener) then
          do i = 1, this%at(n_con)%N
             ener = ener - this%e0(this%at(n_con)%Z(i))
          enddo
       endif

       if(has_local_property) then
          do i = 1, this%at(n_con)%N
             local_property(i) = local_property(i) - this%local_property0(this%at(n_con)%Z(i))
          enddo
       endif

       if( has_ener .and. has_local_property ) then
          RAISE_ERROR("fit_data_from_xyz: energy and local_property both present in configuration, currently not allowed.",error)
       endif

       if( this%at(n_con)%cutoff < my_cutoff ) then
          call set_cutoff(this%at(n_con),my_cutoff)
          call calc_connect(this%at(n_con))
       endif

       if( .not. has_energy_sigma ) then
          if( this%sigma_per_atom ) then
             energy_sigma = this%sigma(1,n_config_type)*sqrt(1.0_dp * this%at(n_con)%N)
          else
             energy_sigma = this%sigma(1,n_config_type)
          endif
       else
          n_energy_sigma = n_energy_sigma + 1
       endif

       if( .not. has_force_sigma ) then
          force_sigma = this%sigma(2,n_config_type)
       else
          n_force_sigma = n_force_sigma + 1
       endif

       if( .not. has_virial_sigma ) then
          if( this%sigma_per_atom ) then
             virial_sigma = this%sigma(3,n_config_type)*sqrt(1.0_dp * this%at(n_con)%N)
          else
             virial_sigma = this%sigma(3,n_config_type)
          endif
       else
          n_virial_sigma = n_virial_sigma + 1
       endif

       if( .not. has_hessian_sigma ) then
          hessian_sigma = this%sigma(4,n_config_type)
       else
          n_hessian_sigma = n_hessian_sigma + 1
       endif

       if( .not. has_local_property_sigma ) then
          local_property_sigma = this%default_local_property_sigma
       else
          n_local_property_sigma = n_local_property_sigma + 1
       endif

       if( has_ener ) then
          if( energy_sigma .feq. 0.0_dp ) then
             RAISE_ERROR("fit_data_from_xyz: too small energy_sigma ("//energy_sigma//"), should be greater than zero",error)
          endif
          ie = gp_addFunctionValue(this%my_gp,ener, energy_sigma)
       elseif( has_local_property ) then
          if( local_property_sigma .feq. 0.0_dp ) then
             RAISE_ERROR("fit_data_from_xyz: too small local_property_sigma ("//local_property_sigma//"), should be greater than zero",error)
          endif
          allocate(local_property_loc(this%at(n_con)%N))
          do i = 1, this%at(n_con)%N
             local_property_loc(i) = gp_addFunctionValue(this%my_gp,local_property(i),local_property_sigma)
          enddo
       endif

       if(has_force) then
          allocate(force_loc(3,this%at(n_con)%N))
          do i = 1, this%at(n_con)%N
             if (has_force_component_sigma) then
                n_force_component_sigma = n_force_component_sigma + 3
                use_force_sigma = huge(1.0_dp) ! Updated later, below
             elseif (has_force_atom_sigma) then
                use_force_sigma = force_atom_sigma(i)
                n_force_atom_sigma = n_force_atom_sigma + 1
             else
                use_force_sigma = force_sigma
             endif

             if( use_force_sigma .feq. 0.0_dp ) then
                RAISE_ERROR("fit_data_from_xyz: too small force_sigma ("//use_force_sigma//"), should be greater than zero",error)
             endif

             exclude_atom = .false.
             if(has_force_mask) exclude_atom = force_mask(i)

             if( exclude_atom ) then
                force_loc(:,i) = EXCLUDE_LOC
             else
                do k = 1, 3
                   if( has_force_component_sigma ) use_force_sigma = force_component_sigma(k,i)
                   force_loc(k,i) = gp_addFunctionDerivative(this%my_gp,-f(k,i),use_force_sigma)
                enddo
             endif
          enddo
       endif
       if(has_virial) then
          ! check if virial is symmetric
          if( sum((virial - transpose(virial))**2) .fne. 0.0_dp ) &
          call print_warning('virial not symmetric, now symmetrised')

          ! Now symmetrise matrix
          virial = ( virial + transpose(virial) ) / 2.0_dp

          if( virial_sigma .feq. 0.0_dp ) then
             RAISE_ERROR("fit_data_from_xyz: too small virial_sigma ("//virial_sigma//"), should be greater than zero",error)
          endif

          do k = 1, 3
             do l = k, 3
                virial_loc(l,k) = gp_addFunctionDerivative(this%my_gp,-virial(l,k),virial_sigma)
             enddo
          enddo
       endif

       if(has_hessian) then
          if( hessian_sigma .feq. 0.0_dp ) then
             RAISE_ERROR("fit_data_from_xyz: too small hessian_sigma ("//hessian_sigma//"), should be greater than zero",error)
          endif

          allocate(hessian_loc(n_hessian))
          do i = 1, n_hessian
             hessian_loc(i) = gp_addFunctionDerivative(this%my_gp,hessian(i),hessian_sigma)
          enddo
       endif

       n_descriptors = 0
       do i_coordinate = 1, this%n_coordinate

          call calc(this%my_descriptor(i_coordinate),this%at(n_con),my_descriptor_data, &
          do_descriptor=.true.,do_grad_descriptor=has_force .or. has_virial)

          allocate(xloc(size(my_descriptor_data%x)))
          n_descriptors = n_descriptors + size(my_descriptor_data%x)

          if( has_ener ) then
             do i = 1, size(my_descriptor_data%x)
                if( .not. my_descriptor_data%x(i)%has_data) cycle
                xloc(i) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(i)%data(:),i_coordinate, &
                cutoff_in=my_descriptor_data%x(i)%covariance_cutoff, current_y=ie,config_type=n_config_type)
             enddo
          elseif( has_local_property ) then
             do i = 1, size(my_descriptor_data%x)
                if( .not. my_descriptor_data%x(i)%has_data) cycle
                xloc(i) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(i)%data(:),i_coordinate, &
                cutoff_in=my_descriptor_data%x(i)%covariance_cutoff, current_y=local_property_loc(my_descriptor_data%x(i)%ci(1)),config_type=n_config_type)
             enddo
          else
             do i = 1, size(my_descriptor_data%x)
                if( .not. my_descriptor_data%x(i)%has_data) cycle
                xloc(i) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(i)%data(:),i_coordinate, &
                cutoff_in=my_descriptor_data%x(i)%covariance_cutoff, config_type=n_config_type)
             enddo
          endif


          if(has_force) then
             do i = 1, size(my_descriptor_data%x)
                do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
                   if( .not. my_descriptor_data%x(i)%has_grad_data(n)) cycle
                   j = my_descriptor_data%x(i)%ii(n)

                   do k = 1, 3
                      if( force_loc(k,j) > EXCLUDE_LOC ) then
                         call gp_addCoordinateDerivatives(this%my_gp,my_descriptor_data%x(i)%grad_data(:,k,n),i_coordinate, &
                         force_loc(k,j), xloc(i), dcutoff_in=my_descriptor_data%x(i)%grad_covariance_cutoff(k,n) )
                      endif
                   enddo
                enddo
             enddo

          endif

          if(has_virial) then
             do k = 1, 3
                do l = k, 3

                   do i = 1, size(my_descriptor_data%x)
                      do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
                         if( .not. my_descriptor_data%x(i)%has_grad_data(n)) cycle
                         j = my_descriptor_data%x(i)%ii(n)
                         pos = my_descriptor_data%x(i)%pos(:,n)
                         call gp_addCoordinateDerivatives(this%my_gp,my_descriptor_data%x(i)%grad_data(:,k,n)*pos(l), i_coordinate, &
                         virial_loc(l,k), xloc(i), dcutoff_in=my_descriptor_data%x(i)%grad_covariance_cutoff(k,n)*pos(l))
                      enddo
                   enddo

                enddo
             enddo
          endif

          if(allocated(xloc)) deallocate(xloc)
       enddo

       if( has_local_property ) then
          if( n_descriptors /= this%at(n_con)%N ) then
             RAISE_ERROR("fit_data_from_xyz: local_propertyes found in configuration, but number of descriptors do not match &
                & the number of atoms. Check your descriptors.",error)
          endif
       endif

       if(allocated(force_loc)) deallocate(force_loc)
       if(allocated(local_property_loc)) deallocate(local_property_loc)

       if( has_hessian ) then
          at = this%at(n_con)
          call set_cutoff( at, my_cutoff )
          do i_coordinate = 1, this%n_coordinate
             allocate( grad_data(descriptor_dimensions(this%my_descriptor(i_coordinate))) )
             
             do i = 1, n_hessian
                if( .not. assign_pointer(this%at(n_con),trim(this%hessian_parameter_name)//i, hessian_eigenvector_i) ) &
                call system_abort("fit_data_from_xyz: could not find "//i//"th of "//n_hessian//" hessian eigenvector.")
                
                do j = -1, 1, 2
                   at%pos = this%at(n_con)%pos + j * this%hessian_delta * hessian_eigenvector_i
                   call calc_connect(at)

                   call calc(this%my_descriptor(i_coordinate),at,my_descriptor_data, &
                   do_descriptor=.true.,do_grad_descriptor=.true.)
                   !hessian_core(i) = hessian_core(i) + j * sum(f_hessian*hessian_eigenvector_i) / 2.0_dp / this%hessian_delta

                   allocate(xloc(size(my_descriptor_data%x)))
                   
                   do k = 1, size(my_descriptor_data%x)
                      if( .not. my_descriptor_data%x(k)%has_data) cycle
                      xloc(k) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(k)%data(:),i_coordinate, &
                      cutoff_in=my_descriptor_data%x(k)%covariance_cutoff,config_type=EXCLUDE_CONFIG_TYPE)
                      !cutoff_in=my_descriptor_data%x(k)%covariance_cutoff,config_type=n_config_type)


                      grad_data = 0.0_dp
                      grad_covariance_cutoff = 0.0_dp
                      do n = lbound(my_descriptor_data%x(k)%ii,1), ubound(my_descriptor_data%x(k)%ii,1)
                         if( .not. my_descriptor_data%x(k)%has_grad_data(n)) cycle
                         l = my_descriptor_data%x(k)%ii(n)
                         grad_data = grad_data + j * matmul(my_descriptor_data%x(k)%grad_data(:,:,n), hessian_eigenvector_i(:,l)) / 2.0_dp / this%hessian_delta
                         grad_covariance_cutoff = grad_covariance_cutoff + &
                         dot_product(my_descriptor_data%x(k)%grad_covariance_cutoff(:,n), hessian_eigenvector_i(:,l)) / 2.0_dp / this%hessian_delta
                      enddo
                      call gp_addCoordinateDerivatives(this%my_gp, grad_data, i_coordinate, &
                      hessian_loc(i), xloc(k), dcutoff_in=grad_covariance_cutoff)

                   enddo !k
                   
                   deallocate(xloc)
                enddo !j = -1, 1, 2
             enddo ! i = 1, n_hessian
             if(allocated(grad_data)) deallocate(grad_data)
          enddo ! i_coordinate = 1, n_coordinate
       endif !has_hessian

       if(allocated(hessian_loc)) deallocate(hessian_loc)
       if(allocated(hessian)) deallocate(hessian)
       call finalise(my_descriptor_data)
    enddo !n_frame

    call print_title("Report on per-configuration/per-atom sigma (error parameter) settings")
    call print("Number of per-configuration setting of energy_"//trim(this%sigma_parameter_name)//" found:     "//sum(this%task_manager%MPI_obj, n_energy_sigma))
    call print("Number of per-configuration setting of force_"//trim(this%sigma_parameter_name)//" found:      "//sum(this%task_manager%MPI_obj, n_force_sigma))
    call print("Number of per-configuration setting of virial_"//trim(this%sigma_parameter_name)//" found:     "//sum(this%task_manager%MPI_obj, n_virial_sigma))
    call print("Number of per-configuration setting of hessian_"//trim(this%sigma_parameter_name)//" found:    "//sum(this%task_manager%MPI_obj, n_hessian_sigma))
    call print("Number of per-configuration setting of local_propery_"//trim(this%sigma_parameter_name)//" found:"//sum(this%task_manager%MPI_obj, n_local_property_sigma))
    call print("Number of per-atom setting of force_atom_"//trim(this%sigma_parameter_name)//" found:          "//sum(this%task_manager%MPI_obj, n_force_atom_sigma))
    call print("Number of per-component setting of force_component_"//trim(this%sigma_parameter_name)//" found:          "//sum(this%task_manager%MPI_obj, n_force_component_sigma))
    call print_title("End of report")

    do i_coordinate = 1, this%n_coordinate
       if( count( (/this%has_theta_file(i_coordinate), this%has_theta_uniform(i_coordinate), &
       this%has_theta_fac(i_coordinate), this%has_zeta(i_coordinate) /) ) /= 1 ) then
          call system_abort("fit_data_from_xyz: only one of theta_file, theta_uniform, theta_fac or zeta may be &
          specified for each GAP.")
       endif
       if( this%covariance_type(i_coordinate) == COVARIANCE_DOT_PRODUCT ) then
          if( .not. this%has_zeta(i_coordinate) ) call system_abort("fit_data_from_xyz: covariance type is DOT_PRODUCT but no zeta was specified.")
       elseif( this%covariance_type(i_coordinate) == COVARIANCE_ARD_SE .or. this%covariance_type(i_coordinate) == COVARIANCE_PP ) then
          if( count( (/this%has_theta_file(i_coordinate), this%has_theta_uniform(i_coordinate), this%has_theta_fac(i_coordinate) /) ) /= 1 ) then
             call system_abort("fit_data_from_xyz: covariance type is Gaussian or PP, so one of theta_file, theta_uniform of theta_fac must be specified")
          endif
       endif

       if( this%has_theta_file(i_coordinate) ) then
          allocate(theta_string_array(this%my_gp%coordinate(i_coordinate)%d))
          allocate(theta(this%my_gp%coordinate(i_coordinate)%d))

          call initialise(theta_inout,trim(this%theta_file(i_coordinate)))
          read(theta_inout%unit,'(a)') theta_string
          call split_string(theta_string,' :;','{}',theta_string_array,d,matching=.true.)
          if(this%my_gp%coordinate(i_coordinate)%d /= d) call system_abort('File '//trim(this%theta_file(i_coordinate))//' does not contain the right number of hyperparameters')
          do i = 1, d
             theta(i) = string_to_real(trim(theta_string_array(i)))
          enddo
          call gp_setTheta(this%my_gp,i_coordinate,theta=theta)
          deallocate(theta_string_array)
          deallocate(theta)
          call finalise(theta_inout)
       elseif(this%has_theta_uniform(i_coordinate)) then
          allocate(theta(this%my_gp%coordinate(i_coordinate)%d))
          theta = this%theta_uniform(i_coordinate)
          call gp_setTheta(this%my_gp,i_coordinate,theta=theta)
          deallocate(theta)
       elseif(this%has_theta_fac(i_coordinate)) then
          allocate(theta_string_array(this%my_gp%coordinate(i_coordinate)%d))
          allocate(theta_fac(this%my_gp%coordinate(i_coordinate)%d))
          call split_string(trim(this%theta_fac_string(i_coordinate))," :;",'{}',theta_string_array,n_theta_fac,matching=.true.)

          if(n_theta_fac == 1) then
             theta_fac = string_to_real(theta_string_array(1))
          elseif(n_theta_fac == this%my_gp%coordinate(i_coordinate)%d) then
             do i = 1, this%my_gp%coordinate(i_coordinate)%d
                theta_fac(i) = string_to_real(theta_string_array(i))
             enddo
          else
             call system_abort("theta_fac can only contain one value or as many as dimensions the descriptor is")
          endif
          call gp_setThetaFactor(this%my_gp,i_coordinate,theta_fac,useSparseX=.false.)
       
          deallocate(theta_fac)
          deallocate(theta_string_array)
       elseif( this%has_zeta(i_coordinate) ) then
          call gp_setTheta(this%my_gp,i_coordinate,zeta=this%zeta(i_coordinate))
       endif
    enddo

    if( this%do_core ) call Finalise(this%core_pot)

    ! @info move this check to gp_sparsify when implementing more methods
    if (this%task_manager%active) then
      if (any(this%sparse_method /= GP_SPARSE_FILE)) then
         call system_abort("Only sparse_method FILE implemented for MPI.")
      end if
    end if

    call gp_sparsify(this%my_gp,n_sparseX=this%config_type_n_sparseX,default_all=(this%n_sparseX/=0), &
       sparseMethod=this%sparse_method, sparse_file=this%sparse_file, &
       use_actual_gpcov=this%sparse_use_actual_gpcov, print_sparse_index = this%print_sparse_index, &
       unique_hash_tolerance=this%unique_hash_tolerance, unique_descriptor_tolerance=this%unique_descriptor_tolerance)

  end subroutine fit_data_from_xyz

  subroutine e0_from_xyz(this)

    type(gap_fit), intent(inout) :: this

    integer :: n_con, n_ener, i, my_n_neighbours
    logical :: has_ener
    real(dp) :: ener, ener_core

    logical, dimension(total_elements) :: found_Z, found_isolated

    if( this%do_core ) call Initialise(this%core_pot, this%core_ip_args, param_str=string(this%quip_string))

    n_ener = 0

    this%e0 = 0.0_dp
    found_isolated = .false.
    found_Z = .false.

    do n_con = 1, this%n_frame

       has_ener = get_value(this%at(n_con)%params,trim(this%energy_parameter_name),ener)

       found_Z(this%at(n_con)%Z) = .true.

       if( has_ener ) then

          ener_core = 0.0_dp
          if( this%do_core ) then
             if( this%at(n_con)%cutoff < cutoff(this%core_pot) ) then
                call set_cutoff(this%at(n_con), cutoff(this%core_pot))
                call calc_connect(this%at(n_con))
             endif
             call calc(this%core_pot,this%at(n_con),energy=ener_core)
          endif

          select case(this%e0_method)
          case(E0_ISOLATED)
             if( this%at(n_con)%N == 1 ) then
                if( this%at(n_con)%cutoff < this%max_cutoff ) then
                   call set_cutoff(this%at(n_con), this%max_cutoff)
                endif
                call calc_connect(this%at(n_con))
                if( n_neighbours(this%at(n_con),1,max_dist = this%max_cutoff) == 0 ) then
                   if( found_isolated(this%at(n_con)%Z(1)) ) then
                      call system_abort("Found more than one isolated atom configuration, which may be ambiguous.")
                   endif
                   this%e0(this%at(n_con)%Z(1)) = ener - ener_core
                   found_isolated(this%at(n_con)%Z(1)) = .true.
                endif
             endif
          case(E0_AVERAGE)
             this%e0 = this%e0 + (ener-ener_core) / this%at(n_con)%N
          case default
             call system_abort("Unknown e0_method")
          endselect

          n_ener = n_ener + 1
       endif
    enddo

    select case(this%e0_method)
    case(E0_ISOLATED)
       if( .not. all(found_isolated .eqv. found_Z) ) then
          do i = 1, size(found_Z)
             if( found_Z(i) .and. .not. found_isolated(i) ) then
                call print("Atom species "//i//" present in teaching XYZ, but not found corresponding isolated &
                   representative")
             endif
          enddo
          call system_abort("Determination of e0 was requested to be based on isolated atom energies, but not all &
             & atom types present in the XYZ had an isolated representative.")
       endif
    case(E0_AVERAGE)
       if( n_ener > 0 ) then
          this%e0 = this%e0 / n_ener
       else
          this%e0 = 0.0_dp
       endif
    case default
       call system_abort("Unknown e0_method")
    endselect

    if( this%do_core ) call Finalise(this%core_pot)

  endsubroutine e0_from_xyz

  subroutine w_Z_from_xyz(this)

    type(gap_fit), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at

    call initialise(xyzfile,this%at_file,mpi=this%mpi_obj)

    call read(xyzfile,at,frame=0)
    !call get_weights(at,this%w_Z)
    call finalise(at)

    call finalise(xyzfile)

  end subroutine w_Z_from_xyz

  subroutine gap_fit_print_xml(this,filename,sparseX_separate_file)

     use iso_c_binding, only : C_NULL_CHAR

     type(gap_fit), intent(in) :: this
     character(len=*), intent(in) :: filename
     logical, intent(in), optional :: sparseX_separate_file

     type(xmlf_t) :: xf
     !type(extendable_str) :: gap_string
     !type(inoutput) :: gp_inout
     character(len=STRING_LENGTH) :: gp_tmp_file, gp_label
     integer :: i
     integer, dimension(8) :: values
     logical :: my_sparseX_separate_file

     call date_and_time(values=values)
     ! Get totally unique label for GAP. This will be used at various places.
     write(gp_label,'("GAP_"7(i0,"_")i0)') values

     ! Unique temporary file
     gp_tmp_file = 'tmp_'//trim(gp_label)//'.xml'

     ! Print GAP part of the potential into the temporary file.
     call xml_OpenFile(gp_tmp_file,xf,addDecl=.false.)

     call xml_NewElement(xf,"GAP_params")
     call xml_AddAttribute(xf,"label",trim(gp_label))
     call xml_AddAttribute(xf,"gap_version",""//gap_version)

     call xml_NewElement(xf,"GAP_data")
     call xml_AddAttribute(xf,"do_core",""//this%do_core)
     
     do i = 1, size(this%e0)
        call xml_NewElement(xf,"e0")
        call xml_AddAttribute(xf,"Z",""//i)
        call xml_AddAttribute(xf,"value",""// (this%e0(i)+this%local_property0(i) ))
        call xml_EndElement(xf,"e0")
     enddo

     call xml_EndElement(xf,"GAP_data")

     my_sparseX_separate_file = optional_default(.false., sparseX_separate_file)

     ! Print GP bit of the potential
     if (my_sparseX_separate_file) then
        call gp_printXML(this%gp_sp,xf,label=gp_label,sparseX_base_filename=trim(filename)//".sparseX")
     else
        call gp_printXML(this%gp_sp,xf,label=gp_label)
     endif

     ! Print the command line used for the fitting
     if(len(trim(this%command_line))> 0 ) then
        call xml_NewElement(xf,"command_line")
        call xml_AddCharacters(xf,trim(this%command_line),parsed=.false.)
        call xml_EndElement(xf,"command_line")
     endif

     if(this%do_copy_at_file) then
        ! Print the fitting configurations used for this potential.
        if(len(trim(this%at_file)) > 0 ) call file_print_xml(this%at_file,xf,ws_significant=.false.)
     endif

     call xml_EndElement(xf,"GAP_params")
     call xml_Close(xf)

     !! Now read back into an extendable string what we have just printed out.
     !call read(gap_string, trim(gp_tmp_file), keep_lf=.true.)

     !! Initialise the final file
     !call initialise(gp_inout,trim(filename),action=OUTPUT)

     ! Open a unique root element for the xml
     !call print('<'//trim(gp_label)//'>',file=gp_inout)
     !!call system_command('echo "<'//trim(gp_label)//'>" >>'//trim(filename))
     call fwrite_line_to_file(trim(filename)//C_NULL_CHAR,'<'//trim(gp_label)//'>'//C_NULL_CHAR,'w'//C_NULL_CHAR)

     if(this%do_core) then
        ! Create the sum potential xml entry (by hand)
        !call print('<Potential label="'//trim(gp_label)//'" init_args="Sum init_args_pot1={'//trim(this%ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}"/>',file=gp_inout)
        !call system_command('echo "<Potential label=\"'//trim(gp_label)//'\" init_args=\"Sum init_args_pot1={'//trim(this%ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}\"/>" >>'//trim(filename))
        call fwrite_line_to_file(trim(filename)//C_NULL_CHAR, &
           '<Potential label="'//trim(gp_label)//'" init_args="Sum init_args_pot1={'//trim(this%core_ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}"/>'//C_NULL_CHAR, &
           'a'//C_NULL_CHAR)

        ! Now add the core potential that was used.
        !call print(string(this%quip_string),file=gp_inout)
        !call system_command('echo "'//string(this%quip_string)//' >>'//trim(filename))
        call fappend_file_to_file(trim(filename)//C_NULL_CHAR,trim(this%core_param_file)//C_NULL_CHAR)
     else
        call fwrite_line_to_file(trim(filename)//C_NULL_CHAR, &
           '<Potential label="'//trim(gp_label)//'" init_args="IP GAP label='//trim(gp_label)//'"/>'//C_NULL_CHAR,'a'//C_NULL_CHAR)
     endif

     ! Add the GAP potential
     !call print(string(gap_string),file=gp_inout)
     !call system_command('cat '//trim(gp_tmp_file)//' >>'//trim(filename))
     call fappend_file_to_file(trim(filename)//C_NULL_CHAR,trim(gp_tmp_file)//C_NULL_CHAR)

     ! Close the root element
     !call print('</'//trim(gp_label)//'>',file=gp_inout)
     !call system_command('echo "</'//trim(gp_label)//'>" >>'//trim(filename))
     call fwrite_line_to_file(trim(filename)//C_NULL_CHAR,'</'//trim(gp_label)//'>'//C_NULL_CHAR,'a'//C_NULL_CHAR)

     !call finalise(gp_inout)
     !call finalise(gap_string)

     ! Delete the temporary file
     !call system_command('rm -f '//trim(gp_tmp_file))
     call frm_file(trim(gp_tmp_file)//C_NULL_CHAR)
     

  endsubroutine gap_fit_print_xml

  subroutine file_print_xml(this,xf,ws_significant)
     character(len=*), intent(in) :: this
     type(xmlf_t), intent(inout) :: xf
     logical, intent(in), optional :: ws_significant

     type(inoutput) :: atfile
     character(len=10240) :: line
     integer :: iostat

     call initialise(atfile,trim(this))
     call xml_NewElement(xf,"XYZ_data")
     call xml_AddNewLine(xf)

     do
        read(atfile%unit,'(a)',iostat=iostat) line
        if(iostat < 0) then
           exit
        elseif(iostat > 0) then
           call system_abort('file_print_xml: unkown error ('//iostat//') while reading '//trim(this))
        endif
        call xml_AddCharacters(xf,trim(line),parsed=.false.,ws_significant=ws_significant)
        call xml_AddNewLine(xf)
     enddo
     call xml_EndElement(xf,"XYZ_data")
     call finalise(atfile)

  endsubroutine file_print_xml

!  subroutine print_sparse(this)
!    type(gap_fit), intent(in) :: this
!    type(cinoutput) :: xyzfile, xyzfile_out
!    type(atoms) :: at, at_out
!
!    integer :: li, ui, n_con
!    logical, dimension(:), allocatable :: x
!    logical, dimension(:), pointer :: sparse
!
!    if(this%do_mark_sparse_atoms) then
!
!       allocate(x(this%n_descriptors))
!       x = .false.
!       x(this%r) = .true.
!
!       call initialise(xyzfile,this%at_file)
!       call initialise(xyzfile_out,this%mark_sparse_atoms,action=OUTPUT)
!
!       li = 0
!       ui = 0
!       do n_con = 1, xyzfile%n_frame
!          call read(xyzfile,at,frame=n_con-1)
!          at_out = at
!
!          call add_property(at_out,'sparse',.false.,ptr=sparse)
!
!          li = ui + 1
!          ui = ui + at%N
!          if(any( x(li:ui) )) sparse(find_indices(x(li:ui))) = .true.
!
!          call write(at_out,xyzfile_out,properties="species:pos:sparse")
!       enddo
!       call finalise(xyzfile)
!       call finalise(xyzfile_out)
!       deallocate(x)
!
!    endif
!
!  endsubroutine print_sparse

  subroutine parse_config_type_sigma(this)
    type(gap_fit), intent(inout) :: this
    character(len=STRING_LENGTH), dimension(200) :: config_type_sigma_fields
    integer :: config_type_sigma_num_fields, i_default, i, n_config_type

    if( this%has_config_type_sigma ) then
       call split_string(this%config_type_sigma_string,' :;','{}',config_type_sigma_fields,config_type_sigma_num_fields,matching=.true.)

       n_config_type = config_type_sigma_num_fields / 5

       ! find "default" if present
       i_default = 0
       do i = 1, config_type_sigma_num_fields, 5
          if( trim(config_type_sigma_fields(i)) == "default" ) i_default = i
       enddo

       if( i_default == 0 ) then
          ! no default present in the string, we add it, and it'll be the last one
          n_config_type = n_config_type + 1
          i_default = n_config_type
          config_type_sigma_fields(config_type_sigma_num_fields+1) = "default"
          config_type_sigma_fields(config_type_sigma_num_fields+2) = ""//this%default_sigma(1)
          config_type_sigma_fields(config_type_sigma_num_fields+3) = ""//this%default_sigma(2)
          config_type_sigma_fields(config_type_sigma_num_fields+4) = ""//this%default_sigma(3)
          config_type_sigma_fields(config_type_sigma_num_fields+5) = ""//this%default_sigma(4)
          config_type_sigma_num_fields = config_type_sigma_num_fields + 5
       endif

       allocate(this%config_type(n_config_type))
       allocate(this%sigma(4,n_config_type))

       do i = 1, n_config_type 
          this%config_type(i) = trim(config_type_sigma_fields(5*(i-1)+1))
          this%sigma(1,i) = string_to_real(config_type_sigma_fields(5*(i-1)+2))
          this%sigma(2,i) = string_to_real(config_type_sigma_fields(5*(i-1)+3))
          this%sigma(3,i) = string_to_real(config_type_sigma_fields(5*(i-1)+4))
          this%sigma(4,i) = string_to_real(config_type_sigma_fields(5*(i-1)+5))
       enddo

       call print('Sparse points and target errors per pre-defined types of configurations')
       do i = 1, n_config_type
          call print(""//trim(this%config_type(i))//"  "//this%sigma(:,i))
       enddo
    else
       allocate(this%config_type(1))
       allocate(this%sigma(4,1))
       this%config_type(1)= "default"
       this%sigma(:,1) = this%default_sigma
    endif

  endsubroutine parse_config_type_sigma

  subroutine parse_config_type_n_sparseX(this)
    type(gap_fit), intent(inout) :: this

    integer :: i, j, i_default, i_coordinate, i_config_type, config_type_n_sparseX_num_fields, n_config_type, new_config_types
    character(len=STRING_LENGTH), dimension(200) :: config_type_n_sparseX_fields
    logical :: config_type_present

    if( .not. allocated(this%config_type) ) call system_abort('config_type not allocated, call parse_config_type_sigma first')

    do i = 1, size(this%config_type)
       if( trim(this%config_type(i)) == "default" ) i_default = i
    enddo

    ! Check first if we have more new config types than we had from config_type_sigma
    do i_coordinate = 1, this%n_coordinate
       if( this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0) then
          call split_string(this%config_type_n_sparseX_string(i_coordinate),' :;','{}',config_type_n_sparseX_fields,config_type_n_sparseX_num_fields,matching=.true.)

          if( mod(config_type_n_sparseX_num_fields,2) /= 0 ) then
             call system_abort("parse_config_type_n_sparseX: config_type_n_sparseX could not be parsed correctly, key/value pairs must always be present")
          endif

          n_config_type = size(this%config_type)
          new_config_types = 0 ! Assume there are no new config_types
          do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
             config_type_present = .false.
             do i = 1, n_config_type ! loop over config_types previously set
                if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) config_type_present = .true. ! Found config_type among old ones
             enddo
             if(.not.config_type_present) new_config_types = new_config_types + 1 ! Increment as it's a genuine new config_type
          enddo
          if( new_config_types > 0 ) then
             call reallocate(this%config_type, n_config_type + new_config_types, copy=.true.)
             call reallocate(this%sigma,4,n_config_type + new_config_types, copy=.true.)

             i_config_type = n_config_type
             do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
                config_type_present = .false.
                do i = 1, n_config_type ! loop over config_types previously set
                   if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) config_type_present = .true. ! Found config_type among old ones
                enddo
                if(.not.config_type_present) then ! it's a genuine new config_type
                   i_config_type = i_config_type + 1
                   this%config_type(i_config_type) = trim(config_type_n_sparseX_fields(j))
                   this%sigma(:,i_config_type) = this%sigma(:,i_default)
                endif
             enddo
          endif

       elseif(this%n_sparseX(i_coordinate) > 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0 .and. len_trim(this%sparse_file(i_coordinate)) ==0 ) then
          call system_abort('Confused: cannot specify both n_sparse and config_type_n_sparse')


       elseif(this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) == 0  .and. len_trim(this%sparse_file(i_coordinate)) == 0) then
          call system_abort('Confused: either n_sparse or config_type_n_sparse has to be specified')
       endif

    enddo

    n_config_type = size(this%config_type)
    allocate(this%config_type_n_sparseX(n_config_type,this%n_coordinate))
    this%config_type_n_sparseX = 0

    do i_coordinate = 1, this%n_coordinate
       if( this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0) then
          call split_string(this%config_type_n_sparseX_string(i_coordinate),' :;','{}',config_type_n_sparseX_fields,config_type_n_sparseX_num_fields,matching=.true.)

          do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
             do i = 1, n_config_type ! loop over config_types previously set
                if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) &
                   this%config_type_n_sparseX(i,i_coordinate) = string_to_int( config_type_n_sparseX_fields(j+1) )
             enddo
          enddo
          !this%n_sparseX(i_coordinate) = sum( this%config_type_n_sparseX(:,i_coordinate) )

       elseif( this%n_sparseX(i_coordinate) > 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) == 0) then
          this%config_type_n_sparseX(i_default,i_coordinate) = this%n_sparseX(i_coordinate)
       endif
    enddo

  endsubroutine parse_config_type_n_sparseX

  subroutine get_species_xyz(this)
    type(gap_fit), intent(inout) :: this

    integer :: n_con, i
    integer, dimension(total_elements) :: species_present

    this%n_species = 0
    species_present = 0

    do n_con = 1, this%n_frame
       do i = 1, this%at(n_con)%N
          if( all(this%at(n_con)%Z(i) /= species_present) ) then
             this%n_species = this%n_species + 1
             species_present(this%n_species) = this%at(n_con)%Z(i)
          endif
       enddo
    enddo

    allocate(this%species_Z(this%n_species))
    this%species_Z = species_present(1:this%n_species)
    
  endsubroutine get_species_xyz

  subroutine add_multispecies_gaps(this)
    type(gap_fit), intent(inout) :: this

    integer :: i_coordinate, i, j, n_gap_str, i_add_species
    character(STRING_LENGTH), dimension(:), allocatable :: gap_str_i, new_gap_str

    ! temporary arrays
    real(dp), dimension(:), allocatable :: delta, f0, theta_uniform, zeta, unique_hash_tolerance, unique_descriptor_tolerance
    integer, dimension(:), allocatable :: n_sparseX, sparse_method, covariance_type
    character(len=STRING_LENGTH), dimension(:), allocatable :: theta_file, sparse_file, theta_fac_string, config_type_n_sparseX_string, print_sparse_index
    logical, dimension(:), allocatable :: mark_sparse_atoms, has_theta_fac, has_theta_uniform, has_theta_file, has_zeta

    n_gap_str = 0
    do i_coordinate = 1, this%n_coordinate
       if( this%add_species(i_coordinate) ) then

          call print('Old GAP: {'//trim(this%gap_str(i_coordinate))//'}')
          call descriptor_str_add_species(this%gap_str(i_coordinate),this%species_Z,gap_str_i)
          call reallocate(new_gap_str, n_gap_str+size(gap_str_i),copy=.true.)

          call reallocate(delta, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(f0, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(n_sparseX, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(config_type_n_sparseX_string, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(theta_fac_string, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(theta_uniform, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(theta_file, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(has_theta_fac, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(has_theta_uniform, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(has_theta_file, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(sparse_file, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(mark_sparse_atoms, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(sparse_method, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(covariance_type, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(zeta, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(has_zeta, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(print_sparse_index, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(unique_hash_tolerance, n_gap_str+size(gap_str_i),copy=.true.)
          call reallocate(unique_descriptor_tolerance, n_gap_str+size(gap_str_i),copy=.true.)

          do i = 1, size(gap_str_i)
             i_add_species = index(gap_str_i(i),'add_species')
             if(i_add_species /= 0) then
                do j = i_add_species, len_trim(gap_str_i(i))
                   if( gap_str_i(i)(j:j) == " " ) exit
                   gap_str_i(i)(j:j) = " "
                   !gap_str_i(i)(i_add_species:i_add_species+len('add_species')-1) = '           '
                enddo
             endif

             new_gap_str(i+n_gap_str) = trim(gap_str_i(i))
             call print('New GAP: {'//trim(gap_str_i(i))//'}')

             delta(i+n_gap_str) = this%delta(i_coordinate)
             f0(i+n_gap_str) = this%f0(i_coordinate)
             n_sparseX(i+n_gap_str) = this%n_sparseX(i_coordinate)
             config_type_n_sparseX_string(i+n_gap_str) = this%config_type_n_sparseX_string(i_coordinate)
             theta_fac_string(i+n_gap_str) = this%theta_fac_string(i_coordinate)
             theta_uniform(i+n_gap_str) = this%theta_uniform(i_coordinate)
             theta_file(i+n_gap_str) = this%theta_file(i_coordinate)

             has_theta_fac(i+n_gap_str) = this%has_theta_fac(i_coordinate)
             has_theta_uniform(i+n_gap_str) = this%has_theta_uniform(i_coordinate)
             has_theta_file(i+n_gap_str) = this%has_theta_file(i_coordinate)

             sparse_file(i+n_gap_str) = this%sparse_file(i_coordinate)
             mark_sparse_atoms(i+n_gap_str) = this%mark_sparse_atoms(i_coordinate)
             sparse_method(i+n_gap_str) = this%sparse_method(i_coordinate)
             covariance_type(i+n_gap_str) = this%covariance_type(i_coordinate)
             zeta(i+n_gap_str) = this%zeta(i_coordinate)
             has_zeta(i+n_gap_str) = this%has_zeta(i_coordinate)
             print_sparse_index(i+n_gap_str) = this%print_sparse_index(i_coordinate)
             unique_hash_tolerance(i+n_gap_str) = this%unique_hash_tolerance(i_coordinate)
             unique_descriptor_tolerance(i+n_gap_str) = this%unique_descriptor_tolerance(i_coordinate)

          enddo
          n_gap_str = n_gap_str + size(gap_str_i)
          deallocate(gap_str_i)

       else
          n_gap_str = n_gap_str + 1

          call reallocate(new_gap_str, n_gap_str,copy=.true.)
          call reallocate(delta, n_gap_str,copy=.true.)
          call reallocate(f0, n_gap_str,copy=.true.)
          call reallocate(n_sparseX, n_gap_str,copy=.true.)
          call reallocate(config_type_n_sparseX_string, n_gap_str,copy=.true.)
          call reallocate(theta_fac_string, n_gap_str,copy=.true.)
          call reallocate(theta_uniform, n_gap_str,copy=.true.)
          call reallocate(theta_file, n_gap_str,copy=.true.)

          call reallocate(has_theta_fac, n_gap_str,copy=.true.)
          call reallocate(has_theta_uniform, n_gap_str,copy=.true.)
          call reallocate(has_theta_file, n_gap_str,copy=.true.)
          
          call reallocate(sparse_file, n_gap_str,copy=.true.)
          call reallocate(mark_sparse_atoms, n_gap_str,copy=.true.)
          call reallocate(sparse_method, n_gap_str,copy=.true.)
          call reallocate(covariance_type, n_gap_str,copy=.true.)
          call reallocate(zeta, n_gap_str,copy=.true.)
          call reallocate(has_zeta, n_gap_str,copy=.true.)
          call reallocate(print_sparse_index, n_gap_str,copy=.true.)

          call reallocate(unique_hash_tolerance, n_gap_str,copy=.true.)
          call reallocate(unique_descriptor_tolerance, n_gap_str,copy=.true.)

          new_gap_str(n_gap_str) = trim(this%gap_str(i_coordinate))
          delta(n_gap_str) = this%delta(i_coordinate)
          f0(n_gap_str) = this%f0(i_coordinate)
          n_sparseX(n_gap_str) = this%n_sparseX(i_coordinate)
          config_type_n_sparseX_string(n_gap_str) = this%config_type_n_sparseX_string(i_coordinate)
          theta_fac_string(n_gap_str) = this%theta_fac_string(i_coordinate)
          theta_uniform(n_gap_str) = this%theta_uniform(i_coordinate)
          theta_file(n_gap_str) = this%theta_file(i_coordinate)

          has_theta_fac(n_gap_str) = this%has_theta_fac(i_coordinate)
          has_theta_uniform(n_gap_str) = this%has_theta_uniform(i_coordinate)
          has_theta_file(n_gap_str) = this%has_theta_file(i_coordinate)
          
          sparse_file(n_gap_str) = this%sparse_file(i_coordinate)
          mark_sparse_atoms(n_gap_str) = this%mark_sparse_atoms(i_coordinate)
          sparse_method(n_gap_str) = this%sparse_method(i_coordinate)
          covariance_type(n_gap_str) = this%covariance_type(i_coordinate)
          zeta(n_gap_str) = this%zeta(i_coordinate)
          has_zeta(n_gap_str) = this%has_zeta(i_coordinate)
          print_sparse_index(n_gap_str) = this%print_sparse_index(i_coordinate)

          unique_hash_tolerance(n_gap_str) = this%unique_hash_tolerance(i_coordinate)
          unique_descriptor_tolerance(n_gap_str) = this%unique_descriptor_tolerance(i_coordinate)

          call print('Unchanged GAP: {'//trim(this%gap_str(i_coordinate))//'}')
       endif

    enddo
    call reallocate(this%delta, n_gap_str)
    call reallocate(this%f0, n_gap_str)
    call reallocate(this%n_sparseX, n_gap_str)
    call reallocate(this%config_type_n_sparseX_string, n_gap_str)
    call reallocate(this%theta_fac_string, n_gap_str)
    call reallocate(this%theta_uniform, n_gap_str)
    call reallocate(this%theta_file, n_gap_str)

    call reallocate(this%has_theta_fac, n_gap_str)
    call reallocate(this%has_theta_uniform, n_gap_str)
    call reallocate(this%has_theta_file, n_gap_str)
    
    call reallocate(this%sparse_file, n_gap_str)
    call reallocate(this%mark_sparse_atoms, n_gap_str)
    call reallocate(this%sparse_method, n_gap_str)
    call reallocate(this%covariance_type, n_gap_str)
    call reallocate(this%zeta, n_gap_str)
    call reallocate(this%has_zeta, n_gap_str)
    call reallocate(this%print_sparse_index, n_gap_str)

    call reallocate(this%unique_hash_tolerance, n_gap_str)
    call reallocate(this%unique_descriptor_tolerance, n_gap_str)

    this%gap_str(1:n_gap_str) = new_gap_str
    this%delta = delta
    this%f0 = f0
    this%n_sparseX = n_sparseX
    this%config_type_n_sparseX_string = config_type_n_sparseX_string
    this%theta_fac_string = theta_fac_string
    this%theta_uniform = theta_uniform
    this%theta_file = theta_file

    this%has_theta_fac = has_theta_fac
    this%has_theta_uniform = has_theta_uniform
    this%has_theta_file = has_theta_file
    
    this%sparse_file = sparse_file
    this%mark_sparse_atoms = mark_sparse_atoms
    this%sparse_method = sparse_method
    this%covariance_type = covariance_type
    this%zeta = zeta
    this%has_zeta = has_zeta
    this%print_sparse_index = print_sparse_index

    this%unique_hash_tolerance = unique_hash_tolerance
    this%unique_descriptor_tolerance = unique_descriptor_tolerance

    this%n_coordinate = n_gap_str

    if(allocated(delta)) deallocate(delta)
    if(allocated(f0)) deallocate(f0)
    if(allocated(n_sparseX)) deallocate(n_sparseX)
    if(allocated(config_type_n_sparseX_string)) deallocate(config_type_n_sparseX_string)
    if(allocated(theta_fac_string)) deallocate(theta_fac_string)
    if(allocated(theta_uniform)) deallocate(theta_uniform)
    if(allocated(theta_file)) deallocate(theta_file)

    if(allocated(has_theta_fac)) deallocate(has_theta_fac)
    if(allocated(has_theta_uniform)) deallocate(has_theta_uniform)
    if(allocated(has_theta_file)) deallocate(has_theta_file)
    
    if(allocated(sparse_file)) deallocate(sparse_file)
    if(allocated(mark_sparse_atoms)) deallocate(mark_sparse_atoms)
    if(allocated(sparse_method)) deallocate(sparse_method)
    if(allocated(covariance_type)) deallocate(covariance_type)
    if(allocated(zeta)) deallocate(zeta)
    if(allocated(has_zeta)) deallocate(has_zeta)
    if(allocated(print_sparse_index)) deallocate(print_sparse_index)

    if(allocated(unique_hash_tolerance)) deallocate(unique_hash_tolerance)
    if(allocated(unique_descriptor_tolerance)) deallocate(unique_descriptor_tolerance)

  endsubroutine add_multispecies_gaps

  subroutine add_template_string(this)
    type(gap_fit), intent(inout) :: this
    character(len=STRING_LENGTH) :: template_string=' '
    character(len=STRING_LENGTH),dimension(:), allocatable :: lines_array    
    type(inoutput) :: tempfile
    integer :: i,n_lines,total_length=0

    if( this%has_template_file ) then
       call print("adding template string, reading from file "//trim(this%template_file))
       call initialise(tempfile,trim(this%template_file))
       call read_file(tempfile,lines_array,n_lines)

       do i=1,n_lines-1
         template_string=trim(template_string)//"{"//trim(lines_array(i))//"};"
         total_length = total_length + len_trim(lines_array(i))
       end do
       template_string=trim(template_string)//"{"//trim(lines_array(n_lines))//"}"
       total_length = total_length + len_trim(lines_array(n_lines))

       if (total_length .ge. STRING_LENGTH) call system_abort("Template atoms object exceeds maximum string size")

       do i=1,len_trim(template_string)
         if(template_string(i:i)==' ') then
           template_string(i:i)='%'
         end if
       end do
       !call print(template_string)

       do i=1,this%n_coordinate
         this%gap_str(i) = trim(this%gap_str(i))//" atoms_template_string={"//trim(template_string)//"}"
       end do
    endif
    
  end subroutine add_template_string

  subroutine gap_fit_read_core_param_file(this)
   type(gap_fit), intent(inout) :: this
   if (this%do_core) then
     call read(this%quip_string, file=trim(this%core_param_file), mpi_comm=this%mpi_obj%communicator, mpi_id=this%mpi_obj%my_proc, keep_lf=.true.)
   end if
 end subroutine gap_fit_read_core_param_file

  subroutine gap_fit_init_mpi_scalapack(this)
    type(gap_fit), intent(inout) :: this

    call initialise(this%mpi_obj)
    call initialise(this%ScaLAPACK_obj, this%mpi_obj, np_r=this%mpi_obj%n_procs, np_c=1)
  end subroutine gap_fit_init_mpi_scalapack

  subroutine gap_fit_init_task_manager(this)
    type(gap_fit), intent(inout) :: this

    this%task_manager%active = this%ScaLAPACK_obj%active
    this%task_manager%MPI_obj = this%MPI_obj
    this%task_manager%ScaLAPACK_obj = this%ScaLAPACK_obj

    call task_manager_init_workers(this%task_manager, this%ScaLAPACK_obj%n_proc_rows)
    call task_manager_init_tasks(this%task_manager, this%n_frame+1) ! mind special task
    this%task_manager%my_worker_id = this%ScaLAPACK_obj%my_proc_row + 1 ! mpi 0-index to tm 1-index

    if (this%task_manager%active) then
      call task_manager_init_idata(this%task_manager, 1)
      this%task_manager%idata(1) = this%mpi_blocksize
    end if
  end subroutine gap_fit_init_task_manager

  subroutine gap_fit_distribute_tasks(this)
    type(gap_fit), intent(inout) :: this

    ! add special task for Cholesky matrix addon to last worker
    call task_manager_add_task(this%task_manager, sum(this%config_type_n_sparseX), n_idata=2, worker_id=SHARED)
    call task_manager_distribute_tasks(this%task_manager)
  end subroutine gap_fit_distribute_tasks

  function gap_fit_is_root(this) result(is_root)
    type(gap_fit), intent(in) :: this
    logical :: is_root
    is_root = (.not. this%MPI_obj%active .or. this%MPI_obj%my_proc == 0)
  end function gap_fit_is_root
  
  subroutine gap_fit_print_linear_system_dump_file(this)
    type(gap_fit), intent(in) :: this
    if (this%has_linear_system_dump_file) then
      call gpFull_print_covariances_lambda(this%my_gp, this%linear_system_dump_file, this%mpi_obj%my_proc)
    end if
  end subroutine gap_fit_print_linear_system_dump_file

  subroutine gap_fit_estimate_memory(this)
    type(gap_fit), intent(in) :: this

    integer(idp), parameter :: rmem = storage_size(1.0_dp, idp) / 8_idp

    integer :: i
    integer(idp) :: s1, s2, entries
    integer(idp) :: mem, memt, memp1  ! scratch, total, peak
    integer(idp) :: sys_total_mem, sys_free_mem

    call print_title("Memory Estimate (per process)")
    
    call print("Descriptors")
    memt = 0
    do i = 1, this%n_coordinate
      s1 = descriptor_dimensions(this%my_descriptor(i))

      entries = s1 * this%n_descriptors(i)
      mem = entries * rmem
      memt = memt + mem
      call print("Descriptor "//i//" :: x "//s1//" "//this%n_descriptors(i)//" memory "//i2si(mem)//"B")

      entries = s1 * this%n_cross(i)
      mem = entries * rmem
      memt = memt + mem
      call print("Descriptor "//i//" :: xPrime "//s1//" "//this%n_cross(i)//" memory "//i2si(mem)//"B")
    end do
    call print("Subtotal "//i2si(memt)//"B")
    call print("")
    memp1 = memt


    call print("Covariances")
    memt = 0
    s1 = sum(this%config_type_n_sparseX)
    s2 = (this%n_ener + this%n_local_property) + (this%n_force + this%n_virial + this%n_hessian)

    entries = s1 * s2
    mem = entries * rmem
    memt = memt + mem * 2
    call print("yY "//s1//" "//s2//" memory "//i2si(mem)//"B * 2")
    memp1 = memp1 + mem

    entries = s1 * s1
    mem = entries * rmem
    memt = memt + mem
    call print("yy "//s1//" "//s1//" memory "//i2si(mem)//"B")

    entries = s1 * (s1 + s2)
    mem = entries * rmem
    memt = memt + mem * 2
    call print("A "//s1//" "//(s1+s2)//" memory "//i2si(mem)//"B * 2")
    call print("Subtotal "//i2si(memt)//"B")
    call print("")

    
    mem = max(memp1, memt)
    call print("Peak1 "//i2si(memp1)//"B")
    call print("Peak2 "//i2si(memt)//"B")
    call print("PEAK  "//i2si(mem)//"B")
    call print("")

    call mem_info(sys_total_mem, sys_free_mem)
    call print("Free system memory  "//i2si(sys_free_mem)//"B")
    call print("Total system memory "//i2si(sys_total_mem)//"B")

    mem = sys_free_mem - mem
    if (mem < 0) then
      call print_warning("Memory estimate exceeds free system memory by "//i2si(-mem)//"B.")
    end if

    call print_title("")
  end subroutine gap_fit_estimate_memory

end module gap_fit_module
