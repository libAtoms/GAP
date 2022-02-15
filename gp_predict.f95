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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Gaussian Process module
!X
!% Module for general GP function interpolations.
!% A gp object contains the training set (fitting points and function values),
!% important temporary matrices, vectors and parameters.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module gp_predict_module

   use iso_c_binding, only : C_NULL_CHAR
   ! use libatoms_module  
   use error_module  
#ifdef _OPENMP  
   use omp_lib  
#endif  
   use system_module, only : idp, dp, qp, optional_default, reallocate, NUMERICAL_ZERO, & 
       system_timer, string_to_numerical, print_warning, progress, progress_timer, &
       current_times, InOutput, OUTPUT, increase_to_multiple, i2si, PRINT_VERBOSE
   use units_module  
   use linearalgebra_module  
   use extendable_str_module  
   use dictionary_module  
   use paramreader_module  
   use descriptors_module
   use fox_wxml
   use FoX_sax, only: xml_t, dictionary_t, haskey, getvalue, parse, &
   open_xml_string, close_xml_t
   use CInOutput_module, only : quip_md5sum
   use task_manager_module
   use matrix_module

   implicit none

   private

   integer, parameter :: besseli_max_n = 20

   real(dp), dimension(besseli_max_n), parameter :: besseli0_c = (/ &
   0.125_dp, &
   7.03125E-002_dp, &
   7.32421875E-002_dp, &
   0.112152099609375_dp, &    
   0.22710800170898438_dp, &    
   0.57250142097473145_dp, &    
   1.7277275025844574_dp, &    
   6.0740420012734830_dp, &    
   24.380529699556064_dp, &    
   110.01714026924674_dp, &    
   551.33589612202059_dp, &    
   3038.0905109223841_dp, &    
   18257.755474293175_dp, &    
   118838.42625678326_dp, &    
   832859.30401628942_dp, &    
   6252951.4934347980_dp, &    
   50069589.531988934_dp, &    
   425939216.50476694_dp, &    
   3836255180.2304339_dp, &    
   36468400807.065559_dp /)    

   real(dp), dimension(besseli_max_n), parameter :: besseli1_c = (/ &
   -0.375_dp, &    
   -0.1171875_dp, &    
   -0.1025390625_dp, &    
   -0.144195556640625_dp, &    
   -0.27757644653320313_dp, &    
   -0.67659258842468262_dp, &    
   -1.9935317337512970_dp, &    
   -6.8839142681099474_dp, &    
   -27.248827311268542_dp, &    
   -121.59789187653587_dp, &    
   -603.84407670507017_dp, &    
   -3302.2722944808525_dp, &    
   -19718.375912236628_dp, &    
   -127641.27264617461_dp, &    
   -890297.87670706783_dp, &    
   -6656367.7188176867_dp, &    
   -53104110.109685220_dp, &    
   -450278600.30503929_dp, &    
   -4043620325.1077542_dp, &    
   -38338575207.427895_dp /)

   real(dp), parameter :: besseli_max_x = 18.0_dp

   real(dp), parameter :: THETA_MIN = 1.0e-8_dp
   integer, parameter, public :: GP_SPARSE_RANDOM = 1
   integer, parameter, public :: GP_SPARSE_PIVOT = 2
   integer, parameter, public :: GP_SPARSE_CLUSTER = 3
   integer, parameter, public :: GP_SPARSE_UNIFORM = 4
   integer, parameter, public :: GP_SPARSE_KMEANS = 5
   integer, parameter, public :: GP_SPARSE_COVARIANCE = 6
   integer, parameter, public :: GP_SPARSE_UNIQ = 7
   integer, parameter, public :: GP_SPARSE_FUZZY = 8
   integer, parameter, public :: GP_SPARSE_FILE = 9
   integer, parameter, public :: GP_SPARSE_INDEX_FILE = 10
   integer, parameter, public :: GP_SPARSE_CUR_COVARIANCE = 11
   integer, parameter, public :: GP_SPARSE_CUR_POINTS = 12
   integer, parameter, public :: GP_SPARSE_NONE = 13

   integer, parameter, public :: GP_COVARIANCE_FITC = 1
   integer, parameter, public :: GP_COVARIANCE_DTC = 2


   integer, parameter, public :: COVARIANCE_NONE             = 0
   integer, parameter, public :: COVARIANCE_ARD_SE           = 1
   integer, parameter, public :: COVARIANCE_DOT_PRODUCT      = 2
   integer, parameter, public :: COVARIANCE_BOND_REAL_SPACE  = 3
   integer, parameter, public :: COVARIANCE_PP               = 4

   integer, parameter, public :: PP_Q = 1

   integer, public :: openmp_chunk_size = 1

   type gpCovariance_bond_real_space

      integer :: n
      real(dp) :: delta
      real(dp) :: atom_sigma

      logical :: initialised = .false.

   endtype gpCovariance_bond_real_space

   type gpCovariance_atom_real_space

      integer :: l_max = 0
      real(dp) :: atom_sigma, delta, zeta
      real(dp) :: cutoff, cutoff_transition_width

      logical :: initialised = .false.

   endtype gpCovariance_atom_real_space

   public :: gpCovariance_bond_real_space  
   public :: gpCovariance_bond_real_space_Calc  
   public :: gpCoordinates_gpCovariance_bond_real_space_Initialise  
  
   public :: gpCovariance_atom_real_space
   public :: gpCovariance_atom_real_space_Calc  

   type gpCoordinates

      integer :: d = 0, n_x, n_xPrime, n_sparseX, n_permutations
      ! dimension of descriptors, number of descriptors, number of derivatives of descriptors

      integer :: current_x, current_xPrime
      ! pointers to the last added values

      real(dp), dimension(:,:), allocatable :: x, xPrime
      ! descriptors (d,n_x), derivatives of descriptors (d, n_xPrime)
      ! for real space covariance descriptors (max(x_size),n_x), derivatives of descriptors (max(x_size),n_xPrime)
      real(dp), dimension(:), allocatable :: cutoff, cutoffPrime
      integer, dimension(:), allocatable :: x_size, xPrime_size
      real(dp), dimension(:), allocatable :: covarianceDiag_x_x, covarianceDiag_xPrime_xPrime

      real(dp), dimension(:,:), allocatable :: sparseX, covarianceDiag_x_xPrime
      real(dp), dimension(:), allocatable :: sparseCutoff
      ! sparse points stored as real array
      ! for real space covariance descriptors
      integer, dimension(:), allocatable :: sparseX_size
      real(dp), dimension(:), allocatable :: covarianceDiag_sparseX_sparseX

      real(dp), dimension(:,:,:), allocatable :: sparseX_permuted  
      real(dp), dimension(:), allocatable :: sparseCovariance  
  
      real(dp), dimension(:), allocatable :: theta
      ! range parameters (d) for descriptors in each directions
      real(dp) :: zeta = 0.0_dp

      real(dp), dimension(:), allocatable :: alpha
      ! 

      real(dp) ::  delta, f0 = 0.0_dp, variance_estimate_regularisation = 0.0_dp  
      ! range of GP (function value) and baseline of function

      integer, dimension(:), allocatable :: map_x_y, map_xPrime_yPrime, map_xPrime_x, config_type
      ! which descriptor is used for a given function value, which derivative descriptor is used for a given derivative function, which descriptor is differentiated

      integer, dimension(:), allocatable :: map_sparseX_globalSparseX
      ! sparse point in this descriptor type -> all sparse points in gpFull

      integer, dimension(:), allocatable :: sparseX_index
      ! sparse points stored as indices of the x array

      integer, dimension(:,:), allocatable :: permutations
      ! Lists the permutations symmetries of the coordinates
      logical, dimension(:,:), allocatable ::  permutation_distance_mask
      ! pairwise distances that may occur given all permutations

      type(gpCovariance_bond_real_space) :: bond_real_space_cov
      integer :: covariance_type = COVARIANCE_NONE

      type(extendable_str) :: descriptor_str

      type(LA_Matrix) :: LA_k_mm

      logical :: initialised = .false.
      logical :: sparsified = .false.
      logical :: variance_estimate_initialised = .false.  
      logical :: sparse_covariance_initialised = .false.  

   endtype gpCoordinates

   public :: gpCoordinates 

   type gpFull

      integer :: n_y, n_yPrime
      ! number of function values, number of derivative function values
      
      integer :: n_globalSparseX
      ! number of all sparse points in every descriptor type

      integer :: n_coordinate
      ! number of different descriptors

      integer :: current_y, current_yPrime

      real(dp) :: sparse_jitter = 1.0e-5_dp

      real(dp), dimension(:), allocatable :: y, yPrime
      ! function values, derivative function values

      real(dp), dimension(:), allocatable :: sigma_y, sigma_yPrime
      ! estimated error of function values, derivatives

      real(dp), dimension(:,:), allocatable :: covariance_subY_y, covariance_subY_subY, covariance_y_y, inverse_sparse_full
      ! covariance matrix

      real(dp), dimension(:), allocatable :: covarianceDiag_y_y, lambda, alpha
      ! covariance matrix

      integer, dimension(:), allocatable :: map_y_globalY, map_yPrime_globalY

      type(gpCoordinates), dimension(:), allocatable :: coordinate

      integer :: covariance_method = GP_COVARIANCE_DTC
      logical :: initialised = .false.

   endtype gpFull

   type gpSparse
      integer :: n_coordinate ! number of different descriptors
      type(gpCoordinates), dimension(:), allocatable :: coordinate
      logical :: initialised = .false.
      logical :: fitted = .false.
   endtype gpSparse

   type cplx_1d_array
      complex(dp), dimension(:), allocatable :: value
   endtype cplx_1d_array

   type cplx_2d_array
      complex(dp), dimension(:,:), allocatable :: value
   endtype cplx_2d_array

   type neighbour_descriptor
      type(cplx_1d_array), dimension(:), allocatable :: spherical_harmonics
      real(dp) :: r
      integer :: n
   endtype neighbour_descriptor

   logical, save :: parse_matched_label, parse_in_gpCoordinates, parse_in_gpFull, parse_in_gpSparse, parse_in_sparseX, parse_sliced, parse_sparseX_separate_file
   integer, save :: parse_i_sparseX, parse_i_x, parse_i_xPrime, parse_i_permutation, parse_slice_start, parse_slice_end
   type(gpCoordinates), pointer :: parse_gpCoordinates
   type(gpFull), pointer :: parse_gpFull
   type(gpSparse), pointer :: parse_gpSparse
   type(extendable_str), save :: parse_cur_data
   integer, dimension(:,:), allocatable :: parse_in_permutations
   character(len=1024), save :: parse_gpCoordinates_label, parse_gpFull_label, parse_gpSparse_label

   public :: gpFull, gpSparse
   public :: gpFull_print_covariances_lambda
   public :: gpSparse_fit

   interface initialise
      module procedure gpSparse_initialise
   endinterface initialise
   public :: initialise

   interface finalise
      module procedure gpFull_Finalise, gpCoordinates_Finalise, gpSparse_finalise, gpNeighbourDescriptor_Finalise
   endinterface finalise
   public :: finalise

   interface gp_setTheta
      module procedure gpCoordinates_setTheta, gpFull_setTheta
   endinterface gp_setTheta
   public :: gp_setTheta

   interface gp_setThetaFactor
      module procedure gpFull_setTheta_thetaFactor !, gpFull_setTheta_thetaFactorArray, gpFull_setTheta_thetaFactorUniform
   endinterface gp_setThetaFactor
   public :: gp_setThetaFactor

   interface gp_setParameters
      module procedure gpFull_setParameters, gpFull_gpCoordinates_setParameters, gpCoordinates_setParameters, &
      gpCoordinates_setParameters_sparse, gpSparse_setParameters
   endinterface gp_setParameters
   public :: gp_setParameters

   interface gp_setPermutations
      module procedure gpCoordinates_setPermutations, gpFull_setPermutations, gpSparse_setPermutations
   endinterface gp_setPermutations
   public :: gp_setPermutations

   interface gp_addFunctionValue
      module procedure gpFull_addFunctionValue
   endinterface gp_addFunctionValue
   public :: gp_addFunctionValue

   interface gp_addFunctionDerivative
      module procedure gpFull_addFunctionDerivative
   endinterface gp_addFunctionDerivative
   public :: gp_addFunctionDerivative

   interface gp_addCoordinates
      module procedure gpFull_addCoordinates_1Darray, gpFull_addCoordinates_2Darray
   endinterface gp_addCoordinates
   public :: gp_addCoordinates

   interface gp_addCoordinateDerivatives
      module procedure gpFull_addCoordinateDerivatives_1Darray, gpFull_addCoordinateDerivatives_2Darray
   endinterface gp_addCoordinateDerivatives
   public :: gp_addCoordinateDerivatives

   interface gp_addDescriptor
      module procedure gpFull_addDescriptor
   endinterface gp_addDescriptor
   public :: gp_addDescriptor

   interface gp_printXML
      module procedure gpCoordinates_printXML, gpFull_printXML, gpSparse_printXML
   endinterface gp_PrintXML
   public :: gp_printXML

   interface gp_readXML
      module procedure gpCoordinates_readXML, gpFull_readXML, gpSparse_readXML, &
      gpCoordinates_readXML_string, gpFull_readXML_string, gpSparse_readXML_string
   endinterface gp_readXML
   public :: gp_readXML

   interface gp_covariance_sparse
      module procedure gpFull_covarianceMatrix_sparse
   endinterface gp_covariance_sparse
   public :: gp_covariance_sparse

   interface gp_covariance_full
      module procedure gpFull_covarianceMatrix
   endinterface gp_covariance_full
   public :: gp_covariance_full

   interface gp_Predict
      module procedure gpCoordinates_Predict
   endinterface gp_Predict
   public :: gp_Predict

   interface gp_log_likelihood
      module procedure gpCoordinates_log_likelihood
   endinterface gp_log_likelihood
   public :: gp_log_likelihood

   public :: gpCoordinates_Covariance  
   public :: gpCoordinates_initialise_variance_estimate  
   public :: covariancePP

   contains

   subroutine gpFull_setParameters(this, n_coordinate, n_y, n_yPrime, sparse_jitter, error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: n_coordinate, n_y, n_yPrime
      real(dp), intent(in) :: sparse_jitter
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      this%n_coordinate = n_coordinate
      this%n_y = n_y
      this%n_yPrime = n_yPrime
      this%current_y = 0
      this%current_yPrime = 0
      this%sparse_jitter = sparse_jitter

      allocate( this%coordinate(n_coordinate) )
      allocate( this%y(n_y), this%yPrime(n_yPrime) )
      allocate( this%map_y_globalY(n_y), this%map_yPrime_globalY(n_yPrime) )
      allocate( this%sigma_y(n_y), this%sigma_yPrime(n_yPrime) )

      this%initialised = .true.

   endsubroutine gpFull_setParameters

   subroutine gpFull_gpCoordinates_setParameters(this, i, d, n_x, n_xPrime, delta, f0, covariance_type, x_size_max, xPrime_size_max, error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i, d, n_x, n_xPrime
      real(dp), intent(in) :: delta, f0
      integer, optional, intent(in) :: covariance_type
      integer, optional, intent(in) :: x_size_max, xPrime_size_max
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_set_gpCoordinates_parameters: object not initialised',error)
      endif

      if( i > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_set_gpCoordinates_parameters: access to descriptor '//i//' is not possible as number of descriptors is '//this%n_coordinate,error )
      endif

      call gpCoordinates_setParameters(this%coordinate(i), d, n_x, n_xPrime, delta, f0, covariance_type = covariance_type, x_size_max=x_size_max, xPrime_size_max=xPrime_size_max, error=error)

   endsubroutine gpFull_gpCoordinates_setParameters

   subroutine gpCoordinates_setParameters(this, d, n_x, n_xPrime, delta, f0, covariance_type, x_size_max, xPrime_size_max, error)

      type(gpCoordinates), intent(inout) :: this
      integer, intent(in) :: d, n_x, n_xPrime
      real(dp), intent(in) :: delta, f0
      integer, optional, intent(in) :: covariance_type
      integer, optional, intent(in) :: x_size_max, xPrime_size_max
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      if( d < 0 ) then
         RAISE_ERROR("gpCoordinates_setParameters: negative value of d = "//d,error)
      else
         this%d = d
      endif

      if( n_x < 0 ) then
         RAISE_ERROR("gpCoordinates_setParameters: negative value of n_x = "//n_x,error)
      else
         this%n_x = n_x
      endif

      if( n_xPrime < 0 ) then
         RAISE_ERROR("gpCoordinates_setParameters: negative value of n_xPrime = "//n_xPrime,error)
      else
         this%n_xPrime = n_xPrime
      endif

      this%delta = delta
      this%f0 = f0

      this%current_x = 0
      this%current_xPrime = 0
      this%n_sparseX = 0
      this%n_permutations = 1

      this%covariance_type = optional_default(COVARIANCE_ARD_SE, covariance_type)

      if(present(x_size_max)) then
         allocate( this%x(x_size_max,n_x) )
      else
         allocate( this%x(d,n_x) )
      endif
      this%x = 0.0_dp

      if(present(xPrime_size_max)) then
         allocate( this%xPrime(xPrime_size_max,n_xPrime) )
      else
         allocate( this%xPrime(d,n_xPrime) )
      endif
      this%xPrime = 0.0_dp

      allocate(this%cutoff(n_x))
      this%cutoff = 1.0_dp
      allocate(this%cutoffPrime(n_xPrime))
      this%cutoffPrime = 0.0_dp

      allocate( this%config_type(n_x) )
      this%config_type = 0

      allocate( this%map_x_y(n_x), this%map_xPrime_yPrime(n_xPrime), this%map_xPrime_x(n_xPrime) )
      this%map_x_y = 0
      this%map_xPrime_yPrime = 0
      this%map_xPrime_x = 0

      allocate(this%covarianceDiag_x_x(n_x), this%covarianceDiag_xPrime_xPrime(n_xPrime))
      this%covarianceDiag_x_x = 1.0_dp
      this%covarianceDiag_xPrime_xPrime = 1.0_dp

      select case(this%covariance_type)
      case(COVARIANCE_BOND_REAL_SPACE) 
         allocate( this%x_size(n_x), this%xPrime_size(n_xPrime) )
         this%x_size = d
         this%xPrime_size = 0
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      case(COVARIANCE_DOT_PRODUCT)
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      case(COVARIANCE_ARD_SE,COVARIANCE_PP)
         allocate( this%theta(d), this%permutations(d,1) )
         this%theta = 0.0_dp
         this%permutations(:,1) = (/ (i, i=1, d) /)

         allocate(this%permutation_distance_mask(this%d,this%d))
         this%permutation_distance_mask = .false.
         forall(i=1:this%d) this%permutation_distance_mask(i,i) = .true.
      endselect

      this%sparsified = .false.
      this%initialised = .true.
   endsubroutine gpCoordinates_setParameters

   subroutine gpCoordinates_setParameters_sparse(this, d, n_sparseX, delta, f0, covariance_type, sparseX_size_max, error)

      type(gpCoordinates), intent(inout) :: this
      integer, intent(in) :: d, n_sparseX
      real(dp), intent(in) :: delta, f0
      integer, optional, intent(in) :: covariance_type
      integer, optional, intent(in) :: sparseX_size_max
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      this%d = d
      this%n_x = 0
      this%n_xPrime = 0
      this%delta = delta
      this%f0 = f0

      this%current_x = 0
      this%current_xPrime = 0
      this%n_sparseX = n_sparseX
      this%n_permutations = 1

      this%covariance_type = optional_default(COVARIANCE_ARD_SE, covariance_type)

      if(present(sparseX_size_max)) then
         allocate( this%sparseX(sparseX_size_max,n_sparseX) )
      else
         allocate( this%sparseX(d,n_sparseX) )
      endif

      allocate( this%alpha(n_sparseX) )
      allocate( this%sparseCutoff(n_sparseX) )

      allocate( this%covarianceDiag_sparseX_sparseX(n_sparseX) )
      this%covarianceDiag_sparseX_sparseX = 1.0_dp

      select case(this%covariance_type)
      case(COVARIANCE_BOND_REAL_SPACE)
         allocate( this%sparseX_size(n_sparseX) )
         this%sparseX_size = d    
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      case(COVARIANCE_DOT_PRODUCT)
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      case(COVARIANCE_ARD_SE,COVARIANCE_PP)
         allocate( this%theta(d), this%permutations(d,1) )
         this%theta = 0.0_dp
         this%permutations(:,1) = (/ (i, i=1, d) /)
         allocate(this%permutation_distance_mask(this%d,this%d))
         this%permutation_distance_mask = .false.
         forall(i=1:this%d) this%permutation_distance_mask(i,i) = .true.
      endselect

      this%sparsified = .true.
      this%initialised = .true.

   endsubroutine gpCoordinates_setParameters_sparse

   subroutine gpSparse_setParameters(this,n_coordinate,error)
      type(gpSparse), intent(inout) :: this
      integer, intent(in) :: n_coordinate
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)
      this%n_coordinate = n_coordinate
      allocate( this%coordinate(this%n_coordinate) )

   endsubroutine gpSparse_setParameters

   subroutine gpCoordinates_setPermutations(this,permutations,error)
      type(gpCoordinates), intent(inout) :: this
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      real(dp), dimension(this%d) :: theta
      integer :: i, d

      INIT_ERROR(error)

      this%n_permutations = size(permutations,2)

      select case(this%covariance_type)
      case(COVARIANCE_ARD_SE,COVARIANCE_PP)
         call reallocate(this%permutations,this%d,this%n_permutations,zero=.true.)
         this%permutations = permutations
         ! Symmetrise theta wrt permutations
         theta = this%theta
         this%theta = 0.0_dp
         do i = 1, this%n_permutations
            this%theta = this%theta + theta(this%permutations(:,i))
         enddo
         this%theta = this%theta / real(this%n_permutations,kind=dp)

         this%permutation_distance_mask = .false.
         do i = 1, this%n_permutations
            do d = 1, this%d
               this%permutation_distance_mask(d,this%permutations(d,i)) = .true.
            enddo
         enddo
      case default

      endselect


   endsubroutine gpCoordinates_setPermutations

   subroutine gpFull_setPermutations(this,i_coordinate,permutations,error)
      type(gpFull), intent(inout) :: this
      integer :: i_coordinate
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_setPermutations: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      call gpCoordinates_setPermutations(this%coordinate(i_coordinate),permutations,error)

   endsubroutine gpFull_setPermutations

   subroutine gpSparse_setPermutations(this,i_coordinate,permutations,error)
      type(gpSparse), intent(inout) :: this
      integer :: i_coordinate
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpSparse_setPermutations: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      call gpCoordinates_setPermutations(this%coordinate(i_coordinate),permutations,error)

   endsubroutine gpSparse_setPermutations

   subroutine gpSparse_initialise(this, from, error)
      type(gpSparse), intent(inout) :: this
      type(gpFull), intent(in) :: from
      integer, optional, intent(out) :: error

      integer :: i

      if( .not. from%initialised ) then
         RAISE_ERROR('gpSparse_initialise: gpFull object not initialised',error)
      endif

      if(this%initialised) call finalise(this,error)

      call gpSparse_setParameters(this, from%n_coordinate)

      do i = 1, this%n_coordinate
         if( from%coordinate(i)%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call gpCoordinates_setParameters_sparse(this%coordinate(i), &
                 from%coordinate(i)%d, from%coordinate(i)%n_sparseX, from%coordinate(i)%delta, from%coordinate(i)%f0, covariance_type = from%coordinate(i)%covariance_type, &
                 sparseX_size_max=maxval(from%coordinate(i)%sparseX_size), error=error)
         else
            call gpCoordinates_setParameters_sparse(this%coordinate(i), &
                 from%coordinate(i)%d, from%coordinate(i)%n_sparseX, from%coordinate(i)%delta, from%coordinate(i)%f0, covariance_type = from%coordinate(i)%covariance_type, &
                 error=error)
         endif

         this%coordinate(i)%alpha = 0.0
         this%coordinate(i)%sparseX = from%coordinate(i)%sparseX
         this%coordinate(i)%covarianceDiag_sparseX_sparseX = from%coordinate(i)%covarianceDiag_sparseX_sparseX

         if(from%coordinate(i)%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            this%coordinate(i)%sparseX_size = from%coordinate(i)%sparseX_size
         endif

         this%coordinate(i)%theta = from%coordinate(i)%theta
         this%coordinate(i)%zeta = from%coordinate(i)%zeta
         this%coordinate(i)%descriptor_str = from%coordinate(i)%descriptor_str
         this%coordinate(i)%sparseCutoff = from%coordinate(i)%sparseCutoff

         call gpSparse_setPermutations(this,i,from%coordinate(i)%permutations,error)
      enddo

      this%initialised = .true.
   end subroutine gpSparse_initialise

   subroutine gpSparse_fit(this, from, task_manager, condition_number_norm, error)
      type(gpSparse), intent(inout) :: this
      type(gpFull), intent(inout) :: from  ! actually input; intent(inout) to free memory early
      type(task_manager_type), intent(in) :: task_manager
      character(len=*), optional, intent(in) :: condition_number_norm
      integer, optional, intent(out) :: error

      character(len=STRING_LENGTH) :: my_condition_number_norm

      integer :: i, j, n, o, t, w, blocksize
      integer :: i_coordinate, i_sparseX, i_global_sparseX, n_globalSparseX, n_globalY, i_y, i_yPrime, &
      i_globalY, i_global_yPrime, nlrows
#ifdef HAVE_QR      
      real(qp) :: rcond
      real(qp), dimension(:,:), allocatable :: c_subYY_sqrtInverseLambda, factor_c_subYsubY, a
      real(qp), dimension(:), allocatable :: globalY, alpha
      type(LA_Matrix) :: LA_c_subYsubY, LA_q_subYsubY
#else
      real(qp), dimension(:,:), allocatable :: c_subYY_inverseLambda, c_subYY_inverseLambda_c_YsubY!, &
!      inverse_q_subYsubY, inverse_c_subYsubY
      real(qp), dimension(:), allocatable :: globalY, alpha
      type(LA_Matrix) :: LA_q_subYsubY
#endif

      INIT_ERROR(error)

      my_condition_number_norm = optional_default(' ', condition_number_norm)
      
      call gpSparse_initialise(this, from, error)

      n_globalSparseX = from%n_globalSparseX
      n_globalY = from%n_y + from%n_yPrime

#ifdef HAVE_QR
      allocate(c_subYY_sqrtInverseLambda(n_globalSparseX,n_globalY))
      call matrix_product_vect_asdiagonal_sub(c_subYY_sqrtInverseLambda,from%covariance_subY_y,sqrt(1.0_qp/from%lambda)) ! O(NM)
      if (allocated(from%covariance_subY_y)) deallocate(from%covariance_subY_y)  ! free input component to save memory
      
      allocate(factor_c_subYsubY(n_globalSparseX,n_globalSparseX))
      call initialise(LA_c_subYsubY,from%covariance_subY_subY)
      call LA_Matrix_Factorise(LA_c_subYsubY,factor_c_subYsubY,error=error)
      call finalise(LA_c_subYsubY)
      if (allocated(from%covariance_subY_subY)) deallocate(from%covariance_subY_subY)  ! free input component to save memory

      do i = 1, n_globalSparseX-1
         do j = i+1, n_globalSparseX
            factor_c_subYsubY(j,i) = 0.0_qp
         end do
      end do

      allocate(alpha(n_globalSparseX))
      if (task_manager%active) then
         blocksize = get_blocksize(task_manager%idata(1), task_manager%unified_workload, n_globalSparseX)
         nlrows = increase_to_multiple(task_manager%unified_workload, blocksize)
         o = nlrows - task_manager%unified_workload
         call print("distA extension: "//o//" "//n_globalSparseX//" memory "//i2si(8_idp * o * n_globalSparseX)//"B", PRINT_VERBOSE)
         allocate(globalY(nlrows))
         allocate(a(nlrows,n_globalSparseX))
         alpha = 0.0_qp
         globalY = 0.0_qp
         a = 0.0_qp
      else
         allocate(globalY(n_globalY+n_globalSparseX))
         allocate(a(n_globalY+n_globalSparseX,n_globalSparseX))
      end if

      a(1:n_globalY,:) = transpose(c_subYY_sqrtInverseLambda)
      if (allocated(c_subYY_sqrtInverseLambda)) deallocate(c_subYY_sqrtInverseLambda)

      if (task_manager%active) then
         ! put L part at the end of local A (take info from last task of each worker)
         w = task_manager%my_worker_id
         t = task_manager%workers(w)%n_tasks
         n = task_manager%workers(w)%tasks(t)%idata(1)
         o = task_manager%workers(w)%tasks(t)%idata(2)
         n = min(o+n, n_globalSparseX) - o ! secure against sparseX reduction since distribution
         a(n_globalY+1:n_globalY+n,:) = factor_c_subYsubY(o+1:o+n,:)
      else
         a(n_globalY+1:,:) = factor_c_subYsubY
      end if
      if (allocated(factor_c_subYsubY)) deallocate(factor_c_subYsubY)

      if (my_condition_number_norm(1:1) /= ' ') then
         if (task_manager%active) then
            call print_warning("Condition number of distributed matrix is not implemented.")
         else
            rcond = matrix_condition_number(a, my_condition_number_norm(1:1))
            call print("Condition number (log10) of matrix A (norm "//my_condition_number_norm(1:1)//"): "//-log10(rcond))
         end if
      end if

      globalY = 0.0_qp
      do i_y = 1, from%n_y
         ! loop over all function values

         i_globalY = from%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         globalY(i_globalY) = from%y(i_y)*sqrt(1.0_qp/from%lambda(i_globalY))
      enddo

      do i_yPrime = 1, from%n_yPrime
         ! loop over all function values

         i_global_yPrime = from%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         globalY(i_global_yPrime) = from%yPrime(i_yPrime)*sqrt(1.0_qp/from%lambda(i_global_yPrime))
      enddo

      if (task_manager%active) then
         call print("Using ScaLAPACK to solve QR")
         call SP_Matrix_QR_Solve(a, globalY, alpha, task_manager%ScaLAPACK_obj, blocksize)
      else
         call print("Using LAPACK to solve QR")
         call initialise(LA_q_subYsubY, a, use_allocate=.false.)
         call LA_Matrix_QR_Solve_Vector(LA_q_subYsubY, globalY, alpha)
         call finalise(LA_q_subYsubY)
      end if

      do i_coordinate = 1, from%n_coordinate
         do i_sparseX = 1, from%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = from%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            this%coordinate(i_coordinate)%alpha(i_sparseX) = real(alpha(i_global_sparseX),kind=dp)
         enddo
      enddo

      if(allocated(a)) deallocate(a)
      if(allocated(globalY)) deallocate(globalY)
      if(allocated(alpha)) deallocate(alpha)
#else
      allocate( c_subYY_inverseLambda(n_globalSparseX,n_globalY), c_subYY_inverseLambda_c_YsubY(n_globalSparseX,n_globalSparseX), &
!      inverse_q_subYsubY(n_globalSparseX,n_globalSparseX), inverse_c_subYsubY(n_globalSparseX,n_globalSparseX), &
      alpha(n_globalSparseX), globalY(n_globalY))

      call matrix_product_vect_asdiagonal_sub(c_subYY_inverseLambda,from%covariance_subY_Y,1.0_qp/from%lambda) ! O(NM)

      c_subYY_inverseLambda_c_YsubY = matmul(c_subYY_inverseLambda,transpose(from%covariance_subY_Y))
      call initialise(LA_q_subYsubY,from%covariance_subY_subY + c_subYY_inverseLambda_c_YsubY)

      globalY = 0.0_qp
      do i_y = 1, from%n_y
         ! loop over all function values

         i_globalY = from%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         globalY(i_globalY) = from%y(i_y) !*sqrt(1.0_qp/from%lambda(i_globalY))
      enddo

      do i_yPrime = 1, from%n_yPrime
         ! loop over all function values

         i_global_yPrime = from%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         globalY(i_global_yPrime) = from%yPrime(i_yPrime) !*sqrt(1.0_qp/from%lambda(i_global_yPrime))
      enddo

      call Matrix_Solve(LA_q_subYsubY,matmul(c_subYY_inverseLambda, globalY),alpha)
      call finalise(LA_q_subYsubY)

      do i_coordinate = 1, from%n_coordinate
         do i_sparseX = 1, from%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = from%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            this%coordinate(i_coordinate)%alpha(i_sparseX) = real(alpha(i_global_sparseX),kind=dp)
         enddo
      enddo

      if(allocated(c_subYY_inverseLambda)) deallocate(c_subYY_inverseLambda)
      if(allocated(c_subYY_inverseLambda_c_YsubY)) deallocate(c_subYY_inverseLambda_c_YsubY)
!      if(allocated(inverse_q_subYsubY)) deallocate(inverse_q_subYsubY)
!      if(allocated(inverse_c_subYsubY)) deallocate(inverse_c_subYsubY)
      if(allocated(alpha)) deallocate(alpha)
      if(allocated(globalY)) deallocate(globalY)
#endif
      this%fitted = .true.

   endsubroutine gpSparse_fit

   function get_blocksize(arg, nrows, ncols) result(blocksize)
      integer, intent(in) :: arg, nrows, ncols
      integer :: blocksize

      integer, parameter :: OVERHEAD_FACTOR = 2  ! worst case: +1/2=50% matrix size

      if (arg > 0) then
         blocksize = arg
      else
         blocksize = merge(nrows, ncols, (nrows < ncols * OVERHEAD_FACTOR))
      end if
   end function get_blocksize

   subroutine gpSparse_finalise(this,error)
      type(gpSparse), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i_coordinate

      INIT_ERROR(error)

      if (allocated(this%coordinate)) then
         do i_coordinate = 1, this%n_coordinate
            call finalise(this%coordinate(i_coordinate), error)
         enddo
         deallocate(this%coordinate)
      end if

      this%n_coordinate = 0
      this%initialised = .false.
      this%fitted = .false.

   endsubroutine gpSparse_finalise

   subroutine gpFull_Finalise(this, error)
      type(gpFull), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) return

      if(allocated(this%coordinate)) then
         do i = 1, this%n_coordinate
            call finalise(this%coordinate(i))
         enddo
         deallocate( this%coordinate )
      endif

      if(allocated(this%y)) deallocate( this%y )
      if(allocated(this%yPrime)) deallocate( this%yPrime )
      if(allocated(this%sigma_y)) deallocate( this%sigma_y )
      if(allocated(this%sigma_yPrime)) deallocate( this%sigma_yPrime )
      if(allocated(this%map_y_globalY)) deallocate( this%map_y_globalY )
      if(allocated(this%map_yPrime_globalY)) deallocate( this%map_yPrime_globalY )
      if(allocated(this%covariance_subY_y)) deallocate( this%covariance_subY_y )
      if(allocated(this%covariance_subY_subY)) deallocate( this%covariance_subY_subY )
      if(allocated(this%covarianceDiag_y_y)) deallocate( this%covarianceDiag_y_y )
      if(allocated(this%lambda)) deallocate( this%lambda )
      if(allocated(this%alpha)) deallocate( this%alpha )
      if(allocated(this%inverse_sparse_full)) deallocate( this%inverse_sparse_full )


      this%n_coordinate = 0
      this%n_y = 0
      this%n_yPrime = 0
      this%current_y = 0
      this%current_yPrime = 0

      this%initialised = .false.

   endsubroutine gpFull_Finalise

   subroutine gpCoordinates_Finalise(this, error)
      type(gpCoordinates), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      if(allocated(this%x)) deallocate( this%x )
      if(allocated(this%xPrime)) deallocate( this%xPrime )
      if(allocated(this%cutoff)) deallocate( this%cutoff )
      if(allocated(this%cutoffPrime)) deallocate( this%cutoffPrime )
      if(allocated(this%theta)) deallocate( this%theta )
      
      if(allocated(this%permutations)) deallocate(this%permutations)
      if(allocated(this%permutation_distance_mask)) deallocate( this%permutation_distance_mask )

      if(allocated(this%map_x_y)) deallocate( this%map_x_y )
      if(allocated(this%map_xPrime_yPrime)) deallocate( this%map_xPrime_yPrime )
      if(allocated(this%map_xPrime_x)) deallocate( this%map_xPrime_x )
      if(allocated(this%map_sparseX_globalSparseX)) deallocate( this%map_sparseX_globalSparseX )
      if(allocated(this%config_type)) deallocate( this%config_type )

      if(allocated(this%sparseX_index)) deallocate(this%sparseX_index)
      if(allocated(this%sparseX)) deallocate(this%sparseX)
      if(allocated(this%alpha)) deallocate(this%alpha)
      if(allocated(this%sparseCutoff)) deallocate(this%sparseCutoff)

      if(allocated(this%x_size)) deallocate( this%x_size )
      if(allocated(this%xPrime_size)) deallocate( this%xPrime_size )
      if(allocated(this%covarianceDiag_x_x)) deallocate( this%covarianceDiag_x_x )
      if(allocated(this%covarianceDiag_x_xPrime)) deallocate( this%covarianceDiag_x_xPrime )
      if(allocated(this%covarianceDiag_xPrime_xPrime)) deallocate( this%covarianceDiag_xPrime_xPrime )

      if(allocated(this%sparseX_size)) deallocate( this%sparseX_size )
      if(allocated(this%covarianceDiag_sparseX_sparseX)) deallocate( this%covarianceDiag_sparseX_sparseX )


      call finalise(this%descriptor_str)
      call gpCoordinates_finalise_variance_estimate(this)  

      if(allocated(this%sparseX_permuted)) deallocate( this%sparseX_permuted )  
      if(allocated(this%sparseCovariance)) deallocate( this%sparseCovariance )  

    
      this%sparse_covariance_initialised = .false.  
  
      this%d = 0
      this%n_x = 0
      this%n_xPrime = 0
      this%delta = 0.0_dp
      this%f0 = 0.0_dp

      this%current_x = 0
      this%current_xPrime = 0

      this%n_sparseX = 0
      this%n_permutations = 0

      this%sparsified = .false.
      this%initialised = .false.

      if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call gpCovariance_bond_real_space_Finalise(this%bond_real_space_cov)

      this%covariance_type = COVARIANCE_NONE

   endsubroutine gpCoordinates_Finalise

   subroutine gpCovariance_bond_real_space_Finalise(this, error)
      type(gpCovariance_bond_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      this%n = 0
      this%delta = 0.0_dp
      this%atom_sigma = 0.0_dp

      this%initialised = .false.

   endsubroutine gpCovariance_bond_real_space_Finalise

   subroutine gpCovariance_atom_real_space_Finalise(this, error)
      type(gpCovariance_atom_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      this%l_max = 0
      this%delta = 0.0_dp

      this%initialised = .false.

   endsubroutine gpCovariance_atom_real_space_Finalise

   function gpFull_addFunctionValue(this,y,sigma_y, error)

      type(gpFull), intent(inout) :: this
      real(dp), intent(in) :: y, sigma_y           ! Function value
      integer :: gpFull_addFunctionValue  ! Which function value we added
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addFunctionValue: object not initialised',error)
      endif

      if( this%current_y == this%n_y ) then
         RAISE_ERROR( 'gpFull_addFunctionValue: object full, no more function values can be added',error)
      endif

      this%current_y = this%current_y + 1
      this%y(this%current_y) = y
      this%sigma_y(this%current_y) = sigma_y

      gpFull_addFunctionValue = this%current_y

   endfunction gpFull_addFunctionValue

   function gpFull_addFunctionDerivative(this, yPrime, sigma_yPrime, error)
      type(gpFull), intent(inout) :: this
      real(dp), intent(in) :: yPrime, sigma_yPrime ! Function value
      integer :: gpFull_addFunctionDerivative  ! Which function value we added
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addFunctionDerivative: object not initialised',error)
      endif

      if( this%current_yPrime == this%n_yPrime ) then
         RAISE_ERROR( 'gpFull_addFunctionDerivative: object full, no more function values can be added',error)
      endif

      this%current_yPrime = this%current_yPrime + 1
      this%yPrime(this%current_yPrime) = yPrime
      this%sigma_yPrime(this%current_yPrime) = sigma_yPrime

      gpFull_addFunctionDerivative = this%current_yPrime

   endfunction gpFull_addFunctionDerivative

   function gpFull_addCoordinates_2Darray(this,x,i_coordinate,cutoff_in, current_y, config_type, error) result(xLocation)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:,:), intent(in) :: x
      integer, intent(in) :: i_coordinate
      integer, optional, intent(in) :: current_y, config_type
      real(dp), dimension(:), intent(in), optional :: cutoff_in
      integer, optional, intent(out) :: error

      integer, dimension(:), pointer :: xLocation

      integer :: previous_x, i
      real(dp), dimension(:,:), allocatable :: new_x

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addCoordinates: object not initialised',error)
      endif

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_addCoordinates: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate ,error)
      endif

      if( .not. this%coordinate(i_coordinate)%initialised ) then
         RAISE_ERROR('gpFull_addCoordinates: '//i_coordinate//'th coordinate object is not initialised',error)
      endif

      if( this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
         if( size(x,1) > size(this%coordinate(i_coordinate)%x,1) ) then
            allocate( new_x(size(x,1),this%coordinate(i_coordinate)%n_x) )
            new_x = 0.0_dp
            new_x(1:size(this%coordinate(i_coordinate)%x,1),:) = this%coordinate(i_coordinate)%x
            deallocate( this%coordinate(i_coordinate)%x )
            allocate( this%coordinate(i_coordinate)%x(size(x,1),this%coordinate(i_coordinate)%n_x) )
            this%coordinate(i_coordinate)%x = new_x
            deallocate( new_x )
            this%coordinate(i_coordinate)%d = size(x,1)
         end if
      else
!         if( size(x,1) /= this%coordinate(i_coordinate)%d ) then
!            RAISE_ERROR('gpFull_addCoordinates: dimensionality of descriptors '//size(x,1)//' does not match what is given in the object '//this%coordinate(i_coordinate)%d,error)
!         endif
      endif

      previous_x = this%coordinate(i_coordinate)%current_x
      this%coordinate(i_coordinate)%current_x = previous_x + size(x,2)

      if( this%coordinate(i_coordinate)%current_x > this%coordinate(i_coordinate)%n_x ) then
         RAISE_ERROR('gpFull_addCoordinates: object full, no more descriptors can be added',error)
      endif

      if( this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
         this%coordinate(i_coordinate)%x(1:size(x,1),previous_x+1:this%coordinate(i_coordinate)%current_x) = x
         this%coordinate(i_coordinate)%x_size(previous_x+1:this%coordinate(i_coordinate)%current_x) = size(x,1)
      else
         this%coordinate(i_coordinate)%x(:,previous_x+1:this%coordinate(i_coordinate)%current_x) = x
      endif

      if(present(cutoff_in)) then
         this%coordinate(i_coordinate)%cutoff(previous_x+1:this%coordinate(i_coordinate)%current_x) = cutoff_in
      endif

      if(present(current_y)) &
         this%coordinate(i_coordinate)%map_x_y(previous_x+1:this%coordinate(i_coordinate)%current_x) = current_y

      if(present(config_type)) &
         this%coordinate(i_coordinate)%config_type(previous_x+1:this%coordinate(i_coordinate)%current_x) = config_type

      allocate(xLocation(size(x,2)))
      xLocation = (/ ( i, i = previous_x+1, this%coordinate(i_coordinate)%current_x ) /)

   endfunction gpFull_addCoordinates_2Darray

   function gpFull_addCoordinates_1Darray(this,x,i_coordinate,cutoff_in,current_y,config_type, error) result(xLocation)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: x
      integer, intent(in) :: i_coordinate
      real(dp), optional, intent(in) :: cutoff_in
      integer, optional, intent(in) :: current_y, config_type
      integer, optional, intent(out) :: error

      integer :: xLocation

      integer, dimension(:), pointer :: xLocation_in

      INIT_ERROR(error)

      xLocation_in => gpFull_addCoordinates_2Darray(this,reshape(x,(/size(x),1/)),i_coordinate,(/cutoff_in/),current_y,config_type,error)

      xLocation =  xLocation_in(1)
      deallocate(xLocation_in)

   endfunction gpFull_addCoordinates_1Darray

   subroutine gpFull_addCoordinateDerivatives_2Darray(this,xPrime,i_coordinate,current_yPrime, xLocation, dcutoff_in, error)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:,:), intent(in) :: xPrime
      integer, intent(in) :: i_coordinate, current_yPrime
      integer, dimension(:), intent(in) :: xLocation
      real(dp), dimension(:), optional, intent(in) :: dcutoff_in
      integer, optional, intent(out) :: error

      integer :: previous_xPrime
      real(dp), dimension(:,:), allocatable :: new_xPrime

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: object not initialised',error)
      endif

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_addCoordinateDerivatives: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      if( .not. this%coordinate(i_coordinate)%initialised ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: '//i_coordinate//'th coordinate object is not initialised',error)
      endif

      if( this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
         if( size(xPrime,1) > size(this%coordinate(i_coordinate)%xPrime,1) ) then
            allocate( new_xPrime(size(xPrime,1),this%coordinate(i_coordinate)%n_xPrime) )
            new_xPrime = 0.0_dp
            new_xPrime(1:size(this%coordinate(i_coordinate)%xPrime,1),:) = this%coordinate(i_coordinate)%xPrime
            deallocate( this%coordinate(i_coordinate)%xPrime )
            allocate( this%coordinate(i_coordinate)%xPrime(size(xPrime,1),this%coordinate(i_coordinate)%n_xPrime) )
            this%coordinate(i_coordinate)%xPrime = new_xPrime
            deallocate( new_xPrime )
         end if
      else
         if( size(xPrime,1) /= this%coordinate(i_coordinate)%d ) then
            RAISE_ERROR('gpFull_addCoordinateDerivatives: dimensionality of descriptors '//size(xPrime,1)//' does not match what is given in the object '//this%coordinate(i_coordinate)%d,error)
         endif
      endif

      if( size(xPrime,2) /= size(xLocation) ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: number of descriptors '//size(xPrime,2)//' has to match the dimensionality of the mapping array '//size(xLocation),error)
      endif

      previous_xPrime = this%coordinate(i_coordinate)%current_xPrime
      this%coordinate(i_coordinate)%current_xPrime = previous_xPrime + size(xPrime,2)

      if( this%coordinate(i_coordinate)%current_xPrime > this%coordinate(i_coordinate)%n_xPrime ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: object full, no more descriptors can be added',error)
      endif

      if( this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
         this%coordinate(i_coordinate)%xPrime(1:size(xPrime,1),previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xPrime
         this%coordinate(i_coordinate)%xPrime_size(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = size(xPrime,1)
      else
         this%coordinate(i_coordinate)%xPrime(:,previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xPrime
      endif

      if(present(dcutoff_in)) then
         this%coordinate(i_coordinate)%cutoffPrime(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = dcutoff_in
      endif

      this%coordinate(i_coordinate)%map_xPrime_yPrime(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = current_yPrime
      this%coordinate(i_coordinate)%map_xPrime_x(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xLocation

   endsubroutine gpFull_addCoordinateDerivatives_2Darray

   subroutine gpFull_addCoordinateDerivatives_1Darray(this,xPrime,i_coordinate,current_yPrime, xLocation, dcutoff_in, error)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: xPrime
      integer, intent(in) :: i_coordinate, current_yPrime
      integer, intent(in) :: xLocation
      real(dp), optional, intent(in) :: dcutoff_in
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call gpFull_addCoordinateDerivatives_2Darray(this, reshape(xPrime,(/size(xPrime),1/)),i_coordinate,current_yPrime,(/xLocation/),(/dcutoff_in/),error)

   endsubroutine gpFull_addCoordinateDerivatives_1Darray

   subroutine gpFull_addDescriptor(this,i_coordinate,descriptor_str,error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate 
      character(len=*), intent(in) :: descriptor_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addDescriptor: object not initialised',error)
      endif

      call gpCoordinates_addDescriptor(this%coordinate(i_coordinate),descriptor_str,error)

   endsubroutine gpFull_addDescriptor

   subroutine gpCoordinates_addDescriptor(this,descriptor_str,error)

      type(gpCoordinates), intent(inout) :: this
      character(len=*), intent(in) :: descriptor_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_addDescriptor: object not initialised',error)
      endif

      call initialise(this%descriptor_str)
      call zero(this%descriptor_str)

      call concat(this%descriptor_str,descriptor_str,keep_lf=.false.,lf_to_whitespace=.true.)

   endsubroutine gpCoordinates_addDescriptor

   subroutine gpCoordinates_setTheta(this, theta, zeta, error)
      type(gpCoordinates), intent(inout), target :: this
      real(dp), dimension(:), intent(in), optional :: theta
      real(dp), intent(in), optional :: zeta
      integer, optional, intent(out) :: error

      integer :: i
      real(dp) :: delta

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_setTheta: object not initialised',error)
      endif

      allocate(this%covarianceDiag_x_xPrime(size(this%x,1),this%n_x))
      this%covarianceDiag_x_xPrime = 0.0_dp

      select case(this%covariance_type)
      case(COVARIANCE_BOND_REAL_SPACE)
         RAISE_ERROR('gpCoordinates_setTheta: this call is not appropriate. Needs to be fixed!!!', error)
         !if(.not. this%bond_real_space_cov%initialised) then
         !   call gpCoordinates_gpCovariance_bond_real_space_Initialise(this)
         !endif

         !delta = this%bond_real_space_cov%delta
         !this%bond_real_space_cov%delta = 1.0_dp

         !do i = 1, this%n_x
         !   this%covarianceDiag_x_x(i) = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i = this%x(:,i), x_i_size = this%x_size(i), x_j = this%x(:,i), x_j_size = this%x_size(i))
         !enddo

!        ! do i = 1, this%n_xPrime
!        ! enddo

         !this%bond_real_space_cov%delta = delta

      case(COVARIANCE_ARD_SE,COVARIANCE_PP)
         call check_size('theta',theta,shape(this%theta),'gpCoordinates_setTheta',error)
         if( .not. present(theta) ) then
            RAISE_ERROR('gpCoordinates_setTheta: no theta present when using ARD_SE or PP for covariance', error)
         endif
         this%theta = theta  
      case(COVARIANCE_DOT_PRODUCT)
         if( .not. present(zeta) ) then
            RAISE_ERROR('gpCoordinates_setTheta: no zeta present when using DOT_PRODUCT for covariance', error)
         endif
         this%zeta = zeta
      endselect

   endsubroutine gpCoordinates_setTheta

   subroutine gpCoordinates_setThetaFactor(this, thetaFactor,useSparseX,error)
      type(gpCoordinates), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: thetaFactor
      logical, optional, intent(in) :: useSparseX
      integer, optional, intent(out) :: error

      integer :: i,p  
      logical :: my_useSparseX
      real(dp), dimension(this%d) :: theta, max_vals, min_vals  
  

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_calculateThetaFactor: object not initialised',error)
      endif

      if( .not. ( (this%covariance_type == COVARIANCE_ARD_SE) .or. (this%covariance_type == COVARIANCE_PP) ) ) then
         RAISE_ERROR('gpCoordinates_calculateThetaFactor: only ARD_SE or PP type covariance may use theta_fac',error)
      endif

      my_useSparseX = .false.
      if( allocated(this%sparseX_index) ) then
         if( sum( this%sparseX_index ) > 0 ) my_useSparseX = optional_default(.true.,useSparseX)
      endif

      if( my_useSparseX ) then
         do i = 1, this%d 
            max_vals(i) = maxval(this%x(i,this%sparseX_index))
            min_vals(i) = minval(this%x(i,this%sparseX_index))
         enddo
         do p=1, this%n_permutations !get max and min value of each dimension including all permutations
            do i=1, this%d
               max_vals(i) = max(max_vals(i),max_vals(this%permutations(i,p)))
               min_vals(i) = min(min_vals(i),min_vals(this%permutations(i,p)))
            enddo  
         enddo
         do i = 1, this%d          
            theta(i) = ( max_vals(i)- min_vals(i) ) * thetaFactor(i)
            if( theta(i) < THETA_MIN ) theta(i) = 1.0_dp
         enddo
      else
         do i = 1, this%d 
            max_vals(i) = maxval(this%x(i,:))
            min_vals(i) = minval(this%x(i,:))
         enddo
         do p=1, this%n_permutations !get max and min value of each dimension including all permutations
            do i=1, this%d
               max_vals(i) = max(max_vals(i),max_vals(this%permutations(i,p)))
               min_vals(i) = min(min_vals(i),min_vals(this%permutations(i,p)))
            enddo  
         enddo
         do i = 1, this%d          
            theta(i) = ( max_vals(i)- min_vals(i) ) * thetaFactor(i)
            if( theta(i) < THETA_MIN ) theta(i) = 1.0_dp
         enddo
      endif

      call gpCoordinates_setTheta(this,theta=theta,error=error)

   endsubroutine gpCoordinates_setThetaFactor

   subroutine gpFull_setTheta(this, i_coordinate, theta, zeta, error)
      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate
      real(dp), dimension(:), intent(in), optional :: theta
      real(dp), intent(in), optional :: zeta
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_setTheta: object not initialised',error)
      endif

      call gpCoordinates_setTheta(this%coordinate(i_coordinate), theta=theta, zeta=zeta, error=error)

   endsubroutine gpFull_setTheta

   subroutine gpFull_setTheta_thetaFactor(this, i_coordinate, thetaFactor, useSparseX, error)
      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate
      real(dp), dimension(:), intent(in) :: thetaFactor
      logical, optional, intent(in) :: useSparseX
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_setTheta_thetaFactor: object not initialised',error)
      endif

      call gpCoordinates_setThetaFactor(this%coordinate(i_coordinate), thetaFactor, useSparseX, error)

   endsubroutine gpFull_setTheta_thetaFactor

!   subroutine gpFull_setTheta_thetaFactorArray(this, thetaFactor, useSparseX, error)
!      type(gpFull), intent(inout) :: this
!      real(dp), dimension(:), intent(in) :: thetaFactor
!      logical, optional, intent(out) :: useSparseX
!      integer, optional, intent(out) :: error
!
!      integer :: i
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpFull_setTheta_thetaFactorArray: object not initialised',error)
!      endif
!
!      call check_size('thetaFactor',thetaFactor,(/this%n_coordinate/),'gpFull_setTheta_thetaFactorArray',error)
!
!      do i = 1, this%n_coordinate
!         call gpCoordinates_setThetaFactor(this%coordinate(i), thetaFactor(i), useSparseX, error)
!      enddo
!
!   endsubroutine gpFull_setTheta_thetaFactorArray
!
!   subroutine gpFull_setTheta_thetaFactorUniform(this, thetaFactor, useSparseX, error)
!      type(gpFull), intent(inout) :: this
!      real(dp), intent(in) :: thetaFactor
!      logical, optional, intent(out) :: useSparseX
!      integer, optional, intent(out) :: error
!
!      integer :: i
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpFull_setTheta_thetaFactorUniform: object not initialised',error)
!      endif
!
!      do i = 1, this%n_coordinate
!         call gpCoordinates_setThetaFactor(this%coordinate(i), thetaFactor, useSparseX, error)
!      enddo
!
!   endsubroutine gpFull_setTheta_thetaFactorUniform

   subroutine gpFull_covarianceMatrix_sparse(this,error)
      type(gpFull), intent(inout), target :: this
      integer, optional, intent(out) :: error
      integer :: i_coordinate, i_global_sparseX, j_global_sparseX, i_sparseX, j_sparseX, &
      n_globalY, i_globalY, i_global_yPrime, i_y, i_yPrime, i_x, j_x, n_x, i_xPrime, j_xPrime, n_xPrime, n, i
      real(dp) :: covariance_xPrime_sparseX, covariance_x_x_single, covariance_xPrime_xPrime, fc_i, fc_j, dfc_i, dfc_j
      real(dp), dimension(:,:,:,:), allocatable :: grad2_Covariance
      real(dp), dimension(:,:,:), allocatable :: grad_Covariance_j
      real(dp), dimension(:,:), allocatable :: grad_Covariance_i, covariance_x_x
      real(dp), dimension(:), allocatable :: covariance_x_sparseX, covariance_subY_currentX_y, covariance_subY_currentX_suby
      real(dp), dimension(:), pointer :: xPrime_i, xPrime_j
      integer, dimension(:), allocatable :: xIndex, xPrime_Index, xPrime_x_Index
      type(LA_matrix) :: LA_covariance_subY_subY
      logical :: found_i_x
      real(dp) :: start_time, cpu_time, wall_time

      
      INIT_ERROR(error)

      call system_timer('gpFull_covarianceMatrix_sparse')

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_covarianceMatrix: object not initialised',error)
      endif

      this%n_globalSparseX = 0
      i_global_sparseX = 0

      do i_coordinate = 1, this%n_coordinate
         if( .not. this%coordinate(i_coordinate)%initialised ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not initialised',error)
         endif

         if( .not. this%coordinate(i_coordinate)%sparsified ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not sparsified',error)
         endif

         if( .not. allocated(this%coordinate(i_coordinate)%x) ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//"th coordinate's x not allocated",error)
         endif

         if( .not. allocated(this%coordinate(i_coordinate)%xPrime) ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//"th coordinate's xPrime not allocated",error)
         endif

         if(this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            if(.not. this%coordinate(i_coordinate)%bond_real_space_cov%initialised) then  
               call gpCoordinates_gpCovariance_bond_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif

         this%n_globalSparseX = this%n_globalSparseX + this%coordinate(i_coordinate)%n_sparseX

         do i_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = i_global_sparseX + 1
            this%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX) = i_global_sparseX
         enddo
      enddo

      n_globalY = this%n_y + this%n_yPrime

      i_globalY = 0
      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = i_globalY + 1
         this%map_y_globalY(i_y) = i_globalY
      enddo
      do i_yPrime = 1, this%n_yPrime
         ! loop over all derivative values

         i_globalY = i_globalY + 1
         this%map_yPrime_globalY(i_yPrime) = i_globalY
      enddo

      call reallocate(this%covariance_subY_y, this%n_globalSparseX, n_globalY, zero = .true.)
      call reallocate(this%covariance_subY_subY, this%n_globalSparseX, this%n_globalSparseX, zero = .true.)
      call reallocate(this%covarianceDiag_y_y, n_globalY, zero = .true.)
      call reallocate(this%lambda, n_globalY, zero = .true.)
      call reallocate(this%inverse_sparse_full, this%n_globalSparseX, n_globalY, zero = .true.)

      allocate( covariance_subY_currentX_y(n_globalY),covariance_subY_currentX_suby(this%n_globalSparseX) )	
      covariance_subY_currentX_y = 0.0_dp	
      covariance_subY_currentX_suby = 0.0_dp

      do i_coordinate = 1, this%n_coordinate
         ! loop over different descriptor types
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate)
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_sparse')
         call print('Started sparse covariance matrix calculation of coordinate '//i_coordinate)
         call current_times(cpu_time, start_time)
         !
         do i_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
            ! loop over sparse points of each descriptor

            i_global_sparseX = this%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            ! find the unique number of the sparse point (to refer to it outside of the context of descriptor)

            allocate( grad_Covariance_i(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%n_x), &
                      covariance_x_sparseX(this%coordinate(i_coordinate)%n_x) )

            covariance_subY_currentX_y = 0.0_dp
            covariance_subY_currentX_suby = 0.0_dp
!$omp parallel do schedule(static,openmp_chunk_size) default(none) shared(openmp_chunk_size,this,i_coordinate,covariance_x_sparseX,grad_Covariance_i,i_sparseX) private(i_x,i_y,i_globalY) reduction(+:covariance_subY_currentX_y)
            do i_x = 1, this%coordinate(i_coordinate)%n_x
               ! loop over all data

               covariance_x_sparseX(i_x) = gpCoordinates_Covariance(this%coordinate(i_coordinate), &
		  i_x = i_x, j_sparseX = i_sparseX, grad_Covariance_i = grad_Covariance_i(:,i_x))

               i_y = this%coordinate(i_coordinate)%map_x_y(i_x)
               ! find which function value depends on the given descriptor

               if( i_y /= 0 ) then
                  i_globalY = this%map_y_globalY(i_y)
                  ! find unique function value/derivative identifier

                  !this%covariance_subY_y(i_global_sparseX, i_globalY) = this%covariance_subY_y(i_global_sparseX, i_globalY) + &
                  covariance_subY_currentX_y(i_globalY) = covariance_subY_currentX_y(i_globalY) + &
		     covariance_x_sparseX(i_x)*this%coordinate(i_coordinate)%cutoff(i_x)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)
               endif
            enddo
!$omp parallel do schedule(static,openmp_chunk_size) default(none) shared(openmp_chunk_size,this,i_coordinate,i_sparseX,grad_Covariance_i,covariance_x_sparseX) private(i_xPrime,i_yPrime,i_x,i_global_yPrime,covariance_xPrime_sparseX) reduction(+:covariance_subY_currentX_y)

	    do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
	       ! loop over all derivative data

	       i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
	       ! find which derivative depends on the given descriptor

	       i_x = this%coordinate(i_coordinate)%map_xPrime_x(i_xPrime)
	       if( i_yPrime /= 0 ) then
		  i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
		  ! find unique function value/derivative identifier

		  ! on Xeon w/ ifort 12, sum is fastest .  ddot is close.  dot_product is terrible
		  ! on Opteron w/ ifort 12 acml 5.2, ddot is 14.95 s, dot_product is 22.5 s, and sum is 13.9 s
		  ! dgemv doesn't seem noticeably faster at Opterons (may be faster on Xeon for 'N' transpose setting)
		  ! covariance_xPrime_sparseX = ddot(size(this%coordinate(i_coordinate)%xPrime,1),grad_Covariance_i(first_nonzero,i_x),1,this%coordinate(i_coordinate)%xPrime(1,i_xPrime),1)*&
		  ! covariance_xPrime_sparseX = dot_product(grad_Covariance_i(first_nonzero:last_nonzero,i_x),this%coordinate(i_coordinate)%xPrime(:,i_xPrime))* &
		  covariance_xPrime_sparseX = sum(grad_Covariance_i(:,i_x)*this%coordinate(i_coordinate)%xPrime(:,i_xPrime))* &
		     this%coordinate(i_coordinate)%cutoff(i_x)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX) + &
		     covariance_x_sparseX(i_x)*this%coordinate(i_coordinate)%cutoffPrime(i_xPrime)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)

		  !this%covariance_subY_y(i_global_sparseX, i_global_yPrime) = this%covariance_subY_y(i_global_sparseX, i_global_yPrime) + covariance_xPrime_sparseX
		  covariance_subY_currentX_y(i_global_yPrime) = covariance_subY_currentX_y(i_global_yPrime) + covariance_xPrime_sparseX
	       endif
	    enddo


            if(allocated(grad_Covariance_i)) deallocate(grad_Covariance_i)
            if(allocated(covariance_x_sparseX)) deallocate(covariance_x_sparseX)

!$omp parallel do schedule(static,openmp_chunk_size) default(none) shared(openmp_chunk_size,this,i_coordinate,covariance_x_sparseX,grad_Covariance_i,i_sparseX,i_global_sparseX) private(j_sparseX,j_global_sparseX) reduction(+:covariance_subY_currentX_suby)
            do j_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
               ! loop over sparse points of each descriptor

               j_global_sparseX = this%coordinate(i_coordinate)%map_sparseX_globalSparseX(j_sparseX)
               ! find the unique number of the sparse point (to refer to it outside of the context of descriptor)
               covariance_subY_currentX_suby(j_global_sparseX) = covariance_subY_currentX_suby(j_global_sparseX) + &
               gpCoordinates_Covariance(this%coordinate(i_coordinate), j_sparseX = j_sparseX, i_sparseX = i_sparseX) * this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)*this%coordinate(i_coordinate)%sparseCutoff(j_sparseX)
            enddo

            this%covariance_subY_subY(i_global_sparseX,i_global_sparseX) = this%covariance_subY_subY(i_global_sparseX,i_global_sparseX) + this%sparse_jitter
            this%covariance_subY_y(i_global_sparseX,:) = this%covariance_subY_y(i_global_sparseX,:) + covariance_subY_currentX_y
            this%covariance_subY_subY(:,i_global_sparseX) = this%covariance_subY_subY(:,i_global_sparseX) + covariance_subY_currentX_suby

            call current_times(cpu_time, wall_time)
            if(mod(i_sparseX,100) == 0) call progress_timer(this%coordinate(i_coordinate)%n_sparseX, i_sparseX, "Covariance matrix", wall_time-start_time)
         enddo
         call current_times(cpu_time, wall_time)
         call progress_timer(this%coordinate(i_coordinate)%n_sparseX, i_sparseX, "Covariance matrix", wall_time-start_time)
         call print('Finished sparse covariance matrix calculation of coordinate '//i_coordinate)
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_sparse')

         select case(this%covariance_method)
         case(GP_COVARIANCE_FITC)
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_diag')
            do i_y = 1, this%n_y
               ! loop over all function values

               i_globalY = this%map_y_globalY(i_y)
               ! find unique function value/derivative identifier

               n_x = count( this%coordinate(i_coordinate)%map_x_y == i_y )
               allocate(xIndex(n_x))
               
               n = 1
               do i = 1, size(this%coordinate(i_coordinate)%map_x_y)
                  if(this%coordinate(i_coordinate)%map_x_y(i) == i_y) then
                     xIndex(n) = i
                     n = n+1
                  endif
               enddo

               ! construct index array of all contributions to a given function value

               do i_x = 1, n_x
                  do j_x = 1, n_x
                     covariance_x_x_single = gpCoordinates_Covariance( this%coordinate(i_coordinate), i_x = xIndex(i_x), j_x = xIndex(j_x) ) * &
                     this%coordinate(i_coordinate)%cutoff(xIndex(i_x))*this%coordinate(i_coordinate)%cutoff(xIndex(j_x))
                     this%covarianceDiag_y_y(i_globalY) = this%covarianceDiag_y_y(i_globalY) + covariance_x_x_single
                  enddo
               enddo
               deallocate(xIndex)
            enddo
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_diag')

            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_grad_diag')

            do i_yPrime = 1, this%n_yPrime
               ! loop over all derivative values
               i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
               ! find unique function value/derivative identifier

               !n_xPrime = count( this%coordinate(i_coordinate)%map_xPrime_yPrime == i_yPrime )
               allocate(xPrime_Index(size(this%coordinate(i_coordinate)%map_xPrime_yPrime)))
               
               ! construct index array of all contributions to a given derivative value
               n = 0
               do i = 1, size(this%coordinate(i_coordinate)%map_xPrime_yPrime)
                  if(this%coordinate(i_coordinate)%map_xPrime_yPrime(i) == i_yPrime) then
                     n = n + 1
                     xPrime_Index(n) = i
                  endif
               enddo
               n_xPrime = n
               call reallocate(xPrime_Index,n_xPrime,copy=.true.)

               ! xIndex: for each xPrime that contributes to the given i_yPrime, this contains the indices 
               ! of the x in the main database. Length is less or equal to the number of xPrimes.
               ! xPrime_x_Index: this is the inverse of xIndex, i.e. for each xPrime this contains the index
               ! of x in xIndex.
               allocate(xPrime_x_Index(n_xPrime),xIndex(n_xPrime))
               xIndex = -1
               n_x = 1
               xIndex(1) = this%coordinate(i_coordinate)%map_xPrime_x(xPrime_Index(1))
               xPrime_x_Index(1) = 1
               ! Search started off.
               do i_xPrime = 2, n_xPrime
                  i_x = this%coordinate(i_coordinate)%map_xPrime_x(xPrime_Index(i_xPrime))
                  found_i_x = .false.
                  do i = 1, n_x
                     if(xIndex(i) == i_x) then
                        xPrime_x_Index(i_xPrime) = i
                        found_i_x = .true.
                        exit
                     endif
                  enddo
                  if(.not. found_i_x) then
                     n_x = n_x + 1
                     xIndex(n_x) = i_x
                     xPrime_x_Index(i_xPrime) = n_x
                  endif

                  !if( any(xIndex(:n_x) == i_x) ) then
                  !   do i = 1, n_x
                  !      if(xIndex(i) == i_x) xPrime_x_Index(i_xPrime) = i
                  !   enddo
                  !   cycle
                  !else
                  !   n_x = n_x + 1
                  !   xIndex(n_x) = i_x
                  !   xPrime_x_Index(i_xPrime) = n_x
                  !endif
               enddo
               call reallocate(xIndex,n_x,copy=.true.)


               allocate( grad2_Covariance(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%d,n_x,n_x), &
              grad_Covariance_j(this%coordinate(i_coordinate)%d,n_x,n_x), covariance_x_x(n_x,n_x) ) !, &
   !            xPrime_i_grad2_Covariance(this%coordinate(i_coordinate)%d) )
               grad2_Covariance = 0.0_dp
               do i_x = 1, n_x
                  do j_x = 1, n_x
                     covariance_x_x(j_x,i_x) = gpCoordinates_Covariance( this%coordinate(i_coordinate), i_x = xIndex(i_x), j_x = xIndex(j_x), &
                     grad_Covariance_j = grad_Covariance_j(:,j_x,i_x), grad2_Covariance=grad2_Covariance(:,:,j_x,i_x) )
                  enddo
               enddo

               do i_xPrime = 1, n_xPrime
                  i_x = xPrime_x_Index(i_xPrime)
                  xPrime_i => this%coordinate(i_coordinate)%xPrime(:,xPrime_Index(i_xPrime))
                  
                  fc_i = this%coordinate(i_coordinate)%cutoff(xIndex(i_x))
                  dfc_i = this%coordinate(i_coordinate)%cutoffPrime(xPrime_Index(i_xPrime))
                  do j_xPrime = 1, n_xPrime
                     j_x = xPrime_x_Index(j_xPrime)
                     xPrime_j => this%coordinate(i_coordinate)%xPrime(:,xPrime_Index(j_xPrime))

                     fc_j = this%coordinate(i_coordinate)%cutoff(xIndex(j_x))
                     dfc_j = this%coordinate(i_coordinate)%cutoffPrime(xPrime_Index(j_xPrime))

                     covariance_xPrime_xPrime = dot_product(matmul(xPrime_i,grad2_Covariance(:,:,j_x,i_x)),xPrime_j)*&
                     fc_i*fc_j + &
                     dot_product(grad_Covariance_j(:,j_x,i_x),xPrime_i)*fc_i*dfc_j + &
                     dot_product(grad_Covariance_j(:,j_x,i_x),xPrime_j)*fc_j*dfc_i + &
                     covariance_x_x(j_x,i_x)*dfc_i*dfc_j
                     !gpCoordinates_Covariance( this%coordinate(i_coordinate), i_xPrime = xPrime_Index(i_xPrime), j_xPrime = xPrime_Index(j_xPrime) )
                     this%covarianceDiag_y_y(i_global_yPrime) = this%covarianceDiag_y_y(i_global_yPrime) + covariance_xPrime_xPrime
                  enddo
               enddo

               if(allocated(xPrime_Index)) deallocate(xPrime_Index)
               if(allocated(xIndex)) deallocate(xIndex)
               if(allocated(xPrime_x_Index)) deallocate(xPrime_x_Index)
               if(allocated(grad2_Covariance)) deallocate(grad2_Covariance)
               if(allocated(grad_Covariance_j)) deallocate(grad_Covariance_j)
               if(allocated(covariance_x_x)) deallocate(covariance_x_x)
            enddo
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_grad_diag')
         case(GP_COVARIANCE_DTC)
            this%covarianceDiag_y_y = 0.0_dp
         case default
            this%covarianceDiag_y_y = 0.0_dp
         endselect

         if (allocated(this%coordinate(i_coordinate)%x)) deallocate(this%coordinate(i_coordinate)%x)
         if (allocated(this%coordinate(i_coordinate)%xPrime)) deallocate(this%coordinate(i_coordinate)%xPrime)
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate)
      enddo

      call system_timer('gpFull_covarianceMatrix_sparse_LinearAlgebra')
      call initialise(LA_covariance_subY_subY,this%covariance_subY_subY)
      call Matrix_Solve(LA_covariance_subY_subY, this%covariance_subY_y, this%inverse_sparse_full,error=error)
      call finalise(LA_covariance_subY_subY)
      call system_timer('gpFull_covarianceMatrix_sparse_LinearAlgebra')

      call system_timer('gpFull_covarianceMatrix_sparse_FunctionValues')
      select case(this%covariance_method)
      case(GP_COVARIANCE_FITC)
         do i_globalY = 1, n_globalY
            ! loop over all function values

            this%lambda(i_globalY) = this%covarianceDiag_y_y(i_globalY) - &
               dot_product( this%inverse_sparse_full(:,i_globalY), this%covariance_subY_y(:,i_globalY) )
         enddo
      case(GP_COVARIANCE_DTC)
         this%lambda = 0.0_dp
      case default
         this%lambda = 0.0_dp
      endselect

      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = this%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         this%lambda(i_globalY) = this%lambda(i_globalY) + &
            this%sigma_y(i_y)**2
      enddo

      do i_yPrime = 1, this%n_yPrime
         ! loop over all function values

         i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         this%lambda(i_global_yPrime) = this%lambda(i_global_yPrime) + &
            this%sigma_yPrime(i_yPrime)**2
      enddo
      call system_timer('gpFull_covarianceMatrix_sparse_FunctionValues')

      call system_timer('gpFull_covarianceMatrix_sparse')

   endsubroutine gpFull_covarianceMatrix_sparse

   subroutine gpFull_covarianceMatrix(this,error)
      type(gpFull), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i_coordinate, n_globalY, i_globalY, j_globalY, i_global_yPrime, j_global_yPrime, i_y, j_y, i_yPrime, j_yPrime, i_x, j_x, i_xPrime, j_xPrime
      real(dp) :: covariance_x_x
      real(dp), dimension(:), allocatable :: globalY
      real(dp), dimension(:,:), allocatable :: grad_Covariance_i
      logical :: is_i_xPrime

      type(LA_matrix) :: LA_covariance_y_y

      INIT_ERROR(error)

      call system_timer('gpFull_covarianceMatrix')

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_covarianceMatrix: object not initialised',error)
      endif

      do i_coordinate = 1, this%n_coordinate
         if( .not. this%coordinate(i_coordinate)%initialised ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not initialised',error)
         endif

         if(this%coordinate(i_coordinate)%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            if(.not. this%coordinate(i_coordinate)%bond_real_space_cov%initialised) then
               call gpCoordinates_gpCovariance_bond_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif
      enddo

      n_globalY = this%n_y + this%n_yPrime

      i_globalY = 0
      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = i_globalY + 1
         this%map_y_globalY(i_y) = i_globalY
      enddo
      do i_yPrime = 1, this%n_yPrime
         ! loop over all derivative values

         i_globalY = i_globalY + 1
         this%map_yPrime_globalY(i_yPrime) = i_globalY
      enddo

      call reallocate(this%covariance_y_y, n_globalY, n_globalY, zero = .true.)

      do i_coordinate = 1, this%n_coordinate
         ! loop over different descriptor types

!!$omp parallel schedule(dynamic) default(none) private(j_x, j_y, j_globalY, i_x, i_y, i_globalY, covariance_x_x, i_xPrime, i_yPrime, i_global_yPrime, j_global_yPrime, j_xPrime, j_yPrime) shared(this,i_coordinate)
!!$omp do
         do j_x = 1, this%coordinate(i_coordinate)%n_x
            ! loop over all data

            j_y = this%coordinate(i_coordinate)%map_x_y(j_x)
            ! find which function value depends on the given descriptor
            if( j_y /= 0 ) then
               j_globalY = this%map_y_globalY(j_y)
               ! find unique function value/derivative identifier

               allocate( grad_Covariance_i(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%n_x) )
               do i_x = 1, this%coordinate(i_coordinate)%n_x
                  ! loop over all data

                  i_y = this%coordinate(i_coordinate)%map_x_y(i_x)
                  ! find which function value depends on the given descriptor

                  is_i_xPrime = any(this%coordinate(i_coordinate)%map_xPrime_x == i_x)

                  if( (i_y /= 0 .and. i_y <= j_y) .or. is_i_xPrime) then

                     if(is_i_xPrime) then
                        covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_x = i_x, j_x = j_x, grad_Covariance_i=grad_Covariance_i(:,i_x))
                     else
                        covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_x = i_x, j_x = j_x)
                     endif

                     if( i_y /= 0 ) then
                        i_globalY = this%map_y_globalY(i_y)
                        ! find unique function value/derivative identifier
                        this%covariance_y_y(i_globalY, j_globalY) = this%covariance_y_y(i_globalY, j_globalY) + covariance_x_x
                     endif
                  endif

               enddo

               do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
                  ! loop over all derivative data

                  i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
                  ! find which derivative depends on the given descriptor

                  i_x = this%coordinate(i_coordinate)%map_xPrime_x(i_xPrime)

                  if( i_yPrime /= 0 ) then
                     i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
                     ! find unique function value/derivative identifier

                     covariance_x_x = dot_product(grad_Covariance_i(:,i_x),this%coordinate(i_coordinate)%xPrime(:,i_xPrime))
                     !gpCoordinates_Covariance(this%coordinate(i_coordinate), i_xPrime = i_xPrime, j_x = j_x)
                     this%covariance_y_y(i_global_yPrime, j_globalY) = this%covariance_y_y(i_global_yPrime, j_globalY) + covariance_x_x
                  endif
               enddo
            endif
            if(allocated(grad_Covariance_i)) deallocate(grad_Covariance_i)
            
         enddo
!!$omp end do

!!$omp do
         do j_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
            ! loop over all derivative data

            j_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(j_xPrime)
            ! find which derivative depends on the given descriptor

            if( j_yPrime /= 0 ) then
               j_global_yPrime = this%map_yPrime_globalY(j_yPrime)
               ! find unique function value/derivative identifier

               do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
                  ! loop over all derivative data

                  i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
                  ! find which derivative depends on the given descriptor

                  if( i_yPrime /= 0 .and. i_yPrime <= j_yPrime) then
                     i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
                     ! find unique function value/derivative identifier

                     call system_abort('not implemented yet')
                     !covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_xPrime = i_xPrime, j_xPrime = j_xPrime)
                     this%covariance_y_y(i_global_yPrime, j_global_yPrime) = this%covariance_y_y(i_global_yPrime, j_global_yPrime) + covariance_x_x
                  endif
               enddo
            endif
         enddo
!!$omp end parallel
      enddo

      do j_y = 1, size(this%covariance_y_y,2)
         do i_y = j_y + 1, size(this%covariance_y_y,1)
            this%covariance_y_y(i_y,j_y) = this%covariance_y_y(j_y,i_y)
         enddo
      enddo

      allocate( globalY(n_globalY) )

      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = this%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         this%covariance_y_y(i_globalY,i_globalY) = this%covariance_y_y(i_globalY,i_globalY) + &
            this%sigma_y(i_y)**2

         globalY(i_globalY) = this%y(i_y)
      enddo

      do i_yPrime = 1, this%n_yPrime
         ! loop over all function values

         i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         this%covariance_y_y(i_global_yPrime,i_global_yPrime) = this%covariance_y_y(i_global_yPrime,i_global_yPrime) + &
            this%sigma_yPrime(i_yPrime)**2

         globalY(i_global_yPrime) = this%y(i_yPrime)
      enddo

      call reallocate(this%alpha, n_globalY, zero = .true.)

      call initialise(LA_covariance_y_y,this%covariance_y_y)
      call Matrix_Solve(LA_covariance_y_y, globalY, this%alpha ,error=error)
      call finalise(LA_covariance_y_y)

      if(allocated(globalY)) deallocate(globalY)

      call system_timer('gpFull_covarianceMatrix')

   endsubroutine gpFull_covarianceMatrix

   function gpCoordinates_Covariance( this, i_x, j_x, i_sparseX, j_sparseX, grad_Covariance_i, grad_Covariance_j, grad2_Covariance, normalise, error )  
      type(gpCoordinates), intent(in), target :: this  
      integer, intent(in), optional :: i_x, j_x, i_sparseX, j_sparseX
      real(dp), dimension(:), optional, intent(out) :: grad_Covariance_i, grad_Covariance_j
      real(dp), dimension(:,:), optional, intent(out) :: grad2_Covariance
      logical, intent(in), optional :: normalise  
      integer, optional, intent(out) :: error

      real(dp) :: gpCoordinates_Covariance

      integer :: i_p, x_i_size, x_j_size, i
      integer :: ii, jj, zeta_int
      real(dp) :: covarianceExp, covarianceDiag_x_x_i, covarianceDiag_x_x_j, covarianceExp_ii, covarianceExp_jj, &
      gpCoordinates_Covariance_ii, gpCoordinates_Covariance_jj, normalisation, &
      covariancePP_ij, covariancePP_ii, covariancePP_jj, grad_covariancePP_ij, r_ij, r_ii, r_jj, grad_covariancePP_ii, grad_covariancePP_jj
      real(dp), dimension(:), pointer :: x_i, x_j, grad_Covariance_Diag_i, grad_Covariance_Diag_j
      !real(dp), dimension(this%d) :: inv_theta2, xI_xJ_theta2, xI_xJ, xI_xI, xI_xI_theta2, xJ_xJ, xJ_xJ_theta2, grad_Covariance_ii, grad_Covariance_jj
      real(dp), dimension(:),allocatable :: inv_theta2, xI_xJ_theta2, xI_xJ, xI_xI, xI_xI_theta2, xJ_xJ, xJ_xJ_theta2, grad_Covariance_ii, grad_Covariance_jj
      real(dp), dimension(:,:), allocatable :: distance_matrix
      logical :: do_normalise  

      INIT_ERROR(error)
      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_Covariance: object not initialised', error)
      endif

      do_normalise = optional_default(.true., normalise)  
      if( count( (/ present(i_x), present(i_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCoordinates_Covariance: exactly one of i_x or i_sparseX can be present', error)
      endif

      if( count( (/ present(j_x), present(j_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCoordinates_Covariance: exactly one of j_x or j_sparseX can be present', error)
      endif

      x_i => null()
      x_j => null()
      grad_Covariance_Diag_i => null()
      grad_Covariance_Diag_j => null()

      x_i_size = 0
      x_j_size = 0

      if(present(i_x)) then
         x_i => this%x(:,i_x)
         covarianceDiag_x_x_i = this%covarianceDiag_x_x(i_x)
         grad_Covariance_Diag_i => this%covarianceDiag_x_xPrime(:,i_x)

         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            x_i_size = this%x_size(i_x)
         endif
      endif

      if(present(j_x)) then
         x_j => this%x(:,j_x)
         covarianceDiag_x_x_j = this%covarianceDiag_x_x(j_x)
         grad_Covariance_Diag_j => this%covarianceDiag_x_xPrime(:,j_x)

         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            x_j_size = this%x_size(j_x)
         endif
      endif

      if(present(i_sparseX)) then
         x_i => this%sparseX(:,i_sparseX)
         covarianceDiag_x_x_i = this%covarianceDiag_sparseX_sparseX(i_sparseX)

         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            x_i_size = this%sparseX_size(i_sparseX)
         endif
      endif

      if(present(j_sparseX)) then
         x_j => this%sparseX(:,j_sparseX)
         covarianceDiag_x_x_j = this%covarianceDiag_sparseX_sparseX(j_sparseX)

         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            x_j_size = this%sparseX_size(j_sparseX)
         endif
      endif

      if( .not. associated(x_i) .or. .not. associated(x_j) ) then  
         RAISE_ERROR('gpCoordinates_Covariance: both i and j indices have to be present', error)
      endif

      gpCoordinates_Covariance = 0.0_dp
      if(present(grad_Covariance_i)) then
         grad_Covariance_i = 0.0_dp
      endif
      if(present(grad_Covariance_j)) then
         grad_Covariance_j = 0.0_dp
      endif
      if(present(grad2_Covariance)) then
         grad2_Covariance = 0.0_dp
      endif

      if(this%covariance_type == COVARIANCE_ARD_SE .or. this%covariance_type == COVARIANCE_PP) then
         allocate(inv_theta2(this%d), &
                  xI_xJ(this%d), &
                  xI_xI(this%d), &
                  xJ_xJ(this%d), &
                  xI_xJ_theta2(this%d), &
                  xI_xI_theta2(this%d), &
                  xJ_xJ_theta2(this%d), &
                  grad_Covariance_ii(this%d), &
                  grad_Covariance_jj(this%d))
      endif

      if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
         if(present(i_x)) then
            if(present(j_x)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=x_i, x_i_size=x_i_size, x_j=x_j, x_j_size=x_j_size) &  
                                          / sqrt(covarianceDiag_x_x_i * covarianceDiag_x_x_j)  
            elseif(present(j_sparseX)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=x_i, x_i_size=x_i_size, x_j=x_j, x_j_size=x_j_size) &  
                                          / sqrt(covarianceDiag_x_x_i * covarianceDiag_x_x_j)  
            elseif(present(grad_Covariance_j)) then  
               RAISE_ERROR('gpCoordinates_Covariance: bond real space derivatives not implemented', error)  
            endif
         elseif(present(i_sparseX)) then
            if(present(j_x)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=x_i, x_i_size=x_i_size, x_j=x_j, x_j_size=x_j_size) &  
                                          / sqrt(covarianceDiag_x_x_i * covarianceDiag_x_x_j)  
            elseif(present(j_sparseX)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=x_i, x_i_size=x_i_size, x_j=x_j, x_j_size=x_j_size) &  
                                          / sqrt(covarianceDiag_x_x_i * covarianceDiag_x_x_j)  
            elseif(present(grad_Covariance_j)) then  
               RAISE_ERROR('gpCoordinates_Covariance: bond real space derivatives not implemented', error)  
            endif
         elseif(present(grad_Covariance_i)) then  
            RAISE_ERROR('gpCoordinates_Covariance: bond real space derivatives not implemented', error)  
         endif
      elseif(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
         gpCoordinates_Covariance = sum(x_i*x_j)

         zeta_int = nint(this%zeta)
         if( zeta_int .feq. this%zeta ) then
            if(present(grad_Covariance_i)) grad_Covariance_i = this%delta**2 * zeta_int * gpCoordinates_Covariance**(zeta_int-1) * x_j
            if(present(grad_Covariance_j)) grad_Covariance_j = this%delta**2 * zeta_int * gpCoordinates_Covariance**(zeta_int-1) * x_i
            gpCoordinates_Covariance = this%delta**2 * gpCoordinates_Covariance**zeta_int
         else
            if(present(grad_Covariance_i)) grad_Covariance_i = this%delta**2 * this%zeta * gpCoordinates_Covariance**(this%zeta-1.0_dp) * x_j
            if(present(grad_Covariance_j)) grad_Covariance_j = this%delta**2 * this%zeta * gpCoordinates_Covariance**(this%zeta-1.0_dp) * x_i
            gpCoordinates_Covariance = this%delta**2 * gpCoordinates_Covariance**this%zeta
         endif
      elseif(this%covariance_type == COVARIANCE_ARD_SE ) then
         inv_theta2 = 1.0_dp / this%theta**2  
         do i_p = 1, this%n_permutations
            ! permute only i. theta should be symmetrised by now.

            do ii = 1, this%d  
               xI_xJ(ii) = x_i(this%permutations(ii,i_p)) - x_j(ii)  
            end do  
            xI_xJ_theta2 = xI_xJ * inv_theta2  

            covarianceExp = this%delta**2 * exp( -0.5_dp * dot_product(xI_xJ_theta2,xI_xJ) )
            gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp

            if(present(grad_Covariance_i)) then
               do ii = 1, this%d  
                  grad_Covariance_i(this%permutations(ii,i_p)) = grad_Covariance_i(this%permutations(ii,i_p)) - covarianceExp * xI_xJ_theta2(ii)  
               end do  
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_j = grad_Covariance_j + covarianceExp * xI_xJ_theta2
            endif

            if(present(grad2_Covariance)) then
               do i = 1, this%d
                  grad2_Covariance(:,this%permutations(i,i_p)) = grad2_Covariance(:,this%permutations(i,i_p)) - covarianceExp * &
                  xI_xJ_theta2*xI_xJ_theta2(i)
                  grad2_Covariance(this%permutations(i,i_p),i) = grad2_Covariance(this%permutations(i,i_p),i) + covarianceExp * inv_theta2(i)  
               enddo
            endif

            !if(present(i_xPrime) .and. .not. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance - covarianceExp * (dot_product(xI_xJ_theta,xPrime_i_theta(this%permutations(:,i_p)))*fc_i - dfc_i)*fc_j

            !elseif(.not. present(i_xPrime) .and. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp * (dot_product(xI_xJ_theta,xPrime_j_theta)*fc_j + dfc_j)*fc_i

            !elseif(present(i_xPrime) .and. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp * ( dot_product( xPrime_i_theta(this%permutations(:,i_p)), xPrime_j_theta )*fc_i*fc_j + &
            !   ( - dot_product( xI_xJ_theta, xPrime_i_theta(this%permutations(:,i_p)) )*fc_i + dfc_i ) * &
            !   ( dot_product( xI_xJ_theta, xPrime_j_theta )*fc_j + dfc_j ) )

            !else
            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp*fc_i*fc_j
            !endif
         enddo

         if( this%n_permutations > 1 .and. do_normalise ) then  
            gpCoordinates_Covariance_ii = 0.0_dp
            gpCoordinates_Covariance_jj = 0.0_dp

            if(present(grad_Covariance_i)) then
               grad_Covariance_ii = 0.0_dp
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_jj = 0.0_dp
            endif

            do i_p = 1, this%n_permutations

               do ii = 1, this%d  
                  xI_xI(ii) = x_i(this%permutations(ii,i_p)) - x_i(ii)  
               enddo
               xI_xI_theta2 = xI_xI * inv_theta2  

               do ii = 1, this%d  
                  xJ_xJ(ii) = x_j(this%permutations(ii,i_p)) - x_j(ii)  
               enddo
               xJ_xJ_theta2 = xJ_xJ * inv_theta2  

               covarianceExp_ii = exp( -0.5_dp * dot_product(xI_xI_theta2,xI_xI) )
               covarianceExp_jj = exp( -0.5_dp * dot_product(xJ_xJ_theta2,xJ_xJ) )

               gpCoordinates_Covariance_ii = gpCoordinates_Covariance_ii + covarianceExp_ii
               gpCoordinates_Covariance_jj = gpCoordinates_Covariance_jj + covarianceExp_jj

               if(present(grad_Covariance_i)) then
                  grad_Covariance_ii = grad_Covariance_ii + covarianceExp_ii * xI_xI_theta2
                  grad_Covariance_ii(this%permutations(:,i_p)) = grad_Covariance_ii(this%permutations(:,i_p)) &
                     - covarianceExp_ii * xI_xI_theta2
               endif

               if(present(grad_Covariance_j)) then
                  grad_Covariance_jj = grad_Covariance_jj + covarianceExp_jj * xJ_xJ_theta2
                  grad_Covariance_jj(this%permutations(:,i_p)) = grad_Covariance_jj(this%permutations(:,i_p)) &
                     - covarianceExp_jj * xJ_xJ_theta2
               endif

               if(present(grad2_Covariance)) then
                  RAISE_ERROR('grad2_Covariance for n_permutations > 1 not implemented yet',error)
               endif
            enddo

            normalisation = sqrt(gpCoordinates_Covariance_ii * gpCoordinates_Covariance_jj)

            if(present(grad_Covariance_i)) then
               grad_Covariance_i = grad_Covariance_i / normalisation - 0.5_dp * grad_Covariance_ii * gpCoordinates_Covariance / normalisation / gpCoordinates_Covariance_ii
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_j = grad_Covariance_j / normalisation - 0.5_dp * grad_Covariance_jj * gpCoordinates_Covariance / normalisation / gpCoordinates_Covariance_jj
            endif
            
            gpCoordinates_Covariance = gpCoordinates_Covariance / normalisation
         else
            normalisation = 1.0_dp
         endif

         gpCoordinates_Covariance = gpCoordinates_Covariance + this%f0**2
      elseif(this%covariance_type == COVARIANCE_PP ) then
         allocate(distance_matrix(this%d, this%d))

         inv_theta2 = 1.0_dp / this%theta**2

         forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) distance_matrix(ii,jj) = ( x_i(ii) - x_j(jj) )**2 / this%theta(ii)**2
         do i_p = 1, this%n_permutations
            if( any( (/ (distance_matrix(this%permutations(ii,i_p),ii) > 1.0_dp, ii=1, this%d) /) ) ) cycle

            r_ij = sqrt( sum( (/ (distance_matrix(this%permutations(ii,i_p),ii), ii=1, this%d) /) ) )
            if( r_ij >= 1.0_dp ) cycle

            covariancePP_ij = this%delta**2 * covariancePP(r_ij,PP_Q, this%d)
            gpCoordinates_Covariance = gpCoordinates_Covariance + covariancePP_ij

            if( ( present(grad_Covariance_i) .or. present(grad_Covariance_j) ) .and. (r_ij > 0.0_dp) ) then
               grad_covariancePP_ij = grad_covariancePP(r_ij,PP_Q, this%d) / r_ij
               !xI_xJ(:) = x_i(this%permutations(:,i_p)) - x_j(:)

               if(present(grad_Covariance_i)) &
                  grad_Covariance_i(this%permutations(:,i_p)) = grad_Covariance_i(this%permutations(:,i_p)) + grad_covariancePP_ij * xI_xJ(:)

               if(present(grad_Covariance_j)) &
                  grad_Covariance_j(:) = grad_Covariance_j(:) - grad_covariancePP_ij * xI_xJ(:)
            endif
         enddo ! i_p

         if(present(grad_Covariance_i)) grad_Covariance_i = grad_Covariance_i * inv_theta2
         if(present(grad_Covariance_j)) grad_Covariance_j = grad_Covariance_j * inv_theta2

         do_normalise = do_normalise .and. ( gpCoordinates_Covariance .fne. 0.0_dp )
         if( this%n_permutations > 1 .and. do_normalise ) then
            gpCoordinates_Covariance_ii = 0.0_dp
            gpCoordinates_Covariance_jj = 0.0_dp

            if(present(grad_Covariance_i)) then
               grad_Covariance_ii = 0.0_dp
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_jj = 0.0_dp
            endif

            forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) distance_matrix(ii,jj) = ( x_i(ii) - x_i(jj) )**2 * inv_theta2(ii)
            do i_p = 1, this%n_permutations
               if( any( (/ (distance_matrix(this%permutations(ii,i_p),ii) > 1.0_dp, ii=1, this%d) /) ) ) cycle

               r_ii = sqrt( sum( (/ (distance_matrix(this%permutations(ii,i_p),ii), ii=1, this%d) /) ) )
               if( r_ii >= 1.0_dp ) cycle

               covariancePP_ii = covariancePP(r_ii,PP_Q, this%d)
               gpCoordinates_Covariance_ii = gpCoordinates_Covariance_ii + covariancePP_ii

               if(present(grad_Covariance_i) .and. (r_ii > 0.0_dp)) then
                  xI_xI(:) = x_i(this%permutations(:,i_p)) - x_i(:)

                  grad_covariancePP_ii = grad_covariancePP(r_ii,PP_Q, this%d) / r_ii
                  grad_Covariance_ii = grad_Covariance_ii + grad_covariancePP_ii * xI_xI(:)
                  grad_Covariance_ii(this%permutations(:,i_p)) = grad_Covariance_ii(this%permutations(:,i_p)) - grad_covariancePP_ii * xI_xI
               endif
            enddo
            if(present(grad_Covariance_i)) grad_Covariance_ii = grad_Covariance_ii * inv_theta2

            forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) distance_matrix(ii,jj) = ( x_j(ii) - x_j(jj) )**2 * inv_theta2(ii)
            do i_p = 1, this%n_permutations
               if( any( (/ (distance_matrix(this%permutations(ii,i_p),ii) > 1.0_dp, ii=1, this%d) /) ) ) cycle

               r_jj = sqrt( sum( (/ (distance_matrix(this%permutations(ii,i_p),ii), ii=1, this%d) /) ) )
               if( r_jj >= 1.0_dp ) cycle

               covariancePP_jj = covariancePP(r_jj,PP_Q, this%d)
               gpCoordinates_Covariance_jj = gpCoordinates_Covariance_jj + covariancePP_jj

               if(present(grad_Covariance_j) .and. (r_jj > 0.0_dp)) then
                  xJ_xJ(:) = x_j(this%permutations(:,i_p)) - x_j(:)

                  grad_covariancePP_jj = grad_covariancePP(r_jj,PP_Q, this%d) / r_jj
                  grad_Covariance_jj = grad_Covariance_jj + grad_covariancePP_jj * xJ_xJ(:)
                  grad_Covariance_jj(this%permutations(:,i_p)) = grad_Covariance_jj(this%permutations(:,i_p)) - grad_covariancePP_jj * xJ_xJ(:)
               endif
            enddo
            if(present(grad_Covariance_j)) grad_Covariance_jj = grad_Covariance_jj * inv_theta2

            normalisation = sqrt(gpCoordinates_Covariance_ii * gpCoordinates_Covariance_jj)

            if(present(grad_Covariance_i)) then
               grad_Covariance_i = grad_Covariance_i / normalisation - 0.5_dp * grad_Covariance_ii * gpCoordinates_Covariance / normalisation / gpCoordinates_Covariance_ii
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_j = grad_Covariance_j / normalisation - 0.5_dp * grad_Covariance_jj * gpCoordinates_Covariance / normalisation / gpCoordinates_Covariance_jj
            endif
            
            gpCoordinates_Covariance = gpCoordinates_Covariance / normalisation
         else
            normalisation = 1.0_dp
         endif

         gpCoordinates_Covariance = gpCoordinates_Covariance + this%f0**2
      endif ! this%covariance_type

      if(allocated(inv_theta2)) deallocate(inv_theta2)
      if(allocated(xI_xJ)) deallocate(xI_xJ)
      if(allocated(xI_xI)) deallocate(xI_xI)
      if(allocated(xJ_xJ)) deallocate(xJ_xJ)
      if(allocated(xI_xJ_theta2)) deallocate(xI_xJ_theta2)
      if(allocated(xI_xI_theta2)) deallocate(xI_xI_theta2)
      if(allocated(xJ_xJ_theta2)) deallocate(xJ_xJ_theta2)
      if(allocated(grad_Covariance_ii)) deallocate(grad_Covariance_ii)
      if(allocated(grad_Covariance_jj)) deallocate(grad_Covariance_jj)
      if( allocated( distance_matrix ) ) deallocate(distance_matrix)

   endfunction gpCoordinates_Covariance

   subroutine gpCoordinates_gpCovariance_bond_real_space_Initialise( this, error )
      type(gpCoordinates), intent(inout) :: this
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      if (.not. this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
         RAISE_ERROR('gpCoordinates_gpCovariance_bond_real_space_Initialise: covariance is not bond_real_space', error)
      endif

      call gpCovariance_bond_real_space_Finalise(this%bond_real_space_cov, error)

      call initialise(params)
      call param_register(params, 'n', '2', this%bond_real_space_cov%n, &
           help_string="Covariance degree for bond_real_space-type descriptors")
      call param_register(params, 'atom_sigma', '0.0', this%bond_real_space_cov%atom_sigma, &
           help_string="Atoms sigma for bond_real_space-type descriptors")

      if (.not. param_read_line(params, string(this%descriptor_str), ignore_unknown=.true., task='gpCoordinates_gpCovariance_bond_real_space_Initialise descriptor_str')) then
         RAISE_ERROR("gpCoordinates_gpCovariance_bond_real_space_Initialise failed to parse descriptor_str='"//trim(string(this%descriptor_str))//"'", error)
      endif
      call finalise(params)

      this%bond_real_space_cov%delta = this%delta

      this%bond_real_space_cov%initialised = .true.

   endsubroutine gpCoordinates_gpCovariance_bond_real_space_Initialise

   function gpCovariance_bond_real_space_Calc( this, x_i, x_i_size, x_j, x_j_size, xPrime_i, xPrime_j, xPrime_ij, error )  
      type(gpCovariance_bond_real_space), intent(in) :: this
      real(dp), intent(in) :: x_i(0:), x_j(0:)  
      integer, intent(in) :: x_i_size, x_j_size  
      real(dp), dimension(:), intent(out), optional, pointer :: xPrime_i, xPrime_j  
      real(dp), dimension(:,:), intent(out), optional, pointer :: xPrime_ij  
      integer, intent(out), optional :: error  

      real(dp) :: gpCovariance_bond_real_space_Calc
      real(dp) :: gpCovariance_bond_real_space_Calc_compensation ! Running compensation for Kahan summation algorithm  

      integer :: i, j, k, m, n
      integer :: x_i_N, x_j_N
      complex(dp), allocatable :: gamma_i(:,:), gamma_j(:,:)  
      real(dp), allocatable :: z_i(:), z_j(:), c_i(:), c_j(:)  
      real(dp) :: x_i_self_overlap, x_j_self_overlap  
      integer, allocatable :: iter_index(:)  
  
      logical :: do_derivative, do_xPrime_i, do_xPrime_j, do_xPrime_ij  

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: object not initialised', error)
      endif

      do_xPrime_i = .false.  
      do_xPrime_j = .false.  
      do_xPrime_ij = .false.  
      if(present(xPrime_i)) do_xPrime_i = associated(xPrime_i)  
      if(present(xPrime_j)) do_xPrime_j = associated(xPrime_j)  
      if(present(xPrime_ij)) do_xPrime_ij = associated(xPrime_ij)  
  
      do_derivative = (do_xPrime_i .or. do_xPrime_j .or. do_xPrime_ij)  
  
      if( do_derivative ) then  
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: derivatives not implemented', error)  
      endif

      !  
      ! x_i: get x_i_N, gamma_i, z_i, c_i, x_i_self_overlap  
      !  
      x_i_N = nint(x_i(0))   

      allocate( gamma_i(x_i_N,x_i_N), z_i(x_i_N), c_i(x_i_N) )  

      gamma_i = cmplx(0.0_dp, 0.0_dp, dp)  
      z_i = 0.0_dp  
      c_i = 0.0_dp  

      do i = 1, x_i_N  
         do j = 1, x_i_N  
            gamma_i(i,j) = cmplx(x_i(1 + x_i_N + (2 * (i - 1) * x_i_N) + (2 * j) - 1), x_i(1 + x_i_N + (2 * (i - 1) * x_i_N) + (2 * j)), dp)  

            if (i == j) then  
               gamma_i(i,j) = cmplx(x_i(1 + x_i_N + (2 * (i - 1) * x_i_N) + (2 * j) - 1), 0.0_dp, dp)  

               z_i(i) = x_i(1 + x_i_N + (2 * (i - 1) * x_i_N) + (2 * j))  
            endif  
         enddo
      enddo  

      x_i_self_overlap = x_i(1)  

      c_i = x_i(2:x_i_N + 1)  

      !  
      ! x_j: get x_j_N, gamma_j, z_j, c_j, x_j_self_overlap  
      !  
      x_j_N = nint(x_j(0))  

      allocate( gamma_j(x_j_N,x_j_N), z_j(x_j_N), c_j(x_j_N) )  

      gamma_j = cmplx(0.0_dp, 0.0_dp, dp)  
      z_j = 0.0_dp  
      c_j = 0.0_dp  

      do i = 1, x_j_N  
         do j = 1, x_j_N  
            gamma_j(i,j) = cmplx(x_j(1 + x_j_N + (2 * (i - 1) * x_j_N) + (2 * j) - 1), x_j(1 + x_j_N + (2 * (i - 1) * x_j_N) + (2 * j)), dp)  

            if (i == j) then  
               gamma_j(i,j) = cmplx(x_j(1 + x_j_N + (2 * (i - 1) * x_j_N) + (2 * j) - 1), 0.0_dp, dp)  

               z_j(i) = x_j(1 + x_j_N + (2 * (i - 1) * x_j_N) + (2 * j))  
            endif  
         enddo
      enddo  

      x_j_self_overlap = x_j(1)  

      c_j = x_j(2:x_j_N + 1)  

      !  
      ! Start with gpCovariance_bond_real_space_Calc = 0  
      !  
      gpCovariance_bond_real_space_Calc = 0.0_dp  
      gpCovariance_bond_real_space_Calc_compensation = 0.0_dp ! Running compensation for Kahan summation algorithm  

      allocate( iter_index(this%n) )  

      call gpCovariance_bond_real_space_sum(1, 1)  

      gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc * (2.0_dp / (x_i_self_overlap + x_j_self_overlap))**this%n  

      contains  

      recursive subroutine gpCovariance_bond_real_space_sum(n, iter)  
         integer :: n, iter  
         integer :: i(n), j(n), k, m, l  
         integer :: iter_index_sorted(n), powers(n)  
         real(dp) :: coefficient, arg_z, arg_r2, arg_gamma  
         complex(dp) :: arg_gamma_complex  
         real(dp) :: x, x_mirror  
         real(dp) :: y, t ! Intermediate variables for Kahan summation algorithm  

         do k = iter, x_i_N*x_j_N  
            iter_index(n) = k  

            if (n < this%n) then  
               call gpCovariance_bond_real_space_sum(n + 1, k)  
            else  
               iter_index_sorted = iter_index  
               powers = 1  
               if (n > 1) then  
                  !  
                  ! For large n could replace with heapsort?  
                  !  
                  call sort_array(iter_index_sorted)  
  
                  do m = 1, (n - 1)  
                     if (iter_index_sorted(m) == iter_index_sorted(m + 1)) then  
                        powers(m + 1) = powers(m + 1) + powers(m)  
                        powers(m) = 0  
                     endif  
                  enddo  
               endif  

               do m = 1, n  
                  i(m) = 1 + ((iter_index(m) - 1) / x_i_N)  
                  j(m) = 1 + mod(iter_index(m) - 1, x_j_N)  
               enddo  

               coefficient = factorial(n)  
               arg_z = 0.0_dp  
               arg_r2 = 0.0_dp  
               arg_gamma_complex = cmplx(0.0_dp, 0.0_dp, dp)  
               do m = 1, n  
                  if (powers(m) /= 0) then  
                     coefficient = coefficient / factorial(powers(m))  
                     coefficient = coefficient * (c_i(i(m)) * c_j(j(m)))**powers(m)  
                  endif  

                  arg_z = arg_z + z_i(i(m)) * z_j(j(m))  
                  arg_r2 = arg_r2 + real(gamma_i(i(m),i(m)), dp) + real(gamma_j(j(m),j(m)), dp)  
                  arg_gamma_complex = arg_gamma_complex + gamma_i(i(m),i(m)) * conjg(gamma_j(j(m),j(m)))  
                  if (m < n) then  
                     do l = (m + 1), n  
                        arg_gamma_complex = arg_gamma_complex + 2.0_dp * gamma_i(i(m),i(l)) * conjg(gamma_j(j(m),j(l)))  
                     enddo  
                  endif  
               enddo  
               arg_z = arg_z / (2.0_dp * this%atom_sigma**2)  
               arg_r2 = - arg_r2 / (4.0_dp * this%atom_sigma**2)  
               arg_gamma = sqrt(real(arg_gamma_complex, dp)) / (2.0_dp * this%atom_sigma**2)  
  
               !  
               ! Use asymptotic expansion of Bessel function  
               !  
               if (arg_gamma < besseli_max_x) then  
                  x = besseli0(arg_gamma)  
  
                  x_mirror = exp(- arg_z + arg_r2) * x  
                  x = exp(arg_z + arg_r2) * x  
               else
                  x = 1.0_dp  
                  do m = 1, besseli_max_n  
                     x = x + besseli0_c(m) / arg_gamma**m  
                  enddo  
                  x = x / sqrt(2.0_dp * pi * arg_gamma)  
  
                  x_mirror = exp(- arg_z + arg_r2 + arg_gamma) * x  
                  x = exp(arg_z + arg_r2 + arg_gamma) * x  
               endif
  
               x_mirror = coefficient * x_mirror  
               x = coefficient * x  
  
               !  
               ! Kahan summation algorithm for x_mirror  
               !  
               y = x_mirror - gpCovariance_bond_real_space_Calc_compensation  
               t = gpCovariance_bond_real_space_Calc + y  
               gpCovariance_bond_real_space_Calc_compensation = (t - gpCovariance_bond_real_space_Calc) - y  
               gpCovariance_bond_real_space_Calc = t  
               !  
               ! Kahan summation algorithm for x  
               !  
               y = x - gpCovariance_bond_real_space_Calc_compensation  
               t = gpCovariance_bond_real_space_Calc + y  
               gpCovariance_bond_real_space_Calc_compensation = (t - gpCovariance_bond_real_space_Calc) - y  
               gpCovariance_bond_real_space_Calc = t  
            endif
         enddo

      endsubroutine gpCovariance_bond_real_space_sum  

   endfunction gpCovariance_bond_real_space_Calc  

   function besseli0(x)  

      real(dp), intent(in) :: x  
      real(dp) :: besseli0  

      real(dp) :: x2, r, k  
      integer :: i  

      x2 = x**2  
  
      if(x == 0.0_dp) then  
         besseli0 = 1.0_dp  
      elseif( x < besseli_max_x ) then  
         besseli0 = 1.0_dp  
         r = 1.0_dp  
         k = 1.0_dp  
         do while ( abs(r/besseli0) > NUMERICAL_ZERO )  
            r = 0.25_dp * r * x2 / k**2  
            besseli0 = besseli0 + r  
            k = k + 1.0_dp  
         enddo  
      else  
         besseli0 = 1.0_dp  
         do i = 1, besseli_max_n  
            besseli0 = besseli0 + besseli0_c(i)/x**i  
         enddo  
         besseli0 = besseli0 * exp(x) / sqrt(2.0_dp*pi*x)  
      endif  

   endfunction besseli0  

   function besseli1(x)  

      real(dp), intent(in) :: x  
      real(dp) :: besseli1  

      real(dp) :: x2, r, k  
      integer :: i  

      x2 = x**2  
  
      if(x == 0.0_dp) then  
         besseli1 = 0.0_dp  
      elseif( x < besseli_max_x ) then  
         besseli1 = 1.0_dp  
         r = 1.0_dp  
         k = 1.0_dp  
         do while ( abs(r/besseli1) > NUMERICAL_ZERO )  
            r = 0.25_dp * r * x2 / (k*(k+1.0_dp))  
            besseli1 = besseli1 + r  
            k = k + 1.0_dp  
         enddo  
         besseli1 = besseli1 * 0.5_dp * x  
      else  
         besseli1 = 1.0_dp  
         do i = 1, besseli_max_n  
            besseli1 = besseli1 + besseli1_c(i)/x**i  
         enddo  
         besseli1 = besseli1 * exp(x) / sqrt(2.0_dp*pi*x)  
      endif  

   endfunction besseli1  

   pure function covariancePP(r,q,d)
      real(dp), intent(in) :: r
      integer, intent(in) :: q, d
      real(dp) :: covariancePP

      real(dp) :: j
      integer :: j_int

      j_int = d/2 + q + 1
      j = float(j_int)

      if( r > 1.0_dp ) then
         covariancePP = 0.0_dp
      elseif( r < 0.0_dp ) then
         covariancePP = 1.0_dp
      else
         if( q == 0 ) then
            covariancePP = (1.0_dp - r)**j_int
         elseif( q == 1 ) then
            covariancePP = (1.0_dp - r)**(j_int+1) * ( (j+1.0_dp)*r + 1.0_dp )
         elseif( q == 2) then
            covariancePP = (1.0_dp - r)**(j_int+2) * ( (j**2 + 4.0_dp*j + 3.0_dp) * r**2 + (3.0_dp*j+6.0_dp)*r + 3.0_dp ) / 3.0_dp
         elseif( q == 3) then
            covariancePP = (1.0_dp - r)**(j_int+3) * &
            ( (j**3 + 9.0_dp*j**2 + 23.0_dp*j + 15.0_dp)*r**3 + &
              (6.0_dp * j**2 + 36.0_dp*j + 45.0_dp) * r**2 + &
              (15.0_dp*j+45.0_dp)*r + 15.0_dp ) / 15.0_dp
         endif
      endif

   endfunction covariancePP

   pure function grad_covariancePP(r,q,d)
      real(dp), intent(in) :: r
      integer, intent(in) :: q, d
      real(dp) :: grad_covariancePP

      real(dp) :: j
      integer :: j_int

      j_int = d/2 + q + 1
      j = float(j_int)

      if( r > 1.0_dp ) then
         grad_covariancePP = 0.0_dp
      elseif( r < 0.0_dp ) then
         grad_covariancePP = 0.0_dp
      else
         if( q == 0 ) then
            grad_covariancePP = - j * (1.0_dp - r)**(j_int-1)
         elseif( q == 1 ) then
            grad_covariancePP = - (j+1.0_dp) * (1.0_dp - r)**j_int * ( (j+1.0_dp)*r + 1.0_dp ) + &
            (1.0_dp - r)**(j_int+1) * (j+1.0_dp)
         elseif( q == 2) then
            grad_covariancePP = - (j+2.0_dp) * (1.0_dp - r)**(j_int+1) * ( (j**2 + 4.0_dp*j + 3.0_dp) * r**2 + (3.0_dp*j+6.0_dp)*r + 3.0_dp ) / 3.0_dp + &
            (1.0_dp - r)**(j_int+2) * ( 2.0_dp * (j**2 + 4.0_dp*j + 3.0_dp) * r + (3.0_dp*j+6.0_dp) ) / 3.0_dp
         elseif( q == 3) then
            grad_covariancePP = -(j+3.0_dp) * (1.0_dp - r)**(j_int+2) * &
            ( (j**3 + 9.0_dp*j**2 + 23.0_dp*j + 15.0_dp)*r**3 + &
              (6.0_dp * j**2 + 36.0_dp*j + 45.0_dp) * r**2 + &
              (15.0_dp*j+45.0_dp)*r + 15.0_dp ) / 15.0_dp + &
            (1.0_dp - r)**(j_int+3) * &
            ( 3.0_dp * (j**3 + 9.0_dp*j**2 + 23.0_dp*j + 15.0_dp)*r**2 + &
              2.0_dp * (6.0_dp * j**2 + 36.0_dp*j + 45.0_dp) * r + &
              (15.0_dp*j+45.0_dp) ) / 15.0_dp
         endif
      endif

   endfunction grad_covariancePP

   function gpCovariance_atom_real_space_Calc( this, x_i, x_i_size, x_j, x_j_size, xPrime_i, xPrime_j, xPrime_ij, error )
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in) :: x_i, x_j
      integer, intent(in) :: x_i_size, x_j_size
      real(dp), dimension(:), intent(out), optional, pointer :: xPrime_i, xPrime_j
      real(dp), dimension(:,:), intent(out), optional, pointer :: xPrime_ij
      integer, intent(out), optional :: error

      real(dp) :: gpCovariance_atom_real_space_Calc

      type(neighbour_descriptor), dimension(:), allocatable :: neighbour_i, neighbour_j, grad_spherical_i, grad_spherical_j
      type(neighbour_descriptor), dimension(:,:), allocatable :: grad_spherical_i_radial_j, grad_spherical_j_radial_i
      
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lm, &
      mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lp, mo_spher_bess_fi_ki_lpp, &
      grad_mo_spher_bess_fi_ki_l, grad2_mo_spher_bess_fi_ki_l, &
      grad_arg_bess1, grad_arg_bess2, grad_fac_exp1, grad_fac_exp2, radial, grad_radial_i, grad_radial_j, grad2_radial_ij, &
      fcut1, fcut2, dfcut1, dfcut2, fac_r1r2, grad_fac_r1r2_1, grad_fac_r1r2_2, grad2_fac_exp, grad2_fac_r1r2

      real(dp), dimension(1) :: real_mould

      integer :: i, j, i_data, j_data, n1, n2, l, l1, l2, m1, m2, n_neighbour_i, n_neighbour_j, real_mould_size

      logical :: do_derivative, do_xPrime_i, do_xPrime_j, do_xPrime_ij

      complex(dp) :: I_lm1m2, tmp_complex
      type(cplx_2d_array), dimension(:), allocatable :: integral_r

      type grad_r_type
         type(cplx_2d_array), dimension(:), allocatable :: integral_r
      endtype grad_r_type

      type real_1d_array
         real(dp), dimension(:), allocatable :: value
      endtype real_1d_array

      type(grad_r_type), dimension(:), allocatable :: grad_ri, grad_rj
      type(grad_r_type), dimension(:,:), allocatable :: grad_rij

      type(real_1d_array), dimension(:,:), allocatable :: grad_spherical_ij

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCovariance_atom_real_space_Calc: object not initialised', error)
      endif

      do_xPrime_i = .false.
      do_xPrime_j = .false.
      do_xPrime_ij = .false.
      if(present(xPrime_i)) do_xPrime_i = associated(xPrime_i)
      if(present(xPrime_j)) do_xPrime_j = associated(xPrime_j)
      if(present(xPrime_ij)) do_xPrime_ij = associated(xPrime_ij)

      do_derivative = (do_xPrime_i .or. do_xPrime_j .or. do_xPrime_ij)

      call gpRealArray_NeighbourDescriptor(this,x_i,x_i_size,neighbour_i,n_neighbour_i)
      call gpRealArray_NeighbourDescriptor(this,x_j,x_j_size,neighbour_j,n_neighbour_j)

      if(do_xPrime_i .or. do_xPrime_ij) then
         allocate(grad_spherical_i(n_neighbour_i))
         allocate(grad_ri(n_neighbour_i))

         do i = 1, n_neighbour_i
            

            allocate(grad_ri(i)%integral_r(0:this%l_max))
            allocate(grad_spherical_i(i)%spherical_harmonics(0:this%l_max))

            do l = 0, this%l_max
               allocate(grad_spherical_i(i)%spherical_harmonics(l)%value(-l:l))
               grad_spherical_i(i)%spherical_harmonics(l)%value = CPLX_ZERO

               allocate(grad_ri(i)%integral_r(l)%value(-l:l,-l:l))
               grad_ri(i)%integral_r(l)%value = CPLX_ZERO
            enddo
         enddo
      endif

      if(do_xPrime_j .or. do_xPrime_ij) then
         allocate(grad_spherical_j(n_neighbour_j))
         allocate(grad_rj(n_neighbour_j))

         do i = 1, n_neighbour_j
            
            allocate(grad_rj(i)%integral_r(0:this%l_max))
            allocate(grad_spherical_j(i)%spherical_harmonics(0:this%l_max))

            do l = 0, this%l_max
               allocate(grad_spherical_j(i)%spherical_harmonics(l)%value(-l:l))
               grad_spherical_j(i)%spherical_harmonics(l)%value = CPLX_ZERO

               allocate(grad_rj(i)%integral_r(l)%value(-l:l,-l:l))
               grad_rj(i)%integral_r(l)%value = CPLX_ZERO
            enddo
         enddo
      endif

      if(do_xPrime_ij) then
         allocate(grad_rij(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_ij(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_i_radial_j(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_j_radial_i(n_neighbour_j,n_neighbour_i))

         do i = 1, n_neighbour_i
            do j = 1, n_neighbour_j

               allocate(grad_spherical_ij(j,i)%value(0:this%l_max))
               grad_spherical_ij(j,i)%value = 0.0_dp
            
               allocate(grad_rij(j,i)%integral_r(0:this%l_max))

               allocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(0:this%l_max))
               allocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(0:this%l_max))

               do l = 0, this%l_max
                  allocate(grad_rij(j,i)%integral_r(l)%value(-l:l,-l:l))
                  grad_rij(j,i)%integral_r(l)%value = CPLX_ZERO

                  allocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value(-l:l))
                  grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value = CPLX_ZERO

                  allocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value(-l:l))
                  grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value = CPLX_ZERO
               enddo
            enddo
         enddo

      endif

      allocate(integral_r(0:this%l_max))

      do l = 0, this%l_max
         allocate(integral_r(l)%value(-l:l,-l:l))
         integral_r(l)%value = CPLX_ZERO
      enddo

      if(do_xPrime_i) xPrime_i = 0.0_dp
      if(do_xPrime_j) xPrime_j = 0.0_dp
      if(do_xPrime_ij) xPrime_ij = 0.0_dp

      ! Overlap of central atoms
      integral_r(0)%value(0,0) = 0.25_dp/PI !/ this%atom_sigma**1.5_dp

      ! Overlaps of central atom of first environment with the other atoms in the second.
      do n1 = 1, n_neighbour_i
         r1 = neighbour_i(n1)%r
         if(r1 > this%cutoff) cycle
         fcut1 = coordination_function(r1,this%cutoff,this%cutoff_transition_width)
         fac_exp = exp(-0.5_dp*this%atom_sigma*r1**2) * 0.25_dp/PI
         fac_r1r2 = fac_exp * fcut1

         integral_r(0)%value(0,0) = integral_r(0)%value(0,0) + fac_r1r2

         if(do_xPrime_i) then
            dfcut1 = dcoordination_function(r1,this%cutoff,this%cutoff_transition_width)
            grad_fac_exp1 = -fac_exp*this%atom_sigma*r1
            grad_fac_r1r2_1 = fac_exp * dfcut1 + grad_fac_exp1 * fcut1
            grad_ri(n1)%integral_r(0)%value(0,0) = grad_ri(n1)%integral_r(0)%value(0,0) + grad_fac_r1r2_1
         endif
      enddo

      ! Overlaps of central atom of second environment with the other atoms in the first.
      do n2 = 1, n_neighbour_j
         r2 = neighbour_j(n2)%r
         if(r2 > this%cutoff) cycle
         fcut2 = coordination_function(r2,this%cutoff,this%cutoff_transition_width)
         fac_exp = exp(-0.5_dp*this%atom_sigma*r2**2) * 0.25_dp/PI
         fac_r1r2 = fac_exp * fcut2

         integral_r(0)%value(0,0) = integral_r(0)%value(0,0) + fac_r1r2

         if(do_xPrime_j) then
            dfcut2 = dcoordination_function(r2,this%cutoff,this%cutoff_transition_width)
            grad_fac_exp2 = -fac_exp*this%atom_sigma*r2
            grad_fac_r1r2_2 = fac_exp * dfcut2 + grad_fac_exp2 * fcut2
            grad_rj(n2)%integral_r(0)%value(0,0) = grad_rj(n2)%integral_r(0)%value(0,0) + grad_fac_r1r2_2
         endif
      enddo

      ! Overlaps of non-central atoms.
      do n1 = 1, n_neighbour_i
         r1 = neighbour_i(n1)%r

         if(r1 > this%cutoff) cycle
         fcut1 = coordination_function(r1,this%cutoff,this%cutoff_transition_width)
         dfcut1 = dcoordination_function(r1,this%cutoff,this%cutoff_transition_width)
         do n2 = 1, n_neighbour_j
            r2 = neighbour_j(n2)%r

            if(r2 > this%cutoff) cycle
            fcut2 = coordination_function(r2,this%cutoff,this%cutoff_transition_width)
            dfcut2 = dcoordination_function(r2,this%cutoff,this%cutoff_transition_width)

            arg_bess = this%atom_sigma*r1*r2
            fac_exp = exp(-0.5_dp*this%atom_sigma*(r1**2+r2**2))
            fac_r1r2 = fac_exp * fcut1 * fcut2

            if(do_xPrime_i .or. do_xPrime_ij) then
               grad_arg_bess1 = this%atom_sigma*r2
               grad_fac_exp1 = -fac_exp*this%atom_sigma*r1
               grad_fac_r1r2_1 = (fac_exp * dfcut1 + grad_fac_exp1 * fcut1) * fcut2
            endif

            if(do_xPrime_j .or. do_xPrime_ij) then
               grad_arg_bess2 = this%atom_sigma*r1
               grad_fac_exp2 = -fac_exp*this%atom_sigma*r2
               grad_fac_r1r2_2 = (fac_exp * dfcut2 + grad_fac_exp2 * fcut2) * fcut1
            endif

            if(do_xPrime_ij) then
               grad2_fac_exp = fac_exp * this%atom_sigma**2 * r1*r2
               grad2_fac_r1r2 = grad2_fac_exp*fcut1*fcut2 + grad_fac_exp1*fcut1*dfcut2 + grad_fac_exp2*dfcut1*fcut2 + fac_exp*dfcut1*dfcut2
            endif

            do l = 0, this%l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                  if(do_derivative) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                  if(do_xPrime_ij) mo_spher_bess_fi_ki_lpp = mo_spher_bess_fi_ki_l - (2*l+3)*mo_spher_bess_fi_ki_lp / arg_bess
               else
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  if(do_derivative) then
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     if(do_xPrime_ij) mo_spher_bess_fi_ki_lpp = mo_spher_bess_fi_ki_l - (2*l+3)*mo_spher_bess_fi_ki_lp / arg_bess
                  else
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                  endif
               endif

               if(do_derivative) grad_mo_spher_bess_fi_ki_l = l * mo_spher_bess_fi_ki_l / arg_bess + mo_spher_bess_fi_ki_lp
               if(do_xPrime_ij) grad2_mo_spher_bess_fi_ki_l = ( l*(2*l+3)*(l-1) + (1.0_dp+2*l)*arg_bess**2 ) * &
               mo_spher_bess_fi_ki_lp / arg_bess**3 + &
               ( 1.0_dp + l*(l-1)/arg_bess**2 ) * mo_spher_bess_fi_ki_lpp

               !radial = mo_spher_bess_fi_ki_l*fac_exp
               radial = mo_spher_bess_fi_ki_l*fac_r1r2

               !if(do_xPrime_i .or. do_xPrime_ij) grad_radial_i = grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp1
               if(do_xPrime_i .or. do_xPrime_ij) grad_radial_i = grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_r1r2 + mo_spher_bess_fi_ki_l * grad_fac_r1r2_1

               !if(do_xPrime_j .or. do_xPrime_ij) grad_radial_j = grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp2
               if(do_xPrime_j .or. do_xPrime_ij) grad_radial_j = grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_r1r2 + mo_spher_bess_fi_ki_l * grad_fac_r1r2_2

               if(do_xPrime_ij) then
                  !grad2_radial_ij = fac_exp * this%atom_sigma**2 * r1 * r2 * mo_spher_bess_fi_ki_l + &
                  !grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * grad_fac_exp2 + &
                  !grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * grad_fac_exp1 + &
                  !fac_exp * ( this%atom_sigma * grad_mo_spher_bess_fi_ki_l + grad_arg_bess1*grad_arg_bess2*grad2_mo_spher_bess_fi_ki_l )
                  grad2_radial_ij = grad2_fac_r1r2 * mo_spher_bess_fi_ki_l + &
                  grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * grad_fac_r1r2_2 + &
                  grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * grad_fac_r1r2_1 + &
                  fac_r1r2 * ( this%atom_sigma * grad_mo_spher_bess_fi_ki_l + grad_arg_bess1*grad_arg_bess2*grad2_mo_spher_bess_fi_ki_l )

                  grad_spherical_ij(n2,n1)%value(l) = grad_spherical_ij(n2,n1)%value(l) + &
                  radial
               endif
                  
               do m1 = -l, l
                  if(do_xPrime_i .or. do_xPrime_ij) then
                     grad_spherical_i(n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_i(n1)%spherical_harmonics(l)%value(m1) + &
                     radial * neighbour_j(n2)%spherical_harmonics(l)%value(m1)
                  endif

                  if(do_xPrime_j .or. do_xPrime_ij) then
                     grad_spherical_j(n2)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_j(n2)%spherical_harmonics(l)%value(m1) + &
                     radial * conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1))
                  endif

                  if(do_xPrime_ij) then
                     grad_spherical_i_radial_j(n2,n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_i_radial_j(n2,n1)%spherical_harmonics(l)%value(m1) + &
                     grad_radial_j * neighbour_j(n2)%spherical_harmonics(l)%value(m1)

                     grad_spherical_j_radial_i(n2,n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_j_radial_i(n2,n1)%spherical_harmonics(l)%value(m1) + &
                     grad_radial_i * conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1))
                  endif

                  do m2 = -l, l
                     I_lm1m2 =  radial * &
                     conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                     neighbour_j(n2)%spherical_harmonics(l)%value(m2)

                     integral_r(l)%value(m2,m1) = integral_r(l)%value(m2,m1) + I_lm1m2

                     if(do_xPrime_i .or. do_xPrime_ij) then
                        grad_ri(n1)%integral_r(l)%value(m2,m1) = grad_ri(n1)%integral_r(l)%value(m2,m1) + &
                        grad_radial_i * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                     if(do_xPrime_j .or. do_xPrime_ij) then
                        grad_rj(n2)%integral_r(l)%value(m2,m1) = grad_rj(n2)%integral_r(l)%value(m2,m1) + &
                        grad_radial_j * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                     if(do_xPrime_ij) then
                        grad_rij(n2,n1)%integral_r(l)%value(m2,m1) = grad_rij(n2,n1)%integral_r(l)%value(m2,m1) + &
                        grad2_radial_ij * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                  enddo
               enddo
            enddo
         enddo
      enddo

      gpCovariance_atom_real_space_Calc = 0.0_dp

      do l = 0, this%l_max
         gpCovariance_atom_real_space_Calc = gpCovariance_atom_real_space_Calc + sum(real(integral_r(l)%value)**2) + sum(aimag(integral_r(l)%value)**2)
      enddo

      if(do_xPrime_i) then
         i_data = 0
         do i = 1, n_neighbour_i
               
            i_data = i_data + 1
            do l = 0, this%l_max
               xPrime_i(i_data) = xPrime_i(i_data) + &
               sum(real(integral_r(l)%value)*real(grad_ri(i)%integral_r(l)%value)) + &
               sum(aimag(integral_r(l)%value)*aimag(grad_ri(i)%integral_r(l)%value)) 
            enddo
            i_data = i_data + 1
            xPrime_i(i_data) = 0.0_dp

            do l = 0, this%l_max
               real_mould_size = size(transfer(grad_spherical_i(i)%spherical_harmonics(l)%value(-l:l),real_mould))
               xPrime_i(i_data+1:i_data+real_mould_size) = transfer(matmul(grad_spherical_i(i)%spherical_harmonics(l)%value,conjg(integral_r(l)%value)),real_mould) 
               i_data = i_data + real_mould_size
            enddo
         enddo

          xPrime_i = xPrime_i * 2.0_dp
      endif

      if(do_xPrime_j) then
         i_data = 0
         do i = 1, n_neighbour_j
               
            i_data = i_data + 1
            do l = 0, this%l_max
               xPrime_j(i_data) = xPrime_j(i_data) + &
               sum(real(integral_r(l)%value)*real(grad_rj(i)%integral_r(l)%value)) + &
               sum(aimag(integral_r(l)%value)*aimag(grad_rj(i)%integral_r(l)%value)) 
            enddo
            i_data = i_data + 1
            xPrime_j(i_data) = 0.0_dp

            do l = 0, this%l_max
               real_mould_size = size(transfer(grad_spherical_j(i)%spherical_harmonics(l)%value(-l:l),real_mould))
               xPrime_j(i_data+1:i_data+real_mould_size) = transfer(matmul(integral_r(l)%value,conjg(grad_spherical_j(i)%spherical_harmonics(l)%value)),real_mould) 
               i_data = i_data + real_mould_size
            enddo
         enddo

         xPrime_j = xPrime_j * 2.0_dp
      endif

      if(do_xPrime_ij) then
         i_data = 0
         do i = 1, n_neighbour_i
            i_data = i_data + 1

            ! i-th neighbour, wrt r
            j_data = 0
            do j = 1, n_neighbour_j
               j_data = j_data + 1

               ! d r_i d r_j
               do l = 0, this%l_max
                  xPrime_ij(i_data,j_data) = xPrime_ij(i_data,j_data) + &
                  sum(real(grad_rj(j)%integral_r(l)%value)*real(grad_ri(i)%integral_r(l)%value)) + &
                  sum(aimag(grad_rj(j)%integral_r(l)%value)*aimag(grad_ri(i)%integral_r(l)%value)) + &
                  sum(real(integral_r(l)%value)*real(grad_rij(j,i)%integral_r(l)%value)) + &
                  sum(aimag(integral_r(l)%value)*aimag(grad_rij(j,i)%integral_r(l)%value)) 
               enddo
               j_data = j_data + 1
               xPrime_ij(i_data,j_data) = 0.0_dp

               ! d r_i d Y^{lm}_j
               do l = 0, this%l_max
                  real_mould_size = size(transfer(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value(-l:l),real_mould))
                  xPrime_ij(i_data,j_data+1:j_data+real_mould_size) = transfer( &
                  matmul(integral_r(l)%value, conjg(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value)) + &
                  matmul(grad_ri(i)%integral_r(l)%value, conjg(grad_spherical_j(j)%spherical_harmonics(l)%value)), real_mould)

                  j_data = j_data + real_mould_size
               enddo
            enddo

            i_data = i_data + 1
            xPrime_ij(i_data,:) = 0.0_dp

            do l1 = 0, this%l_max
               do m1 = -l1, l1
                  j_data = 0
                  i_data = i_data + 1

                  ! d Y^{lm}_i d r_j
                  do j = 1, n_neighbour_j
                     j_data = j_data + 1

                     tmp_complex = dot_product(grad_rj(j)%integral_r(l1)%value(-l1:l1,m1), grad_spherical_i(i)%spherical_harmonics(l1)%value(-l1:l1)) + &
                     dot_product(integral_r(l1)%value(-l1:l1,m1),grad_spherical_i_radial_j(j,i)%spherical_harmonics(l1)%value(-l1:l1))
                     xPrime_ij(i_data:i_data+1,j_data) = (/real(tmp_complex),aimag(tmp_complex)/)

                     j_data = j_data + 1
                     xPrime_ij(i_data:i_data+1,j_data) = 0.0_dp

                     do l2 = 0, this%l_max
                        real_mould_size = size(transfer(grad_spherical_j(j)%spherical_harmonics(l2)%value(-l2:l2),real_mould))
                        if(l1 == l2) then
                           xPrime_ij(i_data,j_data+1:j_data+real_mould_size) = transfer( &
                           grad_spherical_i(i)%spherical_harmonics(l1)%value*conjg(grad_spherical_j(j)%spherical_harmonics(l2)%value(m1)) + &
                           grad_spherical_ij(j,i)%value(l1)*integral_r(l1)%value(-l1:l1,m1), real_mould) 

                           xPrime_ij(i_data+1,j_data+1:j_data+real_mould_size) = transfer( &
                           -CPLX_IMAG*grad_spherical_i(i)%spherical_harmonics(l1)%value*conjg(grad_spherical_j(j)%spherical_harmonics(l2)%value(m1)) + &
                           CPLX_IMAG*grad_spherical_ij(j,i)%value(l1)*integral_r(l1)%value(-l1:l1,m1), real_mould) 
                        endif
                        j_data = j_data + real_mould_size
                     enddo

                  enddo
                  i_data = i_data + 1
               enddo
            enddo
         enddo
         xPrime_ij = xPrime_ij * 2.0_dp
      endif


      if(allocated(integral_r)) then
         do l = 0, this%l_max
            if(allocated(integral_r(l)%value)) deallocate(integral_r(l)%value)
         enddo
         deallocate(integral_r)
      endif
      if(allocated(grad_ri)) then
         do i = 1, size(grad_ri)
            if(allocated(grad_ri(i)%integral_r)) then
               do l = 0, this%l_max
                  if(allocated(grad_ri(i)%integral_r(l)%value)) deallocate(grad_ri(i)%integral_r(l)%value)
               enddo
               deallocate(grad_ri(i)%integral_r)
            endif
         enddo
         deallocate(grad_ri)
      endif

      if(allocated(grad_rj)) then
         do i = 1, size(grad_rj)
            if(allocated(grad_rj(i)%integral_r)) then
               do l = 0, this%l_max
                  if(allocated(grad_rj(i)%integral_r(l)%value)) deallocate(grad_rj(i)%integral_r(l)%value)
               enddo
               deallocate(grad_rj(i)%integral_r)
            endif
         enddo
         deallocate(grad_rj)
      endif

      if(do_xPrime_ij) then

         do i = 1, n_neighbour_i
            do j = 1, n_neighbour_j


               do l = 0, this%l_max
                  deallocate(grad_rij(j,i)%integral_r(l)%value)
                  deallocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value)
                  deallocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value)
               enddo

               deallocate(grad_spherical_ij(j,i)%value)
               deallocate(grad_rij(j,i)%integral_r)
               deallocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics)
               deallocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics)
            enddo
         enddo

         deallocate(grad_rij)
         deallocate(grad_spherical_ij)
         deallocate(grad_spherical_i_radial_j)
         deallocate(grad_spherical_j_radial_i)

      endif

      call finalise(neighbour_i)
      call finalise(neighbour_j)
      call finalise(grad_spherical_i)
      call finalise(grad_spherical_j)

   endfunction gpCovariance_atom_real_space_Calc

!   function gpCovariance_soap_Calc( this, x_i, x_j, xPrime_i, xPrime_j, xPrime_ij, error )
!      type(gpCovariance_soap), intent(in) :: this
!      real(dp), dimension(:), intent(in) :: x_i, x_j
!      real(dp), dimension(:), intent(out), optional, pointer :: xPrime_i, xPrime_j
!      real(dp), dimension(:,:), intent(out), optional, pointer :: xPrime_ij
!      integer, intent(out), optional :: error
!
!      real(dp) :: gpCovariance_soap_Calc
!
!      integer :: l, m, m1, m2, a, i
!      logical :: do_xPrime_i, do_xPrime_j, do_xPrime_ij, do_derivative
!
!      type(cplx_1d_array), dimension(:,:), allocatable :: fourier1_so3, fourier2_so3, dcov_dfourier1, dcov_dfourier2
!      type(cplx_2d_array), dimension(:), allocatable :: int_soap
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpCovariance_soap_Calc: object not initialised', error)
!      endif
!
!      do_xPrime_i = .false.
!      do_xPrime_j = .false.
!      do_xPrime_ij = .false.
!      if(present(xPrime_i)) do_xPrime_i = associated(xPrime_i)
!      if(present(xPrime_j)) do_xPrime_j = associated(xPrime_j)
!      if(present(xPrime_ij)) do_xPrime_ij = associated(xPrime_ij)
!
!      do_derivative = (do_xPrime_i .or. do_xPrime_j .or. do_xPrime_ij)
!
!      allocate( fourier1_so3(0:this%l_max,this%n_max), fourier2_so3(0:this%l_max,this%n_max), int_soap(0:this%l_max) )
!
!      if(do_xPrime_i) allocate( dcov_dfourier1(0:this%l_max,this%n_max) )
!      if(do_xPrime_j) allocate( dcov_dfourier2(0:this%l_max,this%n_max) )
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            allocate(fourier1_so3(l,a)%value(-l:l))
!            allocate(fourier2_so3(l,a)%value(-l:l))
!            if(do_xPrime_i) allocate(dcov_dfourier1(l,a)%value(-l:l))
!            if(do_xPrime_j) allocate(dcov_dfourier2(l,a)%value(-l:l))
!         enddo
!      enddo
!
!      do l = 0, this%l_max
!         allocate(int_soap(l)%value(-l:l,-l:l))
!         int_soap(l)%value = CPLX_ZERO
!      enddo
!
!
!      i = 0
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            do m = -l, l
!               fourier1_so3(l,a)%value(m) = cmplx(x_i(i+1), x_i(i+2))
!               fourier2_so3(l,a)%value(m) = cmplx(x_j(i+1), x_j(i+2))
!               i = i + 2
!            enddo
!         enddo
!      enddo
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            do m1 = -l, l
!               do m2 = -l, l
!                  int_soap(l)%value(m2,m1) = int_soap(l)%value(m2,m1) + &
!                     fourier1_so3(l,a)%value(m1) * conjg(fourier2_so3(l,a)%value(m2))
!               enddo
!            enddo
!         enddo
!      enddo
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            if(do_xPrime_i) dcov_dfourier1(l,a)%value = matmul(fourier2_so3(l,a)%value,int_soap(l)%value)
!            if(do_xPrime_j) dcov_dfourier2(l,a)%value = matmul(conjg(int_soap(l)%value),fourier1_so3(l,a)%value)
!         enddo
!      enddo
!
!      gpCovariance_soap_Calc = 0.0_dp
!      do l = 0, this%l_max
!         gpCovariance_soap_Calc = gpCovariance_soap_Calc + sum(real(int_soap(l)%value)**2+aimag(int_soap(l)%value)**2)
!      enddo
!
!      if(do_derivative) then
!         i = 0
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               do m = -l, l
!                  if(do_xPrime_i) then
!                     xPrime_i(i+1) = real(dcov_dfourier1(l,a)%value(m))
!                     xPrime_i(i+2) = aimag(dcov_dfourier1(l,a)%value(m))
!                  endif
!                  if(do_xPrime_j) then
!                     xPrime_j(i+1) = real(dcov_dfourier2(l,a)%value(m))
!                     xPrime_j(i+2) = aimag(dcov_dfourier2(l,a)%value(m))
!                  endif
!                  i = i + 2
!               enddo
!            enddo
!         enddo
!         if(do_xPrime_i) xPrime_i = xPrime_i*2.0_dp
!         if(do_xPrime_j) xPrime_j = xPrime_j*2.0_dp
!      endif
!
!      if(allocated(fourier1_so3)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(fourier1_so3(l,a)%value)
!            enddo
!         enddo
!         deallocate(fourier1_so3)
!      endif
!
!      if(allocated(fourier2_so3)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(fourier2_so3(l,a)%value)
!            enddo
!         enddo
!         deallocate(fourier2_so3)
!      endif
!
!      if(allocated(int_soap)) then
!         do l = 0, this%l_max
!            deallocate(int_soap(l)%value)
!         enddo
!         deallocate(int_soap)
!      endif
!
!      if(allocated(dcov_dfourier1)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(dcov_dfourier1(l,a)%value)
!            enddo
!         enddo
!         deallocate(dcov_dfourier1)
!      endif
!
!      if(allocated(dcov_dfourier2)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(dcov_dfourier2(l,a)%value)
!            enddo
!         enddo
!         deallocate(dcov_dfourier2)
!      endif
!
!   endfunction gpCovariance_soap_Calc

   subroutine gpRealArray_NeighbourDescriptor(this,x,x_size,neighbour,n_neighbour)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in) :: x
      integer, intent(in) :: x_size
      type(neighbour_descriptor), dimension(:), allocatable, intent(out) :: neighbour
      integer, intent(out) :: n_neighbour

      integer :: l, i_data, i, real_mould_size
      real(dp), dimension(1) :: real_mould
      complex(dp), dimension(1) :: complex_mould

      n_neighbour = x_size / ( 2 * (this%l_max+1)**2 + 2 )
      
      call finalise(neighbour)

      allocate(neighbour(n_neighbour))

      i_data = 0
      do i = 1, n_neighbour
            
         i_data = i_data + 1
         neighbour(i)%r = x(i_data)
         i_data = i_data + 1
         neighbour(i)%n = abs(nint(x(i_data)))

         allocate(neighbour(i)%spherical_harmonics(0:this%l_max))
         do l = 0, this%l_max

            allocate(neighbour(i)%spherical_harmonics(l)%value(-l:l))

            real_mould_size = size(transfer(neighbour(i)%spherical_harmonics(l)%value(-l:l),real_mould))
            neighbour(i)%spherical_harmonics(l)%value = transfer(x(i_data+1:i_data+real_mould_size),complex_mould)
            i_data = i_data + real_mould_size
         enddo
      enddo

   endsubroutine gpRealArray_NeighbourDescriptor

   subroutine gpNeighbourDescriptor_Finalise(this)
      type(neighbour_descriptor), dimension(:), allocatable, intent(inout) :: this

      integer :: i, l

      if(allocated(this)) then
         do i = 1, size(this)
            do l = lbound(this(i)%spherical_harmonics,dim=1), ubound(this(i)%spherical_harmonics,dim=1)
               if(allocated(this(i)%spherical_harmonics(l)%value)) deallocate(this(i)%spherical_harmonics(l)%value)
            enddo
            if(allocated(this(i)%spherical_harmonics)) deallocate(this(i)%spherical_harmonics)
         enddo
         deallocate(this)
      endif

   endsubroutine gpNeighbourDescriptor_Finalise

   subroutine gp_atom_real_space_RealArray_XYZ(this,x_array,x_array_size,xyz_array,xyz_array_size)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in), target :: x_array
      integer, intent(in) :: x_array_size
      real(dp), dimension(:), allocatable, intent(out) :: xyz_array
      integer, intent(out) :: xyz_array_size

      integer :: l, i_data, i, n_neighbour, n, n_data, xyz_start, xyz_end
      real(dp), dimension(:), pointer :: Y1_array
      real(dp), pointer :: Re_Y1m1, Im_Y1m1, Re_Y10
      real(dp) :: r, x, y, z

      real(dp), parameter :: xy_factor = 0.5_dp * sqrt(1.5_dp / PI)
      real(dp), parameter :: z_factor = 0.5_dp * sqrt(3.0_dp / PI)

      n_neighbour = x_array_size / ( 2 * (this%l_max+1)**2 + 2 )
      xyz_array_size = n_neighbour*3
      
      if(allocated(xyz_array)) deallocate(xyz_array)

      allocate(xyz_array(xyz_array_size))

      i_data = 0
      do i = 1, n_neighbour
            
         i_data = i_data + 1
         r = x_array(i_data)
         i_data = i_data + 1
         n = abs(nint(x_array(i_data)))

         do l = 0, this%l_max
            n_data = 2*(2*l + 1)

            if(l == 1) then
               Y1_array => x_array(i_data+1:i_data+n_data)
               Re_Y1m1 => Y1_array(1)
               Im_Y1m1 => Y1_array(2)
               Re_Y10 =>  Y1_array(3)
               !Im_Y10 =>  Y1_value(4)
               !Re_Y1p1 => Y1_value(5)
               !Im_Y1p1 => Y1_value(6)

               z = Re_Y10 * r / z_factor
               x = Re_Y1m1 * r / xy_factor
               y = -Im_Y1m1 * r / xy_factor
            endif

            i_data = i_data + n_data
         enddo

         xyz_start = (i-1)*3+1
         xyz_end = 3*i
         xyz_array(xyz_start:xyz_end) = (/x,y,z/)
      enddo

   endsubroutine gp_atom_real_space_RealArray_XYZ

   subroutine gp_atom_real_space_XYZ_RealArray(this,xyz_array,xyz_array_size,x_array,x_array_size)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in), target :: xyz_array
      integer, intent(in) :: xyz_array_size
      real(dp), dimension(:), allocatable, intent(out) :: x_array
      integer, intent(out) :: x_array_size

      integer :: l, m, i_data, i, n_neighbour, xyz_start, xyz_end
      real(dp), dimension(:), pointer :: xyz
      complex(dp) :: Y_lm

      n_neighbour = xyz_array_size / 3
      x_array_size = n_neighbour * ( 2 * (this%l_max+1)**2 + 2 )
      
      if(allocated(x_array)) deallocate(x_array)

      allocate(x_array(x_array_size))

      i_data = 0
      do i = 1, n_neighbour
            
         xyz_start = (i-1)*3+1
         xyz_end = 3*i
         xyz => xyz_array(xyz_start:xyz_end)

         i_data = i_data + 1
         x_array(i_data) = norm(xyz)
         i_data = i_data + 1
         x_array(i_data) = real(i,dp)

         do l = 0, this%l_max
            do m = -l, l
               Y_lm = SphericalYCartesian(l,m,xyz)
               x_array(i_data+1:i_data+2) = (/real(Y_lm),aimag(Y_lm)/)
               i_data = i_data + 2
            enddo
         enddo

      enddo

   endsubroutine gp_atom_real_space_XYZ_RealArray

   function fast_pow_1d(v, e)
      real(dp), intent(in) :: v(:), e
      real(dp) :: fast_pow_1d(size(v))
      integer :: e_int

      if (e .feq. 0) then
         fast_pow_1d = 1.0_dp
      elseif (e .feq. 1) then
         fast_pow_1d = v
      elseif (e .feq. 2) then
         fast_pow_1d = v*v
      elseif (e .feq. 3) then
         fast_pow_1d = v*v*v
      elseif (e .feq. 4) then
         fast_pow_1d = v*v
         fast_pow_1d = fast_pow_1d*fast_pow_1d
      else
         e_int = nint(e) 
         if (e .feq. e_int) then
            fast_pow_1d = v**e_int
         else
            fast_pow_1d = v**e
         endif
      endif
   end function fast_pow_1d

   function fast_pow_2d(v, e)
       real(dp), intent(in) :: v(:,:), e
       real(dp) :: fast_pow_2d(size(v,1),size(v,2))

      if (e == 0.0) then
        fast_pow_2d = 1.0_dp
      elseif (e == 1.0) then
        fast_pow_2d = v
      elseif (e == 2.0) then
        fast_pow_2d = v*v
      elseif (e == 3.0) then
        fast_pow_2d = v*v*v
      elseif (e == 4.0) then
        fast_pow_2d = v*v
        fast_pow_2d = fast_pow_2d*fast_pow_2d
      else
        fast_pow_2d = v**e
      endif
   end function fast_pow_2d

   function gpCoordinates_Predict( this, xStar, gradPredict, variance_estimate, do_variance_estimate, grad_variance_estimate, error )  
      type(gpCoordinates), intent(inout), target :: this
      real(dp), dimension(:), intent(in) :: xStar
      real(dp), dimension(:), intent(out), optional :: gradPredict
      real(dp), intent(out), optional :: variance_estimate  
      logical, intent(in), optional :: do_variance_estimate  
      real(dp), dimension(:), intent(out), optional :: grad_variance_estimate  
      integer, optional, intent(out) :: error

      real(dp) :: gpCoordinates_Predict

      real(dp) :: covarianceExp, gpCoordinates_Covariance_ii, gpCoordinates_Covariance_jj, covarianceExp_ii, covarianceExp_jj, normalisation, r_ij, r_jj, covariancePP_ij, covariancePP_jj
      real(dp), pointer :: fc_i
      real(dp), dimension(:), pointer :: x_i
      real(dp), dimension(:,:), pointer :: x_i_permuted_theta  
      real(dp), dimension(this%d) :: xI_xJ_theta, xStar_theta  
      !real(dp), dimension(this%d,this%n_permutations) :: xStar_permuted  
      real(dp), dimension(this%n_sparseX) :: k
      integer :: i_sparseX, i_p, ii, jj
      real(dp) :: delta, covariance_x_x, diag_covariance
      real(dp), dimension(:), allocatable :: covariance_x_xStars, alpha_scaled, grad_Covariance_jj  
      real(dp), dimension(:), pointer :: xPrime_i
      real(dp), dimension(:), allocatable, target :: grad_kStar, k_mm_k
      real(dp), dimension(:,:), allocatable, target :: grad_k
      real(dp), dimension(:,:), allocatable :: distance_matrix
      logical :: my_do_variance_estimate  

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_Predict: object not initialised', error)
      endif

      if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
#ifdef _OPENMP  
         if (OMP_IN_PARALLEL()) then  
            RAISE_ERROR('gpCoordinates_Predict: bond_real_space covariance not OpenMP/thread-safe', error)  
         end if  
#endif  
         if(.not. this%bond_real_space_cov%initialised) then  
            call gpCoordinates_gpCovariance_bond_real_space_Initialise(this)
         endif
      endif

      k = 0.0_dp
      xPrime_i => null()
      if(present(gradPredict) .and. this%covariance_type /= COVARIANCE_DOT_PRODUCT) then  
         allocate(grad_k(size(xStar),this%n_sparseX))
         grad_k = 0.0_dp
         if( this%n_permutations > 1 ) allocate(grad_Covariance_jj(this%d))
      endif

      covariance_ard_se_calc_cov_jj: if (this%covariance_type == COVARIANCE_ARD_SE) then  
         if (.not. this%sparse_covariance_initialised) then  
            RAISE_ERROR('gpCoordinates_Predict: gpCoordinates_precalculate_sparse needs to be called first', error)  
         end if  

         xStar_theta = xStar / this%theta  
  
         !do i_p = 1, this%n_permutations  
         !   xStar_permuted(:,i_p) = xStar(this%permutations(:,i_p))  
         !end do  
         if( this%n_permutations > 1 ) then  
            gpCoordinates_Covariance_jj = 0.0_dp  
            if(present(gradPredict)) grad_Covariance_jj = 0.0_dp  
  
            do i_p = 1, this%n_permutations  
               xI_xJ_theta = ( xStar_theta(this%permutations(:,i_p)) - xStar_theta(:) )  
               covarianceExp_jj = exp( -0.5_dp * dot_product(xI_xJ_theta,xI_xJ_theta) )  
               gpCoordinates_Covariance_jj = gpCoordinates_Covariance_jj + covarianceExp_jj  
  
               if(present(gradPredict)) then  
                  !grad_Covariance_jj = grad_Covariance_jj + covarianceExp_jj * xI_xJ_theta / this%theta
                  grad_Covariance_jj = grad_Covariance_jj + covarianceExp_jj * xI_xJ_theta / this%theta  
                  grad_Covariance_jj(this%permutations(:,i_p)) = grad_Covariance_jj(this%permutations(:,i_p)) - covarianceExp_jj * xI_xJ_theta / this%theta
               endif  
            enddo  
         endif  
      end if covariance_ard_se_calc_cov_jj  
  
      covariance_pp_calc_cov_jj: if (this%covariance_type == COVARIANCE_PP) then
         if (.not. this%sparse_covariance_initialised) then
            RAISE_ERROR('gpCoordinates_Predict: gpCoordinates_precalculate_sparse needs to be called first', error)
         end if

         xStar_theta = xStar / this%theta
         allocate(distance_matrix(this%d, this%d))
         forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) &
            distance_matrix(ii,jj) = ( xStar_theta(ii) - xStar_theta(jj) )**2

         if( this%n_permutations > 1 ) then
            gpCoordinates_Covariance_jj = 0.0_dp
            if(present(gradPredict)) grad_Covariance_jj = 0.0_dp

            do i_p = 1, this%n_permutations
               if( any( (/ (distance_matrix(this%permutations(ii,i_p),ii) > 1.0_dp, ii=1, this%d) /) ) ) cycle
            
               r_jj = sqrt( sum( (/ (distance_matrix(this%permutations(ii,i_p),ii), ii=1, this%d) /) ) )
               if( r_jj >= 1.0_dp ) cycle

               covariancePP_jj = covariancePP(r_jj,PP_Q, this%d)
               gpCoordinates_Covariance_jj = gpCoordinates_Covariance_jj + covariancePP_jj

               if(present(gradPredict) .and. ( r_jj .fne. 0.0_dp ) ) then
                  xI_xJ_theta = ( xStar_theta(:) - xStar_theta(this%permutations(:,i_p)) )

                  grad_Covariance_jj = grad_Covariance_jj + grad_covariancePP(r_jj,PP_Q, this%d) * xI_xJ_theta / this%theta / r_jj
                  grad_Covariance_jj(this%permutations(:,i_p)) = grad_Covariance_jj(this%permutations(:,i_p)) - grad_covariancePP(r_jj,PP_Q, this%d) * xI_xJ_theta / this%theta / r_jj
               endif
            enddo
         endif
      end if covariance_pp_calc_cov_jj

      covariance_type_calc_k: if (this%covariance_type == COVARIANCE_DOT_PRODUCT) then  
         allocate(covariance_x_xStars(this%n_sparseX))  
         call dgemv('T', size(this%sparseX,1), size(this%sparseX,2), 1.0_dp, this%sparseX(1,1), size(this%sparseX, 1), &  
                    xStar(1), 1, 0.0_dp, covariance_x_xStars(1), 1)  
! now a single dgemv call outside the loop   
!         do i_sparseX = 1, this%n_sparseX  
!            covariance_x_xStar = dot_product(xStar,this%sparseX(:,i_sparseX))  
!  
!            k(i_sparseX) = this%delta**2 * covariance_x_xStar**this%theta(1)  
!  
!            if(present(gradPredict)) grad_k(:,i_sparseX) = this%delta**2 * this%theta(1) * covariance_x_xStar**(this%theta(1)-1.0_dp) * this%sparseX(:,i_sparseX)  
!         end do  
         
         k(:) = this%delta**2 * fast_pow_1d(covariance_x_xStars(:), this%zeta)

         k = k * this%sparseCutoff
         if(present(gradPredict)) then
            allocate(alpha_scaled(size(this%alpha)))
            alpha_scaled(:) = this%alpha(:) * this%delta**2 * this%zeta * fast_pow_1d(covariance_x_xStars, this%zeta-1.0_dp)
            alpha_scaled = alpha_scaled * this%sparseCutoff
         endif
         deallocate(covariance_x_xStars)

      else if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then covariance_type_calc_k  
         xPrime_i => null()  
         do i_sparseX = 1, this%n_sparseX  
            delta = this%bond_real_space_cov%delta  
            this%bond_real_space_cov%delta = 1.0_dp  
            covariance_x_x = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=xStar, x_i_size=(size(xStar) - 1), x_j=xStar, x_j_size=(size(xStar) - 1))  
            this%bond_real_space_cov%delta = delta  
            k(i_sparseX) = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=xStar, x_i_size=(size(xStar) - 1), x_j=this%sparseX(:,i_sparseX), x_j_size=this%sparseX_size(i_sparseX)) &  
                           / sqrt(covariance_x_x * this%covarianceDiag_sparseX_sparseX(i_sparseX))  
         enddo  

      else if(this%covariance_type == COVARIANCE_ARD_SE) then covariance_type_calc_k  
         xPrime_i => null()  
         do i_sparseX = 1, this%n_sparseX  
            !x_i => this%sparseX(:,i_sparseX)  
            x_i_permuted_theta => this%sparseX_permuted(:,:,i_sparseX)  
            fc_i => this%sparseCutoff(i_sparseX)  
            do i_p = 1, this%n_permutations  

               xI_xJ_theta = (x_i_permuted_theta(:,i_p) - xStar_theta(:))  
               !xI_xJ_theta = (x_i(:) - xStar_permuted(:,i_p)) / this%theta  
               !xI_xJ_theta = (this%sparseX_permuted(:,i_p,i_sparseX) - xStar(:)) / this%theta  

               covarianceExp = this%delta**2 * exp( -0.5_dp * dot_product(xI_xJ_theta,xI_xJ_theta) )  

               if(present(gradPredict)) grad_k(:,i_sparseX) = grad_k(:,i_sparseX) + covarianceExp*xI_xJ_theta / this%theta  
               k(i_sparseX) = k(i_sparseX) + covarianceExp  
            enddo  

            if( this%n_permutations > 1 ) then

               normalisation = sqrt(this%sparseCovariance(i_sparseX) * gpCoordinates_Covariance_jj)
               if(present(gradPredict)) then
                  grad_k(:,i_sparseX) = grad_k(:,i_sparseX) / normalisation - 0.5_dp * grad_Covariance_jj * k(i_sparseX) / normalisation / gpCoordinates_Covariance_jj
               endif

               k(i_sparseX) = k(i_sparseX) / normalisation

            endif
            k(i_sparseX) = ( k(i_sparseX) + this%f0**2 ) * fc_i
            if(present(gradPredict)) grad_k(:,i_sparseX) = grad_k(:,i_sparseX) * fc_i
         enddo
      else if(this%covariance_type == COVARIANCE_PP) then covariance_type_calc_k
         xPrime_i => null()
         do i_sparseX = 1, this%n_sparseX
            x_i => this%sparseX(:,i_sparseX)
            fc_i => this%sparseCutoff(i_sparseX)

            forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) distance_matrix(ii,jj) = ( x_i(ii) - xStar(jj) )**2 / this%theta(ii)**2

            do i_p = 1, this%n_permutations
               if( any( (/ (distance_matrix(this%permutations(ii,i_p),ii) > 1.0_dp, ii=1, this%d) /) ) ) cycle

               r_ij = sqrt( sum( (/ (distance_matrix(this%permutations(ii,i_p),ii), ii=1, this%d) /) ) )
               if( r_ij >= 1.0_dp ) cycle

               covariancePP_ij = this%delta**2 * covariancePP(r_ij,PP_Q, this%d)
               if(present(gradPredict) .and. ( r_ij /= 0.0_dp ) ) grad_k(:,i_sparseX) = grad_k(:,i_sparseX) + &
                  this%delta**2 * grad_covariancePP(r_ij,PP_Q, this%d) * ( xStar(:) - x_i(this%permutations(:,i_p)) ) / r_ij / this%theta(:)**2

               k(i_sparseX) = k(i_sparseX) + covariancePP_ij
            enddo

            if( this%n_permutations > 1 ) then  

               normalisation = sqrt(this%sparseCovariance(i_sparseX) * gpCoordinates_Covariance_jj)  

               if(present(gradPredict)) then  
                  grad_k(:,i_sparseX) = grad_k(:,i_sparseX) / normalisation - 0.5_dp * grad_Covariance_jj * k(i_sparseX) / normalisation / gpCoordinates_Covariance_jj
               endif
  
               k(i_sparseX) = k(i_sparseX) / normalisation  
  
            endif  
            k(i_sparseX) = ( k(i_sparseX) + this%f0**2 ) * fc_i  
            if(present(gradPredict)) grad_k(:,i_sparseX) = grad_k(:,i_sparseX) * fc_i  
         enddo  
      end if covariance_type_calc_k  
      gpCoordinates_Predict = dot_product( k, this%alpha )

      if (this%covariance_type == COVARIANCE_DOT_PRODUCT) then
         if(present(gradPredict)) &  
            call dgemv('N', size(this%sparseX,1), size(this%sparseX,2), 1.0_dp, this%sparseX(1,1), size(this%sparseX,1), &  
               alpha_scaled(1), 1, 0.0_dp, gradPredict(1), 1)  
      else
         if(present(gradPredict)) &  
            call dgemv('N', size(grad_k,1), size(grad_k,2), 1.0_dp, grad_k(1,1), size(grad_k,1), &  
               this%alpha(1), 1, 0.0_dp, gradPredict(1), 1)  
      endif
      my_do_variance_estimate = present(variance_estimate) .and. optional_default(.false.,do_variance_estimate)  

      if(my_do_variance_estimate) then  
         allocate(k_mm_k(this%n_sparseX))

         if(.not.this%variance_estimate_initialised) then  
            RAISE_ERROR('gpCoordinates_Predict: variance_estimate not initialised',error)  
         endif

         call Matrix_Solve(this%LA_k_mm, k, k_mm_k)  
         diag_covariance = this%delta**2 + this%f0**2 + this%variance_estimate_regularisation**2  
  
         variance_estimate = diag_covariance - dot_product(k,k_mm_k)  
         if( variance_estimate < 0.0_dp ) then  
            RAISE_ERROR('gpCoordinates_Predict: variance_estimate: negative variance predicted: '//variance_estimate ,error)  
         endif

         if( present(gradPredict) .and. present(grad_variance_estimate) ) then  
            if (this%covariance_type == COVARIANCE_DOT_PRODUCT) then
               grad_variance_estimate = - 2.0_dp * matmul(this%sparseX, alpha_scaled / this%alpha * k_mm_k)  
            else
               call dgemv('N', size(grad_k,1), size(grad_k,2), 1.0_dp, grad_k(1,1), size(grad_k,1), &
                  k_mm_k(1), 1, 0.0_dp, grad_variance_estimate(1), 1)  
               grad_variance_estimate = - 2.0_dp * grad_variance_estimate  
            endif
         endif
         if(allocated(k_mm_k)) deallocate(k_mm_k)
      endif

      if(allocated(alpha_scaled)) deallocate(alpha_scaled)
      if(allocated(grad_k)) deallocate(grad_k)
      if(allocated(grad_kStar)) deallocate(grad_kStar)
      if(allocated(grad_Covariance_jj)) deallocate(grad_Covariance_jj)
      if( allocated( distance_matrix ) ) deallocate(distance_matrix)

   endfunction gpCoordinates_Predict

   subroutine gpCoordinates_precalculate_sparse(this)  
      type(gpCoordinates), intent(inout), target :: this  
  
      integer :: i_sparseX, i_p, ii, jj
      real(dp), dimension(:), pointer :: x_i  
      real(dp), dimension(:,:), pointer :: x_i_permuted_theta  
      real(dp) :: gpCoordinates_Covariance_ii, covarianceExp_ii, r_ii
      real(dp), dimension(this%d) :: xI_xI_theta  
      real(dp), dimension(:,:), allocatable :: distance_matrix
  
      initialise_sparse_covariance: if( .not. this%sparse_covariance_initialised ) then  
         select case(this%covariance_type)
         case(COVARIANCE_ARD_SE)
            if (allocated(this%sparseX_permuted)) deallocate( this%sparseX_permuted )
            allocate(this%sparseX_permuted(this%d, this%n_permutations, this%n_sparseX))
  
            do i_sparseX = 1, this%n_sparseX  
               x_i => this%sparseX(:,i_sparseX)  
               x_i_permuted_theta => this%sparseX_permuted(:,:,i_sparseX)  
               do i_p = 1, this%n_permutations  
                  x_i_permuted_theta(:,i_p) = x_i(this%permutations(:,i_p)) / this%theta
               end do
            end do  

            if( this%n_permutations > 1 ) then
               call reallocate(this%sparseCovariance,this%n_sparseX)

               do i_sparseX = 1, this%n_sparseX
                  x_i => this%sparseX(:,i_sparseX)
                  x_i_permuted_theta => this%sparseX_permuted(:,:,i_sparseX)
                  gpCoordinates_Covariance_ii = 0.0_dp
                  do i_p = 1, this%n_permutations
                     xI_xI_theta = x_i_permuted_theta(:,i_p) - (x_i / this%theta)
                     covarianceExp_ii = exp( -0.5_dp * dot_product(xI_xI_theta,xI_xI_theta) )
                     gpCoordinates_Covariance_ii = gpCoordinates_Covariance_ii + covarianceExp_ii
                  enddo
                  this%sparseCovariance(i_sparseX) = gpCoordinates_Covariance_ii
               end do
            end if
         case(COVARIANCE_PP)
            if( this%n_permutations > 1 ) then
               call reallocate(this%sparseCovariance,this%n_sparseX)
               allocate(distance_matrix(this%d,this%d))

               do i_sparseX = 1, this%n_sparseX
                  x_i => this%sparseX(:,i_sparseX)
                  forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) &
                    distance_matrix(ii,jj) = ( x_i(ii) - x_i(jj) )**2 / this%theta(ii)**2
                  gpCoordinates_Covariance_ii = 0.0_dp
                  do i_p = 1, this%n_permutations
                     if( any( (/ (distance_matrix(ii,this%permutations(ii,i_p)) > 1.0_dp, ii=1, this%d) /) ) ) cycle
                     r_ii = sqrt( sum( (/ (distance_matrix(ii,this%permutations(ii,i_p)), ii=1, this%d) /) ) )
                     if( r_ii >= 1.0_dp ) cycle

                     gpCoordinates_Covariance_ii = gpCoordinates_Covariance_ii + covariancePP(r_ii,PP_Q, this%d)
                  enddo
                  this%sparseCovariance(i_sparseX) = gpCoordinates_Covariance_ii
               end do

               deallocate(distance_matrix)
            end if

         endselect
  
         this%sparse_covariance_initialised = .true.  
      end if initialise_sparse_covariance  
  
   end subroutine gpCoordinates_precalculate_sparse  
  
   subroutine gpCoordinates_initialise_variance_estimate(this, regularisation, error)  
      type(gpCoordinates), intent(inout), target :: this
      real(dp), intent(in) :: regularisation  
      integer, intent(out), optional :: error

      real(dp) :: r_ij
      real(dp), dimension(:,:), allocatable :: k_mm, distance_matrix
      real(dp), dimension(:), pointer :: x_i, x_j
      real(dp), pointer :: fc_i, fc_j
      real(dp), dimension(this%d) :: xI_xJ_theta

      integer :: i, j, i_p, ii, jj, zeta_int

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_initialise_variance_estimate: object not initialised', error)  
      endif

      if( this%variance_estimate_initialised ) then  
         if( regularisation .feq. this%variance_estimate_regularisation) then  
            return
         else
            call gpCoordinates_finalise_variance_estimate(this,error)  
         endif
      endif

      if( regularisation < 0.0_dp ) then
         RAISE_ERROR("gpCoordinates_initialise_variance_estimate: regularisation ("//regularisation//") is negative.",error)  
      elseif( regularisation == 0.0_dp ) then
         call print_warning("gpCoordinates_initialise_variance_estimate: regularisation = 0.0, proceed with caution")  
      endif

      this%variance_estimate_regularisation = regularisation  

      allocate(k_mm(this%n_sparseX,this%n_sparseX))
      zeta_int = int(this%zeta)

      if( this%covariance_type == COVARIANCE_PP ) allocate(distance_matrix(this%d,this%d))

      if (this%covariance_type == COVARIANCE_DOT_PRODUCT) then
          call dgemm('T', 'N', size(this%sparseX,2), size(this%sparseX,2), size(this%sparseX,1), &
            1.0_dp, this%sparseX(1,1), size(this%sparseX,1), this%sparseX(1,1), size(this%sparseX, 1), & 
            0.0_dp, k_mm(1,1), size(k_mm,1))
         k_mm = fast_pow_2d(k_mm, this%zeta)
      else
          k_mm = 0.0_dp
          do i = 1, this%n_sparseX
             x_i => this%sparseX(:,i)
             fc_i => this%sparseCutoff(i)
             do j = i, this%n_sparseX
                x_j => this%sparseX(:,j)
                fc_j => this%sparseCutoff(j)
                if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
                   if( .not. this%initialised ) then
                      RAISE_ERROR('gpCoordinates_initialise_variance_estimate: bond real space sparse score not implemented', error)
                   endif
                elseif(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
                   if( zeta_int .feq. this%zeta ) then
                      k_mm(j,i) = fc_i*fc_i * sum( x_i * x_j )**zeta_int
                   else
                      k_mm(j,i) = fc_i*fc_i * sum( x_i * x_j )**this%zeta
                   endif
                elseif( this%covariance_type == COVARIANCE_ARD_SE ) then
                   do i_p = 1, this%n_permutations
                      xI_xJ_theta = (x_i(this%permutations(:,i_p)) - x_j) / this%theta
                      !xI_xJ_theta = (x_i - x_j(this%permutations(:,i_p))) / this%theta
                      k_mm(j,i) = k_mm(j,i) + exp( -0.5_dp * dot_product(xI_xJ_theta,xI_xJ_theta) )
                   enddo
                elseif( this%covariance_type == COVARIANCE_PP ) then
                   forall( ii = 1:this%d, jj = 1:this%d, this%permutation_distance_mask(ii,jj) ) distance_matrix(ii,jj) = ( x_i(ii) - x_j(jj) )**2 / this%theta(ii)**2
                   do i_p = 1, this%n_permutations
                      if( any( (/ (distance_matrix(ii,this%permutations(ii,i_p)) > 1.0_dp, ii=1, this%d) /) ) ) cycle
                      r_ij = sqrt( sum( (/ (distance_matrix(ii,this%permutations(ii,i_p)), ii=1, this%d) /) ) )
                      if( r_ij >= 1.0_dp ) cycle

                      k_mm(j,i) = k_mm(j,i) + covariancePP(r_ij,PP_Q, this%d)
                   enddo
                endif
                if( i /= j ) k_mm(i,j) = k_mm(j,i)
             enddo
          enddo
      endif

      if( this%covariance_type == COVARIANCE_ARD_SE .or. this%covariance_type == COVARIANCE_PP ) then
         do i = 1, this%n_sparseX
            fc_i => this%sparseCutoff(i)
            do j = i+1, this%n_sparseX
               fc_j => this%sparseCutoff(j)
               k_mm(j,i) = k_mm(j,i)  * fc_i * fc_j / sqrt(k_mm(j,j)*k_mm(i,i))
               k_mm(i,j) = k_mm(j,i)
            enddo
         enddo
         do i = 1, this%n_sparseX
            fc_i => this%sparseCutoff(i)
            k_mm(i,i) = fc_i**2
         enddo
      endif

      k_mm = k_mm * this%delta**2
      k_mm = k_mm + this%f0**2

      do i = 1, this%n_sparseX
         k_mm(i,i) = k_mm(i,i) + regularisation**2  
      enddo

      call initialise(this%LA_k_mm, k_mm)
      call LA_Matrix_Factorise(this%LA_k_mm,error=error)  
      if(allocated(k_mm)) deallocate(k_mm)
      if(allocated(distance_matrix)) deallocate(distance_matrix)

      this%variance_estimate_initialised = .true.  

   endsubroutine gpCoordinates_initialise_variance_estimate  

   function gpCoordinates_log_likelihood(this,regularisation,error) result(log_likelihood)
      type(gpCoordinates), intent(inout) :: this
      real(dp), intent(in), optional :: regularisation
      integer, intent(out), optional :: error
      real(dp) :: log_likelihood
                     
      real(dp) :: my_regularisation
      logical :: was_initialised

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_log_likelihood: object not initialised', error)
      endif

      was_initialised = this%variance_estimate_initialised
      if( this%variance_estimate_initialised ) then
         my_regularisation = optional_default(this%variance_estimate_regularisation,regularisation)
      else
         my_regularisation = optional_default(0.001_dp,regularisation)
      endif

      call gpCoordinates_initialise_variance_estimate(this,my_regularisation,error)

      log_likelihood = -0.5_dp * sum(matmul(this%LA_k_mm%matrix,this%alpha)*this%alpha) &
         - 0.5_dp*LA_Matrix_LogDet(this%LA_k_mm) - this%n_sparseX * log(2.0_dp*pi)

      if( .not. was_initialised ) call gpCoordinates_finalise_variance_estimate(this,error)

   endfunction gpCoordinates_log_likelihood

   subroutine gpCoordinates_finalise_variance_estimate(this,error)  
      type(gpCoordinates), intent(inout) :: this
      integer, intent(out), optional :: error
                     
      INIT_ERROR(error)

      if( .not. this%variance_estimate_initialised) return  

      call finalise(this%LA_k_mm)

      this%variance_estimate_regularisation = 0.0_dp  
      this%variance_estimate_initialised = .false.  

   endsubroutine gpCoordinates_finalise_variance_estimate  

   subroutine gpCoordinates_print_sparseX_file(this,sparseX_filename,error)
      type(gpCoordinates), intent(in) :: this
      character(len=*), intent(in) :: sparseX_filename
      integer, intent(out), optional :: error

      INIT_ERROR(error)

      call fwrite_array_d(size(this%sparseX), this%sparseX(1,1), trim(sparseX_filename)//C_NULL_CHAR)

   end subroutine gpCoordinates_print_sparseX_file

   ! print covariances and lambda to process-dependent files, one value per line
   subroutine gpFull_print_covariances_lambda(this, file_prefix, my_proc)
      type(gpFull), intent(in) :: this
      character(*), intent(in) :: file_prefix
      integer, intent(in), optional :: my_proc

      integer :: my_proc_opt

      my_proc_opt = optional_default(0, opt_val=my_proc)
      
      if (my_proc_opt == 0) then
         call fwrite_array_d(size(this%covariance_subY_subY), this%covariance_subY_subY, trim(file_prefix)//'_Kmm'//C_NULL_CHAR)
      end if
      call fwrite_array_d(size(this%covariance_subY_y), this%covariance_subY_y, trim(file_prefix)//'_Kmn.'//my_proc_opt//C_NULL_CHAR)
      call fwrite_array_d(size(this%lambda), this%lambda, trim(file_prefix)//'_lambda.'//my_proc_opt//C_NULL_CHAR)
   end subroutine gpFull_print_covariances_lambda

   subroutine gpCoordinates_printXML(this,xf,label,sparseX_base_filename,error)
      type(gpCoordinates), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      character(len=*), intent(in), optional :: sparseX_base_filename
      integer, intent(out), optional :: error

      integer :: i, j, j_end, slash_ind
      type(extendable_str) :: sparseX_filename
      character(len=32) :: sparseX_md5sum
      logical :: have_sparseX_base_filename

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_printXML: object not initialised', error)
      endif

      have_sparseX_base_filename = .false.
      if (present(sparseX_base_filename)) then
	 if (len_trim(sparseX_base_filename) > 0) have_sparseX_base_filename = .true.
      endif

      call xml_NewElement(xf,"gpCoordinates")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"dimensions", ""//this%d)
      call xml_AddAttribute(xf,"signal_variance", ""//this%delta)
      call xml_AddAttribute(xf,"signal_mean", ""//this%f0)
      call xml_AddAttribute(xf,"sparsified", ""//this%sparsified)
      call xml_AddAttribute(xf,"n_permutations", ""//this%n_permutations)
      call xml_AddAttribute(xf,"covariance_type", ""//this%covariance_type)

      if( this%covariance_type == COVARIANCE_DOT_PRODUCT ) &
          call xml_AddAttribute(xf,"zeta", ""//this%zeta)

      if(this%sparsified) then
         call xml_AddAttribute(xf,"n_sparseX",""//this%n_sparseX)
         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"sparseX_size_max", ""//maxval(this%sparseX_size))
	 if (have_sparseX_base_filename) then
            sparseX_filename = trim(sparseX_base_filename)
	    if (present(label)) then
               call concat(sparseX_filename,"."//trim(label))
            endif
            call gpCoordinates_print_sparseX_file(this,trim(string(sparseX_filename)),error=error)
            call quip_md5sum(trim(string(sparseX_filename)),sparseX_md5sum)
            ! remove leading path, since file will be read in from path of
            ! xml file
            slash_ind = index(sparseX_filename, "/")
            do while (slash_ind > 0)
                call substr_replace(sparseX_filename, 1, slash_ind, "")
                slash_ind = index(sparseX_filename, "/")
            end do
            call xml_AddAttribute(xf,"sparseX_filename",trim(string(sparseX_filename)))
            call xml_AddAttribute(xf,"sparseX_md5sum",trim(sparseX_md5sum))
	 endif
      else
         call xml_AddAttribute(xf,"n_x",""//this%n_x)
         call xml_AddAttribute(xf,"n_xPrime",""//this%n_xPrime)
         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"x_size_max", ""//maxval(this%x_size))
         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"xPrime_size_max", ""//maxval(this%xPrime_size))
      endif

      if( this%covariance_type == COVARIANCE_ARD_SE .or. this%covariance_type == COVARIANCE_PP ) then
         call xml_NewElement(xf,"theta")
         call xml_AddCharacters(xf, ""//this%theta//" ")
         call xml_EndElement(xf,"theta")
      endif

      call xml_NewElement(xf,"descriptor")
      call xml_AddCharacters(xf, string(this%descriptor_str))
      call xml_EndElement(xf,"descriptor")

      do i = 1, this%n_permutations
         call xml_NewElement(xf,"permutation")
         call xml_AddAttribute(xf,"i",""//i)
         call xml_AddCharacters(xf,""//this%permutations(:,i)//" ")
         call xml_EndElement(xf,"permutation")
      enddo

      if(this%sparsified) then
         do i = 1, this%n_sparseX
            call xml_NewElement(xf,"sparseX")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"alpha", ""//this%alpha(i))
            call xml_AddAttribute(xf,"sparseCutoff", ""//this%sparseCutoff(i))
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
               call xml_AddAttribute(xf,"covariance_sparseX_sparseX", ""//this%covarianceDiag_sparseX_sparseX(i))
            endif
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
               call xml_AddAttribute(xf,"sparseX_size", ""//this%sparseX_size(i))
               call xml_AddCharacters(xf, ""//this%sparseX(:this%sparseX_size(i),i)//"  ")
            elseif (.not. have_sparseX_base_filename) then
	       if(this%d <= 50) then
		  call xml_AddCharacters(xf, ""//this%sparseX(:,i)//"  ")
	       else
		  call xml_AddAttribute(xf,"sliced", "T")
		  do j = 1, this%d, 50
		     j_end = min(j-1+50,this%d)
		     call xml_NewElement(xf,"sparseX_slice")
		     call xml_AddAttribute(xf,"start", ""//j)
		     call xml_AddAttribute(xf,"end", ""//j_end)
		     call xml_AddCharacters(xf, ""//this%sparseX(j:j_end,i)//"  ")
		     call xml_EndElement(xf,"sparseX_slice")
		  enddo
	       endif
            endif
            call xml_EndElement(xf,"sparseX")
         enddo
      else
         do i = 1, this%n_x
            call xml_NewElement(xf,"x")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"map_x_y", ""//this%map_x_y(i))
            call xml_AddAttribute(xf,"cutoff", ""//this%cutoff(i))
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"x_size", ""//this%x_size(i))
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"covariance_x_x", ""//this%covarianceDiag_x_x(i))
            call xml_AddCharacters(xf, ""//this%x(:,i)//" ")
            call xml_EndElement(xf,"x")
         enddo
         do i = 1, this%n_xPrime
            call xml_NewElement(xf,"xPrime")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"map_xPrime_yPrime", ""//this%map_xPrime_yPrime(i))
            call xml_AddAttribute(xf,"map_xPrime_x", ""//this%map_xPrime_x(i))
            call xml_AddAttribute(xf,"cutoffPrime", ""//this%cutoffPrime(i))
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"xPrime_size", ""//this%xPrime_size(i))
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) call xml_AddAttribute(xf,"covariance_xPrime_xPrime", ""//this%covarianceDiag_xPrime_xPrime(i))
            call xml_AddCharacters(xf, ""//this%xPrime(:,i)//" ")
            call xml_EndElement(xf,"xPrime")
         enddo
      endif

      call xml_EndElement(xf,"gpCoordinates")

   endsubroutine gpCoordinates_printXML

   subroutine gpFull_printXML(this,xf,label,error)
      type(gpFull), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_printXML: object not initialised', error)
      endif

      call xml_NewElement(xf,"gpFull")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"n_y", ""//this%n_y)
      call xml_AddAttribute(xf,"n_yPrime", ""//this%n_yPrime)
      call xml_AddAttribute(xf,"n_globalSparseX", ""//this%n_globalSparseX)
      call xml_AddAttribute(xf,"n_coordinate", ""//this%n_coordinate)
      call xml_AddAttribute(xf,"sparse_jitter", ""//this%sparse_jitter)

      do i = 1, this%n_y
         call xml_NewElement(xf,"y")
         call xml_AddAttribute(xf,"i", ""//i)
         call xml_AddAttribute(xf,"map_y_globalY", ""//this%map_y_globalY(i))
         call xml_AddAttribute(xf,"alpha", ""//this%alpha(this%map_y_globalY(i)) )
         call xml_EndElement(xf,"y")
      enddo

      do i = 1, this%n_yPrime
         call xml_NewElement(xf,"yPrime")
         call xml_AddAttribute(xf,"i", ""//i)
         call xml_AddAttribute(xf,"map_yPrime_globalY", ""//this%map_yPrime_globalY(i))
         call xml_AddAttribute(xf,"alpha", ""//this%alpha(this%map_yPrime_globalY(i)) )
         call xml_EndElement(xf,"yPrime")
      enddo

      do i = 1, this%n_coordinate
         call gpCoordinates_printXML(this%coordinate(i),xf,label=trim(optional_default("",label))//i,error=error)
      enddo

      call xml_EndElement(xf,"gpFull")

   endsubroutine gpFull_printXML

   subroutine gpSparse_printXML(this,xf,label,sparseX_base_filename,error)
      type(gpSparse), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      character(len=*), intent(in), optional :: sparseX_base_filename
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpSparse_printXML: object not initialised', error)
      endif

      call xml_NewElement(xf,"gpSparse")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"n_coordinate", ""//this%n_coordinate)
      call xml_AddAttribute(xf,"fitted", ""//this%fitted)

      do i = 1, this%n_coordinate
         call gpCoordinates_printXML(this%coordinate(i),xf,label=trim(optional_default("",label))//i,&
				     sparseX_base_filename=sparseX_base_filename, error=error)
      enddo

      call xml_EndElement(xf,"gpSparse")

   endsubroutine gpSparse_printXML

   subroutine gpCoordinates_readXML(this,xp,label,error)
      type(gpCoordinates), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpCoordinates = .false.
      parse_matched_label = .false.
      parse_gpCoordinates => this
      parse_gpCoordinates_label = optional_default("",label)

      call initialise(parse_cur_data)
      call parse(xp, &
         characters_handler = gpCoordinates_characters_handler, &
         startElement_handler = gpCoordinates_startElement_handler, &
         endElement_handler = gpCoordinates_endElement_handler)

      call finalise(parse_cur_data)

      this%initialised = .true.

   endsubroutine gpCoordinates_readXML

   subroutine gpFull_readXML(this,xp,label,error)
      type(gpFull), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpFull = .false.
      parse_matched_label = .false.
      parse_gpFull => this
      parse_gpFull_label = optional_default("",label)

      call initialise(parse_cur_data)

      call parse(xp, &
         characters_handler = gpFull_characters_handler, &
         startElement_handler = gpFull_startElement_handler, &
         endElement_handler = gpFull_endElement_handler)

      call finalise(parse_cur_data)

      do i = 1, this%n_coordinate
         call gpCoordinates_readXML(this%coordinate(i),xp,label=trim(parse_gpFull_label)//i,error=error)
      enddo

      this%initialised = .true.

   endsubroutine gpFull_readXML

   subroutine gpSparse_readXML(this,xp,label,error)
      type(gpSparse), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

!      integer :: i

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpSparse = .false.
      parse_gpSparse => this
      parse_matched_label = .false.
      parse_gpSparse_label = optional_default("",label)

      call initialise(parse_cur_data)

      call parse(xp, &
         characters_handler = gpSparse_characters_handler, &
         startElement_handler = gpSparse_startElement_handler, &
         endElement_handler = gpSparse_endElement_handler)

      call finalise(parse_cur_data)

!      do i = 1, this%n_coordinate
!         call gpCoordinates_readXML(this%coordinate(i),xp,label=trim(parse_gpSparse_label)//i,error=error)
!      enddo

      this%initialised = .true.

   endsubroutine gpSparse_readXML

   subroutine gpFull_readXML_string(this,params_str,label,error)
      type(gpFull), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

   endsubroutine gpFull_readXML_string

   subroutine gpCoordinates_readXML_string(this,params_str,label,error)
      type(gpCoordinates), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

   endsubroutine gpCoordinates_readXML_string

   subroutine gpSparse_readXML_string(this,params_str,label,error)
      type(gpSparse), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp
      integer :: i

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

      do i = 1, this%n_coordinate
         call gp_readXML(this%coordinate(i),params_str,label=trim(parse_gpSparse_label)//i,error=error)
         call gpCoordinates_precalculate_sparse(this%coordinate(i))
      enddo

   endsubroutine gpSparse_readXML_string

   subroutine gpCoordinates_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      real(dp) :: delta, f0
      integer :: status, d, n_sparseX, n_x, n_xPrime, n_permutations, i, x_size_max, xPrime_size_max, sparseX_size_max, covariance_type
      logical :: sparsified, exist_sparseX_filename
      character(len=32) :: sparseX_md5sum
      character(len=1024) :: value

      if(name == 'gpCoordinates') then ! new GP_data
         if(parse_in_gpCoordinates) then
            call system_abort("gpCoordinates_startElement_handler entered gpCoordinates with parse_in_gpCoordinates true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpCoordinates_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpCoordinates_label)) then
               parse_matched_label = .true.
               parse_in_gpCoordinates = .true.
            else ! no match
               parse_in_gpCoordinates = .false.
            endif
         else ! no label passed in
            parse_in_gpCoordinates = .true.
         endif

         if(parse_in_gpCoordinates) then
            if(parse_gpCoordinates%initialised) call finalise(parse_gpCoordinates)

            call GP_FoX_get_value(attributes, 'dimensions', value, status)
            if (status == 0) then
               read (value,*) d
            else
               call system_abort("gpCoordinates_startElement_handler did not find the dimensions attribute.")
            endif

            call GP_FoX_get_value(attributes, 'signal_variance', value, status)
            if (status == 0) then
               read (value,*) delta
            else
               call system_abort("gpCoordinates_startElement_handler did not find the signal_variance attribute.")
            endif

            call GP_FoX_get_value(attributes, 'signal_mean', value, status)
            if (status == 0) then
               read (value,*) f0
            else
               call system_abort("gpCoordinates_startElement_handler did not find the signal_variance attribute.")
            endif


            call GP_FoX_get_value(attributes, 'sparsified', value, status)
            if (status == 0) then
               read (value,*) sparsified
            else
               call system_abort("gpCoordinates_startElement_handler did not find the sparsified attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_permutations', value, status)
            if (status == 0) then
               read (value,*) n_permutations
            else
               call system_abort("gpCoordinates_startElement_handler did not find the n_permutations attribute.")
            endif

            call GP_FoX_get_value(attributes, 'covariance_type', value, status)
            if (status == 0) then
               read (value,*) covariance_type
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_type attribute.")
               covariance_type = COVARIANCE_NONE
            endif

            call GP_FoX_get_value(attributes, 'zeta', value, status)
            if (status == 0) then
               if (covariance_type == COVARIANCE_DOT_PRODUCT) then
                  read (value,*) parse_gpCoordinates%zeta
               else
                  call system_abort("gpCoordinates_startElement_handler found zeta attribute but the covariance is not &
                  dot product.")
               endif
            else
               if (covariance_type == COVARIANCE_DOT_PRODUCT) then
                  call print_warning("gpCoordinates_startElement_handler: covariance type is dot product, but no &
                  zeta attribute is present. This may mean an XML generated by an older version. If found, the single &
                  value from the theta element will be used, to ensure backwards compatibility")
               endif
            endif

            call GP_FoX_get_value(attributes, 'x_size_max', value, status)
            if (status == 0) then
               read (value,*) x_size_max
            else
               if ((covariance_type == COVARIANCE_BOND_REAL_SPACE) .and. (.not. sparsified)) call system_abort("gpCoordinates_startElement_handler did not find the x_size_max attribute.")
               x_size_max = 0
            endif

            call GP_FoX_get_value(attributes, 'xPrime_size_max', value, status)
            if (status == 0) then
               read (value,*) xPrime_size_max
            else
               if ((covariance_type == COVARIANCE_BOND_REAL_SPACE) .and. (.not. sparsified)) call system_abort("gpCoordinates_startElement_handler did not find the xPrime_size_max attribute.")
               xPrime_size_max = 0
            endif

            call GP_FoX_get_value(attributes, 'sparseX_size_max', value, status)
            if (status == 0) then
               read (value,*) sparseX_size_max
            else
               if ((covariance_type == COVARIANCE_BOND_REAL_SPACE) .and. sparsified) call system_abort("gpCoordinates_startElement_handler did not find the sparseX_size_max attribute.")
               sparseX_size_max = 0
            endif

            if(sparsified) then
               call GP_FoX_get_value(attributes, 'n_sparseX', value, status)
               if (status == 0) then
                  read (value,*) n_sparseX
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_sparseX attribute.")
               endif

               if (covariance_type == COVARIANCE_BOND_REAL_SPACE) then
                  call gpCoordinates_setParameters_sparse(parse_gpCoordinates,d,n_sparseX,delta,f0,covariance_type=covariance_type,sparseX_size_max=sparseX_size_max)
               else
                  call gpCoordinates_setParameters_sparse(parse_gpCoordinates,d,n_sparseX,delta,f0, covariance_type=covariance_type)
		  call GP_FoX_get_value(attributes, 'sparseX_filename', value, status)
		  if (status == 0) then
                     inquire(file=trim(value),exist=exist_sparseX_filename)
                     if(.not.exist_sparseX_filename) call system_abort("gpCoordinates_startElement_handler: sparseX file "//trim(value)//" does not exist.")

                     call quip_md5sum(trim(value),sparseX_md5sum)
                     if( len_trim(sparseX_md5sum) == 0 ) call print_warning("gpCoordinates_startElement_handler: could not obtain md5 sum of sparse file, will not be able &
                        & to verify consistency with the XML")

		     call fread_array_d(size(parse_gpCoordinates%sparseX), parse_gpCoordinates%sparseX(1,1), trim(value)//C_NULL_CHAR)
		     parse_sparseX_separate_file = .true.
		  else
		     parse_sparseX_separate_file = .false.
		  endif

                  if(parse_sparseX_separate_file) then
                     call GP_FoX_get_value(attributes, 'sparseX_md5sum', value, status)
                     if (status == 0) then
                        if( len_trim(value) /= 32 ) call print_warning("gpCoordinates_startElement_handler: recorded md5 sum in the XML is not 32 characters. &
                           & This could have happened because the md5 tool was not available when the XML was written.")
                        if( len_trim(value) > 0 .and. len_trim(sparseX_md5sum) > 0 .and. trim(sparseX_md5sum) /= trim(value) ) then
                           call system_abort("gpCoordinates_startElement_handler: md5 check sum failed. Sparse file ("//sparseX_md5sum// &
                              ") does not match record in XML ("//trim(value)//")")
                        endif
                     endif
                  endif
               endif
            else
               call GP_FoX_get_value(attributes, 'n_x', value, status)
               if (status == 0) then
                  read (value,*) n_x
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_x attribute.")
               endif

               call GP_FoX_get_value(attributes, 'n_xPrime', value, status)
               if (status == 0) then
                  read (value,*) n_xPrime
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_xPrime attribute.")
               endif

               if (covariance_type == COVARIANCE_BOND_REAL_SPACE) then
                  call gpCoordinates_setParameters(parse_gpCoordinates,d,n_x,n_xPrime,delta,f0,covariance_type=covariance_type,x_size_max=x_size_max,xPrime_size_max=xPrime_size_max)
               else
                  call gpCoordinates_setParameters(parse_gpCoordinates,d,n_x,n_xPrime,delta,f0,covariance_type=covariance_type)
               endif
            endif

            if (covariance_type == COVARIANCE_BOND_REAL_SPACE .or. covariance_type == COVARIANCE_DOT_PRODUCT) then
               allocate(parse_in_permutations(1,n_permutations))
            else
               allocate(parse_in_permutations(d,n_permutations))
            endif

         endif

      elseif(parse_in_gpCoordinates .and. name == 'theta') then
         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'descriptor') then
         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'permutation') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         parse_i_permutation = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. name == 'sparseX') then

         parse_in_sparseX = .true.

         if( .not. parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: not sparsified data and sparseX element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%alpha(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the alpha attribute.")
         endif

         call GP_FoX_get_value(attributes, 'sparseCutoff', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%sparseCutoff(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the cutoff attribute.")
         endif

         call GP_FoX_get_value(attributes, 'sliced', value, status)
         if (status == 0) then
            read (value,*) parse_sliced
         else
            parse_sliced = .false.
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'sparseX_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%sparseX_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the sparseX_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'covariance_sparseX_sparseX', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_sparseX_sparseX(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_sparseX_sparseX attribute.")
            endif
         endif

         parse_i_sparseX = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. parse_in_sparseX .and. name == 'sparseX_slice') then

         call GP_FoX_get_value(attributes, 'start', value, status)
         if (status == 0) then
            read (value,*) parse_slice_start
         else
            call system_abort("gpCoordinates_startElement_handler did not find the start attribute.")
         endif

         call GP_FoX_get_value(attributes, 'end', value, status)
         if (status == 0) then
            read (value,*) parse_slice_end
         else
            call system_abort("gpCoordinates_startElement_handler did not find the end attribute.")
         endif

         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'x') then
         if( parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: sparsified=T but x element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_x_y', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_x_y(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_x_y attribute.")
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'x_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%x_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the x_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'covariance_x_x', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_x_x(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_x_x attribute.")
            endif
         endif

         parse_i_x = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. name == 'xPrime') then
         if( parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: sparsified=T but xPrime element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_xPrime_yPrime', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_xPrime_yPrime(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_xPrime_yPrime attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_xPrime_x', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_xPrime_x(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_xPrime_x attribute.")
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'xPrime_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%xPrime_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the xPrime_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
            call GP_FoX_get_value(attributes, 'covariance_xPrime_xPrime', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_xPrime_xPrime(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_xPrime_xPrime attribute.")
            endif
         endif

         parse_i_xPrime = i

         call zero(parse_cur_data)

      endif

   endsubroutine gpCoordinates_startElement_handler

   subroutine gpCoordinates_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      if(parse_in_gpCoordinates) then
         if(name == 'gpCoordinates') then
            call gpCoordinates_setPermutations(parse_gpCoordinates,parse_in_permutations)
            deallocate(parse_in_permutations)
            parse_in_gpCoordinates = .false.
         elseif(name == 'theta') then
            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%theta
            call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%theta)
            if( parse_gpCoordinates%covariance_type == COVARIANCE_DOT_PRODUCT ) then
               parse_gpCoordinates%zeta = parse_gpCoordinates%theta(1)
               parse_gpCoordinates%theta(1) = 0.0_dp
               call print_warning("gpCoordinates_endElement_handler: dot product covariance is used, but found a theta element &
               in the XML. This may be a sign of an XML generated by an older version. The first and only element of theta will &
               be used as zeta.")
            endif
         elseif(name == 'descriptor') then
            parse_gpCoordinates%descriptor_str = parse_cur_data
         elseif(name == 'permutation') then
            
            if( parse_i_permutation > size(parse_in_permutations,2) ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_permutation ("//parse_i_permutation//") greater than n_permutations ("//size(parse_in_permutations,2)//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_in_permutations(:,parse_i_permutation)
            call string_to_numerical(string(parse_cur_data),parse_in_permutations(:,parse_i_permutation))
         elseif(name == 'sparseX') then
            
            if( .not. allocated(parse_gpCoordinates%sparseX) ) then
               call system_abort("gpCoordinates_endElement_handler: sparseX not allocated")
            endif
            
            if( parse_i_sparseX > parse_gpCoordinates%n_sparseX ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_sparseX ("//parse_i_sparseX//") greater than n_sparseX ("//parse_gpCoordinates%n_sparseX//")")
            endif 

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%sparseX(:,parse_i_sparseX)
            if( parse_gpCoordinates%covariance_type == COVARIANCE_BOND_REAL_SPACE ) then
               parse_gpCoordinates%sparseX(:,parse_i_sparseX) = 0.0_dp
               call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(:parse_gpCoordinates%sparseX_size(parse_i_sparseX),parse_i_sparseX))
            else
               if(.not. parse_sparseX_separate_file .and. .not. parse_sliced) call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(:,parse_i_sparseX))
            endif

            parse_in_sparseX = .false.
         elseif(name == 'sparseX_slice') then
            if(parse_slice_start < 1) then
               call system_abort("gpCoordinates_endElement_handler: slice start less than 1")
            endif

            if(parse_slice_end > parse_gpCoordinates%d) then
               call system_abort("gpCoordinates_endElement_handler: slice start greater than dimension")
            endif

            if(.not. parse_sparseX_separate_file .and. parse_sliced) call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(parse_slice_start:parse_slice_end,parse_i_sparseX))
         elseif(name == 'x') then

            if( .not. allocated(parse_gpCoordinates%x) ) then
               call system_abort("gpCoordinates_endElement_handler: x not allocated")
            endif
            
            if( parse_i_x > parse_gpCoordinates%n_x ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_x ("//parse_i_x//") greater than n_x ("//parse_gpCoordinates%n_x//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%x(:,parse_i_x)
            call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%x(:,parse_i_x))
         elseif(name == 'xPrime') then

            if( .not. allocated(parse_gpCoordinates%xPrime) ) then
               call system_abort("gpCoordinates_endElement_handler: xPrime not allocated")
            endif
            
            if( parse_i_xPrime > parse_gpCoordinates%n_xPrime ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_xPrime ("//parse_i_xPrime//") greater than n_xPrime ("//parse_gpCoordinates%n_xPrime//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%xPrime(:,parse_i_xPrime)
            call string_to_numerical(string(parse_cur_data), parse_gpCoordinates%xPrime(:,parse_i_xPrime))
         endif
      endif

   endsubroutine gpCoordinates_endElement_handler

   subroutine gpCoordinates_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpCoordinates) then
         call concat(parse_cur_data, in, keep_lf=.false.,lf_to_whitespace=.true.)
      endif
   endsubroutine gpCoordinates_characters_handler

   subroutine gpFull_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      integer :: status, n_y, n_yPrime, n_coordinate, i
      real(dp) :: sparse_jitter
      character(len=1024) :: value

      if(name == 'gpFull') then ! new GP_data
         if(parse_in_gpFull) then
            call system_abort("gpFull_startElement_handler entered gpFull with parse_in_gpFull true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpFull_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpFull_label)) then
               parse_matched_label = .true.
               parse_in_gpFull = .true.
            else ! no match
               parse_in_gpFull = .false.
            endif
         else ! no label passed in
            parse_in_gpFull = .true.
         endif

         if(parse_in_gpFull) then
            if(parse_gpFull%initialised) call finalise(parse_gpFull)

            call GP_FoX_get_value(attributes, 'n_y', value, status)
            if (status == 0) then
               read (value,*) n_y
            else
               call system_abort("gpFull_startElement_handler did not find the n_y attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_yPrime', value, status)
            if (status == 0) then
               read (value,*) n_yPrime
            else
               call system_abort("gpFull_startElement_handler did not find the n_yPrime attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_coordinate', value, status)
            if (status == 0) then
               read (value,*) n_coordinate
            else
               call system_abort("gpFull_startElement_handler did not find the n_coordinate attribute.")
            endif

            call GP_FoX_get_value(attributes, 'sparse_jitter', value, status)
            if (status == 0) then
               read (value,*) sparse_jitter
            else
               call print_warning("gpFull_startElement_handler did not find the sparse_jitter attribute, using default value 1.0e-5.")
               sparse_jitter = 1.0e-5_dp
            endif
            call gpFull_setParameters(parse_gpFull,n_coordinate, n_y, n_yPrime, sparse_jitter)

         endif

      elseif(parse_in_gpFull .and. name == 'y') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpFull_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_y_globalY', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%map_y_globalY(i)
         else
            call system_abort("gpFull_startElement_handler did not find the map_y_globalY attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%alpha(parse_gpFull%map_y_globalY(i))
         else
            call system_abort("gpFull_startElement_handler did not find the alpha attribute.")
         endif

      elseif(parse_in_gpFull .and. name == 'yPrime') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpFull_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_yPrime_globalY', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%map_yPrime_globalY(i)
         else
            call system_abort("gpFull_startElement_handler did not find the map_yPrime_globalY attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%alpha(parse_gpFull%map_yPrime_globalY(i))
         else
            call system_abort("gpFull_startElement_handler did not find the alpha attribute.")
         endif

      endif

   endsubroutine gpFull_startElement_handler

   subroutine gpFull_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      if(parse_in_gpFull) then
         if(name == 'gpFull') then
            parse_in_gpFull = .false.
         endif
      elseif(name == 'y') then

      elseif(name == 'yPrime') then

      endif

   endsubroutine gpFull_endElement_handler

   subroutine gpFull_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpFull) then
         call concat(parse_cur_data, in, keep_lf=.false.,lf_to_whitespace=.true.)
      endif
   endsubroutine gpFull_characters_handler

   subroutine gpSparse_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      integer :: status, n_coordinate
      character(len=1024) :: value

      if(name == 'gpSparse') then ! new GP_data
         if(parse_in_gpSparse) then
            call system_abort("gpSparse_startElement_handler entered gpSparse with parse_in_gpSparse true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpSparse_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpSparse_label)) then
               parse_matched_label = .true.
               parse_in_gpSparse = .true.
            else ! no match
               parse_in_gpSparse = .false.
            endif
         else ! no label passed in
            parse_in_gpSparse = .true.
         endif

         if(parse_in_gpSparse) then
            if(parse_gpSparse%initialised) call finalise(parse_gpSparse)

            call GP_FoX_get_value(attributes, 'n_coordinate', value, status)
            if (status == 0) then
               read (value,*) n_coordinate
            else
               call system_abort("gpSparse_startElement_handler did not find the n_coordinate attribute.")
            endif
            call gpSparse_setParameters(parse_gpSparse,n_coordinate)

            call GP_FoX_get_value(attributes, 'fitted', value, status)
            if (status == 0) then
               read (value,*) parse_gpSparse%fitted
            else
               parse_gpSparse%fitted = .true.  ! for backward compatibility
            endif

         endif

      endif

   endsubroutine gpSparse_startElement_handler

   subroutine gpSparse_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      if(parse_in_gpSparse) then
         if(name == 'gpSparse') then
            parse_in_gpSparse = .false.
         endif
      endif

   endsubroutine gpSparse_endElement_handler

   subroutine gpSparse_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpSparse) then
         call concat(parse_cur_data, in, keep_lf=.false.,lf_to_whitespace=.true.)
      endif
   endsubroutine gpSparse_characters_handler

   subroutine gp_FoX_get_value(attributes, key, val, status)
     type(dictionary_t), intent(in) :: attributes
     character(len=*), intent(in) :: key
     character(len=*), intent(inout) :: val
     integer, intent(out), optional :: status
            
     if (HasKey(attributes,key)) then
       val = GetValue(attributes, trim(key))
       if (present(status)) status = 0
     else
       val = "" 
       if (present(status)) status = 1
     endif
   end subroutine gp_FoX_get_value

end module gp_predict_module
