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

module gp_fit_module

   use iso_c_binding, only : C_NULL_CHAR
   ! use libatoms_module
   use error_module
   use system_module
   use extendable_str_module
   use linearalgebra_module
   use dictionary_module, only : STRING_LENGTH
   use gp_predict_module
   use clustering_module
   implicit none
   private

   integer, parameter, public :: EXCLUDE_CONFIG_TYPE = -10

   interface gp_sparsify
      module procedure gpFull_sparsify_array_config_type
   endinterface gp_sparsify
   public :: gp_sparsify

   public :: count_entries_in_sparse_file

   contains

   subroutine gpCoordinates_sparsify_config_type(this, n_sparseX, default_all, sparseMethod, sparse_file, use_actual_gpcov, print_sparse_index, &
         unique_hash_tolerance, unique_descriptor_tolerance, error)
      type(gpCoordinates), intent(inout) :: this
      integer, dimension(:), intent(in) :: n_sparseX
      logical, intent(in) :: default_all
      integer, optional, intent(in) :: sparseMethod
      character(len=STRING_LENGTH), intent(in), optional :: sparse_file, print_sparse_index
      logical, intent(in), optional :: use_actual_gpcov
      real(dp), intent(in), optional :: unique_descriptor_tolerance, unique_hash_tolerance
      integer, optional, intent(out) :: error

      integer :: my_sparseMethod, i, j, li, ui, i_config_type, n_config_type, d, n_x, n_sparse_file
      integer, dimension(:), allocatable :: config_type_index, sparseX_index, my_n_sparseX
      real(dp), dimension(:,:), allocatable :: x, sparseX_array
      real(dp), dimension(:), allocatable :: x_hash
      real(dp), pointer, dimension(:,:) :: dm
      integer, dimension(:), allocatable :: x_index
      logical, dimension(:), allocatable :: x_unique

      character(len=STRING_LENGTH) :: my_sparse_file
      type(Inoutput) :: inout_sparse_index
      logical :: exist_sparse_file

      real(dp) :: my_unique_hash_tolerance, my_unique_descriptor_tolerance

      INIT_ERROR(error)

      my_sparseMethod = optional_default(GP_SPARSE_RANDOM,sparseMethod)
      my_sparse_file = optional_default("",sparse_file)

      my_unique_hash_tolerance = optional_default(1.0e-10_dp,unique_hash_tolerance)
      my_unique_descriptor_tolerance = optional_default(1.0e-10_dp,unique_descriptor_tolerance)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_sparsify: : object not initialised',error)
      endif

      d = size(this%x,1)
      !n_x = size(this%x,2)
      n_x = count(EXCLUDE_CONFIG_TYPE /= this%config_type)

      allocate(my_n_sparseX(size(n_sparseX)))
      my_n_sparseX = 0

      ! Remove duplicates
      allocate(x_hash(n_x))
      allocate(x_index(n_x))

      ! Compute 1-norm hash on all descriptors that we want to include, and the mapping to the full vector
      j = 0
      do i = 1, size(this%x,2)
         if( this%config_type(i) /= EXCLUDE_CONFIG_TYPE ) then
            j = j + 1
            x_hash(j) = sum(abs(this%x(:,i)))
            x_index(j) = i
         endif
      enddo

      ! Sort hashes
      call heap_sort(x_hash,i_data=x_index)

      ! Compare neighbouring hashes. If they're within tolerance, compare the corresponding descriptors using the eucledian norm.
      ! Update the config type if they're equivalent.
      do j = 2, n_x
         if( abs( x_hash(j-1) - x_hash(j) ) < my_unique_hash_tolerance ) then
            !if ( sum( ( this%x(:,x_index(j))-this%x(:,x_index(j-1)) )**2 ) < my_unique_descriptor_tolerance ) then
            if ( maxval( abs( this%x(:,x_index(j))-this%x(:,x_index(j-1)) ) ) < my_unique_descriptor_tolerance ) then
               this%config_type(x_index(j-1)) = EXCLUDE_CONFIG_TYPE
            endif
         endif
      enddo

      deallocate(x_hash)
      deallocate(x_index)

      ! need to update n_x once duplicates are removed
      n_x = count(EXCLUDE_CONFIG_TYPE /= this%config_type)

      if(my_sparseMethod == GP_SPARSE_UNIQ) then
         RAISE_ERROR('gpCoordinates_sparsify: UNIQ is no longer in use, please use NONE instead.',error)

      !   allocate(x(d,n_x))
      !   allocate(x_index(n_x))
      !   allocate(x_unique(n_x))

      !   !x = this%x
      !   !x_index = (/(i,i=1,n_x)/)
      !   j = 0
      !   do i = 1, size(this%x,2)
      !      if( this%config_type(i) /= EXCLUDE_CONFIG_TYPE ) then
      !         j = j + 1
      !         x(:,j) = this%x(:,i)
      !         x_index(j) = i
      !      endif
      !   enddo

      !   call heap_sort(x,i_data=x_index)
      !   call uniq(x,unique=x_unique)
      !   this%n_sparseX = count(x_unique)

      !   call print('UNIQ type sparsification specified. The number of sparse points was changed to '//this%n_sparseX//' from '//n_sparseX//'.')

      elseif(my_sparseMethod == GP_SPARSE_NONE) then

         allocate(x_index(n_x))

         j = 0
         do i = 1, size(this%x,2)
            if( this%config_type(i) /= EXCLUDE_CONFIG_TYPE ) then
               j = j + 1
               x_index(j) = i
            endif
         enddo

         this%n_sparseX = n_x

         call print('NONE type sparsification specified. The number of sparse points was changed to '//this%n_sparseX//' from '//n_sparseX//'.')

      elseif(my_sparseMethod == GP_SPARSE_FILE .or. my_sparseMethod == GP_SPARSE_INDEX_FILE) then
         this%n_sparseX = count_entries_in_sparse_file(my_sparse_file, my_sparseMethod, d, error)
      else
         do i_config_type = 1, size(n_sparseX)
            if(default_all) then

               if( n_x < sum(n_sparseX) ) then
                  call print_warning('gpCoordinates_sparsify: number of data points ('//n_x//') less than the number of sparse points ('//sum(n_sparseX)//'), &
                  number of sparse points changed to '//n_x)
                  call print_warning('gpCoordinates_sparsify: affected descriptor : '//this%descriptor_str)
                  my_n_sparseX(1) = n_x
               else
                  my_n_sparseX(1) = sum(n_sparseX)
               endif

            else
               if( n_sparseX(i_config_type) == 0 ) cycle

               n_config_type = count(i_config_type == this%config_type)

               if( n_config_type < n_sparseX(i_config_type) ) then
                  call print_warning('gpCoordinates_sparsify: number of data points ('//n_config_type//') less than the number of sparse points ('//n_sparseX(i_config_type)//'), &
                  number of sparse points changed to '//n_config_type)
                  call print_warning('gpCoordinates_sparsify: affected descriptor : '//this%descriptor_str)
                  my_n_sparseX(i_config_type) = n_config_type
               else
                  my_n_sparseX(i_config_type) = n_sparseX(i_config_type)
               endif
            endif
         
            if(default_all) exit
         enddo
         this%n_sparseX = sum(my_n_sparseX)
      endif

      call reallocate(this%sparseX, this%d,this%n_sparseX, zero = .true.)

      call reallocate(this%sparseX_index, this%n_sparseX, zero = .true.)
      call reallocate(this%map_sparseX_globalSparseX, this%n_sparseX, zero = .true.)
      call reallocate(this%alpha, this%n_sparseX, zero = .true.)
      call reallocate(this%sparseCutoff, this%n_sparseX, zero = .true.)
      this%sparseCutoff = 1.0_dp

      if( my_sparseMethod /= GP_SPARSE_FILE .and. my_sparseMethod /= GP_SPARSE_INDEX_FILE) then
         ui = 0
         do i_config_type = 1, size(my_n_sparseX)
            
            if( my_sparseMethod == GP_SPARSE_NONE) exit

            if(default_all) then

               allocate(config_type_index(n_x), sparseX_index(this%n_sparseX))
               !config_type_index = (/(i,i=1,n_x)/)
               j = 0
               do i = 1, size(this%x,2)
                  if( this%config_type(i) /= EXCLUDE_CONFIG_TYPE ) then
                     j = j + 1
                     config_type_index(j) = i
                  endif
               enddo

               li = 1
               ui = this%n_sparseX
               n_config_type = n_x
            else
               if( my_n_sparseX(i_config_type) == 0 ) cycle

               n_config_type = count(i_config_type == this%config_type)

               allocate(config_type_index(n_config_type),sparseX_index(my_n_sparseX(i_config_type)))
               config_type_index = find(i_config_type == this%config_type)

               li = ui + 1
               ui = ui + my_n_sparseX(i_config_type)
            endif
            
            select case(my_sparseMethod)
            case(GP_SPARSE_RANDOM)
               call fill_random_integer(sparseX_index, n_config_type)
            case(GP_SPARSE_PIVOT)
               if(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
                  call pivot(this%x(:,config_type_index), sparseX_index)
               else
                  call pivot(this%x(:,config_type_index), sparseX_index, theta = this%theta)
               endif
            case(GP_SPARSE_CLUSTER)
               if(use_actual_gpcov) then
                  call print('Started kernel distance matrix calculation')
                  dm => kernel_distance_matrix(this, config_type_index = config_type_index)
                  call print('Finished kernel distance matrix calculation')
               endif
               call print('Started kmedoids clustering')
               if(use_actual_gpcov) then
                  call bisect_kmedoids(dm, my_n_sparseX(i_config_type), med = sparseX_index)
                  !call bisect_kmedoids(dm, my_n_sparseX(i_config_type), c = c, med = sparseX_index)
               else
                  if(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
                     call bisect_kmedoids(this%x(:,config_type_index), my_n_sparseX(i_config_type), med = sparseX_index, is_distance_matrix = .false.)
                  else
                     call bisect_kmedoids(this%x(:,config_type_index), my_n_sparseX(i_config_type), med = sparseX_index, theta = this%theta, is_distance_matrix = .false.)
                  endif
               endif
               call print('Finished kmedoids clustering')
               if(use_actual_gpcov) deallocate(dm)
            case(GP_SPARSE_UNIFORM)
               call select_uniform(this%x(:,config_type_index), sparseX_index)
            case(GP_SPARSE_KMEANS)
               call print('Started kmeans clustering')
               if(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
                  call cluster_kmeans(this%x(:,config_type_index), sparseX_index)
               else
                  call cluster_kmeans(this%x(:,config_type_index), sparseX_index, theta = this%theta)
               endif
               call print('Finished kmeans clustering')
            case(GP_SPARSE_COVARIANCE)
               call sparse_covariance(this,sparseX_index,config_type_index,use_actual_gpcov)
            case(GP_SPARSE_FUZZY)
               call print('Started fuzzy cmeans clustering')
               if(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
                  call cluster_fuzzy_cmeans(this%x(:,config_type_index), sparseX_index, fuzziness=2.0_dp)
               else
                  call cluster_fuzzy_cmeans(this%x(:,config_type_index), sparseX_index, theta=this%theta,fuzziness=2.0_dp)
               endif
               call print('Finished fuzzy cmeans clustering')
            case(GP_SPARSE_CUR_COVARIANCE)
               call print("Started covariance matrix calculation")
               dm => kernel_distance_matrix(this, config_type_index=config_type_index, covariance_only = .true.)
               call print("Finished covariance matrix calculation")
               call print("Started CUR decomposition")
               call cur_decomposition(dm, sparseX_index)
               call print("Finished CUR decomposition")
               deallocate(dm)
            case(GP_SPARSE_CUR_POINTS)
               call print("Started CUR decomposition")
               call cur_decomposition(this%x(:,config_type_index), sparseX_index)
               call print("Finished CUR decomposition")
            case default
               RAISE_ERROR('gpCoordinates_sparsify: '//my_sparseMethod//' method is unknown', error)
            endselect
            this%sparseX_index(li:ui) = config_type_index(sparseX_index)
            deallocate(config_type_index,sparseX_index)

            if(default_all) exit
         enddo

      elseif(my_sparseMethod == GP_SPARSE_INDEX_FILE) then
         call print('Started reading sparse indices from file '//trim(my_sparse_file))
         call fread_array_i(size(this%sparseX_index),this%sparseX_index(1),trim(my_sparse_file)//C_NULL_CHAR)
         call print('Finished reading sparse indices from file, '//size(this%sparseX_index)//' of them.')
      endif

      if(allocated(this%covarianceDiag_sparseX_sparseX)) deallocate(this%covarianceDiag_sparseX_sparseX)
      allocate(this%covarianceDiag_sparseX_sparseX(this%n_sparseX))

      if(my_sparseMethod == GP_SPARSE_FILE) then
         call print('Started reading sparse descriptors from file '//trim(my_sparse_file))
         allocate(sparseX_array(d+1,this%n_sparseX))
         call fread_array_d(size(sparseX_array),sparseX_array(1,1),trim(my_sparse_file)//C_NULL_CHAR)
         this%sparseCutoff = sparseX_array(1,:)
         this%sparseX = sparseX_array(2:,:)
         this%covarianceDiag_sparseX_sparseX = 1.0_dp ! only used for COVARIANCE_BOND_REAL_SPACE
         deallocate(sparseX_array)
         call print('Finished reading sparse descriptors from file, '//size(this%sparseCutoff)//'  of them.')
      else
         if(my_sparseMethod == GP_SPARSE_NONE) this%sparseX_index = x_index

         call sort_array(this%sparseX_index)
         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            if(allocated(this%sparseX)) deallocate(this%sparseX)
            allocate(this%sparseX(maxval(this%x_size(this%sparseX_index)),this%n_sparseX))
            if(allocated(this%sparseX_size)) deallocate(this%sparseX_size)
            allocate(this%sparseX_size(this%n_sparseX))
            this%sparseX(:,:) = this%x(1:maxval(this%x_size(this%sparseX_index)),this%sparseX_index)
            this%sparseX_size = this%x_size(this%sparseX_index)
         else
            this%sparseX(:,:) = this%x(:,this%sparseX_index)
         endif

         this%covarianceDiag_sparseX_sparseX = this%covarianceDiag_x_x(this%sparseX_index)

         this%sparseCutoff = this%cutoff(this%sparseX_index)

         if(present(print_sparse_index)) then
            if(len_trim(print_sparse_index) > 0) then
               call initialise(inout_sparse_index, trim(print_sparse_index), action=OUTPUT, append=.true.)
               call print(""//this%sparseX_index,file=inout_sparse_index)
               call finalise(inout_sparse_index)
            endif
         endif
      endif


      if(allocated(x)) deallocate(x)
      if(allocated(x_index)) deallocate(x_index)
      if(allocated(x_unique)) deallocate(x_unique)
      if(allocated(my_n_sparseX)) deallocate(my_n_sparseX)

      this%sparsified = .true.

   endsubroutine gpCoordinates_sparsify_config_type

   function count_entries_in_sparse_file(sparse_file, sparse_method, d, error) result(res)
      character(len=*), intent(in) :: sparse_file
      integer, intent(in) :: sparse_method
      integer, intent(in) :: d ! coordinate_length
      integer, intent(out), optional :: error
      integer :: res

      logical :: exist_sparse_file
      integer :: n_sparse_file

      INIT_ERROR(error)

      inquire(file=trim(sparse_file), exist=exist_sparse_file)
      if (.not. exist_sparse_file) then
         RAISE_ERROR('count_entries_in_sparse_file: '//trim(sparse_file)//' does not exist', error)
      end if

      call fwc_l(trim(sparse_file)//C_NULL_CHAR, n_sparse_file)

      select case (sparse_method)
         case (GP_SPARSE_INDEX_FILE)
            res = n_sparse_file
         case (GP_SPARSE_FILE)
            if (mod(n_sparse_file, d+1) /= 0) then
               RAISE_ERROR('count_entries_in_sparse_file: file '//trim(sparse_file)//' contains '//n_sparse_file//" lines, not conforming with descriptor size "//d, error)
            end if
            res = n_sparse_file / (d + 1)
         case default
            RAISE_ERROR('count_entries_in_sparse_file: given sparse_method is not implemented: '//sparse_method, error)
      end select
   end function count_entries_in_sparse_file

   subroutine gpFull_sparsify_array_config_type(this, n_sparseX, default_all, sparseMethod, sparse_file, use_actual_gpcov, print_sparse_index, &
         unique_hash_tolerance, unique_descriptor_tolerance, error)
      type(gpFull), intent(inout) :: this
      integer, dimension(:,:), intent(in) :: n_sparseX
      logical, dimension(:), intent(in) :: default_all
      integer, dimension(:), optional, intent(in) :: sparseMethod
      character(len=STRING_LENGTH), dimension(:), intent(in), optional :: sparse_file, print_sparse_index
      logical, intent(in), optional :: use_actual_gpcov
      real(dp), dimension(:), optional, intent(in) :: unique_hash_tolerance, unique_descriptor_tolerance
      integer, optional, intent(out) :: error

      integer :: i
      integer, dimension(:), allocatable :: my_sparseMethod
      character(len=STRING_LENGTH), dimension(:), allocatable :: my_sparse_file

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_sparsify_array: object not initialised',error)
      endif

      allocate(my_sparseMethod(this%n_coordinate))
      allocate(my_sparse_file(this%n_coordinate))
      my_sparseMethod = optional_default((/ (GP_SPARSE_RANDOM, i=1,this%n_coordinate) /),sparseMethod)
      my_sparse_file = optional_default((/ ("", i=1,this%n_coordinate) /),sparse_file)

      do i = 1, this%n_coordinate
         call gpCoordinates_sparsify_config_type(this%coordinate(i),n_sparseX(:,i), default_all(i), &
            sparseMethod=my_sparseMethod(i), sparse_file=my_sparse_file(i), use_actual_gpcov=use_actual_gpcov, &
            print_sparse_index = print_sparse_index(i), &
            unique_hash_tolerance=unique_hash_tolerance(i), unique_descriptor_tolerance=unique_descriptor_tolerance(i), &
            error = error)
      enddo

      if(allocated(my_sparseMethod)) deallocate(my_sparseMethod)
      if(allocated(my_sparse_file)) deallocate(my_sparse_file)

   endsubroutine gpFull_sparsify_array_config_type

   function kernel_distance_matrix(this, config_type_index, covariance_only) result(k_nn)
      type(gpCoordinates), intent(in) :: this
      integer, dimension(:), intent(in), optional :: config_type_index
      logical, intent(in), optional :: covariance_only

      real(dp), pointer, dimension(:,:) :: k_nn ! actually the kernel distance matrix

      !real(dp), dimension(:,:), allocatable :: k_nn
      real(dp), dimension(:), allocatable :: k_self
      logical :: do_kernel_distance
      integer :: i, j, n, ii, jj
      integer :: stat

      call system_timer('kernel_distance_matrix')

      do_kernel_distance = .not. optional_default(.false., covariance_only)

      if(present(config_type_index)) then
         n = size(config_type_index)
      else
         n = size(this%x,2)
      endif
 
      allocate(k_self(n))

      allocate(k_nn(n,n), stat=stat)
      if(stat /= 0) call system_abort('kernel_distance_matrix: could not allocate matrix.')
 
!$omp parallel do default(none) shared(this,n,config_type_index,k_self) private(i,ii)
      do i = 1, n
         if(present(config_type_index)) then
            ii = config_type_index(i)
         else
            ii = i
         endif

         k_self(i) = gpCoordinates_Covariance(this, i_x = ii, j_x = ii, normalise = .false.)
      enddo

      do j = 1, n
         if(present(config_type_index)) then
            jj = config_type_index(j)
         else
            jj = j
         endif

         !k_nn(j,j) = 1.0_dp ! normalised kernel self-covariance
         k_nn(j,j) = 0.0_dp ! distance to itself = 0

!$omp parallel do default(none) shared(n,this,k_nn,jj,j,k_self,config_type_index,do_kernel_distance) private(i,ii)
         do i = j+1, n
            if(present(config_type_index)) then
               ii = config_type_index(i)
            else
               ii = i
            endif

            ! kernel covariance
            k_nn(j,i) = gpCoordinates_Covariance(this, i_x = ii, j_x = jj, normalise = .false.)
            ! then normalise
            k_nn(j,i) = k_nn(j,i) / sqrt(k_self(i)*k_self(j))

            if (do_kernel_distance) then
              ! now convert to distance
              k_nn(j,i) = sqrt(2.0_dp * (1.0_dp - k_nn(j,i)))
            endif

            ! finally, symmetrise
            k_nn(i,j) = k_nn(j,i)
         enddo ! i
      enddo ! j

      !dm = sqrt(2.0_dp * (1.0_dp - k_nn))
      !do i = 1, n
      !   do j = i+1, n
      !      dm(i,j) = sqrt(2.0_dp*(1.0_dp - kij))
      !      dm(j,i) = dm(i,j)
      !   end do
      !end do

      !deallocate(k_nn, k_self)
      deallocate(k_self)
      call system_timer('kernel_distance_matrix')
   end function kernel_distance_matrix

   subroutine sparse_covariance(this, index_out, config_type_index, use_actual_gpcov)
      type(gpCoordinates), intent(in) :: this
      integer, dimension(:), intent(out) :: index_out
      integer, dimension(:), intent(in), optional :: config_type_index
      logical, intent(in), optional :: use_actual_gpcov
 
      real(dp), dimension(:), allocatable :: score, k_self !, xI_xJ
      real(dp), dimension(:,:), allocatable :: k_mn, k_mm_k_m
      real(dp), dimension(1,1) :: k_mm
      integer :: m, n, i, ii, j, jj, i_p, zeta_int
      integer, dimension(1) :: j_loc
      logical, dimension(:), allocatable :: not_yet_added
      logical :: do_use_actual_gpcov
 
      type(LA_Matrix) :: LA_k_mm
 
      call system_timer('sparse_covariance')
      if(present(config_type_index)) then
         n = size(config_type_index)
      else
         n = size(this%x,2)
      endif
      m = size(index_out)

      do_use_actual_gpcov = optional_default(.false., use_actual_gpcov)
      if(do_use_actual_gpcov) then
         call print("sparse_covariance using actual gpCoordinates_Covariance")
      else
         call print("sparse_covariance using manual 'covariance'")
      endif
 
      allocate(k_mn(m,n), score(n), k_mm_k_m(m,n), k_self(n), not_yet_added(n))
      k_mn = 0.0_dp
      not_yet_added = .true.
 
      !allocate(xI_xJ(this%d))
 
      j = 1
      index_out(j) = 1 !ceiling(ran_uniform() * n)
      not_yet_added(index_out(j)) = .false.
 
      k_mm = 1.0_dp+1.0e-6_dp
      zeta_int = nint(this%zeta)
      call initialise(LA_k_mm,k_mm)
 
!$omp parallel do default(none) shared(this,n,config_type_index,k_self,do_use_actual_gpcov,zeta_int) private(i,ii,i_p)
      do i = 1, n
         if(present(config_type_index)) then
            ii = config_type_index(i)
         else
            ii = i
         endif

         if(do_use_actual_gpcov) then
            k_self(i) = gpCoordinates_Covariance(this, i_x = ii, j_x = ii, normalise = .false.)
         else
         if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
         elseif(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
            if( zeta_int .feq. this%zeta ) then
               k_self(i) = dot_product( this%x(:,ii), this%x(:,ii) )**zeta_int
            else
               k_self(i) = dot_product( this%x(:,ii), this%x(:,ii) )**this%zeta
            endif
         elseif( this%covariance_type == COVARIANCE_ARD_SE ) then
            k_self(i) = 0.0_dp
            do i_p = 1, this%n_permutations
               !xI_xJ = (this%x(this%permutations(:,i_p),i) - this%x(:,j)) / 4.0_dp
               k_self(i) = k_self(i) + exp( -0.5_dp * sum((this%x(this%permutations(:,i_p),ii) - this%x(:,ii))**2) / 16.0_dp )
            enddo
         elseif( this%covariance_type == COVARIANCE_PP ) then
            k_self(i) = 0.0_dp
            do i_p = 1, this%n_permutations
               !xI_xJ = (this%x(this%permutations(:,i_p),i) - this%x(:,j)) / 4.0_dp
               k_self(i) = k_self(i) + covariancePP( sqrt( sum((this%x(this%permutations(:,i_p),ii) - this%x(:,ii))**2) ) / 4.0_dp, PP_Q, this%d)
            enddo
         endif
         endif
      enddo

      do j = 1, m-1
 
         if(present(config_type_index)) then
            jj = config_type_index(index_out(j))
         else
            jj = index_out(j)
         endif

!$omp parallel do default(none) shared(n,this,k_mn,jj,j,LA_k_mm,k_mm_k_m,score,k_self,config_type_index,index_out,do_use_actual_gpcov,zeta_int) private(i,i_p,ii)
         do i = 1, n

            if(present(config_type_index)) then
               ii = config_type_index(i)
            else
               ii = i
            endif
            if(do_use_actual_gpcov) then
               k_mn(j,i) = gpCoordinates_Covariance(this, i_x = ii, j_x = jj, normalise = .false.)
            else
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            elseif(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
               if( zeta_int .feq. this%zeta ) then
                  k_mn(j,i) = dot_product( this%x(:,ii), this%x(:,jj) )**zeta_int
               else
                  k_mn(j,i) = dot_product( this%x(:,ii), this%x(:,jj) )**this%zeta
               endif
            elseif( this%covariance_type == COVARIANCE_ARD_SE ) then
               k_mn(j,i) = 0.0_dp
               do i_p = 1, this%n_permutations
                  !xI_xJ = (this%x(this%permutations(:,i_p),i) - this%x(:,j)) / 4.0_dp
                  k_mn(j,i) = k_mn(j,i) + exp( -0.5_dp * sum((this%x(this%permutations(:,i_p),ii) - this%x(:,jj))**2) / 16.0_dp )
               enddo
            elseif( this%covariance_type == COVARIANCE_PP ) then
               k_mn(j,i) = 0.0_dp
               do i_p = 1, this%n_permutations
                  !xI_xJ = (this%x(this%permutations(:,i_p),i) - this%x(:,j)) / 4.0_dp
                  k_mn(j,i) = k_mn(j,i) + covariancePP( sqrt( sum((this%x(this%permutations(:,i_p),ii) - this%x(:,jj))**2) ) / 4.0_dp, PP_Q, this%d)
               enddo
            endif
            endif
            k_mn(j,i) = k_mn(j,i) / sqrt(k_self(i)*k_self(index_out(j)))
 
            call Matrix_Solve(LA_k_mm,k_mn(1:j,i),k_mm_k_m(1:j,i))
            score(i) = sum( k_mn(1:j,i) * k_mm_k_m(1:j,i) )
         enddo
 
         j_loc = minloc(score, mask=not_yet_added)
         jj = j_loc(1)
         index_out(j+1) = jj
         not_yet_added(jj) = .false.

         if(j == 1) then
            call print('Initial score: '//score)
         endif
         call print('Min score: '//minval(score))
 
         !k_mm(1:j_i,j_i+1) = k_mn(1:j_i,j)
         !k_mm(j_i+1,1:j_i) = k_mn(1:j_i,j)
         !k_mm(j_i+1,j_i+1) = 1.0_dp
         call LA_Matrix_Expand_Symmetrically(LA_k_mm,(/k_mn(1:j,jj),1.0_dp+1.0e-6_dp/))
         !call initialise(LA_k_mm,k_mm(1:j_i+1,1:j_i+1))
 
      enddo
      call print('Final score: '//score)
      call print('Min score: '//minval(score))

      deallocate(k_mn, score, k_mm_k_m, k_self, not_yet_added)
      !if(allocated(xI_xJ)) deallocate(xI_xJ)
      call finalise(LA_k_mm)
      call system_timer('sparse_covariance')
 
   endsubroutine sparse_covariance

end module gp_fit_module
