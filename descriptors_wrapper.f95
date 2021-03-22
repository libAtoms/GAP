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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X descriptors_wrapper subroutine
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine descriptors_wrapper_distances(N,lattice,symbol,coord,descriptor_str,descriptor_str_len, &
      calc_args_str,calc_args_str_len,i,fractional,previous_accepted,distances)

  use system_module
  use linearalgebra_module
  use dictionary_module
  use periodictable_module
  use atoms_types_module
  use connection_module
  use atoms_module

  use descriptors_module


  implicit none

  integer, intent(in) :: N
  real(dp), dimension(3,3), intent(inout) :: lattice
  character(len=3), dimension(N), intent(in) :: symbol
  integer, intent(in) :: descriptor_str_len
  character(len=descriptor_str_len) :: descriptor_str
  integer, intent(in) :: calc_args_str_len
  character(len=calc_args_str_len) :: calc_args_str
  real(dp), dimension(3,N), intent(in) :: coord
  integer, intent(in) :: i
  logical, intent(in) :: fractional, previous_accepted
  real(dp), dimension(N,N), intent(out) :: distances
  
  type(atoms), save           :: at
  type(Connection), save      :: at_connect_last_accepted, at_connect_previous
  type(descriptor), save      :: desc
  type(descriptor_data)       :: desc_data
  real(dp), dimension(:,:), allocatable, save :: desc_array_last_accepted, distances_in_last_accepted, desc_array_previous, distances_in_previous
  logical, dimension(:), pointer :: desc_mask

  integer, save :: d
  integer :: j, k, l, n_i

  logical, save :: first_run = .true.
  logical :: recalculate

  recalculate = .false.

  if( first_run ) then
     call system_initialise(verbosity=PRINT_SILENT)
     call initialise(desc,trim(descriptor_str))
     call initialise(at,N,lattice)
     call add_property(at,'desc_mask',.true.,ptr=desc_mask)

     d = descriptor_dimensions(desc)
     allocate(desc_array_previous(d,N), desc_array_last_accepted(d,N))
     allocate(distances_in_previous(N,N), distances_in_last_accepted(N,N))

     recalculate = .true.
  endif
  
  if( .not. first_run .and. (N /= at%N) ) then
     call finalise(at)
     call initialise(at,N,lattice)
     call add_property(at,'desc_mask',.true.,ptr=desc_mask)

     if(allocated(desc_array_previous)) deallocate(desc_array_previous)
     allocate(desc_array_previous(d,N))
     if(allocated(desc_array_last_accepted)) deallocate(desc_array_last_accepted)
     allocate(desc_array_last_accepted(d,N))

     if(allocated(distances_in_previous)) deallocate(distances_in_previous)
     allocate(distances_in_previous(N,N))
     if(allocated(distances_in_last_accepted)) deallocate(distances_in_last_accepted)
     allocate(distances_in_last_accepted(N,N))

     recalculate = .true.
  endif

  if( .not. first_run ) then
     if( previous_accepted ) then
        at_connect_last_accepted = at_connect_previous
        desc_array_last_accepted = desc_array_previous
        distances_in_last_accepted = distances_in_previous
     else
        at_connect_previous = at_connect_last_accepted
        desc_array_previous = desc_array_last_accepted
        distances_in_previous = distances_in_last_accepted
     endif
  endif

  if( at%lattice .fne. lattice ) then 
     call set_lattice(at,lattice, scale_positions=.false.)
     recalculate = .true.
  endif
  
  do k = 1, at%N
     at%Z(k) = atomic_number_from_symbol(symbol(k))
  enddo 

  if( i > 0 .and. previous_accepted .and. .not. recalculate ) then
     if( fractional ) then
        at%pos(:,i) = matmul(at%lattice,coord(:,i))
     else
        at%pos(:,i) = coord(:,i)
     endif
  else
     if( fractional ) then
        at%pos = matmul(at%lattice,coord)
     else
        at%pos = coord
     endif
  endif

  call set_cutoff(at,cutoff(desc)+0.5_dp)
  call calc_connect(at)

  if( .not. assign_pointer(at,'desc_mask',desc_mask) ) call system_abort("descriptors_wrapper: could not assign pointer desc_mask")

  if( i > 0 .and. .not. recalculate ) then

     if( i > at%N ) call system_abort("descriptors_wrapper: argument i = "//i//" greater than number of atoms "//at%N)

     desc_mask = .false.

     desc_mask(i) = .true.

     if( at_connect_previous%initialised ) then
        do n_i = 1, n_neighbours(at,i,alt_connect=at_connect_previous)
           desc_mask(neighbour(at,i,n_i,alt_connect=at_connect_previous)) = .true.
        enddo
     endif

     do n_i = 1, n_neighbours(at,i)
        desc_mask(neighbour(at,i,n_i)) = .true.
     enddo

     call calc(desc,at,desc_data,do_descriptor=.true.,do_grad_descriptor=.false.,args_str="atom_mask_name=desc_mask "//trim(calc_args_str))

     do k = 1, count(desc_mask)
        j = desc_data%x(k)%ci(1)
        desc_array_previous(:,j) = desc_data%x(k)%data(:)
     enddo

     do k = 1, count(desc_mask)
        j = desc_data%x(k)%ci(1)
        do l = 1, at%N
           distances_in_previous(l,j) = sum( desc_array_previous(:,l) * desc_array_previous(:,j) )
           distances_in_previous(j,l) = distances_in_previous(l,j)
        enddo
     enddo
     desc_mask = .true.
     at_connect_previous = at%connect
  else
     call calc(desc,at,desc_data,do_descriptor=.true.,do_grad_descriptor=.false.,args_str=trim(calc_args_str))
     do j = 1, at%N
        desc_array_previous(:,j) = desc_data%x(j)%data
     enddo

     distances_in_previous = matmul(transpose(desc_array_previous),desc_array_previous)

     at_connect_previous = at%connect
  endif

  distances = -log(distances_in_previous)

  call finalise(desc_data)

  if( first_run ) then
     at_connect_last_accepted = at_connect_previous
     desc_array_last_accepted = desc_array_previous
     distances_in_last_accepted = distances_in_previous
  endif
  first_run = .false.

endsubroutine descriptors_wrapper_distances

module descriptors_wrapper_module

use system_module
use periodictable_module, only : atomic_number_from_symbol
use atoms_module
use linearalgebra_module
use descriptors_module, only : descriptor, initialise, finalise, cutoff, calc, descriptor_sizes, descriptor_dimensions

implicit none

#ifdef HAVE_GAP
type(descriptor), save :: desc
#endif

logical :: first_run = .true.

contains

   subroutine descriptors_wrapper_initialise(descriptor_str)

      character(len=*) :: descriptor_str

#ifdef HAVE_GAP
      if( first_run ) then
         call system_initialise(verbosity=PRINT_SILENT)
         call initialise(desc,trim(descriptor_str))
      else
         call finalise(desc)
         call initialise(desc,trim(descriptor_str))
      endif
      first_run = .false.
#endif
   endsubroutine descriptors_wrapper_initialise

   subroutine descriptors_wrapper_initialise_C(descriptor_str,n_descriptor_str) bind(c)

      use iso_c_binding, only: c_int, c_char

      character(kind=c_char), intent(in) :: descriptor_str(n_descriptor_str)
      integer(kind=c_int), intent(in) :: n_descriptor_str

      call descriptors_wrapper_initialise(trim(a2s(descriptor_str)))

   endsubroutine descriptors_wrapper_initialise_C

   function descriptors_wrapper_dimensions()
      integer :: descriptors_wrapper_dimensions
      
      if(.not. first_run) then
         descriptors_wrapper_dimensions = descriptor_dimensions(desc)
      else
         call system_abort("descriptors_wrapper_dimensions: initialise with calling descriptors_wrapper_initialise() first.")
      endif

   endfunction descriptors_wrapper_dimensions

   function descriptors_wrapper_dimensions_C() bind(c)
      use iso_c_binding, only: c_int
      integer(kind=c_int) :: descriptors_wrapper_dimensions_C
      
      descriptors_wrapper_dimensions_C = descriptors_wrapper_dimensions()
   endfunction descriptors_wrapper_dimensions_C

   function descriptors_wrapper_size(N,lattice,symbol,coord,fractional)
      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional

      integer :: descriptors_wrapper_size

      type(atoms), save           :: at
      integer :: n_descriptors,n_cross

      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call descriptor_sizes(desc,at,n_descriptors,n_cross)

      descriptors_wrapper_size = n_descriptors

   endfunction descriptors_wrapper_size

   function descriptors_wrapper_size_C(N,lattice,symbol,coord,fractional) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional

      integer(kind=c_int) :: descriptors_wrapper_size_C

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      descriptors_wrapper_size_C = descriptors_wrapper_size(N,lattice,my_symbol,coord,my_fractional)

   endfunction descriptors_wrapper_size_C

   function descriptors_wrapper_gradient_size(N,lattice,symbol,coord,fractional)
      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional

      integer :: descriptors_wrapper_gradient_size

      type(atoms), save           :: at
      integer :: n_descriptors,n_cross

      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call descriptor_sizes(desc,at,n_descriptors,n_cross)

      descriptors_wrapper_gradient_size = n_cross

   endfunction descriptors_wrapper_gradient_size

   function descriptors_wrapper_gradient_size_C(N,lattice,symbol,coord,fractional) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional

      integer(kind=c_int) :: descriptors_wrapper_gradient_size_C

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      descriptors_wrapper_gradient_size_C = descriptors_wrapper_gradient_size(N,lattice,my_symbol,coord,my_fractional)

   endfunction descriptors_wrapper_gradient_size_C

   subroutine descriptors_wrapper_both_sizes(N,lattice,symbol,coord,fractional,n_descriptors,n_cross)
      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional
      integer, intent(out) :: n_descriptors, n_cross

      type(atoms), save           :: at

      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call descriptor_sizes(desc,at,n_descriptors,n_cross)

   endsubroutine descriptors_wrapper_both_sizes

   subroutine descriptors_wrapper_both_sizes_C(N,lattice,symbol,coord,fractional,n_descriptors,n_cross) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional
      integer(kind=c_int), intent(out) :: n_descriptors, n_cross

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      call descriptors_wrapper_both_sizes(N,lattice,my_symbol,coord,my_fractional,n_descriptors,n_cross)

   endsubroutine descriptors_wrapper_both_sizes_C

   subroutine descriptors_wrapper_array(N,lattice,symbol,coord,fractional,descriptor_array,d_descriptor, n_descriptor)

      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional
      integer, intent(in) :: d_descriptor, n_descriptor
      real(dp), dimension(d_descriptor,n_descriptor), intent(out):: descriptor_array
      
      type(atoms), save           :: at

      if( first_run ) then
          call system_abort("descriptors_wrapper_array: initialise with calling descriptors_wrapper_initialise() first.")
      endif
      
      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call calc(desc,at,descriptor_array)

   endsubroutine descriptors_wrapper_array

   subroutine descriptors_wrapper_array_C(N,lattice,symbol,coord,fractional,descriptor_array,d_descriptor, n_descriptor) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional
      integer(kind=c_int), intent(in) :: d_descriptor, n_descriptor
      real(kind=c_double), dimension(d_descriptor,n_descriptor), intent(out):: descriptor_array

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      call descriptors_wrapper_array(N,lattice,my_symbol,coord,my_fractional,descriptor_array,d_descriptor,n_descriptor)

   endsubroutine descriptors_wrapper_array_C

   subroutine descriptors_wrapper_gradient_array(N,lattice,symbol,coord,fractional,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_cross)

      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional
      integer, intent(in) :: d_descriptor, n_cross
      real(dp), dimension(d_descriptor,3,n_cross), intent(out):: grad_descriptor_array
      integer, dimension(2,n_cross), intent(out):: grad_descriptor_index
      real(dp), dimension(3,n_cross), intent(out):: grad_descriptor_pos
      
      type(atoms), save           :: at

      if( first_run ) then
          call system_abort("descriptors_wrapper_gradient_array: initialise with calling descriptors_wrapper_initialise() first.")
      endif
      
      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call calc(desc,at,grad_descriptor_out=grad_descriptor_array,grad_descriptor_index=grad_descriptor_index,grad_descriptor_pos=grad_descriptor_pos)

   endsubroutine descriptors_wrapper_gradient_array

   subroutine descriptors_wrapper_gradient_array_C(N,lattice,symbol,coord,fractional,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_cross) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional
      integer(kind=c_int), intent(in) :: d_descriptor, n_cross
      real(kind=c_double), dimension(d_descriptor,3,n_cross), intent(out):: grad_descriptor_array
      integer(kind=c_int), dimension(2,n_cross), intent(out):: grad_descriptor_index
      real(kind=c_double), dimension(3,n_cross), intent(out):: grad_descriptor_pos

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      call descriptors_wrapper_gradient_array(N,lattice,my_symbol,coord,my_fractional,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_cross)

   endsubroutine descriptors_wrapper_gradient_array_C

   subroutine descriptors_wrapper_both_arrays(N,lattice,symbol,coord,fractional,descriptor_array,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_descriptor,n_cross)

      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional
      integer, intent(in) :: d_descriptor, n_descriptor, n_cross
      real(dp), dimension(d_descriptor,n_descriptor), intent(out):: descriptor_array
      real(dp), dimension(d_descriptor,3,n_cross), intent(out):: grad_descriptor_array
      integer, dimension(2,n_cross), intent(out):: grad_descriptor_index
      real(dp), dimension(3,n_cross), intent(out):: grad_descriptor_pos
      
      type(atoms), save           :: at

      if( first_run ) then
          call system_abort("descriptors_wrapper_both_arrays: initialise with calling descriptors_wrapper_initialise() first.")
      endif
      
      call copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      call calc(desc,at,descriptor_out=descriptor_array,grad_descriptor_out=grad_descriptor_array,grad_descriptor_index=grad_descriptor_index,grad_descriptor_pos=grad_descriptor_pos)

   endsubroutine descriptors_wrapper_both_arrays

   subroutine descriptors_wrapper_both_arrays_C(N,lattice,symbol,coord,fractional,descriptor_array,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_descriptor,n_cross) bind(c)

      use iso_c_binding, only: c_double, c_int, c_bool, c_char

      integer(kind=c_int), intent(in) :: N
      real(kind=c_double), dimension(3,3), intent(inout) :: lattice
      character(kind=c_char), dimension(3,N), intent(in) :: symbol
      real(kind=c_double), dimension(3,N), intent(in) :: coord
      logical(kind=c_bool), intent(in) :: fractional
      integer(kind=c_int), intent(in) :: d_descriptor, n_descriptor, n_cross
      real(kind=c_double), dimension(d_descriptor,n_descriptor), intent(out):: descriptor_array
      real(kind=c_double), dimension(d_descriptor,3,n_cross), intent(out):: grad_descriptor_array
      integer(kind=c_int), dimension(2,n_cross), intent(out):: grad_descriptor_index
      real(kind=c_double), dimension(3,n_cross), intent(out):: grad_descriptor_pos

      character(len=3), dimension(N) :: my_symbol
      integer :: i
      logical :: my_fractional

      do i = 1, N
         my_symbol(i) = a2s(symbol(:,i))
      enddo

      my_fractional = logical(fractional,kind=kind(my_fractional))
      call descriptors_wrapper_both_arrays(N,lattice,my_symbol,coord,my_fractional,descriptor_array,grad_descriptor_array,grad_descriptor_index,grad_descriptor_pos,d_descriptor,n_descriptor,n_cross)

   endsubroutine descriptors_wrapper_both_arrays_C

   subroutine copy_data_to_atoms(at,N,lattice,symbol,coord,fractional)
      type(atoms), intent(inout) :: at
      integer, intent(in) :: N
      real(dp), dimension(3,3), intent(inout) :: lattice
      character(len=3), dimension(N), intent(in) :: symbol
      real(dp), dimension(3,N), intent(in) :: coord
      logical, intent(in) :: fractional

      integer :: k

      if( N /= at%N ) then
         call finalise(at)
         call initialise(at,N,lattice)
      endif

      if( at%lattice .fne. lattice ) then 
         call set_lattice(at,lattice, scale_positions=.false.)
      endif
      
      do k = 1, at%N
         at%Z(k) = atomic_number_from_symbol(symbol(k))
      enddo 

      if( fractional ) then
         at%pos = matmul(at%lattice,coord)
      else
         at%pos = coord
      endif

      call set_cutoff(at,cutoff(desc)+0.5_dp)
      call calc_connect(at)

   endsubroutine copy_data_to_atoms

endmodule descriptors_wrapper_module
