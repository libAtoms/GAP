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

! HND X
! HND X
! HND X This file contains descriptors written by others, not A. P. Bartok and Gabor Csanyi
! HND X and the code here is owned by their respective authors, as indicated below. 
! HND X



!!!!!!!!
!!! Work of Wojciech Szlachta BEGINS here
!!!!!!!!
   subroutine bond_real_space_initialise(this,args_str,error)
      type(bond_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'bond_cutoff', '0.00', this%bond_cutoff, help_string="Bond cutoff for bond_real_space-type descriptors")
      call param_register(params, 'bond_transition_width', '0.00', this%bond_transition_width, help_string="Bond transition width for bond_real_space-type descriptors")
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for bond_real_space-type descriptors")
      call param_register(params, 'transition_width', '0.00', this%transition_width, help_string="Space transition width for bond_real_space-type descriptors")
      call param_register(params, 'atom_gaussian_width', '0.00', this%atom_sigma, help_string="Atom sigma for bond_real_space-type descriptors", altkey='atom_sigma')
      call param_register(params, 'max_neighbours', '0', this%max_neighbours, help_string="Maximum number of neighbours")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bond_real_space_initialise args_str')) then
         RAISE_ERROR("bond_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine bond_real_space_initialise

   subroutine bond_real_space_finalise(this,error)
      type(bond_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%bond_cutoff = 0.0_dp
      this%bond_transition_width = 0.0_dp
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%atom_sigma = 0.0_dp
      this%max_neighbours = 0

      this%initialised = .false.

   endsubroutine bond_real_space_finalise

   subroutine bond_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(atoms) :: at_copy
      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: n_descriptors, n_cross, i_desc, i, j, n, k, m, m_index, l, &
         ij_neighbours, n_index
      integer, dimension(3) :: shift_j, shift_k
      real(dp) :: r_ij, r_ijk
      real(dp) :: atom_i(3), atom_j(3), atom_k(3), bond(3), bond_len
      real(dp) :: atom_i_cross_atom_j(3), atom_i_normsq_min_atom_j_normsq
      real(dp), allocatable :: r(:,:), z(:), c(:)
      real(dp) :: self_overlap
      real(dp), allocatable :: dr(:,:,:,:), dz(:,:,:), dc(:,:,:)
      real(dp), allocatable :: dself_overlap(:,:)
      integer, allocatable :: ii(:)
      real(dp), allocatable :: pos(:,:)
      real(dp) :: r_m_cross_r_l(3)

      INIT_ERROR(error)

      call system_timer('bond_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_calc: descriptor object not initialised", error)
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

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='bond_real_space_calc args_str')) then
            RAISE_ERROR("bond_real_space_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("bond_real_space_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("bond_real_space_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0

      do i = 1, at%N
         do n = 1, n_neighbours(at, i)
            j = neighbour(at, i, n, shift=shift_j, distance=r_ij, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            i_desc = i_desc + 1

            atom_i = at%pos(:,i)
            atom_j = at%pos(:,j) + matmul(at%lattice, shift_j)

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (atom_i + atom_j), 1)
            call calc_connect(at_copy)

            ij_neighbours = 0

            do m = 1, n_neighbours(at_copy, at%N + 1)
               k = neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               ij_neighbours = ij_neighbours + 1
            enddo

            if(ij_neighbours > this%max_neighbours) then
               RAISE_ERROR("bond_real_space_calc: number of neighbours exceeds max_neighbours", error)
            endif

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               allocate(r(3,ij_neighbours), z(ij_neighbours), c(ij_neighbours))
               allocate(ii(ij_neighbours), pos(3,ij_neighbours))

               r = 0.0_dp
               z = 0.0_dp
               c = 0.0_dp
               self_overlap = 0.0_dp
               bond = atom_i - atom_j
               bond_len = norm(bond)
               atom_i_cross_atom_j = atom_i .cross. atom_j
               atom_i_normsq_min_atom_j_normsq = normsq(atom_i) - normsq(atom_j)
               ii = 0
               pos = 0.0_dp

               if(my_do_grad_descriptor) then
                  allocate(dr(3,ij_neighbours,3,ij_neighbours), dz(ij_neighbours,3,ij_neighbours), dc(ij_neighbours,3,ij_neighbours), dself_overlap(3,ij_neighbours))

                  dr = 0.0_dp
                  dz = 0.0_dp
                  dc = 0.0_dp
                  dself_overlap = 0.0_dp
               endif

               m_index = 2

               do m = 1, n_neighbours(at_copy, at%N + 1)
                  k = neighbour(at_copy, at%N + 1, m, shift=shift_k, distance=r_ijk, max_dist=this%cutoff)

                  if(k == 0) cycle

                  if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

                  atom_k = at_copy%pos(:,k) + matmul(at_copy%lattice, shift_k)

                  if(atom_k .feq. atom_i) then
                     ! r remains zero
                     z(1) = 0.5_dp * bond_len
                     c(1) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(1) = k
                     pos(:,1) = atom_k

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(1,:,1) = 0.5_dp * bond / bond_len
                        dz(1,:,2) = - dz(1,:,1)
                        dc(1,:,1) = 0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(1,:,2) = - dc(1,:,1)
                     endif
                  elseif(atom_k .feq. atom_j) then
                     ! r remain zero
                     z(2) = -0.5_dp * bond_len
                     c(2) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(2) = k
                     pos(:,2) = atom_k

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(2,:,1) = -0.5_dp * bond / bond_len
                        dz(2,:,2) = - dz(2,:,1)
                        dc(2,:,1) = -0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(2,:,2) = - dc(2,:,1)
                     endif
                  else
                     m_index = m_index + 1

                     r(:,m_index) = ((atom_k .cross. bond) + atom_i_cross_atom_j) / bond_len
                     z(m_index) = ((atom_k .dot. bond) - 0.5_dp * atom_i_normsq_min_atom_j_normsq) / bond_len
                     c(m_index) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     ii(m_index) = k
                     pos(:,m_index) = atom_k

                     if(my_do_grad_descriptor) then
                        dr(:,m_index,1,1) = ((/ 0.0_dp, atom_k(3) - atom_j(3), atom_j(2) - atom_k(2) /) / bond_len) - (r(:,m_index) * bond(1) / bond_len**2)
                        dr(:,m_index,2,1) = ((/ atom_j(3) - atom_k(3), 0.0_dp, atom_k(1) - atom_j(1) /) / bond_len) - (r(:,m_index) * bond(2) / bond_len**2)
                        dr(:,m_index,3,1) = ((/ atom_k(2) - atom_j(2), atom_j(1) - atom_k(1), 0.0_dp /) / bond_len) - (r(:,m_index) * bond(3) / bond_len**2)
                        dz(m_index,:,1) = ((atom_k - atom_i) / bond_len) - (z(m_index) * bond / bond_len**2)
                        dc(m_index,:,1) = -0.5_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * (atom_k - at_copy%pos(:,at%N + 1)) / r_ijk

                        dr(:,m_index,1,2) = - dr(:,m_index,1,1) + ((/ 0.0_dp, bond(3), - bond(2) /) / bond_len)
                        dr(:,m_index,2,2) = - dr(:,m_index,2,1) + ((/ - bond(3), 0.0_dp, bond(1) /) / bond_len)
                        dr(:,m_index,3,2) = - dr(:,m_index,3,1) + ((/ bond(2), - bond(1), 0.0_dp /) / bond_len)
                        dz(m_index,:,2) = - dz(m_index,:,1) - (bond / bond_len)
                        dc(m_index,:,2) = dc(m_index,:,1)

                        dr(:,m_index,1,m_index) = (/ 0.0_dp, - bond(3), bond(2) /) / bond_len
                        dr(:,m_index,2,m_index) = (/ bond(3), 0.0_dp, - bond(1) /) / bond_len
                        dr(:,m_index,3,m_index) = (/ - bond(2), bond(1), 0.0_dp /) / bond_len
                        dz(m_index,:,m_index) = bond / bond_len
                        dc(m_index,:,m_index) = -2.0_dp * dc(m_index,:,1)
                     endif
                  endif
               enddo
            endif

            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(2 + (1 + 2 * this%max_neighbours) * this%max_neighbours))
               allocate(descriptor_out%x(i_desc)%ci(n_index))

               descriptor_out%x(i_desc)%data = 0.0_dp

               do m = 1, ij_neighbours
                  self_overlap = self_overlap + c(m)**2

                  if(m == 1) then
                     descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m)) = z(m)
                  elseif(m == 2) then
                     descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m)) = z(m)

                     self_overlap = self_overlap + 2.0_dp * c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 )
                  else
                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = normsq(r(:,m))
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = z(m)
                        else
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = r(:,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = ((r(:,m) .cross. r(:,l)) .dot. bond) / bond_len
                        endif

                        if(l < m) then
                           self_overlap = self_overlap + 2.0_dp * c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                        endif
                     enddo
                  endif
               enddo

               descriptor_out%x(i_desc)%data(1) = real(ij_neighbours, dp)

               descriptor_out%x(i_desc)%data(2) = self_overlap

               descriptor_out%x(i_desc)%data(3:ij_neighbours + 2) = c

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij, this%bond_cutoff, this%bond_transition_width)

               descriptor_out%x(i_desc)%ci(:) = (/ i, j /)
               descriptor_out%x(i_desc)%has_data = .true.
            endif

            if(my_do_grad_descriptor) then
               allocate(descriptor_out%x(i_desc)%grad_data(2 + (1 + 2 * ij_neighbours) * ij_neighbours,3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%ii(ij_neighbours))
               allocate(descriptor_out%x(i_desc)%pos(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%has_grad_data(ij_neighbours))

               descriptor_out%x(i_desc)%grad_data = 0.0_dp

               do m = 1, ij_neighbours
                  dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * c(m) * dc(m,:,1)
                  dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * c(m) * dc(m,:,2)

                  if(m == 1) then
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,1) = dz(m,:,1)
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,2) = dz(m,:,2)
                  elseif(m == 2) then
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,1) = dz(m,:,1)
                     descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,2) = dz(m,:,2)

                     dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + 2.0_dp * c(m) * dc(m - 1,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             * (pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2
                     dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + 2.0_dp * c(m) * dc(m - 1,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             + c(m) * c(m - 1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,m - 1)) / this%atom_sigma**2 ) \
                                                             * (pos(:,m - 1) - pos(:,m)) / this%atom_sigma**2
                  else
                     dself_overlap(:,m) = dself_overlap(:,m) + 2.0_dp * c(m) * dc(m,:,m)

                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,1) = dz(m,:,1)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,2) = dz(m,:,2)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,m))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,m) = dz(m,:,m)
                        else
                           r_m_cross_r_l = r(:,m) .cross. r(:,l)

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = (dr(:,m,1,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = (dr(:,m,2,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = (dr(:,m,3,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,1))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,1) = ((((dr(:,m,1,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 1.0_dp, 0.0_dp, 0.0_dp /) - (bond * bond(1) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,1) = ((((dr(:,m,2,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 1.0_dp, 0.0_dp /) - (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,1) = ((((dr(:,m,3,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, 1.0_dp /) - (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = (dr(:,m,1,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = (dr(:,m,2,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = (dr(:,m,3,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,2))
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,2) = ((((dr(:,m,1,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ -1.0_dp, 0.0_dp, 0.0_dp /) + (bond * bond(1) / bond_len**2)))) / bond_len 
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,2) = ((((dr(:,m,2,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, -1.0_dp, 0.0_dp /) + (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,2) = ((((dr(:,m,3,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, -1.0_dp /) + (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = dr(:,m,1,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = dr(:,m,2,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = dr(:,m,3,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,m) = ((dr(:,m,1,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,m) = ((dr(:,m,2,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,m) = ((dr(:,m,3,m) .cross. r(:,l)) .dot. bond) / bond_len

                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,l) = r(:,m) .dot. dr(:,l,1,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,l) = r(:,m) .dot. dr(:,l,2,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,l) = r(:,m) .dot. dr(:,l,3,l)
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,l) = ((r(:,m) .cross. dr(:,l,1,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,l) = ((r(:,m) .cross. dr(:,l,2,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(2 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,l) = ((r(:,m) .cross. dr(:,l,3,l)) .dot. bond) / bond_len
                        endif

                        if(l < m) then
                           dself_overlap(:,m) = dself_overlap(:,m) + 2.0_dp * dc(m,:,m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                   + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                   * (pos(:,l) - pos(:,m)) / this%atom_sigma**2

                           if(l == 1) then
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                           elseif(l == 2) then
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                           else
                              dself_overlap(:,1) = dself_overlap(:,1) + 2.0_dp * dc(m,:,1) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,1) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,2) = dself_overlap(:,2) + 2.0_dp * dc(m,:,2) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + 2.0_dp * c(m) * dc(l,:,2) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 )
                              dself_overlap(:,l) = dself_overlap(:,l) + 2.0_dp * c(m) * dc(l,:,l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      + c(m) * c(l) * exp( -0.25_dp * normsq(pos(:,m) - pos(:,l)) / this%atom_sigma**2 ) \
                                                                      * (pos(:,m) - pos(:,l)) / this%atom_sigma**2
                           endif
                        endif
                     enddo
                  endif
               enddo

               !descriptor_out%x(i_desc)%grad_data(1,:,:) = 0.0_dp

               descriptor_out%x(i_desc)%grad_data(2,:,:) = dself_overlap

               descriptor_out%x(i_desc)%grad_data(3:ij_neighbours + 2,:,:) = dc

               descriptor_out%x(i_desc)%ii = ii
               descriptor_out%x(i_desc)%pos = pos

               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dcoordination_function(r_ij, this%bond_cutoff, this%bond_transition_width) * bond / r_ij
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)

               descriptor_out%x(i_desc)%has_grad_data = .true.
            endif

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               deallocate(r, z, c)
               deallocate(ii, pos)

               if(my_do_grad_descriptor) then
                  deallocate(dr, dz, dc, dself_overlap)
               endif
            endif

            call finalise(at_copy)
         enddo
      enddo

      call system_timer('bond_real_space_calc')

   endsubroutine bond_real_space_calc

   function bond_real_space_dimensions(this,error) result(i)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 2 + (1 + 2 * this%max_neighbours) * this%max_neighbours

   endfunction bond_real_space_dimensions

  function bond_real_space_cutoff(this,error)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bond_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_cutoff: descriptor object not initialised", error)
      endif

      bond_real_space_cutoff = max(this%cutoff, this%bond_cutoff)

   endfunction bond_real_space_cutoff

   subroutine bond_real_space_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      type(atoms) :: at_copy
      integer :: i, j, k, n, m, shift_j(3)

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + n_neighbours(at, i, max_dist=this%bond_cutoff)

         do n = 1, n_neighbours(at, i)
            j = neighbour(at, i, n, shift=shift_j, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (at%pos(:,i) + at%pos(:,j) + matmul(at%lattice,shift_j)), 1)
            call calc_connect(at_copy)

            do m = 1, n_neighbours(at_copy, at%N + 1)
               k = neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               n_cross = n_cross + 1
            enddo

            call finalise(at_copy)
         enddo
      enddo

      if( present(n_index) ) n_index = 2

   endsubroutine bond_real_space_sizes

!!!!!!!!
!!! Work of Wojciech Szlachta ENDS here
!!!!!!!!


!!!!!!!!
!!! Work of Alan Nichol and S. T. John BEGINS here
!!!!!!!!

      subroutine descriptor_general_monomer_nmer_MPI_setup(this,at,mpi,mpi_mask,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(MPI_Context), intent(in) :: mpi
      logical, dimension(:), intent(out) :: mpi_mask
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_index
      logical, dimension(at%N) :: associated_to_monomer
      integer :: n_monomer

      integer :: i

      INIT_ERROR(error)

      associated_to_monomer = .false.
      select case(this%descriptor_type)
      case(DT_GENERAL_MONOMER)
         call find_general_monomer(at, monomer_index, &
              this%descriptor_general_monomer%signature, associated_to_monomer, &
              this%descriptor_general_monomer%cutoff, &
              this%descriptor_general_monomer%atom_ordercheck, error)
      case(DT_GENERAL_DIMER)
         call find_general_monomer(at, monomer_index, &
              this%descriptor_general_dimer%signature_one, associated_to_monomer, &
              this%descriptor_general_dimer%monomer_one_cutoff, &
              this%descriptor_general_dimer%atom_ordercheck, error)
      case(DT_GENERAL_TRIMER)
         call find_general_monomer(at, monomer_index, &
              this%descriptor_general_trimer%signature_one, associated_to_monomer, &
              this%descriptor_general_trimer%monomer_one_cutoff,&
              this%descriptor_general_trimer%atom_ordercheck,error)
      case(DT_COM_DIMER)
         call find_general_monomer(at, monomer_index, &
              this%descriptor_com_dimer%signature_one, associated_to_monomer,&
              this%descriptor_com_dimer%monomer_one_cutoff,&
              this%descriptor_com_dimer%atom_ordercheck,error)
      case default
         RAISE_ERROR("descriptor_general_monomer_nmer_MPI_setup: descriptor type "//this%descriptor_type//" not recognised.",error)
      endselect

      n_monomer = size(monomer_index,2)

      mpi_mask = .false.
      do i = 1, n_monomer ! for dimer, trimer this is the first monomer (signature_one)
         if( mod(i-1, mpi%n_procs) == mpi%my_proc ) then
            mpi_mask(monomer_index(:,i)) = .true.
         endif
      enddo

      deallocate(monomer_index)

   endsubroutine descriptor_general_monomer_nmer_MPI_setup

      subroutine AN_monomer_initialise(this,args_str,error)
      type(AN_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for AN_monomer-type descriptors")
      call param_register(params, 'atomic_number', '1', this%atomic_number, help_string="Atomic number in AN_monomer-type descriptors")
      call param_register(params, 'N', '4', this%N, help_string="Number of atoms in cluster")
      call param_register(params, 'do_atomic', 'T', this%do_atomic, help_string="Descriptors are cluster based or atom-based")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='AN_monomer_initialise args_str')) then
         RAISE_ERROR("AN_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine AN_monomer_initialise

   subroutine AN_monomer_finalise(this,error)
      type(AN_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%atomic_number = 0
      this%N = 0

      this%do_atomic = .false.
      this%initialised = .false.

   endsubroutine AN_monomer_finalise

   subroutine general_monomer_initialise(this,args_str,error)
      type(general_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_string
      character(len=STRING_LENGTH_SHORT), dimension(99) :: signature_fields
      integer, optional, intent(out) :: error
      integer :: i,n_atoms,j


      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for general_monomer-type descriptors")
      call param_register(params, 'signature', PARAM_MANDATORY, signature_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer")
      call param_register(params, 'power', '1.0', this%power, help_string="Power of distances to be used in the kernel")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_monomer_initialise args_str')) then
         RAISE_ERROR("general_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call split_string(signature_string,' ','{}',signature_fields(:),n_atoms,matching=.true.)
      allocate(this%signature(n_atoms))

      do i=1,n_atoms
        this%signature(i) = string_to_int(signature_fields(i))
      end do

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)

      this%initialised = .true.

   endsubroutine general_monomer_initialise

   subroutine general_monomer_finalise(this,error)
      type(general_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%power = 1.0_dp
      if(allocated(this%signature)) deallocate(this%signature)

      this%initialised = .false.

   endsubroutine general_monomer_finalise

   subroutine com_dimer_initialise(this,args_str,error)
      type(com_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string
      character(len=STRING_LENGTH_SHORT), dimension(99) :: signature_one_fields, signature_two_fields
      integer, optional, intent(out) :: error
      integer :: i, n_atoms_one, n_atoms_two

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff(intermolecular) for com_dimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, help_string="Cutoff(mono1) for com_dimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, help_string="Cutoff(mono2) for com_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for com_dimer-type descriptors")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'mpifind', 'false', this%mpifind, help_string="Use find_monomer_pairs_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='com_dimer_initialise args_str')) then
         RAISE_ERROR("com_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%transfer_parameters, args_str, error)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do

      this%monomers_identical=.False.
      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%monomers_identical = .True.
         end if
      end if

      this%initialised = .true.
   endsubroutine com_dimer_initialise

   subroutine com_dimer_finalise(this,error)
      type(com_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%use_smooth_cutoff = .false.
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)

      this%initialised = .false.

   endsubroutine com_dimer_finalise

   subroutine general_dimer_initialise(this,args_str,error)
      type(general_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string
      character(len=STRING_LENGTH_SHORT), dimension(99) :: signature_one_fields, signature_two_fields
      integer, optional, intent(out) :: error
      integer :: i,j, n_atoms_one, n_atoms_two, dimer_size, start, finish, d
      logical, dimension(:,:), allocatable :: intermolecular
      integer, dimension(:), allocatable :: signature

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, &
           help_string="Cutoff(intermolecular) for general_dimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, &
           help_string="Cutoff(mono1) for general_dimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, &
           help_string="Cutoff(mono2) for general_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, &
           help_string="Width of smooth cutoff region for general_dimer-type descriptors")
      call param_register(params, 'internal_swaps_only', 'true', this%internal_swaps_only, &
           help_string="F: energies will be symmetrised over swaps of nuclei between monomers")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, &
           help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'double_count', 'false', this%double_count, &
           help_string="T: double count when constructing the dimers, for compatibility with water dimer descriptor, default False")
      call param_register(params, 'strict', 'true', this%strict, &
           help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'strict_mask', 'true', this%strict_mask, &
           help_string="Raise error if atom mask includes only part of a monomer")
      call param_register(params, 'use_com', 'false', this%use_com, &
           help_string="Use COM instead of COG")
      call param_register(params, 'mpifind', 'false', this%mpifind, &
           help_string="Use find_monomer_pairs_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, &
           help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, &
           help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'power', '1.00', this%power, &
           help_string="Power of interatomic distances to be used in the kernel.")
      call param_register(params, 'dist_shift', '0.0', this%dist_shift, &
           help_string="Distance shift for inverse distance descriptors.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_dimer_initialise args_str')) then
         RAISE_ERROR("general_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%transfer_parameters, args_str, error)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do

      this%monomers_identical=.False.
      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%monomers_identical = .True.
         end if
      end if

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature_one,signature_two=this%signature_two,internal_swaps_only=this%internal_swaps_only,error=error)

      dimer_size=n_atoms_one + n_atoms_two
      d=dimer_size*(dimer_size-1)/2

      allocate(signature(dimer_size))
      allocate(intermolecular(dimer_size,dimer_size))
      allocate(this%is_intermolecular(d))
      allocate(this%cutoff_contributor(d))
      allocate(this%component_atoms(d,2))

      signature(1:n_atoms_one) = this%signature_one
      signature(1+n_atoms_one:dimer_size) = this%signature_two
      intermolecular = .false.
      this%cutoff_contributor=.false.

      do i=1,n_atoms_one
        do j=1+n_atoms_one,dimer_size
          intermolecular(i,j)=.true.
        end do
      end do

      start = 0
      finish=dimer_size-1
      do i=1,dimer_size
        do j=1,finish-start
          this%is_intermolecular(start+j) = intermolecular(i,i+j)
          this%component_atoms(start+j,:) = (/ i, i+j /)
        end do
        start = finish
        finish=finish + dimer_size-i-1
      end do


      do i=1,d
        if (this%is_intermolecular(i)) then
          if (.not. signature(this%component_atoms(i,1))==1 ) then
            if (.not. signature(this%component_atoms(i,2))==1 ) then
              this%cutoff_contributor(i)=.true.
            end if
          end if
        end if
      end do

      this%initialised = .true.

      deallocate(signature)
      deallocate(intermolecular)
   endsubroutine general_dimer_initialise

   subroutine general_dimer_finalise(this,error)
      type(general_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%internal_swaps_only = .true.
      this%use_smooth_cutoff = .false.
      this%power = 1.0_dp
      this%dist_shift = 0.0_dp
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)
      if(allocated(this%is_intermolecular)) deallocate(this%is_intermolecular)
      if(allocated(this%component_atoms)) deallocate(this%component_atoms)
      if(allocated(this%cutoff_contributor)) deallocate(this%cutoff_contributor)

      this%initialised = .false.

   endsubroutine general_dimer_finalise

   subroutine general_trimer_initialise(this,args_str,error)
      type(general_trimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_one_string, signature_two_string, signature_three_string
      character(len=STRING_LENGTH_SHORT), dimension(99) :: signature_one_fields, signature_two_fields, signature_three_fields
      integer, optional, intent(out) :: error
      integer :: i,j, n_atoms_one, n_atoms_two, n_atoms_three, trimer_size, start, finish,d
      logical, dimension(:,:), allocatable :: intermolecular
      integer, dimension(:), allocatable :: signature

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, &
           help_string="Cutoff(intermolecular) for general_trimer-type descriptors")
      call param_register(params, 'monomer_one_cutoff', '0.00', this%monomer_one_cutoff, &
           help_string="Cutoff(mono1) for general_trimer-type descriptors")
      call param_register(params, 'monomer_two_cutoff', '0.00', this%monomer_two_cutoff, &
           help_string="Cutoff(mono2) for general_trimer-type descriptors")
      call param_register(params, 'monomer_three_cutoff', '0.00', this%monomer_three_cutoff, &
           help_string="Cutoff(mono3) for general_trimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, &
           help_string="Width of smooth cutoff region for general_trimer-type descriptors")
      call param_register(params, 'internal_swaps_only', 'true', this%internal_swaps_only, &
           help_string="F: energies will be symmetrised over swaps of nuclei between monomers")
      call param_register(params, 'atom_ordercheck', 'true', this%atom_ordercheck, &
           help_string="T: find molecules. F: go by order of atoms")
      call param_register(params, 'strict', 'true', this%strict, &
           help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
      call param_register(params, 'use_com', 'false', this%use_com, &
           help_string="Use COM instead of COG")
      call param_register(params, 'mpifind', 'false', this%mpifind, &
           help_string="Use find_monomer_triplets_MPI")
      call param_register(params, 'signature_one', PARAM_MANDATORY, signature_one_string, &
           help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_two', PARAM_MANDATORY, signature_two_string, &
           help_string="Atomic numbers of monomer two, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'signature_three', PARAM_MANDATORY, signature_three_string, &
           help_string="Atomic numbers of monomer three, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'power', '1.0', this%power, &
           help_string="Power of distances to be used in the kernel")
      call param_register(params, 'dist_shift', '0.0', this%dist_shift, &
           help_string="Distance shift for inverse distance descriptors.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='general_trimer_initialise args_str')) then
         RAISE_ERROR("general_trimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call split_string(signature_one_string,' ','{}',signature_one_fields(:),n_atoms_one,matching=.true.)
      call split_string(signature_two_string,' ','{}',signature_two_fields(:),n_atoms_two,matching=.true.)
      call split_string(signature_three_string,' ','{}',signature_three_fields(:),n_atoms_three,matching=.true.)
      allocate(this%signature_one(n_atoms_one))
      allocate(this%signature_two(n_atoms_two))
      allocate(this%signature_three(n_atoms_three))

      do i=1,n_atoms_one
        this%signature_one(i) = string_to_int(signature_one_fields(i))
      end do
      do i=1,n_atoms_two
        this%signature_two(i) = string_to_int(signature_two_fields(i))
      end do
      do i=1,n_atoms_three
        this%signature_three(i) = string_to_int(signature_three_fields(i))
      end do

      this%one_two_identical = .false.
      this%one_three_identical = .false.
      this%two_three_identical = .false.

      if (size(this%signature_one) == size(this%signature_two)) then
         if (all(this%signature_one == this%signature_two)) then
            this%one_two_identical = .True.
         end if
      end if

      if (size(this%signature_one) == size(this%signature_three)) then
         if (all(this%signature_one == this%signature_three)) then
            this%one_three_identical = .True.
         end if
      end if

      if (size(this%signature_two) == size(this%signature_three)) then
         if (all(this%signature_two == this%signature_three)) then
            this%two_three_identical = .True.
         end if
      end if

      call permutation_data_initialise(this%permutation_data,signature_one=this%signature_one,signature_two=this%signature_two,signature_three=this%signature_three,internal_swaps_only=this%internal_swaps_only,error=error)

      trimer_size=n_atoms_one + n_atoms_two + n_atoms_three
      d=trimer_size*(trimer_size-1)/2

      allocate(signature(trimer_size))
      allocate(intermolecular(trimer_size,trimer_size))
      allocate(this%is_intermolecular(d))
      allocate(this%cutoff_contributor(d))
      allocate(this%component_atoms(d,2))

      signature(1:n_atoms_one) = this%signature_one
      signature(1+n_atoms_one:n_atoms_one+n_atoms_two) = this%signature_two
      signature(1+n_atoms_one+n_atoms_two:trimer_size) = this%signature_three
      intermolecular = .false.
      this%cutoff_contributor=.false.

      do i=1,n_atoms_one
        do j=1+n_atoms_one,trimer_size
          intermolecular(i,j)=.true.
          intermolecular(j,i)=.true.
        end do
      end do
      do i=1+n_atoms_one,n_atoms_one+n_atoms_two
        do j=1+n_atoms_one+n_atoms_two,trimer_size
          intermolecular(i,j)=.true.
          intermolecular(j,i)=.true.
        end do
      end do

      start = 0
      finish=trimer_size-1
      do i=1,trimer_size
        do j=1,finish-start
          this%is_intermolecular(start+j) = intermolecular(i,i+j)
          this%component_atoms(start+j,:) = (/ i, i+j /)
        end do
        start = finish
        finish=finish + trimer_size-i-1
      end do

      do i=1,d
        if (this%is_intermolecular(i)) then
          if (.not. signature(this%component_atoms(i,1))==1 ) then
            if (.not. signature(this%component_atoms(i,2))==1 ) then
              this%cutoff_contributor(i)=.true.
            end if
          end if
        end if
      end do

      this%initialised = .true.
      deallocate(signature)
      deallocate(intermolecular)
   endsubroutine general_trimer_initialise


   subroutine water_trimer_initialise(this,args_str,error)
    type(water_trimer), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    integer, optional, intent(out) :: error
    integer :: i,j, start, finish
    integer, parameter :: trimer_size = 9, d = 36, at_per_mol = 3
    integer, parameter :: signature(*) = [8,1,1,8,1,1,8,1,1]
    logical :: intermolecular(9,9), ordercheck
    type(Dictionary) :: params
    
    print*, "water_trimer_initialize!"


    INIT_ERROR(error)

    call finalise(this)

    call initialise(params)
    call param_register(params, 'cutoff', '4.00', this%cutoff, help_string="Cutoff(intermolecular) for water_trimer-type descriptors")
    call param_register(params, 'monomer_cutoff', '1.90', this%monomer_cutoff, help_string="Monomer cutoff for water_trimer-type descriptors")
    call param_register(params, 'cutoff_transition_width', '1.00', this%cutoff_transition_width, help_string="Width of smooth cutoff region for water_trimer-type descriptors")
    call param_register(params, 'atom_ordercheck', 'false', ordercheck, help_string="T: find molecules. F: go by order of atoms")
    call param_register(params, 'strict', 'true', this%strict, help_string="Raise error if not all atoms assigned to monomer or if no monomer pairs found")
    call param_register(params, 'power', '1.0', this%power, help_string="Power of distances to be used in the kernel")
    call param_register(params, 'dist_shift', '0.0', this%dist_shift, help_string="Distance shift for inverse distance descriptors.")
    
        
    if (ordercheck)then
        RAISE_ERROR("This rouine doesn't do ordercheck, OHHOHHOHH is necessary",error)
    endif
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_trimer_initialise args_str')) then
       RAISE_ERROR("water_trimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)

    call permutation_data_initialise(this%permutation_data,signature_one=this%signature,signature_two=this%signature,signature_three=this%signature,internal_swaps_only=this%internal_swaps_only,error=error)


    ! allocate(signature(trimer_size))
    ! allocate(intermolecular(trimer_size,trimer_size))
    
    ! allocate(this%is_intermolecular(d))
    ! allocate(this%cutoff_contributor(d))
    ! allocate(this%component_atoms(d,2))

    ! signature(1:n_atoms_one) = this%signature_one
    ! signature(1+n_atoms_one:n_atoms_one+n_atoms_two) = this%signature_two
    ! signature(1+n_atoms_one+n_atoms_two:trimer_size) = this%signature_three
    intermolecular = .false.
    this%cutoff_contributor=.false.
    
    ! Set the intermolecular interactions
    do i=1,at_per_mol
      do j=1+at_per_mol,trimer_size
        intermolecular(i,j)=.true.
        intermolecular(j,i)=.true.
      end do
    end do
    do i=1+at_per_mol,at_per_mol+at_per_mol
      do j=1+at_per_mol+at_per_mol,trimer_size
        intermolecular(i,j)=.true.
        intermolecular(j,i)=.true.
      end do
    end do
    
    ! Linearize intermolecular matrix, and component_atoms
    start = 0
    finish=trimer_size-1
    do i=1,trimer_size
      do j=1,finish-start
        this%is_intermolecular(start+j) = intermolecular(i,i+j)
        this%component_atoms(start+j,:) = (/ i, i+j /)
      end do
      start = finish
      finish=finish + trimer_size-i-1
    end do
    
    ! Only non-hydrogens contribute to general cutoff
    do i=1,d
      if (this%is_intermolecular(i)) then
        if (.not. signature(this%component_atoms(i,1))==1 ) then
          if (.not. signature(this%component_atoms(i,2))==1 ) then
            this%cutoff_contributor(i)=.true.
          end if
        end if
      end if
    end do

    this%initialised = .true.
    print*, "water trimer initialisED!"
 endsubroutine water_trimer_initialise

   
   
   
   subroutine general_trimer_finalise(this,error)
      type(general_trimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)
      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_one_cutoff = 0.0_dp
      this%monomer_two_cutoff = 0.0_dp
      this%monomer_three_cutoff = 0.0_dp
      this%atom_ordercheck = .true.
      this%internal_swaps_only = .true.
      this%use_smooth_cutoff = .false.
      this%power = 1.0_dp
      this%dist_shift = 0.0_dp
      if(allocated(this%signature_one)) deallocate(this%signature_one)
      if(allocated(this%signature_two)) deallocate(this%signature_two)
      if(allocated(this%signature_three)) deallocate(this%signature_three)

      this%initialised = .false.

   endsubroutine general_trimer_finalise

   subroutine water_trimer_finalise(this,error)
    type(water_trimer), intent(inout) :: this
    integer, optional, intent(out) :: error
    print*, "  wt finilize!"
    INIT_ERROR(error)
    if(.not. this%initialised) return
    this%cutoff = 0.0_dp
    this%cutoff_transition_width = 0.0_dp
    this%atom_ordercheck = .false.
    this%internal_swaps_only = .true.
    this%use_smooth_cutoff = .false.
    this%power = 1.0_dp
    this%dist_shift = 0.0_dp
    this%initialised = .false.
    print*, "  wt finalizED!"
 endsubroutine water_trimer_finalise


   subroutine molecule_lo_d_initialise(this,args_str,error)
      type(molecule_lo_d), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      character(len=STRING_LENGTH) :: signature_string, atoms_template_string, symmetry_string, symmetry_property_name, append_file, append_string
      character(len=STRING_LENGTH), dimension(99) ::  append_rows,template_rows !symmetry_rows
      character(len=STRING_LENGTH_SHORT), dimension(99) :: signature_fields, row_fields
      integer, optional, intent(out) :: error
      integer :: i,n_atoms,j,n_symm_rows, current_depth, start, finish, i_component, atom_j, N_atom_pairs, atom_k,n_append_rows, old_size, n_perms,n_template_rows
      integer, dimension(:,:), allocatable :: equivalents_input, bonds_to_append, tmp_permutations
      logical :: signature_given, symmetries_given, add_bond, j_k_present, k_j_present
      integer, dimension(:,:), pointer :: symm_2d
      integer, dimension(:), pointer :: symm_1d

      type(Table) :: atom_a, atom_b
      type(CInOutput) :: tempatoms
      type(inoutput) :: tempfile
      type(inoutput) :: symmetry_inout
      type(inoutput) :: append_inout
      type(Connection) :: at_connect

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for molecule_lo_d-type descriptors")
      call param_register(params, 'atoms_template_string', PARAM_MANDATORY, atoms_template_string , help_string="Atoms object which serves as a template - written to a string")
      call param_register(params, 'neighbour_graph_depth', '2', this%neighbour_graph_depth, help_string="Ignore distances between atoms separated by more than this number of bonds")
      call param_register(params, 'signature', '', signature_string, help_string="Atomic numbers of monomer one, format {Z1 Z2 Z3 ...}")
      call param_register(params, 'symmetry_property_name', 'symm', symmetry_property_name, help_string="Integer arrays specifying symmetries - see header of make_permutations_v2.f95 for format")
      call param_register(params, 'append_file', '', append_file, help_string="Pairs of atoms for which we want the distance to be additionally included in the descriptor")
      call param_register(params, 'atom_ordercheck', 'T', this%atom_ordercheck, help_string= &
                                                                    "T: basic check that atoms in same order as in template F: assume all xyz frames have atoms in same order")
      call param_register(params, 'distance_transform', '0', this%distance_transform, help_string="0: distance matrix (no transform), 1: inverse distance matrix, 2: Coulomb matrix, 3: exponential, < 0: integer power of distance")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='molecule_lo_d_initialise args_str')) then
         RAISE_ERROR("molecule_lo_d_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      do i=1,len_trim(atoms_template_string)
        if(atoms_template_string(i:i)=='%') then
          atoms_template_string(i:i)=' '
        end if
      end do

      ! read in the atoms object in the sample geometry file
      call initialise(tempfile, filename="temp.xyz",action=OUTPUT)
      call split_string(atoms_template_string,';','{}',template_rows(:),n_template_rows,matching=.true.)
      do i=1,n_template_rows
        call print(template_rows(i),file=tempfile)
      end do
      call finalise(tempfile)
      call initialise(tempatoms,"temp.xyz")
      call read(this%template_atoms, tempatoms, error=error)

      this%n_atoms = this%template_atoms%N
      ! make a table of bonds - this copied from topology module private function create_bond_list
      call calc_connect(at_connect,this%template_atoms,error=error)
      if (this%neighbour_graph_depth > 0) then
        do i=1,this%template_atoms%N

           call initialise(atom_a,4,0,0,0,0)
           call append(atom_a,(/i,0,0,0/))
           call bfs_step(this%template_atoms,atom_a,atom_b,nneighb_only=.false.,min_images_only=.true.,alt_connect=at_connect)

           do j = 1,atom_b%N
              atom_j = atom_b%int(1,j)
              if (atom_j.gt.i) then
                 add_bond = .true.

                 if (add_bond) then
                    call append(this%bonds,(/i,atom_j/))
                 endif

              else
!                 call print('not added '//i//' -- '//atom_j)
              endif
           enddo
           call finalise(atom_a)
           call finalise(atom_b)
        enddo
        ! add atom pairs separated by up to neighbour_graph_depth bonds
        current_depth = 1
        this%atom_pairs=this%bonds
        do while (current_depth < this%neighbour_graph_depth)
           current_depth = current_depth + 1
           call bond_list_next_layer(this%bonds,this%atom_pairs)
        end do
      endif

     ! append any manually specified bonds, if not already present
      if (.not. append_file .eq. '') then
         call initialise(append_inout,trim(append_file))
         read(append_inout%unit,'(a)') append_string
         call split_string(append_string,',','{}',append_rows(:),n_append_rows,matching=.true.)
         
         do i = 1, n_append_rows
            call split_string(append_rows(i),' ','{}',row_fields(:),n_atoms,matching=.true.)
            if ( n_atoms .ne. 2) then
              RAISE_ERROR("append_file incorrectly formatted, expected atoms in pairs",error)
            end if
            if (.not. allocated(bonds_to_append)) allocate(bonds_to_append(n_append_rows,n_atoms))
            do j=1,n_atoms
              bonds_to_append(i,j) = string_to_int(row_fields(j))
            end do
         end do


         do i=1,size(bonds_to_append,1)
           atom_j = bonds_to_append(i,1)
           atom_k = bonds_to_append(i,2)
           if ( this%neighbour_graph_depth > 0) then
             N_atom_pairs = this%atom_pairs%N
             j_k_present = any(this%atom_pairs%int(1,:N_atom_pairs) .eq. atom_j .and. (this%atom_pairs%int(2,:N_atom_pairs) .eq. atom_k))
             k_j_present = any(this%atom_pairs%int(1,:N_atom_pairs) .eq. atom_k .and. (this%atom_pairs%int(2,:N_atom_pairs) .eq. atom_j))
           else
             j_k_present = .false.
             k_j_present = .false.
           end if
           if (.not. j_k_present .and. .not. k_j_present) then
             call append(this%atom_pairs,(/atom_j,atom_k/))
           end if
         end do
      end if

      !call print('table of atom pairs included in descriptor')
      !call print(this%atom_pairs)
      signature_given = .False.
      symmetries_given = .False.
      call permutation_data_finalise(this%permutation_data)
      ! parse signature, if given
      if (.not. signature_string .eq. '') then
         signature_given=.True.
         call split_string(signature_string,' ','{}',signature_fields(:),n_atoms,matching=.true.)
         if (.not. n_atoms .eq. this%n_atoms) then
           RAISE_ERROR('signature does not have correct number of entries', error)
         end if
         allocate(this%signature(this%n_atoms))

         do i=1,this%n_atoms
           this%signature(i) = string_to_int(signature_fields(i))
         end do
         call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)
      end if

      ! parse symmetries, if given
      symmetries_given = has_property(this%template_atoms,symmetry_property_name)
      if (symmetries_given) then
        if (.not. assign_pointer(this%template_atoms, symmetry_property_name, symm_1d)) then
           if (.not. assign_pointer(this%template_atoms, symmetry_property_name, symm_2d)) then
             RAISE_ERROR('IPModel_Coulomb_Calc failed to assign pointer to "'//trim(symmetry_property_name)//'" property', error)
           else
             if (.not. allocated(equivalents_input)) allocate(equivalents_input(size(symm_2d,1),size(symm_2d,2)))
             equivalents_input=symm_2d
           end if
        else
          if (.not. allocated(equivalents_input)) allocate(equivalents_input(1,size(symm_1d)))
             equivalents_input(1,:)=symm_1d
        endif
        call permutation_data_initialise(this%permutation_data,equivalents_input=equivalents_input,error=error)
      end if

      ! If no signature is given and no symmetries are manually specified we just make a dummy signature to initialise the permutation data
      if (.not. (signature_given) .and. .not. (symmetries_given)) then
         allocate(this%signature(this%n_atoms))
         do i=1,this%n_atoms
           this%signature(i) = i
         end do
         call permutation_data_initialise(this%permutation_data,signature_one=this%signature,error=error)
      end if

      ! make a mapping from i,j pairs to dist_vec components, and compile the list of these which actually make up the descriptor

      this%max_dimension = this%n_atoms * (this%n_atoms -1 ) / 2
      N_atom_pairs = this%atom_pairs%N
      allocate(this%component_atoms(this%max_dimension,2))
      allocate(this%included_components(N_atom_pairs))
      this%included_components=0

      i_component=1
      start = 0
      finish=this%n_atoms-1
      do i=1,this%n_atoms
        do j=1,finish-start
          this%component_atoms(start+j,:) = (/ i, i+j /)
          if (any(this%atom_pairs%int(1,:n_atom_pairs) .eq. i .and. this%atom_pairs%int(2,:n_atom_pairs) .eq. i+j) .or. &
              any(this%atom_pairs%int(2,:n_atom_pairs) .eq. i .and. this%atom_pairs%int(1,:n_atom_pairs) .eq. i+j))  then
               this%included_components(i_component) = start+j
               i_component = i_component + 1
          end if
        end do
        start = finish
        finish=finish + this%n_atoms-i-1
      end do

      if (any(this%included_components .eq. 0)) then
        RAISE_ERROR('molecule_lo_d_initialise : something went wrong picking out the correct interatomic distances',error)
      end if

      this%initialised = .true.

   endsubroutine molecule_lo_d_initialise

   subroutine molecule_lo_d_finalise(this,error)
      type(molecule_lo_d), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      if (allocated(this%signature)) deallocate(this%signature)
      if (allocated(this%component_atoms)) deallocate(this%component_atoms)
      if (allocated(this%included_components)) deallocate(this%included_components)



      this%initialised = .false.

   endsubroutine molecule_lo_d_finalise


   subroutine AN_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(AN_monomer), intent(in) :: this
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
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, n_index
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk
      real(dp), dimension(3) :: d_ij, d_ik, d_jk, u_ij, u_jk

      INIT_ERROR(error)

      call system_timer('AN_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_calc: descriptor object not initialised", error)
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

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='AN_monomer_calc args_str')) then
            RAISE_ERROR("AN_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("AN_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("AN_monomer_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      d = AN_monomer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)

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
            allocate(descriptor_out%x(i)%grad_data(d,3,0:this%N-1))
            allocate(descriptor_out%x(i)%ii(0:this%N-1))
            allocate(descriptor_out%x(i)%pos(3,0:this%N-1))
            allocate(descriptor_out%x(i)%has_grad_data(0:this%N-1))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:this%N-1))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if(at%N /= this%N) then
         RAISE_ERROR("AN_monomer_calc: number of atoms is "//at%N//" instead of "//this%N,error)
      endif

      do i = 1, at%N

         i_desc = 0

         if(my_do_descriptor) then
            if(this%do_atomic) then
               descriptor_out%x(i)%ci(1) = i
            else
               descriptor_out%x(i)%ci(:) = (/(m,m=1,this%N)/)
            endif
         endif

         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(0) = i
            descriptor_out%x(i)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i)%has_grad_data(:) = .true.
         endif

         do n = 1, n_neighbours(at,i)
            j = neighbour(at,i,n,distance=r_ij, cosines=u_ij, shift=shift_ij)

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i)%has_data = .true.
               descriptor_out%x(i)%data(i_desc) = r_ij
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i)%ii(n) = j
               descriptor_out%x(i)%pos(:,n) = at%pos(:,j) + matmul(at%lattice,shift_ij)

               descriptor_out%x(i)%grad_data(i_desc,:,n) =  u_ij
               descriptor_out%x(i)%grad_data(i_desc,:,0) = -u_ij
            endif

            do m = 1, n_neighbours(at,i)
               if(n >= m) cycle

               k = neighbour(at,i,m,distance=r_ik, shift=shift_ik)

               d_jk = ( at%pos(:,j) + matmul(at%lattice,shift_ij) ) - ( at%pos(:,k) + matmul(at%lattice,shift_ik) )
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               i_desc = i_desc + 1
               if(my_do_descriptor) then
                  descriptor_out%x(i)%has_data = .true.
                  descriptor_out%x(i)%data(i_desc) = r_jk
               endif

               if(my_do_grad_descriptor) then
                  descriptor_out%x(i)%grad_data(i_desc,:,n) =  u_jk 
                  descriptor_out%x(i)%grad_data(i_desc,:,m) = -u_jk 
               endif

            enddo
         enddo

         if(.not. this%do_atomic) exit

      enddo

      call system_timer('AN_monomer_calc')

   endsubroutine AN_monomer_calc

   subroutine general_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(general_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, monomer_size, i, &
         i_atomic, j_atomic, k, start, finish, n_index
      integer, dimension(3) :: temp_shift
      real(dp), dimension(:), allocatable :: dist_vec
      real(dp), dimension(:,:), allocatable :: interatomic_distances
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(:), allocatable :: atomic_index
      integer, dimension(:,:), allocatable :: monomer_index, shifts
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)

      call system_timer('general_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_calc: descriptor object not initialised", error)
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

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_monomer_calc args_str')) then
            RAISE_ERROR("general_monomer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("general_monomer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      monomer_size=size(this%signature)
      d = general_monomer_dimensions(this,error)

      allocate(shifts(monomer_size,3))
      allocate(dist_vec(d))
      allocate(atomic_index(monomer_size))
      allocate(associated_to_monomer(at%N))
      allocate(interatomic_vectors(monomer_size,monomer_size,3))
      allocate(interatomic_distances(monomer_size,monomer_size))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp
      associated_to_monomer=.False.

      call find_general_monomer(at, monomer_index, &
           this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)
      if(.not. all(associated_to_monomer)) then
         !RAISE_ERROR("general_monomer_calc: not all atoms assigned to a monomer", error)
         call print("Not all atoms can be assigned to a monomer with atomic numbers "//this%signature)
      endif
      n_descriptors = size(monomer_index,2)
      call print("found "//n_descriptors//" monomers", PRINT_VERBOSE)
      n_index = size(this%signature)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,monomer_size))
            allocate(descriptor_out%x(i)%ii(monomer_size))
            allocate(descriptor_out%x(i)%pos(3,monomer_size))
            allocate(descriptor_out%x(i)%has_grad_data(monomer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,monomer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      do i = 1, n_descriptors

         atomic_index = monomer_index(:,i) !stores the indices of atoms in this monomer
!write(*,*) "THE ATOMS IN THE MONOMER ARE : "// atomic_index

         if(associated(atom_mask_pointer)) then
            if(.not. any(atom_mask_pointer(atomic_index))) then
               cycle
            else
               if(.not. all(atom_mask_pointer(atomic_index))) then
                  RAISE_ERROR("general_monomer_calc: atom mask has to encompass either all or none of the atoms of a monomer",error)
               endif
            endif
         endif

         !calc all positions relative to atom 1
         do i_atomic=2,monomer_size
           temp_shift=0
           interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index(1),atomic_index(i_atomic),shift=temp_shift)
           shifts(i_atomic,:) = temp_shift
         end do

         !find other relative positions through vector addition
         do j_atomic=2,monomer_size
           do i_atomic=2,j_atomic-1
             interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
           end do
         end do

         !Now convert vectors to scalar distances
         do i_atomic=1,monomer_size
           do j_atomic=i_atomic+1,monomer_size
             interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
           end do
         end do
!!$do i_atomic=1,size(interatomic_distances,1)
!!$  write(*,'(6F12.8)') interatomic_distances(i_atomic,:)
!!$end do
         !and convert this NxN matrix into the required vector length N(N-1)/2
         start = 1
         finish = monomer_size-1
         do i_atomic=1,monomer_size-1
           dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:monomer_size)
           start = finish+1
           finish=finish + monomer_size-i_atomic-1
         end do

         if(my_do_descriptor) then
            descriptor_out%x(i)%ci(:) = atomic_index
            descriptor_out%x(i)%has_data = .true.
            descriptor_out%x(i)%data = dist_vec**this%power
         endif
         call print("distances: "//dist_vec, PRINT_VERBOSE)
         if(my_do_grad_descriptor) then
!!$write(*,*) "doing grad descriptor"
            descriptor_out%x(i)%ii(:) = atomic_index
!!$do i_atomic=1,at%N
!!$write(*,*) at%pos(:,atomic_index(i_atomic))
!!$end do
            descriptor_out%x(i)%pos(:,1) = at%pos(:,atomic_index(1))
            do i_atomic =2,monomer_size
              descriptor_out%x(i)%pos(:,i_atomic) = at%pos(:,atomic_index(i_atomic)) + matmul(at%lattice,shifts(i_atomic,:))
            end do

            !build the grad_data matrix
            descriptor_out%x(i)%has_grad_data(:) = .true.
            do k=1,d
             !find the pair of atoms contributing to this descriptor
             do i_atomic=1,monomer_size
               do j_atomic=i_atomic+1,monomer_size
                 if (interatomic_distances(i_atomic,j_atomic)==dist_vec(k)) then
                   descriptor_out%x(i)%grad_data(k,:,i_atomic) = -this%power * dist_vec(k)**(this%power-1.0_dp) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! kth descriptor wrt atom i_atomic
                   descriptor_out%x(i)%grad_data(k,:,j_atomic) = -descriptor_out%x(i)%grad_data(k,:,i_atomic)        ! kth descriptor wrt j_atomic
!write(*,*) "descriptor dimension "//k//" wrt atoms "//atomic_index(i_atomic)//" and "//atomic_index(j_atomic)
                 end if
               end do
             end do
            end do

       
         endif

      enddo

      deallocate(shifts)
      deallocate(dist_vec)
      deallocate(atomic_index)
      deallocate(associated_to_monomer)
      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)
      deallocate(monomer_index)
      call system_timer('general_monomer_calc')

   endsubroutine general_monomer_calc

   subroutine com_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

      type(com_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor, use_smooth_cutoff, double_count
      integer :: n_descriptors, dimer_size, i, j, i_atomic, j_atomic, i_desc, n_index
      integer :: monomer_one_size, monomer_two_size, n_monomer_one, n_monomer_two, this_pair, diff_loc
      integer, dimension(1) :: unit_array
      real(dp), dimension(3) :: diff_one_two, com_pos_one, com_pos_two, transdirvec
      real(dp) :: dist, primitive_cutoff, primitive_cutoff_grad
      real(dp), dimension(:,:), allocatable :: diffs_one, diffs_two, com_pos_diffs
      real(dp), dimension(:), allocatable :: weight_one, weight_two
      integer, dimension(:), allocatable :: atomic_index, atomic_index_one, atomic_index_two, pairs_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index,  monomer_pairs
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)
      use_smooth_cutoff = .false.
      double_count = .false.
      call system_timer('com_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_calc: descriptor object not initialised", error)
      endif

      if (.not. has_property(at, 'mass')) then
         RAISE_ERROR('com_dimer_calc: Atoms has no mass property', error)
      end if
      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      monomer_one_size =size(this%signature_one)
      monomer_two_size =size(this%signature_two)
      dimer_size = monomer_one_size + monomer_two_size

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
         help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='com_dimer_calc args_str')) then
            RAISE_ERROR("com_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("com_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      allocate(atomic_index(dimer_size))
      allocate(atomic_index_one(monomer_one_size))
      allocate(atomic_index_two(monomer_two_size))
      allocate(diffs_one(3,monomer_one_size))
      allocate(diffs_two(3,monomer_two_size))
      allocate(weight_one(monomer_one_size))
      allocate(weight_two(monomer_two_size))

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,&
           this%signature_one,associated_to_monomer,this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at,monomer_two_index,&
              this%signature_two,associated_to_monomer,this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         call print("WARNING: com_dimer_calc: not all atoms assigned to a monomer, if you have molecules present other than the following, this is OK")
         call print("signature of molecule 1 ")
         call print(this%signature_one)
         call print("signature of molecule 2 ")
         call print(this%signature_two)
      endif

      n_monomer_one = size(monomer_one_index,2)
      n_monomer_two = size(monomer_two_index,2)
      if (n_monomer_one < 1 .or. n_monomer_two < 1) then
        if ( this%strict ) then
          RAISE_ERROR("com_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs", error)
        else
          call print("WARNING: com_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs")
        end if
      end if

      call system_timer('com_dimer_calc: find_monomer_pairs')
      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         if(associated(atom_mask_pointer)) then
            call find_monomer_pairs_MPI(at,monomer_pairs,com_pos_diffs,pairs_diffs_map,&
                 monomer_one_index,monomer_two_index,this%monomers_identical,double_count,&
                 this%cutoff,error=error,use_com=.true.,atom_mask=atom_mask_pointer)
         else
            call find_monomer_pairs_MPI(at,monomer_pairs,com_pos_diffs,pairs_diffs_map,&
                 monomer_one_index,monomer_two_index,this%monomers_identical,double_count,&
                 this%cutoff,error=error,use_com=.true.)
         end if
      else
         call find_monomer_pairs    (at,monomer_pairs,com_pos_diffs,pairs_diffs_map,&
              monomer_one_index,monomer_two_index,this%monomers_identical,double_count,&
              this%cutoff,use_com=.true.,error=error)
      end if
      call system_timer('com_dimer_calc: find_monomer_pairs')

      if ( size(pairs_diffs_map) < 1) then
        if ( this%strict ) then
          RAISE_ERROR("com_dimer_calc did not find any monomer pairs to make a dimer", error)
        else
          call print("WARNING: com_dimer_calc did not find any monomer pairs to make a dimer")
        end if
      end if

      n_descriptors = size(pairs_diffs_map)
      n_index = size(this%signature_one) + size(this%signature_two)

      call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
      allocate(descriptor_out%x(n_descriptors))
      loop_descriptor_init: do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(1))
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(1,3,dimer_size))
            allocate(descriptor_out%x(i)%ii(dimer_size))
            allocate(descriptor_out%x(i)%pos(3,dimer_size))
            allocate(descriptor_out%x(i)%has_grad_data(dimer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,dimer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo loop_descriptor_init

      if (n_descriptors > 0) then ! only loop over monomers if we actually found any dimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_pairs(1,:) .eq. i)) cycle

            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i)

            if(associated(atom_mask_pointer)) then
               if(.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if(.not. all(atom_mask_pointer(atomic_index_one))) then
                     RAISE_ERROR("com_dimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                  endif
               endif
            endif

            do i_atomic=1,monomer_one_size
               weight_one(i_atomic) = at%mass(atomic_index_one(i_atomic))
               call print("weight_one("//i_atomic//") = "//weight_one(i_atomic), PRINT_NERD)
            end do
            call print("weight_one = "//weight_one, PRINT_NERD)
            call print("sum(weight_one) = "//sum(weight_one), PRINT_NERD)
            weight_one = weight_one / sum(weight_one)
            call print("weight_one = "//weight_one, PRINT_NERD)

            com_pos_one = centre_of_mass(at, index_list=atomic_index_one)
            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),com_pos_one)
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_pairs(2,:), monomer_pairs(1,:) .eq. i) ! find a monomer paired with i
               this_pair = unit_array(1)

               if (this_pair == 0) exit

               !get indices of monomer two
               j = monomer_pairs(2,this_pair)
               monomer_pairs(:,this_pair) = 0 ! make sure this pair isn't found again
               atomic_index_two = monomer_two_index(:,j)
               atomic_index=(/atomic_index_one,atomic_index_two/)

               do i_atomic=1,monomer_two_size
                  weight_two(i_atomic) = at%mass(atomic_index_two(i_atomic))
               end do
               weight_two = weight_two / sum(weight_two)
               call print("weight_two="//weight_two, PRINT_NERD)

               com_pos_two = centre_of_mass(at, index_list=atomic_index_two)

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do j_atomic=1,monomer_two_size
                  diffs_two(:,j_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(j_atomic)),com_pos_two)
               end do

               loop_different_shifts: do
                  unit_array = maxloc(pairs_diffs_map, pairs_diffs_map .eq. this_pair) ! find repeats of this pair with different shifts
                  diff_loc = unit_array(1)
                  if (diff_loc == 0) exit

                  i_desc = i_desc + 1
                  
                  diff_one_two = com_pos_diffs(:,diff_loc) ! shift between mean positions of these two monomers
                  pairs_diffs_map(diff_loc) = 0 ! make sure this shifted pair isn't found again

                  ! calculate distance
                  dist = norm(diff_one_two)
                  call print("COM distance = "//dist, PRINT_NERD)

                  calc_descriptor: if(my_do_descriptor) then
                     descriptor_out%x(i_desc)%has_data = .true.
                     if(this%transfer_parameters%do_transfer) then
                        descriptor_out%x(i_desc)%data = transferfunction(dist, this%transfer_parameters)
                     else
                        descriptor_out%x(i_desc)%data = dist
                     end if
                     descriptor_out%x(i_desc)%ci(:) = atomic_index
                     descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp
                     primitive_cutoff = coordination_function(dist,this%cutoff,this%cutoff_transition_width)
                     descriptor_out%x(i_desc)%covariance_cutoff = primitive_cutoff
                  end if calc_descriptor

                  calc_grad_descriptor: if(my_do_grad_descriptor) then !calc grads and update

                     descriptor_out%x(i_desc)%ii(:) = atomic_index
                     descriptor_out%x(i_desc)%has_grad_data(:) = .true.

                     primitive_cutoff_grad = dcoordination_function(dist,this%cutoff,this%cutoff_transition_width)

                     if(this%transfer_parameters%do_transfer) then
                        transdirvec = transferfunction_grad(dist, this%transfer_parameters) * diff_one_two / dist
                     else
                        transdirvec = diff_one_two / dist
                     endif

                     do i_atomic=1,monomer_one_size
                        descriptor_out%x(i_desc)%pos(:,i_atomic) = com_pos_one - diffs_one(:,i_atomic)
                        descriptor_out%x(i_desc)%grad_data(1,:,i_atomic) = - weight_one(i_atomic) * transdirvec ! descriptor wrt atom i_atomic
                        descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = primitive_cutoff_grad * descriptor_out%x(i_desc)%grad_data(1,:,i_atomic)
                     end do

                     do j_atomic=1,monomer_two_size
                        descriptor_out%x(i_desc)%pos(:,monomer_one_size+j_atomic) = com_pos_one + diff_one_two - diffs_two(:,j_atomic)
                        descriptor_out%x(i_desc)%grad_data(1,:,monomer_one_size+j_atomic) = weight_two(j_atomic) * transdirvec ! descriptor wrt atom j_atomic
                        descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+j_atomic) = primitive_cutoff_grad * descriptor_out%x(i_desc)%grad_data(1,:,monomer_one_size+j_atomic)
                     end do

                  endif calc_grad_descriptor
               enddo loop_different_shifts
            enddo loop_monomer_pairs
         enddo loop_monomer_one
      endif ! (n_descriptors > 0) ... still need to deallocate if no dimers were found:


      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_pairs)
      deallocate(com_pos_diffs)
      deallocate(pairs_diffs_map)

      deallocate(atomic_index)
      deallocate(atomic_index_one)
      deallocate(atomic_index_two)
      deallocate(weight_one)
      deallocate(weight_two)
      deallocate(diffs_one)
      deallocate(diffs_two)

      deallocate(associated_to_monomer)

      call system_timer('com_dimer_calc')

   endsubroutine com_dimer_calc

   subroutine general_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

      type(general_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional      :: do_descriptor, do_grad_descriptor !, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str
      integer, optional, intent(out) :: error
      real(dp) :: r_one_two
      type(Dictionary) :: params

      character(STRING_LENGTH)       :: atom_mask_name
      logical                        :: has_atom_mask_name = .false.
      logical, dimension(:), pointer :: atom_mask_pointer

      character(STRING_LENGTH)       :: persistent_indices_name
      logical                        :: has_persistent_indices_name = .false.
      integer, dimension(:), pointer :: persistent_indices_pointer

      ! lammps stuff: if we are called from lammps, atoms will be ordered
      logical               :: called_from_lammps = .false.
      type(atoms), pointer  :: at_ordered

      real(dp) :: com_dist, cutoff_grad
      logical :: my_do_descriptor, my_do_grad_descriptor, monomers_identical, use_smooth_cutoff, compound_cutoff
      integer :: d, n_descriptors, n_cross, dimer_size, i, j, k, n, m, &
         i_atomic, j_atomic, start, finish, i_desc, cutoff_pos, n_index
      integer :: monomer_one_size, monomer_two_size, n_monomer_one, n_monomer_two, n_products, this_pair, this_shift,diff_loc
      integer, dimension(1) :: unit_array
      real(dp), dimension(3) :: diff_one_two, temp_diff, mean_pos_one, mean_pos_two, transdirvec
      real(dp), dimension(:), allocatable :: dist_vec, primitive_cutoffs, primitive_cutoff_grads
      real(dp), dimension(:,:), allocatable :: interatomic_distances, diffs_one, diffs_two, mean_pos_diffs
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(3) :: temp_shift, shift_one_two
      real(dp), dimension(:), allocatable :: weight_one, weight_two
      integer, dimension(:), allocatable :: atomic_index, atomic_index_one, atomic_index_two, pairs_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index,  monomer_pairs
      logical, dimension(:), allocatable :: associated_to_monomer


      INIT_ERROR(error)
      use_smooth_cutoff = .false.
      call system_timer('general_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      monomer_one_size =size(this%signature_one)
      monomer_two_size =size(this%signature_two)
      dimer_size = monomer_one_size + monomer_two_size

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      atom_mask_pointer => null()
      if(present(args_str)) then
         call initialise(params)
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
              help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
         "calculated.")
         call param_register(params, 'persistent_indices_name', 'NONE', persistent_indices_name, has_value_target=has_persistent_indices_name, &
              help_string="Name of a property in the atoms object giving persistent indices, i.e. atom IDs that do not change from frame to frame")
         call param_register(params, 'lammps', 'F', called_from_lammps, &
              help_string="True if the potential was called from LAMMPS")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_dimer_calc args_str')) then
            RAISE_ERROR("general_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         call finalise(params)

         if (called_from_lammps .and. .not. (has_atom_mask_name .and.  has_persistent_indices_name)) then
            RAISE_ERROR("general_dimer_calc needs BOTH atom_mask and persistent_indices if called from LAMMPS", error)
         endif

      endif

      if (this%use_com .and. .not. has_property(at, 'mass')) then
         RAISE_ERROR('general_dimer_calc: Atoms has no mass property', error)
      end if

      d = general_dimer_dimensions(this,error)

      compound_cutoff=.True.
      if (count(this%cutoff_contributor) == 1) then
        compound_cutoff=.False.
      else if (count(this%cutoff_contributor) == 0) then
        RAISE_ERROR("general_dimer_calc, initialisation of general dimer did not find a pair of heavy atoms to use for calculating cutoff", error)
      end if

      allocate(dist_vec(d))
      allocate(primitive_cutoffs(d))
      allocate(primitive_cutoff_grads(d))

      allocate(atomic_index(dimer_size))
      allocate(atomic_index_one(monomer_one_size))
      allocate(atomic_index_two(monomer_two_size))
      allocate(diffs_one(3,monomer_one_size))
      allocate(diffs_two(3,monomer_two_size))
      allocate(weight_one(monomer_one_size))
      allocate(weight_two(monomer_two_size))

      allocate(interatomic_vectors(dimer_size,dimer_size,3))
      allocate(interatomic_distances(dimer_size,dimer_size))
      allocate(associated_to_monomer(at%N))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp
      associated_to_monomer=.false.

      ! if we are called from lammps, we need to sort atoms, so that they respect the persistent indices
      ! Warning: Potential performance hogm, but descriptors expect the atoms object not to be modified !
      ! Should we waste time sorting if we're not finding monomers based on index alone (i.e. if atom_ordercheck=T)?
      if (called_from_lammps) then
         allocate(at_ordered)
         at_ordered = at
         at_ordered%own_this = .true.
         call sort(at_ordered, persistent_indices_name, error=error)
      else
         call shallowcopy(at_ordered, at)
      endif

      if( has_atom_mask_name ) then
         call assign_property_pointer(at_ordered, trim(atom_mask_name), atom_mask_pointer, error)
      else
         atom_mask_pointer => null()
      endif

      call find_general_monomer(at_ordered, monomer_one_index,&
           this%signature_one,associated_to_monomer,&
           this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at_ordered, monomer_two_index,&
              this%signature_two,associated_to_monomer,&
              this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         call print("WARNING: general_dimer_calc: not all atoms assigned to a monomer, if you have molecules present other than the following, this is OK")
         call print("signature of molecule 1 ")
         call print(this%signature_one)
         call print("signature of molecule 2 ")
         call print(this%signature_two)
      endif

      n_monomer_one = size(monomer_one_index,2)
      n_monomer_two = size(monomer_two_index,2)
      if (n_monomer_one < 1 .or. n_monomer_two < 1) then
        if ( this%strict ) then
          RAISE_ERROR("general_dimer_calc,failed to find at least one of the monomer types, try increasing monomer cutoffs", error)
        else
          call print("WARNING: general_dimer_calc failed to find at least one of the monomer types, try increasing monomer cutoffs")
        end if
      end if

      call system_timer('general_dimer_calc: find_monomer_pairs')
      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         if(associated(atom_mask_pointer)) then
            call find_monomer_pairs_MPI(at_ordered, monomer_pairs, mean_pos_diffs, &
                 pairs_diffs_map,monomer_one_index,monomer_two_index,&
                 this%monomers_identical,this%double_count,this%cutoff,&
                 error=error,use_com=this%use_com,atom_mask=atom_mask_pointer)
         else
            call find_monomer_pairs_MPI(at_ordered, monomer_pairs, mean_pos_diffs, &
                 pairs_diffs_map,monomer_one_index,monomer_two_index,&
                 this%monomers_identical,this%double_count,this%cutoff,&
                 error=error,use_com=this%use_com)
         end if
      else
         call find_monomer_pairs(at_ordered,monomer_pairs,mean_pos_diffs,&
              pairs_diffs_map,monomer_one_index,monomer_two_index,&
              this%monomers_identical,this%double_count,this%cutoff,&
              error=error,use_com=this%use_com)
      end if
      call system_timer('general_dimer_calc: find_monomer_pairs')

      if ( size(pairs_diffs_map) < 1) then
        if ( this%strict ) then
          RAISE_ERROR("general_dimer_calc did not find any monomer pairs to make a dimer", error)
        else
          call print("WARNING: general_dimer_calc did not find any monomer pairs to make a dimer")
        end if
      end if

      n_descriptors = size(pairs_diffs_map)
      call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)

      n_index = size(this%signature_one) + size(this%signature_two)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if (my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%ci = 0
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         end if
         if (my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,dimer_size))
            allocate(descriptor_out%x(i)%ii(dimer_size))
            allocate(descriptor_out%x(i)%pos(3,dimer_size))
            allocate(descriptor_out%x(i)%has_grad_data(dimer_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,dimer_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         end if
      end do

      has_descriptors: if (n_descriptors > 0) then ! only loop over monomers if we actually found any dimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_pairs(1,:) .eq. i)) cycle

            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i)

            if (associated(atom_mask_pointer)) then
               if (.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if (.not. all(atom_mask_pointer(atomic_index_one))) then
                     if (this%strict_mask) then
                        RAISE_ERROR("general_dimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                     else
                        call print("general_dimer_calc: atom mask encompasses only part of a monomer; deciding based on first atom.", PRINT_NERD)
                        if (.not. atom_mask_pointer(atomic_index_one(1))) then
                            cycle
                        end if
                     end if
                  end if
               end if
            end if

            if (this%use_com) then
               do i_atomic=1,monomer_one_size
                  weight_one(i_atomic) = at_ordered%mass(atomic_index_one(i_atomic))
               end do
               weight_one = weight_one / sum(weight_one)
               mean_pos_one = centre_of_mass(at_ordered, index_list=atomic_index_one)
            else
               weight_one = 1.0_dp / monomer_one_size
               mean_pos_one = calc_mean_pos(at_ordered,atomic_index_one)
            end if


            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at_ordered,at_ordered%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
               interatomic_vectors(1,i_atomic,:) = diff_min_image(at_ordered,atomic_index_one(1),atomic_index_one(i_atomic))
            end do

            !find other relative positions through vector addition
            do j_atomic=2,monomer_one_size
               do i_atomic=2,j_atomic-1
                  interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
               end do
            end do

            !And convert vectors to scalar distances
            do i_atomic=1,monomer_one_size
               do j_atomic=i_atomic+1,monomer_one_size
                  interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
               end do
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_pairs(2,:), monomer_pairs(1,:) .eq. i) ! find a monomer paired with i
               this_pair = unit_array(1)

               if (this_pair == 0) exit

               !get indices of monomer two
               j = monomer_pairs(2,this_pair)
               monomer_pairs(:,this_pair) = 0 ! make sure this pair isn't found again
               atomic_index_two = monomer_two_index(:,j)
               atomic_index=(/atomic_index_one,atomic_index_two/)

               if (this%use_com) then
                  do j_atomic=1,monomer_two_size
                     weight_two(j_atomic) = at_ordered%mass(atomic_index_two(j_atomic))
                  end do
                  weight_two = weight_two / sum(weight_two)
                  mean_pos_two = centre_of_mass(at_ordered, index_list=atomic_index_two)
               else
                  weight_two = 1.0_dp / monomer_two_size
                  mean_pos_two = calc_mean_pos(at_ordered,atomic_index_two)
               end if

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do i_atomic=1,monomer_two_size
                  diffs_two(:,i_atomic) = diff_min_image(at_ordered,at_ordered%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
                  interatomic_vectors(monomer_one_size+1,monomer_one_size+i_atomic,:) = diff_min_image(at_ordered,atomic_index_two(1),atomic_index_two(i_atomic))
               end do

               !find other relative positions through vector addition
               do j_atomic=monomer_one_size+2,dimer_size
                  do i_atomic=monomer_one_size+2,j_atomic-1
                     interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(monomer_one_size+1,j_atomic,:) - interatomic_vectors(monomer_one_size+1,i_atomic,:)
                  end do
               end do

               !And convert vectors to scalar distances
               do i_atomic=monomer_one_size+1,dimer_size
                  do j_atomic=i_atomic+1,dimer_size
                     interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                  end do
               end do

               loop_pair_shifts: do
                  unit_array = maxloc(pairs_diffs_map, pairs_diffs_map .eq. this_pair) ! find repeats of this pair with different shifts
                  diff_loc = unit_array(1)
                  if (diff_loc == 0) exit

                  i_desc = i_desc + 1
                  
                  diff_one_two = mean_pos_diffs(:,diff_loc) ! shift between mean positions of these two monomers
                  pairs_diffs_map(diff_loc) = 0 ! make sure this shifted pair isn't found again

                  com_dist = norm(diff_one_two)

                  ! calculate intermolecular distances, also by vector addition
                  do i_atomic=1,monomer_one_size
                     do j_atomic = monomer_one_size+1,dimer_size
                        interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_two - diffs_two(:,j_atomic-monomer_one_size)
                        interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                     end do
                  end do

                  !Now take the whole matrix of scalar distances and combine into 1D array
                  start = 1
                  finish = dimer_size-1
                  do i_atomic=1,dimer_size-1
                     dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:dimer_size)
                     start = finish+1
                     finish=finish + dimer_size-i_atomic-1
                  end do
                  call print("dist vec "//dist_vec,PRINT_NERD)

                  primitive_cutoffs=1.0_dp
                  do k=1,d
                     if (this%cutoff_contributor(k)) then
                        primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                     end if
                  end do

                  calc_descriptor: if (my_do_descriptor) then
                     descriptor_out%x(i_desc)%has_data = .true.

                     do_transfer: if (this%transfer_parameters%do_transfer) then
                        do k=1,d
                           if (this%is_intermolecular(k)) then
                              descriptor_out%x(i_desc)%data(k) = transferfunction(dist_vec(k), this%transfer_parameters)
                           else
                              ! don't apply transfer function to intra-molecular (bond) distances
                              descriptor_out%x(i_desc)%data(k) = dist_vec(k)
                           end if
                        end do
                     else
                        descriptor_out%x(i_desc)%data = (dist_vec+this%dist_shift)**this%power
                     end if do_transfer

                     descriptor_out%x(i_desc)%ci(:) = atomic_index
                     descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp

                     is_com_cutoff: if (this%mpifind) then
                        descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(com_dist, this%cutoff, this%cutoff_transition_width)
                     else

                        is_compound_cutoff: if (compound_cutoff) then
                        ! Covariance cutoff is sum of pairwise products of primitive cutoffs for all *inter*molecular distances excluding H atoms

                           n_products=0
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                      n_products = n_products + 1
                                      descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff + primitive_cutoffs(k)*primitive_cutoffs(m)
                                    end if
                                 end do
                              end if
                           end do
                           ! normalise
                           descriptor_out%x(i_desc)%covariance_cutoff = descriptor_out%x(i_desc)%covariance_cutoff / n_products

                        else ! Covariance cutoff is primitive cutoff, i.e. there is only one pair of heavy atoms

                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 cutoff_pos=k
                                 descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                                 exit
                              end if
                           end do

                        end if is_compound_cutoff

                     end if is_com_cutoff
                  end if calc_descriptor

                  calc_grad_descriptor: if (my_do_grad_descriptor) then !calc grads and update

                     descriptor_out%x(i_desc)%ii(:) = atomic_index

                     do i_atomic=1,monomer_one_size
                       descriptor_out%x(i_desc)%pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
                     end do
                     do i_atomic=1,monomer_two_size
                       descriptor_out%x(i_desc)%pos(:,monomer_one_size+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
                     end do

                     call print("DIMER  "//dimer_size,PRINT_NERD)
                     call print('DIMER  Lattice="10.0000000       0.00000000       0.00000000       0.00000000      10.0000000       0.00000000       0.00000000       0.00000000      10.00000" Properties=Z:I:1:pos:R:3',PRINT_NERD)
                     do i_atomic=1,dimer_size
                       call print("DIMER  "//at_ordered%Z(atomic_index(i_atomic))//"    "//descriptor_out%x(i_desc)%pos(:,i_atomic),PRINT_NERD)
                     end do

                     !build the grad_data matrix
                     descriptor_out%x(i_desc)%has_grad_data(:) = .true.

                     do k=1,d
                        !get pair of atoms contributing to this component
                        i_atomic = this%component_atoms(k,1)
                        j_atomic = this%component_atoms(k,2)
                        do_transfer_grad: if (this%transfer_parameters%do_transfer .and. this%is_intermolecular(k)) then
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -transferfunction_grad(dist_vec(k), this%transfer_parameters) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic, using a transfer function
                        else
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = - this%power * (dist_vec(k)+this%dist_shift)**(this%power-1.0_dp) * &
                              interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
                        end if do_transfer_grad
                        descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic
                     end do

                     is_com_cutoff_grad: if (this%mpifind) then
                        transdirvec = diff_one_two / com_dist

                        cutoff_grad = dcoordination_function(com_dist, this%cutoff, this%cutoff_transition_width)
                        do i_atomic=1,monomer_one_size
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = - weight_one(i_atomic) * transdirvec * cutoff_grad
                        end do
                        do j_atomic=1,monomer_two_size
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+j_atomic) = weight_two(j_atomic) * transdirvec * cutoff_grad
                        end do
                     else

                        primitive_cutoff_grads=0.0_dp
                        do k=1,d
                           if (this%cutoff_contributor(k)) then
                             primitive_cutoff_grads(k) = dcoordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                           end if
                        end do

                        is_compound_cutoff_grad: if (compound_cutoff) then

                           do i_atomic=1,dimer_size ! for each atom in the dimer...
                              do k=1,d ! ...iterate over all distance...
                                 if (this%cutoff_contributor(k)) then
                                    do m=k+1,d ! ...pairs involved...
                                       if (this%cutoff_contributor(m)) then
                                          if (any(this%component_atoms(k,:) == i_atomic) .or. any(this%component_atoms(m,:) == i_atomic)) then
                                             descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                             = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                                + primitive_cutoffs(m)*primitive_cutoff_grads(k)*descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
                                                + primitive_cutoffs(k)*primitive_cutoff_grads(m)*descriptor_out%x(i_desc)%grad_data(m,:,i_atomic)
                                          end if
                                       end if
                                    end do
                                 end if
                              end do
                           end do
                           !normalisation factor
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) / n_products

                        else

                           i_atomic = this%component_atoms(cutoff_pos,1)
                           j_atomic = this%component_atoms(cutoff_pos,2)
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = primitive_cutoff_grads(cutoff_pos) *descriptor_out%x(i_desc)%grad_data(cutoff_pos,:,i_atomic)
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,j_atomic) = primitive_cutoff_grads(cutoff_pos) *descriptor_out%x(i_desc)%grad_data(cutoff_pos,:,j_atomic)

                        end if is_compound_cutoff_grad

                     end if is_com_cutoff_grad
                  end if calc_grad_descriptor

               end do loop_pair_shifts
            end do loop_monomer_pairs
         end do loop_monomer_one
      end if has_descriptors ! (n_descriptors > 0) ... still need to deallocate if no dimers were found:

      call finalise_ptr(at_ordered)

      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_pairs)
      deallocate(mean_pos_diffs)
      deallocate(pairs_diffs_map)

      deallocate(dist_vec)
      deallocate(primitive_cutoffs)
      deallocate(primitive_cutoff_grads)

      deallocate(atomic_index)
      deallocate(atomic_index_one)
      deallocate(atomic_index_two)
      deallocate(weight_one)
      deallocate(weight_two)
      deallocate(diffs_one)
      deallocate(diffs_two)

      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)
      deallocate(associated_to_monomer)

      call system_timer('general_dimer_calc')

   end subroutine general_dimer_calc

   subroutine general_trimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)

     type(general_trimer), intent(in) :: this
     type(atoms), intent(in) :: at
     type(descriptor_data), intent(out) :: descriptor_out
     logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
     character(len=*), intent(in), optional :: args_str
     integer, optional, intent(out) :: error
     real(dp) :: r_one_two

     type(Dictionary) :: params
     character(STRING_LENGTH) :: atom_mask_name
     logical :: has_atom_mask_name
     logical, dimension(:), pointer :: atom_mask_pointer

     logical :: my_do_descriptor, my_do_grad_descriptor, one_two_identical, one_three_identical, two_three_identical, &
          use_smooth_cutoff, done_this_monomer, done_this_dimer

     integer :: d, n_descriptors, n_cross, dimer_size, trimer_size, &
        i, j, k, n, m, i_atomic, j_atomic, start, finish, i_desc,this_diff,triplet_pos, &
        n_index
     integer :: monomer_one_size, monomer_two_size, monomer_three_size, n_monomer_one, n_monomer_two, n_monomer_three, n_products!, contributor_loc
     integer, dimension(1) :: unit_array
     real(dp) :: temp_dist
     real(dp), dimension(3) :: diff_one_two, diff_one_three,diff_two_three,mean_pos_one,mean_pos_two,mean_pos_three
     real(dp), dimension(3) :: grad_cut_one, grad_cut_two, grad_cut_three
     real(dp) :: dist_one_two, dist_one_three, dist_two_three, cut12, cut13, cut23, dcut12, dcut13, dcut23
     real(dp), dimension(:), allocatable :: dist_vec, primitive_cutoffs, primitive_cutoff_grads
     real(dp), dimension(:,:), allocatable :: interatomic_distances,diffs_one,diffs_two,diffs_three,triplets_diffs
     real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
     integer, dimension(3) :: temp_shift, shift_one_two,shift_one_three
     integer, dimension(:), allocatable :: atomic_index_dimer, atomic_index_trimer, atomic_index_one, atomic_index_two, atomic_index_three,triplets_diffs_map
     integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_three_index, monomer_triplets
     logical, dimension(:), allocatable :: associated_to_monomer
     logical :: double_count

     INIT_ERROR(error)
     use_smooth_cutoff = .false.
     double_count=.false.
     call system_timer('general_trimer_calc')

     if(.not. this%initialised) then
        RAISE_ERROR("general_trimer_calc: descriptor object not initialised", error)
     endif

     my_do_descriptor = optional_default(.false., do_descriptor)
     my_do_grad_descriptor = .true.!optional_default(.false., do_grad_descriptor)

     monomer_one_size =size(this%signature_one)
     monomer_two_size =size(this%signature_two)
     monomer_three_size =size(this%signature_three)
     dimer_size = monomer_one_size + monomer_two_size
     trimer_size = monomer_one_size + monomer_two_size + monomer_three_size

     if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

     call finalise(descriptor_out)

     atom_mask_pointer => null()
     if(present(args_str)) then
        call initialise(params)
        
        call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
        help_string="Name of a logical property in the atoms object. For atoms where this property is true descriptors are " // &
        "calculated.")

        if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_trimer_calc args_str')) then
           RAISE_ERROR("general_trimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
        endif
        
        call finalise(params)

        if( has_atom_mask_name ) then
           if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
              RAISE_ERROR("general_trimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
           endif
        else
           atom_mask_pointer => null()
        endif

     endif

     d = general_trimer_dimensions(this,error)

     allocate(dist_vec(d))
     allocate(primitive_cutoffs(d))
     allocate(primitive_cutoff_grads(d))

     allocate(atomic_index_one(monomer_one_size))
     allocate(atomic_index_two(monomer_two_size))
     allocate(atomic_index_three(monomer_three_size))
     allocate(atomic_index_dimer(dimer_size))
     allocate(atomic_index_trimer(trimer_size))

     allocate(interatomic_vectors(trimer_size,trimer_size,3))
     allocate(interatomic_distances(trimer_size,trimer_size))
     allocate(associated_to_monomer(at%N))
     allocate(diffs_one(3,monomer_one_size))
     allocate(diffs_two(3,monomer_two_size))
     allocate(diffs_three(3,monomer_three_size))

     interatomic_vectors = 0.0_dp
     interatomic_distances = 0.0_dp
     associated_to_monomer=.false.

     call find_general_monomer(at,monomer_one_index,&
          this%signature_one,associated_to_monomer,&
          this%monomer_one_cutoff,this%atom_ordercheck,error)
     if (this%one_two_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
     else
        call find_general_monomer(at,monomer_two_index,&
             this%signature_two,associated_to_monomer,&
             this%monomer_two_cutoff,this%atom_ordercheck,error)
     end if
     if (this%one_three_identical) then
        allocate(monomer_three_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_three_index = monomer_one_index
     else if (this%two_three_identical) then
        allocate(monomer_three_index(size(monomer_two_index,1),size(monomer_two_index,2)))
        monomer_three_index = monomer_two_index
     else
        call find_general_monomer(at,monomer_three_index,&
             this%signature_three,associated_to_monomer,&
             this%monomer_three_cutoff,this%atom_ordercheck,error)
     end if

     if(.not. all(associated_to_monomer)) then
        if(this%strict) then
           RAISE_ERROR("general_trimer_calc: not all atoms assigned to a monomer", error)
        else
           call print("WARNING: general_trimer_calc: not all atoms assigned to a monomer")
        endif
     endif

     n_monomer_one = size(monomer_one_index,2)
     n_monomer_two = size(monomer_two_index,2)
     n_monomer_three = size(monomer_three_index,2)

     if (this%use_com) then
        RAISE_ERROR("general_trimer_calc: use_com=T not implemented yet", error)
     end if
     call system_timer('general_trimer_calc: find_monomer_triplets')
     if(this%mpifind) then
        call print("Using find_monomer_triplets_MPI", PRINT_NERD)
        if(associated(atom_mask_pointer)) then
           call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
                monomer_one_index,monomer_two_index,monomer_three_index,&
                this%one_two_identical,this%one_three_identical,this%two_three_identical,&
                this%cutoff,error,use_com=.false.,atom_mask=atom_mask_pointer)
        else
           call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
                monomer_one_index,monomer_two_index,monomer_three_index,&
                this%one_two_identical,this%one_three_identical,this%two_three_identical,&
                this%cutoff,error,use_com=.false.)
        end if
     else
        call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
             monomer_one_index,monomer_two_index,monomer_three_index,&
             this%one_two_identical,this%one_three_identical,this%two_three_identical,&
             this%cutoff,error)
     end if
     call system_timer('general_trimer_calc: find_monomer_triplets')
     call print("monomer_triplets("//size(monomer_triplets,1)//", "//size(monomer_triplets,2)//")", PRINT_NERD)
     call print("triplets_diffs("//size(triplets_diffs,1)//", "//size(triplets_diffs,2)//")", PRINT_NERD)
     call print("triplets_diffs_map("//size(triplets_diffs_map)//")", PRINT_NERD)

     !call print("monomer_triplets",PRINT_NERD)
     !call print(monomer_triplets,PRINT_NERD)

     n_descriptors = size(triplets_diffs_map)
     n_index = size(this%signature_one) + size(this%signature_two) + size(this%signature_three)

     call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
     allocate(descriptor_out%x(n_descriptors))
     do i = 1, n_descriptors
        if(my_do_descriptor) then
           allocate(descriptor_out%x(i)%data(d))
           allocate(descriptor_out%x(i)%ci(n_index))
           descriptor_out%x(i)%data = 0.0_dp
           descriptor_out%x(i)%ci = 0
           descriptor_out%x(i)%has_data = .false.
           descriptor_out%x(i)%covariance_cutoff = 1.0_dp
        endif
        if(my_do_grad_descriptor) then
           allocate(descriptor_out%x(i)%grad_data(d,3,trimer_size))
           allocate(descriptor_out%x(i)%ii(trimer_size))
           allocate(descriptor_out%x(i)%pos(3,trimer_size))
           allocate(descriptor_out%x(i)%has_grad_data(trimer_size))
           descriptor_out%x(i)%grad_data = 0.0_dp
           descriptor_out%x(i)%ii = 0
           descriptor_out%x(i)%pos = 0.0_dp
           descriptor_out%x(i)%has_grad_data = .false.

           allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,trimer_size))
           descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
        endif
     enddo

      has_descriptors: if (n_descriptors > 0) then ! only loop over monomers if we actually found any trimers
         i_desc = 0
         loop_monomer_one: do i = 1, n_monomer_one
            if (.not. any(monomer_triplets(1,:) .eq. i)) cycle
            done_this_monomer = .false.
            !get indices of monomer and calc internal distances
            atomic_index_one = monomer_one_index(:,i) !store the indices of atoms in this monomer

            if(associated(atom_mask_pointer)) then
               if(.not. any(atom_mask_pointer(atomic_index_one))) then
                  cycle
               else
                  if(.not. all(atom_mask_pointer(atomic_index_one))) then
                     RAISE_ERROR("general_trimer_calc: atom mask has to encompass either all or none of the atoms of monomer one",error)
                  endif
               endif
            endif

            mean_pos_one = calc_mean_pos(at,atomic_index_one)

            !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
            do i_atomic=1,monomer_one_size
               diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
               interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index_one(1),atomic_index_one(i_atomic))
            end do

            !find other relative positions through vector addition
            do j_atomic=2,monomer_one_size
               do i_atomic=2,j_atomic-1
                  interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
               end do
            end do

            !Now convert vectors to scalar distances
            do i_atomic=1,monomer_one_size
               do j_atomic=i_atomic+1,monomer_one_size
                  interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
               end do
            end do

            ! Loop through monomers paired with this one to make dimers
            loop_monomer_pairs: do
               unit_array = maxloc(monomer_triplets(2,:), monomer_triplets(1,:) .eq. i) ! find a monomer paired with i
               if (all(unit_array .eq. 0)) exit

               !get indices of monomer two
               j = monomer_triplets(2,unit_array(1))
               atomic_index_two = monomer_two_index(:,j)
               atomic_index_dimer=(/atomic_index_one,atomic_index_two/)

               mean_pos_two = calc_mean_pos(at,atomic_index_two)

               ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
               do i_atomic=1,monomer_two_size
                  diffs_two(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
                  interatomic_vectors(monomer_one_size+1,monomer_one_size+i_atomic,:) = diff_min_image(at,atomic_index_two(1),atomic_index_two(i_atomic))
               end do

               !find other relative positions through vector addition
               do j_atomic=monomer_one_size+2,dimer_size
                  do i_atomic=monomer_one_size+2,j_atomic-1
                     interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(monomer_one_size+1,j_atomic,:) - interatomic_vectors(monomer_one_size+1,i_atomic,:)
                  end do
               end do

               !And convert vectors to scalar distances
               do i_atomic=monomer_one_size+1,dimer_size
                  do j_atomic=i_atomic+1,dimer_size
                     interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                  end do
               end do

               !! Now make trimers based on this dimer
               loop_triplets: do
                  unit_array = maxloc(monomer_triplets(3,:), monomer_triplets(1,:) .eq. i .and. monomer_triplets(2,:) .eq. j ) ! look for  a triplet
                  if (all(unit_array .eq. 0)) exit
                  triplet_pos=unit_array(1)

                  !get indices of monomer three
                  k = monomer_triplets(3,triplet_pos)
                  monomer_triplets(:,triplet_pos) = 0 ! make sure this triplet isn't found again
                  atomic_index_three = monomer_three_index(:,k)
                  atomic_index_trimer=(/atomic_index_dimer,atomic_index_three/)
                  mean_pos_three = calc_mean_pos(at,atomic_index_three)
                  ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
                  do i_atomic=1,monomer_three_size
                     diffs_three(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_three(i_atomic)),mean_pos_three)
                     interatomic_vectors(dimer_size+1,dimer_size+i_atomic,:) = diff_min_image(at,atomic_index_three(1),atomic_index_three(i_atomic))
                  end do

                  !find other relative positions through vector addition
                  do j_atomic=dimer_size+2,trimer_size
                     do i_atomic=dimer_size+2,j_atomic-1
                        interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(dimer_size+1,j_atomic,:) -interatomic_vectors(dimer_size+1,i_atomic,:)
                     end do
                  end do

                  !Now convert vectors to scalar distances
                  do i_atomic=dimer_size+1,trimer_size
                     do j_atomic=i_atomic+1,trimer_size
                        interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                     end do
                  end do

                  ! Loop over all trimers which can be created from these three monomers, i.e. with different periodic images
                  loop_trimers: do
                     unit_array = maxloc(triplets_diffs_map, triplets_diffs_map .eq. triplet_pos)
                     this_diff = unit_array(1)
                     if (this_diff == 0) exit

                     i_desc = i_desc + 1

                     diff_one_two = triplets_diffs(1:3,this_diff)
                     diff_one_three = triplets_diffs(4:6,this_diff)
                     diff_two_three = diff_one_three - diff_one_two
                     triplets_diffs_map(this_diff)=0 ! make sure this shifted triplet isn't found again
                     dist_one_two = norm(diff_one_two)
                     dist_one_three = norm(diff_one_three)
                     dist_two_three = norm(diff_two_three)

                     ! calculate intermolecular distances, monomers one and two
                     do i_atomic=1,monomer_one_size
                        do j_atomic = monomer_one_size+1,dimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_two - diffs_two(:,j_atomic-monomer_one_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     ! calculate intermolecular distances, monomers one and three
                     do i_atomic=1,monomer_one_size
                        do j_atomic = dimer_size+1,trimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_three - diffs_three(:,j_atomic-dimer_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     ! calculate intermolecular distances, monomers two and three
                     do i_atomic=monomer_one_size+1,dimer_size
                        do j_atomic = dimer_size+1,trimer_size
                           interatomic_vectors(i_atomic,j_atomic,:) = diffs_two(:,i_atomic-monomer_one_size) + diff_two_three - diffs_three(:,j_atomic-dimer_size)
                           interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
                        end do
                     end do

                     !Now take the whole matrix of scalar distances and combine into 1D array
                     start = 1
                     finish=trimer_size-1
                     do i_atomic=1,trimer_size-1
                        dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:trimer_size)
                        start = finish+1
                        finish=finish + trimer_size-i_atomic-1
                     end do
                     call print( "list of distances: "//dist_vec,PRINT_NERD)

                     calc_descriptor: if(my_do_descriptor) then
                        descriptor_out%x(i_desc)%has_data = .true.
                        descriptor_out%x(i_desc)%data = (dist_vec+this%dist_shift)**this%power
                        descriptor_out%x(i_desc)%ci(:) = atomic_index_trimer

                        cutoff_type: if (this%mpifind) then
                           cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
                        else
                           descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp

                           primitive_cutoffs=1.0_dp

                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                              end if
                           end do
                           n_products=0
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                 do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                       n_products=n_products+1
                                       descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff + primitive_cutoffs(k)*primitive_cutoffs(m)
                                    end if
                                 end do
                              end if
                           end do
                           ! normalise
                           descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff / n_products
                        end if cutoff_type
                     end if calc_descriptor

                     calc_grad_descriptor: if(my_do_grad_descriptor) then !calc grads and update

                        descriptor_out%x(i_desc)%ii(:) = atomic_index_trimer
                        do i_atomic=1,monomer_one_size
                           descriptor_out%x(i_desc)%pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
                        end do
                        do i_atomic=1,monomer_two_size
                           descriptor_out%x(i_desc)%pos(:,monomer_one_size+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
                        end do
                        do i_atomic=1,monomer_three_size
                           descriptor_out%x(i_desc)%pos(:,dimer_size+i_atomic) = mean_pos_one + diff_one_three - diffs_three(:,i_atomic)
                        end do

                        call print("TRIM  9",PRINT_NERD)
                        call print('TRIM  Lattice="10.0000000       0.00000000       0.00000000       0.00000000      10.0000000       0.00000000       0.00000000       0.00000000      10.00000" Properties=Z:I:1:pos:R:3',PRINT_NERD)
                        do i_atomic=1,trimer_size
                        call print("TRIM  "//at%Z(atomic_index_trimer(i_atomic))//"    "//descriptor_out%x(i_desc)%pos(:,i_atomic),PRINT_NERD)
                        end do

                        !build the grad_data matrix
                        descriptor_out%x(i_desc)%has_grad_data(:) = .true.
                        do k=1,d
                           i_atomic = this%component_atoms(k,1)
                           j_atomic = this%component_atoms(k,2)
                           descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = - this%power * (dist_vec(k)+this%dist_shift)**(this%power-1.0_dp) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
                           descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic
                        end do

                        cutoff_type_grad: if (this%mpifind) then
                           cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           dcut12 = dcoordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
                           dcut13 = dcoordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
                           dcut23 = dcoordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
                           !descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
                           grad_cut_one = (- diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
                                           - diff_one_three / dist_one_three * dcut13 * (cut12+cut23)) / 3.0_dp
                           grad_cut_two = (  diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
                                           - diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
                           grad_cut_three = (diff_one_three / dist_one_three * dcut13 * (cut12+cut23) &
                                           + diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
                           do i_atomic=1,monomer_one_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = grad_cut_one / real(monomer_one_size,dp)
                           end do
                           do i_atomic=1,monomer_two_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+i_atomic) = grad_cut_two / real(monomer_two_size,dp)
                           end do
                           do i_atomic=1,monomer_three_size
                              descriptor_out%x(i_desc)%grad_covariance_cutoff(:,dimer_size+i_atomic) = grad_cut_three / real(monomer_three_size,dp)
                           end do
                        else
                           primitive_cutoff_grads=0.0_dp
                           do k=1,d
                              if (this%cutoff_contributor(k)) then
                                primitive_cutoff_grads(k) = dcoordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
                              end if
                           end do

                           do i_atomic=1,trimer_size
                             do k=1,d
                                if (this%cutoff_contributor(k) ) then
                                  do m=k+1,d
                                    if (this%cutoff_contributor(m)) then
                                      if (any(this%component_atoms(k,:) == i_atomic) .or. any(this%component_atoms(m,:) == i_atomic)) then
                                      descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
                                        + primitive_cutoffs(m)*primitive_cutoff_grads(k)*descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
                                        + primitive_cutoffs(k)*primitive_cutoff_grads(m)*descriptor_out%x(i_desc)%grad_data(m,:,i_atomic)
                                      end if
                                    end if
                                  end do
                                end if
                             end do
                           end do

                           !normalisation factor
                           descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) / n_products
                        end if cutoff_type_grad

                     end if calc_grad_descriptor
                     write(13,'(a,*(f10.5))') "  desc = ", descriptor_out%x(i_desc)%data
                     write(13,'(a,*(f10.5))') "g desc = ", descriptor_out%x(i_desc)%grad_data
                     write(13,'(a,*(f10.5))') "   cut = ", descriptor_out%x(i_desc)%covariance_cutoff
                     write(13,'(a,*(f10.5))') "g  cut = ", descriptor_out%x(i_desc)%grad_covariance_cutoff
                     stop "testing"
                

                  end do loop_trimers
               end do loop_triplets
            end do loop_monomer_pairs
         end do loop_monomer_one
      end if has_descriptors ! (n_descriptors > 0) ... still need to deallocate if no trimers were found:

     deallocate(monomer_one_index)
     deallocate(monomer_two_index)
     deallocate(monomer_three_index)
     if (allocated(monomer_triplets)) deallocate(monomer_triplets)
     if (allocated(triplets_diffs)) deallocate(triplets_diffs)
     if (allocated(triplets_diffs_map)) deallocate(triplets_diffs_map)

     deallocate(dist_vec)
     deallocate(primitive_cutoffs)
     deallocate(primitive_cutoff_grads)

     deallocate(atomic_index_dimer)
     deallocate(atomic_index_trimer)
     deallocate(atomic_index_one)
     deallocate(atomic_index_two)
     deallocate(atomic_index_three)
     deallocate(diffs_one)
     deallocate(diffs_two)
     deallocate(diffs_three)

     deallocate(interatomic_vectors)
     deallocate(interatomic_distances)
     deallocate(associated_to_monomer)


     call system_timer('general_trimer_calc')

   endsubroutine general_trimer_calc





#include "fast_water_trimer.inc"








   
!    subroutine water_trimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error) !//////////////////////////////////////////////////////////////////////////////////
!      type(water_trimer), intent(in) :: this
!      type(atoms), intent(in) :: at
!      type(descriptor_data), intent(out) :: descriptor_out
!      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
!      character(len=*), intent(in), optional :: args_str
!      integer, optional, intent(out) :: error
!      integer, parameter :: monomer_one_size = 3, monomer_two_size = 3, monomer_three_size = 3, dimer_size = 6, trimer_size = 9, d = 36 , n_index = 9 
!     !  real(dp) :: r_one_two
!     !  type(Dictionary) :: params
!      logical :: my_do_descriptor, my_do_grad_descriptor, use_smooth_cutoff, done_this_dimer 
!      !done_this_monomer, one_two_identical, one_three_identical, two_three_identical, 
!      integer :: n_descriptors, i, j, k, m !n, n_cross, 
!      integer :: i_atomic, j_atomic, start, finish, i_desc,this_diff!, n_index !triplet_pos,
!      integer :: n_monomer_one, n_monomer_two, n_monomer_three, n_products!, contributor_loc
!     !  integer, dimension(1) :: unit_array
!     !  real(dp) :: temp_dist
!      real(dp), dimension(3) :: diff_one_two, diff_one_three,diff_two_three,mean_pos_one,mean_pos_two,mean_pos_three
!     !  real(dp), dimension(3) :: grad_cut_one, grad_cut_two, grad_cut_three
!     !  real(dp) :: dist_one_two, dist_one_three, dist_two_three!, cut12, cut13, cut23, dcut12, dcut13, dcut23
!     !  real(dp), dimension(:), allocatable :: dist_vec, primitive_cutoffs, primitive_cutoff_grads
!      real(dp), dimension(:,:), allocatable :: triplets_diffs !interatomic_distances,diffs_one,diffs_two,diffs_three,
!     !  real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
!     !  integer, dimension(3) :: temp_shift!, shift_one_two,shift_one_three
!      integer, dimension(:), allocatable :: triplets_diffs_map !atomic_index_dimer, atomic_index_trimer, atomic_index_one, atomic_index_two, atomic_index_three, 
!      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_triplets , monomer_two_index, monomer_three_index
!     !  logical, dimension(:), allocatable :: associated_to_monomer
!      logical :: double_count
!      real(dp) :: dist_vec(d)
!      real(dp) :: primitive_cutoffs(d)
!      real(dp) :: primitive_cutoff_grads(d)
!      integer :: atomic_index_one(monomer_one_size)
!      integer :: atomic_index_two(monomer_two_size)
!      integer :: atomic_index_three(monomer_three_size)
!      integer :: atomic_index_dimer(dimer_size)
!      integer :: atomic_index_trimer(trimer_size)
!      real(dp) :: interatomic_vectors(trimer_size,trimer_size,3)
!      real(dp) :: interatomic_distances(trimer_size,trimer_size)
!      logical :: associated_to_monomer(at%N)
!      real(dp) :: diffs_one(3,monomer_one_size)
!      real(dp) :: diffs_two(3,monomer_two_size)
!      real(dp) :: diffs_three(3,monomer_three_size)
!      integer n_mon_trip
!      INIT_ERROR(error)
     
!     !  print*, "Entered water_trimer_calc"
!      use_smooth_cutoff = .false.
!      double_count=.false.
!      call system_timer('water_trimer_calc')

!      if(.not. this%initialised) then
!         RAISE_ERROR("water_trimer_calc: descriptor object not initialised", error)
!      endif

!      my_do_descriptor = optional_default(.false., do_descriptor)
!      my_do_grad_descriptor = .true. !optional_default(.false., do_grad_descriptor) ! fix!!!


!      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

!      call finalise(descriptor_out)

!      !water_trimer_dimensions(this,error)
!     !  print*, "water_trimer_dimensions(this,error) =",d


!      interatomic_vectors = 0
!      interatomic_distances = 0
!      associated_to_monomer=.false. 

!      call find_general_monomer(at,monomer_one_index,&
!           this%signature,associated_to_monomer,&
!           this%monomer_cutoff,.false.,error) !this%atom_ordercheck=false
!     ! Allocation on assignment
!     block 
!         integer s1,s2
!         s1 = size(monomer_one_index,1)
!         s2 = size(monomer_one_index,2)
!         allocate(monomer_two_index(s1,s2))
!         allocate(monomer_three_index(s1,s2))
!         monomer_two_index   = monomer_one_index
!         monomer_three_index = monomer_one_index
!     endblock
        
!      n_monomer_one = size(monomer_one_index,2) !, monomer_two_index, monomer_three_index
!      n_monomer_two = size(monomer_two_index,2)
!      n_monomer_three = size(monomer_three_index,2)

!     call system_timer('water_trimer_calc: find_monomer_triplets')
!     call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
!             monomer_one_index,monomer_two_index,monomer_three_index,&
!             ! this%one_two_identical,this%one_three_identical,this%two_three_identical,&
!             .true.,.true.,.true.,this%cutoff,error)
    
!      n_descriptors = size(triplets_diffs_map)
!      n_mon_trip = size(monomer_triplets,2)
     
!      if (n_descriptors /= n_mon_trip) then 
!         RAISE_ERROR("water_triplets_calc: #diffs="//n_descriptors//" and #tiplets="//n_mon_trip//", make sure box is orthogonal and > 2*cutoff ",error)
!      endif
     
!      print*, "  wt calc: found #mono="//n_monomer_one//" #trip="//n_mon_trip
     
!      call system_timer('water_trimer_calc: find_monomer_triplets')
!      call print("monomer_triplets("//size(monomer_triplets,1)//", "//size(monomer_triplets,2)//")", PRINT_NERD)
!      call print("triplets_diffs("//size(triplets_diffs,1)//", "//size(triplets_diffs,2)//")", PRINT_NERD)
!      call print("triplets_diffs_map("//size(triplets_diffs_map)//")", PRINT_NERD)

!      !call print("monomer_triplets",PRINT_NERD)
!      !call print(monomer_triplets,PRINT_NERD)

!      !size(this%signature_one) + size(this%signature_two) + size(this%signature_three)

!      call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
!      allocate(descriptor_out%x(n_descriptors))
!      do i = 1, n_descriptors
!         if(my_do_descriptor) then
!            allocate(descriptor_out%x(i)%data(d))
!            allocate(descriptor_out%x(i)%ci(n_index))
!            descriptor_out%x(i)%data = 0.0_dp
!            descriptor_out%x(i)%ci = 0
!            descriptor_out%x(i)%has_data = .false.
!            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
!         endif
!         if(my_do_grad_descriptor) then
!            allocate(descriptor_out%x(i)%grad_data(d,3,trimer_size))
!            allocate(descriptor_out%x(i)%ii(trimer_size))
!            allocate(descriptor_out%x(i)%pos(3,trimer_size))
!            allocate(descriptor_out%x(i)%has_grad_data(trimer_size))
!            descriptor_out%x(i)%grad_data = 0.0_dp
!            descriptor_out%x(i)%ii = 0
!            descriptor_out%x(i)%pos = 0.0_dp
!            descriptor_out%x(i)%has_grad_data = .false.

!            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,trimer_size))
!            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
!         endif
!      enddo

!     !   has_descriptors: if (n_descriptors > 0) then ! only loop over monomers if we actually found any trimers
        
!     the_loop : do i_desc = 1, n_descriptors
!             i = monomer_triplets(1,i_desc)
!             j = monomer_triplets(2,i_desc)
!             k = monomer_triplets(3,i_desc)
!         !  loop_monomer_one: do i = 1, n_monomer_one
!             ! if (.not. any(monomer_triplets(1,:) .eq. i)) cycle
!             ! done_this_monomer = .false.
!             !get indices of monomer and calc internal distances
!             atomic_index_one = monomer_one_index(:,i) !store the indices of atoms in this monomer
!             atomic_index_two = monomer_two_index(:,j) 
!             atomic_index_three = monomer_three_index(:,k)
            
!             atomic_index_dimer=(/atomic_index_one,atomic_index_two/)
!             atomic_index_trimer=(/atomic_index_dimer,atomic_index_three/)
            
!             mean_pos_one = calc_mean_pos(at,atomic_index_one)
!             mean_pos_two = calc_mean_pos(at,atomic_index_two)
!             mean_pos_three = calc_mean_pos(at,atomic_index_three)
            
            
!             this_diff = i_desc
                     
!             diff_one_two = triplets_diffs(1:3,this_diff)
!             diff_one_three = triplets_diffs(4:6,this_diff)
!             diff_two_three = diff_one_three - diff_one_two
!            !  triplets_diffs_map(this_diff)=0 ! make sure this shifted triplet isn't found again
!             ! dist_one_two = norm(diff_one_two)
!             ! dist_one_three = norm(diff_one_three)
!             ! dist_two_three = norm(diff_two_three)
            
            
!             do i_atomic=1,3
!                diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
!                diffs_two(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
!                diffs_three(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_three(i_atomic)),mean_pos_three)
!             end do
            
            
!             block 
!                 real(dp) a0, shft, pow, cut
!                 real(dp), dimension(3,9) :: pos, gcuts
!                 real(dp), dimension(3) :: a3, b3 , cuts
!                 real(dp), dimension(36) :: desc0, rlen
!                 real(dp), dimension(3,36) :: rvec2, rvec
!                 ! real(dp), dimension(3,9,36) :: gdesc0
!                 real(dp), dimension(36,3,9) :: gdesc0
!                 integer ii, jj, ij
!                 do i_atomic=1,3
!                     pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
!                     pos(:,3+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
!                     pos(:,6+i_atomic) = mean_pos_one + diff_one_three - diffs_three(:,i_atomic)
!                 end do
!                 ! ij = 0
!                 ! do jj = 1,9
!                 !     do ii = jj+1, 9
!                 !         ij = ij + 1
!                 !         a3 = pos(:,ii)-pos(:,jj)
!                 !         rvec(:,ij) = a3
!                 !         rlen(ij) = norm(a3)
!                 !     enddo
!                 ! enddo
!                 shft  = this%dist_shift
!                 pow   = this%power
!                 gdesc0 = 0
!                 gcuts = 0
!                 cuts = 0
!                 do ij = 1, 36
!                     jj = this%component_atoms(ij,1)
!                     ii = this%component_atoms(ij,2)
!                     a3 = pos(:,ii)-pos(:,jj)
!                     a0 = norm(a3)
!                     rvec(:,ij) = a3
!                     rlen(ij) = a0
!                     desc0(ij) = (a0 + shft) ** pow
                    
!                     b3 = - this%power * (a0+shft)**(pow-1.0_dp)  * a3 / a0
!                     gdesc0(ij,:,jj) = b3
!                     gdesc0(ij,:,ii) = -b3
                    
!                 enddo
                
!                 block 
!                     real(dp), dimension(3) :: oo12, oo13, oo23
!                     real(dp), dimension(3) :: gcut1, gcut2, gcut3
!                     real(dp) :: o12, o13, o23
!                     real(dp) :: cut12, cut13, cut23
!                     real(dp) :: dcut12, dcut13, dcut23
!                     oo12 = pos(:,4)-pos(:,1)
!                     oo13 = pos(:,7)-pos(:,1)
!                     oo23 = pos(:,7)-pos(:,4)
!                     o12 = norm(oo12)
!                     o13 = norm(oo13)
!                     o23 = norm(oo23)
!                     cut12 = coordination_function(o12,this%cutoff,this%cutoff_transition_width)
!                     cut13 = coordination_function(o13,this%cutoff,this%cutoff_transition_width)
!                     cut23 = coordination_function(o23,this%cutoff,this%cutoff_transition_width)
!                     cut = (cut12*cut13 + cut12*cut23 + cut13*cut23) / 3
                    
                
!                     dcut12 = dcoordination_function(o12, this%cutoff, this%cutoff_transition_width)
!                     dcut13 = dcoordination_function(o13, this%cutoff, this%cutoff_transition_width)
!                     dcut23 = dcoordination_function(o23, this%cutoff, this%cutoff_transition_width)
!                     gcut1 = (- oo12 / o12 * dcut12 * (cut13+cut23) - oo13 / o13 * dcut13 * (cut12+cut23)) / 3.0_dp
!                     gcut2 = (  oo12 / o12 * dcut12 * (cut13+cut23) - oo23 / o23 * dcut23 * (cut12+cut13)) / 3.0_dp
!                     gcut3 = (  oo13 / o13 * dcut13 * (cut12+cut23) + oo23 / o23 * dcut23 * (cut12+cut13)) / 3.0_dp
                    
!                     gcuts = 0 
!                     gcuts(:,1) = gcut1
!                     gcuts(:,4) = gcut2
!                     gcuts(:,7) = gcut3
!                 endblock
                
!                         !    cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
!                         !    cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
!                         !    cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
!                         !    dcut12 = dcoordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
!                         !    dcut13 = dcoordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
!                         !    dcut23 = dcoordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
!                         !    !descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
!                         !    grad_cut_one = (- diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
!                         !                    - diff_one_three / dist_one_three * dcut13 * (cut12+cut23)) / 3.0_dp
!                         !    grad_cut_two = (  diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
!                         !                    - diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
!                         !    grad_cut_three = (diff_one_three / dist_one_three * dcut13 * (cut12+cut23) &
!                         !                    + diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
                
                
                
!                 ! roo = rlen(this%cutoff_contributor)
!                 ! cuts = 
                
!                 ! do k=1,d
!                 !     if (this%cutoff_contributor(k)) then
!                 !         primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
!                 !     end if
!                 ! end do
                
!                 ! do k=1,d
!                 !     descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
!                 !          = - this%power * (dist_vec(k)+this%dist_shift)**(this%power-1.0_dp) &
!                 !              * interatomic_vectors(i_atomic,j_atomic,:) &
!                 !              / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
!                 !  end do
         
                
!                 !##########################################################################################################################################################
!                 !##########################################################################################################################################################
!                 !##########################################################################################################################################################
                
                
                
                
                
!             ! print'(a,*(f10.5))', "rlen0 = ",rlen 
!                 ! print'(a,*(f10.5))', "rvec0 = ",rvec 
            
            
!             !calc atomic positions and shifts relative to mean pos for monomer one, and also distances wrt atom 1
!             do i_atomic=1,monomer_one_size
!                diffs_one(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_one(i_atomic)),mean_pos_one)
!                interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index_one(1),atomic_index_one(i_atomic))
!             end do

!             !find other relative positions through vector addition
!             do j_atomic=2,monomer_one_size
!                do i_atomic=2,j_atomic-1
!                   interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:)
!                end do
!             end do

!             !Now convert vectors to scalar distances
!             do i_atomic=1,monomer_one_size
!                do j_atomic=i_atomic+1,monomer_one_size
!                   interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                end do
!             end do

!             ! Loop through monomers paired with this one to make dimers
!             ! loop_monomer_pairs: do
!             !    unit_array = maxloc(monomer_triplets(2,:), monomer_triplets(1,:) .eq. i) ! find a monomer paired with i
!             !    if (all(unit_array .eq. 0)) exit

!                !get indices of monomer two
!             !    j = monomer_triplets(2,unit_array(1))


!                ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
!                do i_atomic=1,monomer_two_size
!                   diffs_two(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_two(i_atomic)),mean_pos_two)
!                   interatomic_vectors(monomer_one_size+1,monomer_one_size+i_atomic,:) = diff_min_image(at,atomic_index_two(1),atomic_index_two(i_atomic))
!                end do

!                !find other relative positions through vector addition
!                do j_atomic=monomer_one_size+2,dimer_size
!                   do i_atomic=monomer_one_size+2,j_atomic-1
!                      interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(monomer_one_size+1,j_atomic,:) - interatomic_vectors(monomer_one_size+1,i_atomic,:)
!                   end do
!                end do

!                !And convert vectors to scalar distances
!                do i_atomic=monomer_one_size+1,dimer_size
!                   do j_atomic=i_atomic+1,dimer_size
!                      interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                   end do
!                end do

!                !! Now make trimers based on this dimer
!             !    loop_triplets: do
!             !       unit_array = maxloc(monomer_triplets(3,:), monomer_triplets(1,:) .eq. i .and. monomer_triplets(2,:) .eq. j ) ! look for  a triplet
!             !       if (all(unit_array .eq. 0)) exit
!             !       print*, "   wt calc: found triplet!"
!             !       triplet_pos=unit_array(1)

!             !       !get indices of monomer three
!             !       k = monomer_triplets(3,triplet_pos)
!             !       monomer_triplets(:,triplet_pos) = 0 ! make sure this triplet isn't found again
!                   ! calc distances and shifts wrt to mean pos this monomer, and distances wrt its first atom
!                   do i_atomic=1,monomer_three_size
!                      diffs_three(:,i_atomic) = diff_min_image(at,at%pos(:,atomic_index_three(i_atomic)),mean_pos_three)
!                      interatomic_vectors(dimer_size+1,dimer_size+i_atomic,:) = diff_min_image(at,atomic_index_three(1),atomic_index_three(i_atomic))
!                   end do

!                   !find other relative positions through vector addition
!                   do j_atomic=dimer_size+2,trimer_size
!                      do i_atomic=dimer_size+2,j_atomic-1
!                         interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(dimer_size+1,j_atomic,:) -interatomic_vectors(dimer_size+1,i_atomic,:)
!                      end do
!                   end do

!                   !Now convert vectors to scalar distances
!                   do i_atomic=dimer_size+1,trimer_size
!                      do j_atomic=i_atomic+1,trimer_size
!                         interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                      end do
!                   end do

!                   ! Loop over all trimers which can be created from these three monomers, i.e. with different periodic images
!                 !   loop_trimers: do
!                 !      unit_array = maxloc(triplets_diffs_map, triplets_diffs_map .eq. triplet_pos)
!                 !      this_diff = unit_array(1)
!                 !      if (this_diff == 0) exit
                       
!                 !      i_desc = i_desc + 1
                     

!                      ! calculate intermolecular distances, monomers one and two
!                      do i_atomic=1,monomer_one_size
!                         do j_atomic = monomer_one_size+1,dimer_size
!                            interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_two - diffs_two(:,j_atomic-monomer_one_size)
!                            interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                         end do
!                      end do

!                      ! calculate intermolecular distances, monomers one and three
!                      do i_atomic=1,monomer_one_size
!                         do j_atomic = dimer_size+1,trimer_size
!                            interatomic_vectors(i_atomic,j_atomic,:) = diffs_one(:,i_atomic) + diff_one_three - diffs_three(:,j_atomic-dimer_size)
!                            interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                         end do
!                      end do

!                      ! calculate intermolecular distances, monomers two and three
!                      do i_atomic=monomer_one_size+1,dimer_size
!                         do j_atomic = dimer_size+1,trimer_size
!                            interatomic_vectors(i_atomic,j_atomic,:) = diffs_two(:,i_atomic-monomer_one_size) + diff_two_three - diffs_three(:,j_atomic-dimer_size)
!                            interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
!                         end do
!                      end do

!                      !Now take the whole matrix of scalar distances and combine into 1D array
!                      start = 1
!                      finish=trimer_size-1
!                      do i_atomic=1,trimer_size-1
!                         dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:trimer_size)
!                         start = finish+1
!                         finish=finish + trimer_size-i_atomic-1
!                      end do
!                      call print( "list of distances: "//dist_vec,PRINT_NERD)
                     
                     
!                         ij = 0
!                         do ii = 1,9
!                             do jj = ii+1, 9
!                                 ij = ij + 1
!                                 rvec2(:,ij) = interatomic_vectors(ii,jj,:)
!                             enddo
!                         enddo
!                         print'(a,*(f10.5))', "rlen d = ",dist_vec-rlen
!                         print'(a,*(f10.5))', "rvec d = ",rvec2-rvec 
                                  
                    
                     
                     
!                      calc_descriptor: if(my_do_descriptor) then
!                         descriptor_out%x(i_desc)%has_data = .true.
!                         descriptor_out%x(i_desc)%data = (dist_vec+this%dist_shift)**this%power
!                         descriptor_out%x(i_desc)%ci(:) = atomic_index_trimer
                        
!                         print'(a,*(f10.5))', "desc d = ",desc0-descriptor_out%x(i_desc)%data
                        
                        
                        
                        
                        
!                         ! cutoff_type: if (this%mpifind) then
!                         !    ! Using average position for cutoff
!                         !    cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
!                         !    cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
!                         !    cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
!                         !    descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
!                         ! else
!                            ! Using heavy atoms for cutoff
!                            descriptor_out%x(i_desc)%covariance_cutoff = 0.0_dp

!                            primitive_cutoffs=1.0_dp

!                            do k=1,d
!                               if (this%cutoff_contributor(k)) then
!                                  primitive_cutoffs(k) = coordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
!                               end if
!                            end do
!                            n_products=0
!                            do k=1,d
!                               if (this%cutoff_contributor(k)) then
!                                  do m=k+1,d
!                                     if (this%cutoff_contributor(m)) then
!                                        n_products=n_products+1
!                                        descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff + primitive_cutoffs(k)*primitive_cutoffs(m)
!                                     end if
!                                  end do
!                               end if
!                            end do
!                            ! normalise
!                            descriptor_out%x(i_desc)%covariance_cutoff =  descriptor_out%x(i_desc)%covariance_cutoff / n_products
!                         ! end if cutoff_type
                           
                           
!                            print'(a,*(f10.5))', "cut d = ", cut, descriptor_out%x(i_desc)%covariance_cutoff, cut - descriptor_out%x(i_desc)%covariance_cutoff
                           
!                      end if calc_descriptor
                     
                     

!                      calc_grad_descriptor: if(my_do_grad_descriptor) then !calc grads and update

!                         descriptor_out%x(i_desc)%ii(:) = atomic_index_trimer
!                         do i_atomic=1,monomer_one_size
!                            descriptor_out%x(i_desc)%pos(:,i_atomic) = mean_pos_one - diffs_one(:,i_atomic)
!                         end do
!                         do i_atomic=1,monomer_two_size
!                            descriptor_out%x(i_desc)%pos(:,monomer_one_size+i_atomic) = mean_pos_one + diff_one_two - diffs_two(:,i_atomic)
!                         end do
!                         do i_atomic=1,monomer_three_size
!                            descriptor_out%x(i_desc)%pos(:,dimer_size+i_atomic) = mean_pos_one + diff_one_three - diffs_three(:,i_atomic)
!                         end do

!                         call print("TRIM  9",PRINT_NERD)
!                         call print('TRIM  Lattice="10.0000000       0.00000000       0.00000000       0.00000000      10.0000000       0.00000000       0.00000000       0.00000000      10.00000" Properties=Z:I:1:pos:R:3',PRINT_NERD)
!                         do i_atomic=1,trimer_size
!                             call print("TRIM  "//at%Z(atomic_index_trimer(i_atomic))//"    "//descriptor_out%x(i_desc)%pos(:,i_atomic),PRINT_NERD)
!                         end do

!                         !build the grad_data matrix
!                         descriptor_out%x(i_desc)%has_grad_data(:) = .true.
!                         do k=1,d
!                            i_atomic = this%component_atoms(k,1)
!                            j_atomic = this%component_atoms(k,2)
!                            descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
!                                 = - this%power * (dist_vec(k)+this%dist_shift)**(this%power-1.0_dp) &
!                                     * interatomic_vectors(i_atomic,j_atomic,:) &
!                                     / interatomic_distances(i_atomic,j_atomic)  ! descriptor wrt atom i_atomic
!                            descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! descriptor wrt j_atomic
!                         end do
                        
                        
!                         print'(a,*(f10.5))', "gdesc d = ",gdesc0(:,1,1)-descriptor_out%x(i_desc)%grad_data(:,1,1)
                        
                        
                        
                        
                        
                        
!                         ! cutoff_type_grad: if (this%mpifind) then
!                         !    !! Using the average position of all atoms in molecule
!                         !    cut12 = coordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
!                         !    cut13 = coordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
!                         !    cut23 = coordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
!                         !    dcut12 = dcoordination_function(dist_one_two, this%cutoff, this%cutoff_transition_width)
!                         !    dcut13 = dcoordination_function(dist_one_three, this%cutoff, this%cutoff_transition_width)
!                         !    dcut23 = dcoordination_function(dist_two_three, this%cutoff, this%cutoff_transition_width)
!                         !    !descriptor_out%x(i_desc)%covariance_cutoff = (cut12*cut13 + cut12*cut23 + cut13*cut23)/3.0_dp
!                         !    grad_cut_one = (- diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
!                         !                    - diff_one_three / dist_one_three * dcut13 * (cut12+cut23)) / 3.0_dp
!                         !    grad_cut_two = (  diff_one_two / dist_one_two * dcut12 * (cut13+cut23) &
!                         !                    - diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
!                         !    grad_cut_three = (diff_one_three / dist_one_three * dcut13 * (cut12+cut23) &
!                         !                    + diff_two_three / dist_two_three * dcut23 * (cut12+cut13)) / 3.0_dp
!                         !    do i_atomic=1,monomer_one_size
!                         !       descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = grad_cut_one / real(monomer_one_size,dp)
!                         !    end do
!                         !    do i_atomic=1,monomer_two_size
!                         !       descriptor_out%x(i_desc)%grad_covariance_cutoff(:,monomer_one_size+i_atomic) = grad_cut_two / real(monomer_two_size,dp)
!                         !    end do
!                         !    do i_atomic=1,monomer_three_size
!                         !       descriptor_out%x(i_desc)%grad_covariance_cutoff(:,dimer_size+i_atomic) = grad_cut_three / real(monomer_three_size,dp)
!                         !    end do
!                         ! else
                        
!                            ! Using only the heavy atoms for the cutoff. 
!                            primitive_cutoff_grads=0.0_dp
!                            do k=1,d
!                               if (this%cutoff_contributor(k)) then
!                                 primitive_cutoff_grads(k) = dcoordination_function(dist_vec(k),this%cutoff,this%cutoff_transition_width)
!                               end if
!                            end do
                           
!                            ! For each atom, for each heavy-atomc pairs it is part of, accumulate cutoff gradients
!                            do i_atomic=1,trimer_size
!                             ! if (signatures(i_atom)/=1)then !? Jonatan
!                              do k=1,d
!                                 if (this%cutoff_contributor(k) ) then
!                                   do m=k+1,d
!                                     if (this%cutoff_contributor(m)) then
!                                       if (any(this%component_atoms(k,:) == i_atomic) .or. any(this%component_atoms(m,:) == i_atomic)) then
!                                       descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,i_atomic) &
!                                         + primitive_cutoffs(m)*primitive_cutoff_grads(k)*descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) &
!                                         + primitive_cutoffs(k)*primitive_cutoff_grads(m)*descriptor_out%x(i_desc)%grad_data(m,:,i_atomic)
!                                       end if
!                                     end if
!                                   end do
!                                 end if
!                              end do
!                            end do

!                            !normalisation factor
!                            descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) / n_products
!                         ! end if cutoff_type_grad
                           
                           
                           
!                         !    print'(a,*(f10.5))', "gcuts 0 = ",gcuts
!                         !    print'(a,*(f10.5))', "gcuts 1 = ",descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:)
!                            print'(a,*(f10.5))', "gcuts d = ",gcuts - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:)
                           
                        
                        
!                            !##########################################################################################################################################################
!                            !##########################################################################################################################################################
!                            !##########################################################################################################################################################
!                            stop "testing"
                              
                           

!                      end if calc_grad_descriptor
!                      endblock
                     
!                      stop "testing 2"


!                 !   end do loop_trimers
!             !    end do loop_triplets
!             ! end do loop_monomer_pairs
!         !  end do loop_monomer_one
!     enddo the_loop
!     !   end if has_descriptors ! (n_descriptors > 0) ... still need to deallocate if no trimers were found:

!      deallocate(monomer_one_index)
!      deallocate(monomer_two_index)
!      deallocate(monomer_three_index)
!      if (allocated(monomer_triplets)) deallocate(monomer_triplets)
!      if (allocated(triplets_diffs)) deallocate(triplets_diffs)
!      if (allocated(triplets_diffs_map)) deallocate(triplets_diffs_map)



!      call system_timer('water_trimer_calc')

!    endsubroutine water_trimer_calc !////////////////////////////////////////////////////////////////////////////////////////////////////////////














   subroutine molecule_lo_d_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(molecule_lo_d), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor!, use_smooth_cutoff
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i, j, i_atomic, j_atomic, k, &
         start, finish, molecule_size, i_component, i_desc, n_index
      integer, dimension(3) :: temp_shift
      real(dp), dimension(:), allocatable :: dist_vec
      real(dp), dimension(:,:), allocatable :: interatomic_distances
      real(dp), dimension(:,:,:), allocatable :: interatomic_vectors
      integer, dimension(:), allocatable :: atomic_index
      integer, dimension(:,:), allocatable :: shifts
      logical, dimension(:), allocatable :: associated_to_molecule


      INIT_ERROR(error)

      call system_timer('molecule_lo_d_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_calc: descriptor object not initialised", error)
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

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='molecule_lo_d_calc args_str')) then
            RAISE_ERROR("molecule_lo_d_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("molecule_lo_d_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
            RAISE_ERROR("molecule_lo_d_calc cannot use atom masks yet.",error)
         else
            atom_mask_pointer => null()
         endif

      endif

      molecule_size=this%n_atoms
      d = molecule_lo_d_dimensions(this,error)

      if (this%atom_ordercheck) then
        do i=1,molecule_size
          if (this%template_atoms%Z(i) /= at%Z(i)) then
            call print("atoms in all input frams should be in the same order as this template : ")
            call print(this%template_atoms)
            RAISE_ERROR("molecule_lo_d_calc: atoms not in same order as in template used to teach", error)
          end if
        end do
      end if

      allocate(shifts(molecule_size,3))
      allocate(dist_vec(this%max_dimension))
      allocate(atomic_index(molecule_size))

      allocate(associated_to_molecule(at%N))
      associated_to_molecule=.True.
      do i=1,molecule_size
        atomic_index(i) = i
      end do


      allocate(interatomic_vectors(molecule_size,molecule_size,3))
      allocate(interatomic_distances(molecule_size,molecule_size))
      interatomic_vectors = 0.0_dp
      interatomic_distances = 0.0_dp

      !      call find_general_monomer(at,monomer_index,&
      !            this%signature,associated_to_monomer,&
      !            this%cutoff,this%atom_ordercheck,error)

      if(.not. all(associated_to_molecule)) then
         RAISE_ERROR("molecule_lo_d_calc: not all atoms assigned to a monomer", error)
      endif

      n_descriptors = 1 ! currently no support for finding multiple nontrivial molecules, otherwise would be size(monomer_index,2)
      n_index = this%n_atoms

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            allocate(descriptor_out%x(i)%ci(n_index))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,molecule_size))
            allocate(descriptor_out%x(i)%ii(molecule_size))
            allocate(descriptor_out%x(i)%pos(3,molecule_size))
            allocate(descriptor_out%x(i)%has_grad_data(molecule_size))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,molecule_size))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo


      do i_desc = 1, n_descriptors

         !atomic_index = monomer_index(:,i) !for fixed ordering don't need this
         ! for now calculate all O(N^2) distances and just pick out the ones we want

         !calc all positions relative to atom 1
         do i_atomic=2,molecule_size
           temp_shift=0
           interatomic_vectors(1,i_atomic,:) = diff_min_image(at,atomic_index(1),atomic_index(i_atomic),shift=temp_shift)
           shifts(i_atomic,:) = temp_shift
         end do

         !find other relative positions through vector addition
         do j_atomic=2,molecule_size
           do i_atomic=2,j_atomic-1
             interatomic_vectors(i_atomic,j_atomic,:) = interatomic_vectors(1,j_atomic,:) -interatomic_vectors(1,i_atomic,:) 
           end do
         end do

         !Now convert vectors to scalar distances
         do i_atomic=1,molecule_size
           do j_atomic=i_atomic+1,molecule_size
             interatomic_distances(i_atomic,j_atomic) = norm(interatomic_vectors(i_atomic,j_atomic,:))
           end do
         end do       

         !and convert this NxN matrix into the required vector length N(N-1)/2
         start = 1
         finish = molecule_size-1
         do i_atomic=1,molecule_size-1
           dist_vec(start:finish) = interatomic_distances(i_atomic,i_atomic+1:molecule_size)  
           start = finish+1
           finish=finish + molecule_size-i_atomic-1
         end do



         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%ci(:) = atomic_index
            descriptor_out%x(i_desc)%has_data = .true.

            if(this%distance_transform == 0) then                                  ! no transform
               descriptor_out%x(i_desc)%data = dist_vec(this%included_components)  ! no transform
            else if(this%distance_transform == 1) then                                    ! inverse of distance
               descriptor_out%x(i_desc)%data = 1.0_dp / dist_vec(this%included_components)! inverse of distance
            else if(this%distance_transform == 2) then                             ! Coulomb matrix
               do k=1,d                                                            ! Coulomb matrix
                  i_component = this%included_components(k)                        ! Coulomb matrix
                  i_atomic = this%component_atoms(i_component,1)                   ! Coulomb matrix
                  j_atomic = this%component_atoms(i_component,2)                   ! Coulomb matrix
                  descriptor_out%x(i_desc)%data(k) = real(at%Z(i_atomic)*at%Z(j_atomic)) / dist_vec(i_component) ! Coulomb matrix
               end do                                                              ! Coulomb matrix
            else if(this%distance_transform == 3) then                                 ! exponential
               descriptor_out%x(i_desc)%data = exp(-dist_vec(this%included_components))! exponential
            else if(this%distance_transform < 0) then                              ! negative power
               descriptor_out%x(i_desc)%data = dist_vec(this%included_components)**this%distance_transform
            else ! unknown
               RAISE_ERROR("molecule_lo_d_calc: not implemented distance transform",error)
            end if ! distance_transform

         endif

!call print(descriptor_out%x(i_desc)%data)
!
!        do i_component=1,d
!            write(*,*) 'component : ',i_component
!            write(*,*) 'value : ',descriptor_out%x(i_desc)%data(i_component)
!            write(*,*) 'atoms : ',this%component_atoms(this%included_components(i_component),:)
!            write(*,*) 'dist_mat_entry : ',interatomic_distances(this%component_atoms(this%included_components(i_component),1),this%component_atoms(this%included_components(i_component),2))
!        end do

         if(my_do_grad_descriptor) then

            descriptor_out%x(i_desc)%ii(:) = atomic_index
            descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,atomic_index(1))
            do i_atomic =2,molecule_size
              descriptor_out%x(i_desc)%pos(:,i_atomic) = at%pos(:,atomic_index(i_atomic)) + matmul(at%lattice,shifts(i_atomic,:))
            end do

            !build the grad_data matrix
            descriptor_out%x(i_desc)%has_grad_data(:) = .true.


            do k=1,d
              !get pair of atoms contributing to this component
              i_component = this%included_components(k)
              i_atomic = this%component_atoms(i_component,1)
              j_atomic = this%component_atoms(i_component,2)

              if(this%distance_transform == 0) then ! no transform
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)    ! derivative of descriptor wrt atom i_atomic 
              else if(this%distance_transform == 1) then ! inverse transform
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)**3  ! derivative of descriptor wrt atom i_atomic 
              else if(this%distance_transform == 2) then ! coulomb matrix
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = real(at%Z(i_atomic)*at%Z(j_atomic)) * interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic)**3  ! derivative of descriptor wrt atom i_atomic
              else if(this%distance_transform == 3) then ! exponential
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = interatomic_vectors(i_atomic,j_atomic,:) / interatomic_distances(i_atomic,j_atomic) * exp(-interatomic_distances(i_atomic,j_atomic)) ! derivative of descriptor wrt atom i_atomic
              else if(this%distance_transform < 0) then ! inverse negative power
                 descriptor_out%x(i_desc)%grad_data(k,:,i_atomic) = -real(this%distance_transform) * interatomic_vectors(i_atomic,j_atomic,:) * interatomic_distances(i_atomic,j_atomic)**(this%distance_transform-2) ! derivative of descriptor wrt atom i_atomic
              else  ! no transform
                 RAISE_ERROR("molecule_lo_d_calc: not implemented distance transform",error)
              end if ! distance_transform
              
              descriptor_out%x(i_desc)%grad_data(k,:,j_atomic) = -descriptor_out%x(i_desc)%grad_data(k,:,i_atomic)        ! derivative of descriptor wrt atom j_atomic is just the negative of that wrt atom i_atomic
              
!   call print(descriptor_out%x(i_desc)%grad_data(k,:,:))
            end do
       
         endif

      enddo
      deallocate(shifts)
      deallocate(dist_vec)
      deallocate(atomic_index)
      deallocate(interatomic_vectors)
      deallocate(interatomic_distances)

      call system_timer('molecule_lo_d_calc')

   endsubroutine molecule_lo_d_calc




   function AN_monomer_dimensions(this,error) result(i)
      type(AN_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_dimensions: descriptor object not initialised", error)
      endif

      i = this%N * (this%N - 1) / 2

   endfunction AN_monomer_dimensions

   function general_monomer_dimensions(this,error) result(i)
      type(general_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_monomer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_monomer_dimensions

   function com_dimer_dimensions(this,error) result(i)
      type(com_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction com_dimer_dimensions

   function general_dimer_dimensions(this,error) result(i)
      type(general_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_dimer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_dimer_dimensions


   function general_trimer_dimensions(this,error) result(i)
      type(general_trimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("general_trimer_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%permutation_data%dist_vec)

   endfunction general_trimer_dimensions


   function water_trimer_dimensions(this,error) result(i)
    type(water_trimer), intent(in) :: this
    integer, optional, intent(out) :: error
    integer :: i

    INIT_ERROR(error)

    if(.not. this%initialised) then
       RAISE_ERROR("water_trimer_dimensions: descriptor object not initialised", error)
    endif
    if(.not. this%permutation_data%initialised) then
       RAISE_ERROR("water_trimer_dimensions: descriptor object's permutation data not initialised", error)
    endif

    i = size(this%permutation_data%dist_vec)

 endfunction water_trimer_dimensions



   function molecule_lo_d_dimensions(this,error) result(i)
      type(molecule_lo_d), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_dimensions: descriptor object not initialised", error)
      endif
      if(.not. this%permutation_data%initialised) then
         RAISE_ERROR("molecule_lo_d_dimensions: descriptor object's permutation data not initialised", error)
      endif

      i = size(this%included_components)

   endfunction molecule_lo_d_dimensions



   function AN_monomer_cutoff(this,error) 
      type(AN_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: AN_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer_cutoff: descriptor object not initialised", error)
      endif

      AN_monomer_cutoff = this%cutoff

   endfunction AN_monomer_cutoff

   function general_monomer_cutoff(this,error) 
      type(general_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_cutoff: descriptor object not initialised", error)
      endif

      general_monomer_cutoff = this%cutoff

   endfunction general_monomer_cutoff

   function com_dimer_cutoff(this,error) 
      type(com_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: com_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_cutoff: descriptor object not initialised", error)
      endif

      com_dimer_cutoff = this%cutoff

   endfunction com_dimer_cutoff

   function general_dimer_cutoff(this,error) 
      type(general_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_cutoff: descriptor object not initialised", error)
      endif

      general_dimer_cutoff = this%cutoff

   endfunction general_dimer_cutoff

   function general_trimer_cutoff(this,error) 
      type(general_trimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: general_trimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_cutoff: descriptor object not initialised", error)
      endif

      general_trimer_cutoff = this%cutoff

   endfunction general_trimer_cutoff

   function water_trimer_cutoff(this,error) 
      type(water_trimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_trimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_trimer_cutoff: descriptor object not initialised", error)
      endif

      water_trimer_cutoff = this%cutoff

   endfunction water_trimer_cutoff

   
   
   function molecule_lo_d_cutoff(this,error) 
      type(molecule_lo_d), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: molecule_lo_d_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_cutoff: descriptor object not initialised", error)
      endif

      molecule_lo_d_cutoff = this%cutoff

   endfunction molecule_lo_d_cutoff


   subroutine AN_monomer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(AN_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AN_monomer: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + this%N
         if(.not.this%do_atomic) exit
      enddo

      if(this%do_atomic) then
         if( present(n_index) ) n_index = 1
      else
         if( present(n_index) ) n_index = this%N
      endif

   endsubroutine AN_monomer_sizes

   subroutine general_monomer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(general_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_index
      integer :: i
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_monomer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_index,&
           this%signature,associated_to_monomer,&
           this%cutoff,this%atom_ordercheck,error)
      n_descriptors = size(monomer_index,2)
      n_cross=size(monomer_index)
      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_monomer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_monomer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      deallocate(monomer_index)
      deallocate(associated_to_monomer)

      if( present(n_index) ) n_index = size(this%signature)

   endsubroutine general_monomer_sizes

   subroutine com_dimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)

      type(com_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_pairs
      integer, dimension(:), allocatable :: pairs_diffs_map
      real(dp), dimension(:,:), allocatable :: mean_pos_diffs
      logical, dimension(:), allocatable :: associated_to_monomer
      logical :: double_count

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("com_dimer_sizes: descriptor object not initialised", error)
      endif

      double_count = .false.

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,&
           this%signature_one,associated_to_monomer,&
           this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at,monomer_two_index,&
              this%signature_two,associated_to_monomer,&
              this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("com_dimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: com_dimer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,&
              monomer_one_index,monomer_two_index,this%monomers_identical,double_count,&
              this%cutoff,error=error,use_com=.true.)
      else
         call find_monomer_pairs(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,&
              monomer_one_index,monomer_two_index,this%monomers_identical,double_count,&
              this%cutoff,error=error,use_com=.true.)
      end if
      n_descriptors = size(pairs_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two))
   
      deallocate(associated_to_monomer)
      deallocate(pairs_diffs_map)
      deallocate(monomer_pairs)
      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(mean_pos_diffs)

      if( present(n_index) ) n_index = size(this%signature_one) + size(this%signature_two)

   endsubroutine com_dimer_sizes

   subroutine general_dimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)

      type(general_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_pairs
      integer, dimension(:), allocatable :: pairs_diffs_map
      real(dp), dimension(:,:), allocatable :: mean_pos_diffs
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_dimer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,&
           this%signature_one,associated_to_monomer,&
           this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%monomers_identical) then
        allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
        monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at,monomer_two_index,&
              this%signature_two,associated_to_monomer,&
              this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_dimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_dimer_sizes: not all atoms assigned to a monomer")
         endif
      endif

      if (this%mpifind) then
         call print("Using find_monomer_pairs_MPI", PRINT_NERD)
         call find_monomer_pairs_MPI(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,&
              monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,&
              this%cutoff,error=error,use_com=this%use_com)
      else
         call find_monomer_pairs(at,monomer_pairs,mean_pos_diffs,pairs_diffs_map,&
              monomer_one_index,monomer_two_index,this%monomers_identical,this%double_count,&
              this%cutoff,error=error,use_com=this%use_com)
      end if
      n_descriptors = size(pairs_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two))
   
      deallocate(associated_to_monomer)
      deallocate(pairs_diffs_map)
      deallocate(monomer_pairs)
      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(mean_pos_diffs)

      if( present(n_index) ) n_index = size(this%signature_one) + size(this%signature_two)
   endsubroutine general_dimer_sizes

   subroutine general_trimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
 
      type(general_trimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:), allocatable ::  triplets_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_two_index, monomer_three_index, monomer_triplets
      real(dp), dimension(:,:), allocatable :: triplets_diffs
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("general_trimer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,&
           this%signature_one,associated_to_monomer,&
           this%monomer_one_cutoff,this%atom_ordercheck,error)
      if (this%one_two_identical) then
         allocate(monomer_two_index(size(monomer_one_index,1),size(monomer_one_index,2)))
         monomer_two_index = monomer_one_index
      else
         call find_general_monomer(at,monomer_two_index,&
              this%signature_two,associated_to_monomer,&
              this%monomer_two_cutoff,this%atom_ordercheck,error)
      end if
      if (this%one_three_identical) then
         allocate(monomer_three_index(size(monomer_one_index,1),size(monomer_one_index,2)))
         monomer_three_index = monomer_one_index
      else if (this%two_three_identical) then
         allocate(monomer_three_index(size(monomer_two_index,1),size(monomer_two_index,2)))
         monomer_three_index = monomer_two_index
      else
         call find_general_monomer(at,monomer_three_index,&
              this%signature_three,associated_to_monomer,&
              this%monomer_three_cutoff,this%atom_ordercheck,error)
      end if

      if(.not. all(associated_to_monomer)) then
         if(this%strict) then
            RAISE_ERROR("general_trimer_sizes: not all atoms assigned to a monomer", error)
         else
            call print("WARNING: general_trimer_sizes: not all atoms assigned to a monomer")
         endif
      endif
    
      if (this%use_com) then
         RAISE_ERROR("general_trimer_calc: use_com=T not implemented yet", error)
      end if
      if(this%mpifind) then
         call print("Using find_monomer_triplets_MPI", PRINT_NERD)
         call find_monomer_triplets_MPI(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
              monomer_one_index,monomer_two_index,monomer_three_index,&
              this%one_two_identical,this%one_three_identical,this%two_three_identical,&
              this%cutoff,error,use_com=.false.)
      else
         call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
              monomer_one_index,monomer_two_index,monomer_three_index,&
              this%one_two_identical,this%one_three_identical,this%two_three_identical,&
              this%cutoff,error)
      end if
      n_descriptors = size(triplets_diffs_map)
      n_cross=n_descriptors*(size(this%signature_one)+size(this%signature_two)+size(this%signature_three))

      deallocate(monomer_one_index)
      deallocate(monomer_two_index)
      deallocate(monomer_three_index)
      if(allocated(monomer_triplets)) deallocate(monomer_triplets)
      if(allocated(triplets_diffs)) deallocate(triplets_diffs)
      deallocate(associated_to_monomer)

      if( present(n_index) ) n_index = size(this%signature_one) + &
         size(this%signature_two) + size(this%signature_three)

   endsubroutine general_trimer_sizes



   subroutine water_trimer_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
 
      type(water_trimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:), allocatable ::  triplets_diffs_map
      integer, dimension(:,:), allocatable :: monomer_one_index, monomer_triplets
      real(dp), dimension(:,:), allocatable :: triplets_diffs
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_trimer_sizes: descriptor object not initialised", error)
      endif

      allocate(associated_to_monomer(at%N))
      associated_to_monomer=.false.

      call find_general_monomer(at,monomer_one_index,&
           this%signature,associated_to_monomer,&
           this%monomer_cutoff,.false.,error) !this%atom_ordercheck

      call find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,&
          monomer_one_index,monomer_one_index,monomer_one_index,&
        !   this%one_two_identical,this%one_three_identical,this%two_three_identical,&
          .true.,.true.,.true., this%cutoff,error)
      n_descriptors = size(triplets_diffs_map)
      n_cross=n_descriptors * 9

      deallocate(monomer_one_index)
      if(allocated(monomer_triplets)) deallocate(monomer_triplets)
      if(allocated(triplets_diffs)) deallocate(triplets_diffs)
      deallocate(associated_to_monomer)

      if( present(n_index) ) n_index = 9

   endsubroutine water_trimer_sizes

      subroutine molecule_lo_d_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(molecule_lo_d), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer, dimension(:,:), allocatable :: monomer_index
      integer :: i
      logical, dimension(:), allocatable :: associated_to_monomer

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("molecule_lo_d_sizes: descriptor object not initialised", error)
      endif

      n_descriptors=1
      n_cross = this%n_atoms

!      allocate(associated_to_monomer(at%N))
!      associated_to_monomer=.false.
!
!      call find_general_monomer(at,monomer_index,this%signature,associated_to_monomer,this%cutoff,this%atom_ordercheck,error)
!      n_descriptors = size(monomer_index,2)
!      n_cross=size(monomer_index)

      if( present(n_index) ) n_index = this%n_atoms

   endsubroutine molecule_lo_d_sizes


   subroutine generate_AN_permutations(this,list,error)
      integer, dimension(:,:), intent(out) :: this
      integer, dimension(:), intent(in), optional :: list
      integer, optional, intent(out) :: error

      integer, dimension(:), allocatable :: my_list, my_list_uniq
      integer :: i, n, m, p, np, min_tail, min_tail_i, tmp_i

      INIT_ERROR(error)

      if(present(list)) then
         n = size(list)
         allocate(my_list(n))
         my_list = list
      else
         n = size(this,1)
         allocate(my_list(n))
         my_list = (/(i, i = 1, n)/)
      endif

      call uniq(my_list, my_list_uniq)

      np = factorial(size(my_list_uniq))

      call check_size('this', this, (/n,np/), 'generate_permutations',error)

      call sort_array(my_list)

      this(:,1) = my_list

      do p = 2, np
         ! Find longest tail that is ordered in decreasing order.
         do m = n - 1, 1, -1
            if(my_list(m) < my_list(m+1)) exit
         enddo

         min_tail = my_list(m+1)
         min_tail_i = m+1
         ! Find the smallest number bigger than my_list(m) in the tail
         do i = m + 1, n
            if(min_tail > my_list(i) .and. my_list(m) < my_list(i)) then
               min_tail = my_list(i)
               min_tail_i = i
            endif
         enddo

         ! swap
         tmp_i = my_list(m)
         my_list(m) = my_list(min_tail_i)
         my_list(min_tail_i) = tmp_i

         
         ! reverse tail
         my_list(m+1:n) = my_list(n:m+1:-1)

         this(:,p) = my_list
      enddo

   endsubroutine generate_AN_permutations



   subroutine bond_list_next_layer(shallow_list,deep_list,error)
      ! shallow list contains all pairs of bonded atoms, whereas deep_list will include pairs of atoms
      ! separated by two bonds, three bonds, etc., depending on how many times this function has been called
      type(Table), intent(in):: shallow_list
      type(Table), intent(inout):: deep_list
      type(Table) :: connected, pairs_containing_atom_i, pairs_containing_atom_j, deep_list_input, bonded_to_atom_i, bonded_to_atom_j
      integer :: i, atom_i, atom_j, j, atom_k, k, N_input,N_deep,N_shallow
      integer, intent(inout), optional :: error
      logical, dimension(:), allocatable :: mask_deep, mask_shallow
      logical :: i_k_present, k_i_present, j_k_present, k_j_present

      allocate(mask_deep(deep_list%N))
      allocate(mask_shallow(shallow_list%N))
      mask_deep = .False.
      mask_shallow = .False.

      ! make a copy of the deep list input
      deep_list_input = deep_list
      N_input=deep_list_input%N
      N_shallow=shallow_list%N

      ! loop over pairs in deep_list_input
      do i=1,N_input


        atom_i = deep_list_input%int(1,i)
        atom_j = deep_list_input%int(2,i)

        ! select the 1st neighbours of atom_i and atom_j
        mask_shallow = shallow_list%int(1,:N_shallow) .eq. atom_i .or. shallow_list%int(2,:N_shallow) .eq. atom_i
        call select(bonded_to_atom_i, shallow_list, row_mask=mask_shallow)

        mask_shallow = shallow_list%int(1,:N_shallow) .eq. atom_j .or. shallow_list%int(2,:N_shallow) .eq. atom_j
        call select(bonded_to_atom_j, shallow_list, row_mask=mask_shallow)


        do k=1,bonded_to_atom_j%N ! and append distances between atom i and all atoms bonded to atom j
          N_deep = deep_list%N
          atom_k = bonded_to_atom_j%int(1,k)
          if (atom_k .eq. atom_j) then
            atom_k = bonded_to_atom_j%int(2,k) 
          end if
          if (atom_k .eq. atom_i) cycle

          i_k_present = any(deep_list%int(1,:N_deep) .eq. atom_i .and. (deep_list%int(2,:N_deep) .eq. atom_k))
          k_i_present = any(deep_list%int(1,:N_deep) .eq. atom_k .and. (deep_list%int(2,:N_deep) .eq. atom_i))

          if (.not. k_i_present .and. .not. i_k_present) then
            call append(deep_list,(/atom_i,atom_k/))
          end if
        end do

        do k=1,bonded_to_atom_i%N ! and append distances between atom i and all atoms bonded to atom j
          N_deep = deep_list%N
          atom_k = bonded_to_atom_i%int(1,k)
          if (atom_k .eq. atom_i) then
            atom_k = bonded_to_atom_i%int(2,k) 
          end if
          if (atom_k .eq. atom_j) cycle

          j_k_present = any(deep_list%int(1,:N_deep) .eq. atom_j .and. (deep_list%int(2,:N_deep) .eq. atom_k))
          k_j_present = any(deep_list%int(1,:N_deep) .eq. atom_k .and. (deep_list%int(2,:N_deep) .eq. atom_j))

          if (.not. k_j_present .and. .not. j_k_present) then
            call append(deep_list,(/atom_j,atom_k/))
          end if
        end do

      end do


   end subroutine bond_list_next_layer



   subroutine transfer_initialise(this, args_str, error)
      type(transfer_parameters_type), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)
      call initialise(params)
      call param_register(params, 'do_transfer', 'false', this%do_transfer, help_string="Enable transfer function")
      call param_register(params, 'transfer_factor', '5.0', this%factor, help_string="Transfer function: stretch factor")
      call param_register(params, 'transfer_width', '1.0', this%width, help_string="Transfer function: transition width")
      call param_register(params, 'transfer_r0', '3.0', this%r0, help_string="Transfer function: transition distance")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true., task='transfer_initialise args_str')) then
         RAISE_ERROR("transfer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (this%do_transfer) then
         call print("Using transfer function with factor="//this%factor//", r0="//this%r0//", width="//this%width)
      endif
   end subroutine transfer_initialise

   function transferfunction_grad(x, params) result(td)
     real(dp), intent(in) :: x
     type(transfer_parameters_type), intent(in) :: params
     real(dp) :: td
     ! for x << r0 - width: params%factor
     ! for x >> r0 + width: 1
     td = (params%factor - 1.0_dp) * 0.5_dp*(tanh((params%r0 - x)/params%width) + 1.0_dp) + 1.0_dp
   end function transferfunction_grad

   function transferfunction(x, params) result(t)
     real(dp), intent(in) :: x
     type(transfer_parameters_type), intent(in) :: params
     real(dp) :: t
     ! for x >> r0 - width: identity (slope=1)
     ! for x << r0 + width: linear slope = factor
     t = (1.0_dp - params%factor) * 0.5_dp*(params%width*(log(2.0_dp*cosh((x - params%r0)/params%width))) + x + params%r0) + params%factor * x

   end function transferfunction

!!!!!!!!
!!! Work of Alan Nichol and S. T. John ENDS here
!!!!!!!!


!!!!!!!!
!!! Joint work of A. P. Bartok and Miguel Caro BEGINS here
!!!!!!!!
   
  subroutine descriptor_soap_express_get_eimphi_conjg(this, phi, r_ij_sigma_angular_scaled, eimphi)
     type(soap_express), intent(in) :: this
     real(dp), intent(in) :: phi, r_ij_sigma_angular_scaled
     complex(dp), dimension(:), intent(out) :: eimphi
     
    real(dp) ::  cosm2, sinm2, cosm1, sinm1, cos0, sin0, cosphi2
    integer :: l, m, k

    real(dp), dimension(0:this%l_max) :: prefactor_l
    complex(dp), dimension(0:this%l_max) :: prefactor_m

    call descriptor_soap_express_get_ilexp(this,r_ij_sigma_angular_scaled,prefactor_l)
!   Complex exponential using Euler's formula and Chebyshev recursion
    cosm2 = cos(phi)
    cosphi2 = 2.0_dp * cosm2
    sinm2 = -sin(phi)
    cosm1 = 1.0_dp
    sinm1 = 0.0_dp
    prefactor_m(0) = CPLX_ONE
    do l = 1, this%l_max
      cos0 = cosphi2 * cosm1 - cosm2
      sin0 = cosphi2 * sinm1 - sinm2
      cosm2 = cosm1
      sinm2 = sinm1
      cosm1 = cos0
      sinm1 = sin0
      prefactor_m(l) = cmplx(cos0,-sin0,kind=dp)
    end do
    k = 1
    do l = 0, this%l_max
      eimphi(k:k+l) = prefactor_l(l) * prefactor_m(0:l)
      k = k + l + 1
    end do
  end subroutine

  subroutine descriptor_soap_express_get_ilexp(this, x, prefactor_l)
     type(soap_express), intent(in) :: this
     real(dp), intent(in) :: x
     real(dp), dimension(0:this%l_max), intent(out) :: prefactor_l

     real(dp) :: x2, x4, fl, flm1, flm2
     integer :: l

     real(dp), parameter :: small_number = 1.0e-7_dp, small_number2 = 1.0e-4_dp

     x2 = x**2
     x4 = x**4

!    Full calculation. This is numerically unstable for small x, that's why
!    we have cases below
     flm2 = abs( (1.0_dp - exp(-2.0_dp * x2)) / 2.0_dp / x2 )
     flm1 = abs( (x2 - 1.0_dp + exp(-2.0_dp * x2)*(x2+1.0_dp)) / 2.0_dp / x4 )

     do l = 0, this%l_max
        if( l == 0 ) then
           if( x < small_number )then
              prefactor_l(0) = 1.0_dp - x2
           else
              prefactor_l(0) = flm2
           endif
        elseif( l == 1 )then
           if( x2 < small_number2 ) then
              prefactor_l(1) = (x2 - x4)/ this%semifactorial_table(1)
           else
              prefactor_l(1) = flm1
           endif
        else
           if( x2**l / this%semifactorial_table(l) * l < small_number )then
              fl = x2**l / this%semifactorial_table(l)
           else
              fl = abs( flm2 - (2.0_dp*l - 1.0_dp)/x2 * flm1 )
           end if
           flm2 = flm1
           flm1 = fl
           prefactor_l(l) = fl
        endif
     end do

  endsubroutine

  subroutine descriptor_soap_express_get_plm_array(this,x,plm_array)
!   Returns an array with the Associated Legendre polynomials Plm(x), where Pll(x)
!   is the highest order in the series, together with, all the lower order ones.
!   l is > 0, m can take values from 0 to l and x is within the [-1, 1] domain
!
!   plm_array is a 1D array with modified index lm -> k, where
!   k = 1 + l*(l+1)/2 + m
     type(soap_express), intent(in) :: this
     real(dp), intent(in) :: x
     real(dp), dimension(this%angular_array_size), intent(out) :: plm_array

    integer :: l, m, i, lmax, k

    if( abs(x) > 1.0_dp ) then
        call system_abort("Bad argument for associated Legendre polynomial")
    end if

!   We need these 6 polynomials to initialize the recursion series
!   P_00(x), k = 1
    plm_array(1) = 1.0_dp
    if( this%l_max > 0 )then
!     P_10(x), k = 2
      plm_array(2) = x
!     P_11(x), k = 3
      plm_array(3) = -sqrt(1.0_dp-x**2)
      if( this%l_max > 1 )then
!       P_20(x), k = 4
        plm_array(4) = 1.5_dp*x**2 - 0.5_dp
!       P_21(x), k = 5
        plm_array(5) = -3.0_dp * x * sqrt(1.0_dp-x**2)
!       P_22(x), k = 6
        plm_array(6) = 3.0_dp - 3.0_dp*x**2
      else
        return
      end if
    else
      return
    end if

    if( this%l_max == 2 )then
      return
    else
      do l = 3, this%l_max
!       First we need to obtain Pl0, Pl1, up to m=l-2 with the recursion formula on l:
!       Plm(x) = ( (2l-1)*x*Pl-1m(x) - (l-1+m)*Pl-2m(x) ) / (l-m)
        do m = 0, l-2
          k = 1 + l*(l+1)/2 + m
          plm_array(k) = (real(2*l-1,kind=dp)*x*plm_array(k-l) - real(l-1+m,kind=dp)*plm_array(k-2*l+1)) &
                         / real(l-m,kind=dp)
        end do
!       Now we get Pll-1 and Pll with the recursion formulas
!       Pll-1(x) = x*(2l-1)*Pl-1l-1(x)
!       Pll(x) = -(2l-1)*sqrt(1-x^2)*Pl-1l-1(x)
        k = k + 1
        plm_array(k) = x * real(2*l-1,kind=dp) * plm_array(k-l)
        k = k + 1
        plm_array(k) = - real(2*l-1,kind=dp) * sqrt(1.0_dp - x**2) * plm_array(k-l-1)
      end do
    end if





  end subroutine

  subroutine descriptor_soap_express_get_radial_array(this,r_ij,coeff)
     type(soap_express), intent(in) :: this
     real(dp), intent(in) :: r_ij
     real(dp), dimension(this%n_max), intent(out) :: coeff

     real(dp) :: r_ij_radial_scaled, cutoff_hard_scaled, cutoff_soft_scaled, atom_sigma_radial_scaled, &
        atom_sigma_radial_scaled_2, I_n, N_n, N_np1, I_np1, C1, C2, atom_sigma_radial_filtered, &
        r_ij_radial_scaled_filtered, atom_sigma_radial_filtered_2, N_np2, I_np2, &
        cutoff_decay_scaled, atom_sigma_radial0_scaled, n_gauss, pref_f

     integer :: a

     cutoff_decay_scaled = this%cutoff_transition_width / this%cutoff
     atom_sigma_radial0_scaled = this%atom_sigma_radial/this%cutoff
     n_gauss = sqrt(2.0_dp / atom_sigma_radial0_scaled) / PI**0.25_dp

     r_ij_radial_scaled = r_ij / this%cutoff
     cutoff_hard_scaled = 1.0_dp
     cutoff_soft_scaled = (this%cutoff - this%cutoff_transition_width) / this%cutoff

     atom_sigma_radial_scaled = ( this%atom_sigma_radial + this%atom_sigma_scaling_radial*r_ij ) / this%cutoff
     atom_sigma_radial_scaled_2 = atom_sigma_radial_scaled**2

     coeff = 0.0_dp

!    We have the recursion series starting at n = 0, which means alpha = -2
!    However, we only need to save the expansion coefficients for alpha >= 1
!    This is I_-1
     I_n = 0.0_dp
     N_n = 1.0_dp
!    This is I_0
     N_np1 = descriptor_soap_express_N_a(cutoff_hard_scaled, -2)
     I_np1 = sqrt(PI/2.0_dp) * atom_sigma_radial_scaled * ( &
        erf( (cutoff_soft_scaled-r_ij_radial_scaled)/SQRT_TWO/atom_sigma_radial_scaled ) - &
        erf( (-r_ij_radial_scaled)/SQRT_TWO/atom_sigma_radial_scaled ) ) / N_np1
!    Speed up the computation of these coefficients
     if( this%cutoff_transition_width .feq. 0.0_dp )then
       C1 = 0.0_dp
     else
       C1 = atom_sigma_radial_scaled_2 / cutoff_decay_scaled * &
          exp(-0.5_dp * (cutoff_soft_scaled - r_ij_radial_scaled)**2 / atom_sigma_radial_scaled_2)
     endif
     C2 = atom_sigma_radial_scaled_2 / cutoff_hard_scaled * exp(-0.5_dp * r_ij_radial_scaled**2 / atom_sigma_radial_scaled_2)
!        This is different wrt the regular polynomial basis, we only go up to alpha_max-1
     do a = -1, this%n_max-1
        C1 = C1 * cutoff_decay_scaled
        C2 = C2 * cutoff_hard_scaled
        N_np2 = descriptor_soap_express_N_a(cutoff_hard_scaled, a)
!       This is I_alpha
        I_np2 = atom_sigma_radial_scaled_2 * (a+1) * N_n/ N_np2 * I_n &
                - N_np1 * (r_ij_radial_scaled - cutoff_hard_scaled) / N_np2 * I_np1 &
                + C1 / N_np2 &
                - C2 / N_np2
        if(a > 0) coeff(a) = I_np2
        N_n = N_np1
        N_np1 = N_np2
        I_n = I_np1
        I_np1 = I_np2
     enddo
     
!    If the atom is less than 4 standard deviations away from the soft cutoff, we add
!    also this correction to the integrals. This corresponds to a Gaussian filter. We
!    integrate between rcut_soft and rcut_hard in this case
     if( (cutoff_soft_scaled - r_ij_radial_scaled) < 4.0_dp*atom_sigma_radial_scaled )then
        atom_sigma_radial_filtered = atom_sigma_radial_scaled * cutoff_decay_scaled / &
           this%cutoff_decay_rate / sqrt(atom_sigma_radial_scaled_2 + &
           cutoff_decay_scaled**2/this%cutoff_decay_rate**2)

        r_ij_radial_scaled_filtered = (atom_sigma_radial_scaled_2 * cutoff_soft_scaled + &
           cutoff_decay_scaled**2/this%cutoff_decay_rate**2*r_ij_radial_scaled) / &
           (atom_sigma_radial_scaled_2 + cutoff_decay_scaled**2/this%cutoff_decay_rate**2)

!       We leave this here because in the future atom_sigma will be rj dependent
        atom_sigma_radial_filtered_2 = atom_sigma_radial_filtered**2
!       The products of two Gaussians is a Gaussian, but we need to add a prefactor
        pref_f = exp( -0.5_dp * (cutoff_soft_scaled - r_ij_radial_scaled)**2 / ( atom_sigma_radial_scaled_2 + &
           cutoff_decay_scaled**2 / this%cutoff_decay_rate**2 ) )
!       We have the recursion series starting at n = 0, which means alpha = -2
!       However, we only need to save the expansion coefficients for alpha >= 1
!       This is I_-1
        I_n = 0.0_dp
        N_n = 1.0_dp
!       This is I_0
        N_np1 = descriptor_soap_express_N_a(cutoff_hard_scaled, -2)
        I_np1 = sqrt(PI/2.0_dp) * atom_sigma_radial_filtered * &
           ( erf( (cutoff_hard_scaled-r_ij_radial_scaled_filtered)/SQRT_TWO/atom_sigma_radial_filtered ) - &
           erf( (cutoff_soft_scaled-r_ij_radial_scaled_filtered)/SQRT_TWO/atom_sigma_radial_filtered ) ) / N_np1
!       Speed up the computation of these coefficients
        C2 = atom_sigma_radial_filtered_2 / cutoff_decay_scaled * &
           exp(-0.5_dp * (cutoff_soft_scaled - r_ij_radial_scaled_filtered)**2 / atom_sigma_radial_filtered_2)
!       This is different wrt the regular polynomial basis, we only go up to this%n_max-1
        do a = -1, this%n_max-1
           C2 = C2 * cutoff_decay_scaled
           N_np2 = descriptor_soap_express_N_a(cutoff_hard_scaled, a)
!          This is I_alpha
           I_np2 = atom_sigma_radial_filtered_2 * (a+1) * N_n/ N_np2 * I_n &
                   - N_np1 * (r_ij_radial_scaled_filtered - cutoff_hard_scaled) / N_np2 * I_np1 &
                   - C2 / N_np2
           if(a > 0)then
             coeff(a) = coeff(a) + pref_f * I_np2
           end if
           N_n = N_np1
           N_np1 = N_np2
           I_n = I_np1
           I_np1 = I_np2
        end do
     endif

!    Now we obtain the overlap integral between the atomic density and the Gaussian centered at the origin
!    We assume that at the soft cutoff the Gaussian basis function is approx. zero, and we only
!    compute the overlap coefficient if the atom is close to the origin
!    Note that the sigma for the Gaussian basis function and the atom's sigma are not the same if there is
!    sigma scaling
     if( r_ij_radial_scaled < 4.0_dp*(atom_sigma_radial0_scaled + atom_sigma_radial_scaled) ) then
       coeff(this%n_max) = exp(-0.5_dp * r_ij_radial_scaled**2 / &
          (atom_sigma_radial_scaled_2 + atom_sigma_radial0_scaled**2) ) * sqrt(PI/2.0_dp) * &
          atom_sigma_radial_scaled*atom_sigma_radial0_scaled / &
          sqrt( atom_sigma_radial_scaled_2 + atom_sigma_radial0_scaled**2 ) * &
          ( 1.0_dp + erf(atom_sigma_radial0_scaled/atom_sigma_radial_scaled*r_ij_radial_scaled/SQRT_TWO/&
          sqrt( atom_sigma_radial_scaled_2 + atom_sigma_radial0_scaled**2 )) ) * n_gauss
     endif
!    Transform from g_alpha to g_n (the orthonormal basis)
     coeff(1:this%n_max) = matmul( this%basis_transformation_coefficients, coeff(1:this%n_max) )
!    Add filter if necessary (only for inverse scaling mode, which applies a global 1/(1+a*rj) filter
!    plus this one inside the buffer zone; the polynomial scaling mode applies a global filter only
!    which appropriately goes to zero at rcut_hard)
     if( r_ij_radial_scaled > cutoff_soft_scaled .and. .false.)then
        coeff(1:this%n_max) = coeff(1:this%n_max) * exp(-0.5_dp*(r_ij_radial_scaled - &
        cutoff_soft_scaled)**2 / cutoff_decay_scaled**2 * this%cutoff_decay_rate**2)
     endif

  endsubroutine descriptor_soap_express_get_radial_array

  function descriptor_soap_express_N_a(cutoff, a)

    integer, intent(in) :: a
    real(dp), intent(in) :: cutoff
    real(dp) :: descriptor_soap_express_N_a

    integer :: b

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    descriptor_soap_express_N_a = sqrt( cutoff / real(b,kind=dp) )
!   *********************************************

  endfunction descriptor_soap_express_N_a

  subroutine soap_express_sizes(this,at,n_descriptors,n_cross,mask,n_index,error)
      type(soap_express), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      logical, dimension(:), intent(in), optional :: mask
      integer, intent(out), optional :: n_index
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_express_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(present(mask)) then
            if(.not. mask(i)) cycle
         endif
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

      if( present(n_index) ) n_index = 1

   endsubroutine soap_express_sizes

   
   function soap_express_cutoff(this,error) 
      type(soap_express), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: soap_express_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_express_cutoff: descriptor object not initialised", error)
      endif

      soap_express_cutoff = this%cutoff

   endfunction soap_express_cutoff

   function soap_express_dimensions(this,error) result(i)
      type(soap_express), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_express_dimensions: descriptor object not initialised", error)
      endif

      i = ( this%l_max+1 ) * ( this%n_max*(this%n_max+1) ) / 2 + 1

   endfunction soap_express_dimensions

   subroutine soap_express_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type(soap_express), intent(in) :: this
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
      integer :: a, b, d, i, j, k, l, m, n, i_n, l_n_neighbours, &
         i_desc, n_descriptors, n_cross, n_index, i_pow
      integer, dimension(3) :: shift
      real(dp) :: r_ij, phi, atom_sigma_angular, atom_sigma_radial, r_ij_sigma_angular_scaled, &
         amplitude, multiplicity, atom_sigma_radial_scaled
      real(dp), dimension(3) :: u_ij
      real(dp), dimension(this%angular_array_size) :: plm_array
      complex(dp), dimension(this%angular_array_size) :: eimphi
      complex(dp), dimension(:,:), allocatable :: exp_coeff_angular
      real(dp), dimension(:,:), allocatable :: exp_coeff_radial
      real(dp), dimension(:), allocatable :: descriptor_i
      complex(dp), dimension(:,:), allocatable :: fourier_so3_array

      INIT_ERROR(error)

      call system_timer('soap_express_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_express_calc: descriptor object not initialised", error)
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
            RAISE_ERROR("soap_express_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("soap_express_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      call finalise(descriptor_out)

      d = soap_express_dimensions(this,error)

      allocate(fourier_so3_array(this%angular_array_size,this%n_max))
      allocate(descriptor_i(d))

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

         allocate(exp_coeff_angular(this%angular_array_size, n_neighbours(at,i, max_dist = this%cutoff)))
         allocate(exp_coeff_radial(this%n_max, n_neighbours(at,i, max_dist = this%cutoff)))

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            phi = atan2( u_ij(2), u_ij(1) )

            atom_sigma_angular = this%atom_sigma_angular + this%atom_sigma_scaling_angular * r_ij

            r_ij_sigma_angular_scaled = r_ij / atom_sigma_angular
            
            call descriptor_soap_express_get_eimphi_conjg( this, phi, r_ij_sigma_angular_scaled, eimphi )
            call descriptor_soap_express_get_plm_array(this, u_ij(3), plm_array)

            exp_coeff_angular(:,i_n) = this%Y_lm_prefactor * plm_array * eimphi / atom_sigma_angular**2 * this%cutoff**2
            call descriptor_soap_express_get_radial_array( this, r_ij, exp_coeff_radial(:,i_n) )

            if( .false. )then
!              This is the "inverse power" scaling originally proposed
               amplitude = 1.0_dp / (1.0_dp + this%amplitude_scaling*r_ij)
            else
!              This is the polynomial scaling that made it to the paper
               if( this%amplitude_scaling == 0 )then
                  amplitude = 1.0_dp
               else if( this%amplitude_scaling == 1 )then
                  amplitude = 1.0_dp * ( 1.d0 + 2.d0*(r_ij/this%cutoff)**3 - 3.d0*(r_ij/this%cutoff)**2 )
               else
                  amplitude = 1.0_dp * ( 1.d0 + 2.d0*(r_ij/this%cutoff)**3 - 3.d0*(r_ij/this%cutoff)**2 )**this%amplitude_scaling
               end if
            endif
            atom_sigma_radial_scaled = ( this%atom_sigma_radial + this%atom_sigma_scaling_radial*r_ij ) / this%cutoff
            exp_coeff_radial(1:this%n_max, i_n) = exp_coeff_radial(1:this%n_max, i_n) * amplitude / atom_sigma_radial_scaled
         enddo

!        This results from the change of variable in the
!        overlap integrals. We only need this if we want to
!        know the actual value of the expansion coefficients.
!        Since this is a global factor, once we have normalized
!        the SOAP vectors it does not have an effect anymore.
         exp_coeff_radial = exp_coeff_radial * sqrt(this%cutoff)

         fourier_so3_array = CPLX_ZERO

         i_n = 0
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance = r_ij)

            if( r_ij >= this%cutoff ) cycle
            i_n = i_n + 1

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = 0, l

                     k = 1 + (l*(l+1)) / 2 + m
!                    It is messy with the prefactor in spherical harmonics but we need to be sure because of the central atom below
                     fourier_so3_array(k, a) = fourier_so3_array(k, a) + 4.0_dp*PI * exp_coeff_radial(a, i_n) * &
                        exp_coeff_angular(k,i_n)

                  enddo
               enddo
            enddo
         enddo

         if( this%central_weight > 0.0_dp )then
             fourier_so3_array(1, 1:this%n_max) = fourier_so3_array(1, 1:this%n_max) + this%central_weight * sqrt(4.0_dp*PI) * &
             PI**0.25_dp * sqrt(this%atom_sigma_radial / 2.0_dp) * &
             this%cutoff**3 / this%atom_sigma_angular**2 / this%atom_sigma_radial * &
             matmul(this%basis_transformation_coefficients, this%overlap(1:this%n_max, this%n_max))
         end if

         if(my_do_descriptor) then
            descriptor_i = 0.0_dp
            i_pow = 0
            do a = 1, this%n_max
               do b = a, this%n_max
                  do l = 0, this%l_max
                     i_pow = i_pow+1
                     do m = 0, l
                        k = 1 + (l*(l+1))/2 + m
                        multiplicity = 1.0_dp
                        if( a /= b ) multiplicity = multiplicity * SQRT_TWO
                        if( m > 0 ) multiplicity = multiplicity * 2.0_dp
                        descriptor_i(i_pow) = descriptor_i(i_pow) + multiplicity * real(fourier_so3_array(k, a) * &
                           conjg(fourier_so3_array(k, b)))
                     enddo
                  enddo
               enddo
            enddo
            descriptor_out%x(i_desc)%data = descriptor_i / sqrt(dot_product(descriptor_i, descriptor_i))
         endif

         if(my_do_grad_descriptor) then
            RAISE_ERROR("soap_express_calc: derivatives not implemented yet.", error)
         endif

         deallocate(exp_coeff_angular)
         deallocate(exp_coeff_radial)

      enddo

      deallocate(fourier_so3_array)
      deallocate(descriptor_i)

      call system_timer('soap_express_calc')

   endsubroutine soap_express_calc

   subroutine soap_express_initialise(this,args_str,error)
      type(soap_express), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      logical :: has_atom_sigma_angular

      integer :: l, k, i, j, m, n
      real(dp) :: fact, fact1, fact2, ppi, atom_sigma_radial_normalised, cutoff_hard,&
         s2, I_n, N_n, N_np1, N_np2, I_np1, I_np2, C2

      type(LA_Matrix) :: LA_overlap
      real(dp), dimension(:), allocatable :: s
      real(dp), dimension(:,:), allocatable :: sqrt_overlap, u, v
      real(dp), parameter :: sqrt_two = sqrt(2.0_dp)

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', PARAM_MANDATORY, this%cutoff, &
           help_string="Cutoff for soap-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, &
           help_string="Cutoff transition width for soap-type descriptors")
      call param_register(params, 'cutoff_decay_rate', '4.', this%cutoff_decay_rate, &
           help_string="Cutoff decay rate for soap-type descriptors")
      call param_register(params, 'atom_sigma_radial', PARAM_MANDATORY, this%atom_sigma_radial, &
           help_string="Width of atomic Gaussians for soap-type descriptors in the radial direction")
      call param_register(params, 'atom_sigma_angular', '0.0', this%atom_sigma_angular, &
           help_string="Width of atomic Gaussians for soap-type descriptors in the angular direction",has_value_target = has_atom_sigma_angular)
      call param_register(params, 'central_weight', '1.0', this%central_weight, &
           help_string="Weight of central atom in environment")
      call param_register(params, 'covariance_sigma0', '0.0', this%covariance_sigma0, &
           help_string="sigma_0 parameter in polynomial covariance function")
      call param_register(params, 'atom_sigma_scaling_radial', '0.0', this%atom_sigma_scaling_radial, &
           help_string="Scaling rate of radial sigma: scaled as a function of neighbour distance")
      call param_register(params, 'atom_sigma_scaling_angular', '0.0', this%atom_sigma_scaling_angular, &
           help_string="Scaling rate of angular sigma: scaled as a function of neighbour distance")
      call param_register(params, 'amplitude_scaling', '0.0', this%amplitude_scaling, &
           help_string="Scaling rate of amplitude: scaled as an inverse function of neighbour distance")
      call param_register(params, 'l_max', PARAM_MANDATORY, this%l_max, &
           help_string="L_max (spherical harmonics basis band limit) for soap-type descriptors")
      call param_register(params, 'n_max', PARAM_MANDATORY, this%n_max, &
           help_string="N_max (number of radial basis functions) for soap-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_express_initialise args_str')) then
         RAISE_ERROR("soap_express_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)
      if( this%n_max > 12 ) then
         RAISE_ERROR("soap_express_initialise: radial expansion unstable for n_max > 12",error)
      endif

      if( this%l_max > 12 ) then
         RAISE_ERROR("soap_express_initialise: radial expansion unstable for l_max > 12",error)
      endif

      if( .not. has_atom_sigma_angular ) this%atom_sigma_angular = this%atom_sigma_radial

      allocate(this%semifactorial_table(0:this%l_max))
      If( this%l_max >= 0 )then
         fact = 1.0_dp
         do l = 0, this%l_max
            fact = fact * (2.0_dp*l + 1.0_dp)
            this%semifactorial_table(l) = fact
         enddo
      else
         RAISE_ERROR("soap_express_initialise: l_max must be greater or equal to zero", error)
      endif
      this%angular_array_size = 1 + ( this%l_max*(this%l_max+1) ) / 2 + this%l_max

      allocate(this%Y_lm_prefactor(this%angular_array_size))

      k = 1
      do l = 0, this%l_max
         ppi = sqrt((2*l+1)/4.0_dp/PI)
         fact2 = 1.0_dp
         do i = 1, l
           fact2 = fact2 * i
         enddo
         fact1 = fact2
         do m = 0, l
           if( m > 0 )then
   !         This is the factorial of (l-m)
             fact1 = fact1 / (l+1-m)
   !         This is the factorial of (l+m)
             fact2 = fact2 * (l+m)
           end if
           this%Y_lm_prefactor(k) = ppi * sqrt(fact1/fact2)
           k = k + 1
         end do
       end do

      allocate( this%overlap(this%n_max, this%n_max) )
      allocate( this%basis_transformation_coefficients(this%n_max, this%n_max) )
      allocate( sqrt_overlap(this%n_max, this%n_max) )

  !   These are the overlap integrals for the polynomial functions
      do i = 1, this%n_max-1
         this%overlap(i,i) = 1.0_dp
         do j = i+1, this%n_max-1
            this%overlap(i,j) = sqrt( (5.0_dp+2.0_dp*i) * (5.0_dp+2.0_dp*j) ) / real(5+i+j,kind=dp)
            this%overlap(j,i) = this%overlap(i,j)
         end do
      end do

!   These are the overlap integrals between the Gaussian and the polynomials
!   See derivation of radial expansion coefficients to understand this code
!   **** New basis ****
      atom_sigma_radial_normalised = this%atom_sigma_radial/this%cutoff
      cutoff_hard = 1.0_dp
!   *******************
      s2 = atom_sigma_radial_normalised**2
      I_n = 0.0_dp
      N_n = 1.0_dp
      N_np1 = descriptor_soap_express_N_a(cutoff_hard, -2)
      I_np1 = sqrt(PI/2.0_dp) * atom_sigma_radial_normalised * &
         erf( cutoff_hard/SQRT_TWO/atom_sigma_radial_normalised ) / N_np1
      C2 = s2 / cutoff_hard
      do n = -1, this%n_max-1
         C2 = C2 * cutoff_hard
         N_np2 = descriptor_soap_express_N_a(cutoff_hard, n)
         I_np2 = s2 * (n+1) * N_n/ N_np2 * I_n &
                 + N_np1 * cutoff_hard / N_np2 * I_np1 &
                 - C2 / N_np2
         if(n > 0)then
   !       Include the normalization factor of the Gaussian
           this%overlap(this%n_max, n) = I_np2 * SQRT_TWO / &
              sqrt(atom_sigma_radial_normalised) / PI**0.25_dp
           this%overlap(n, this%n_max) = this%overlap(this%n_max, n)
         end if
         N_n = N_np1
         N_np1 = N_np2
         I_n = I_np1
         I_np1 = I_np2
       end do
       this%overlap(this%n_max,this%n_max) = 1.0_dp

       call initialise(LA_overlap,this%overlap)
       call LA_Matrix_SVD_Allocate(LA_overlap,s,u,v,error)
       call LA_Matrix_SVD(LA_overlap,s,u,v,error)

       sqrt_overlap = 0.0_dp
       do i = 1, this%n_max
          sqrt_overlap(i,i) = sqrt(s(i))
       end do
       sqrt_overlap = matmul(U,sqrt_overlap)
       sqrt_overlap = matmul(sqrt_overlap,transpose(V))

       LA_overlap = sqrt_overlap
       call LA_Matrix_Factorise(LA_overlap,error=error)
       call LA_Matrix_Inverse(LA_overlap,this%basis_transformation_coefficients)
       call Finalise(LA_overlap)
       if(allocated(sqrt_overlap)) deallocate(sqrt_overlap)
       if(allocated(s)) deallocate(s)
       if(allocated(u)) deallocate(u)
       if(allocated(v)) deallocate(v)
      
       this%initialised = .true.

   endsubroutine soap_express_initialise

   subroutine soap_express_finalise(this,error)
      type(soap_express), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%cutoff_decay_rate = 0.0_dp
      this%l_max = 0
      this%n_max = 0

      this%atom_sigma_radial = 0.0_dp
      this%atom_sigma_angular = 0.0_dp
      this%central_weight = 0.0_dp
      this%covariance_sigma0 = 0.0_dp
      this%atom_sigma_scaling_radial = 0.0_dp
      this%atom_sigma_scaling_angular = 0.0_dp
      this%amplitude_scaling = 0.0_dp

      if(allocated(this%basis_transformation_coefficients)) deallocate(this%basis_transformation_coefficients)
      if(allocated(this%overlap)) deallocate(this%overlap)
      if(allocated(this%semifactorial_table)) deallocate(this%semifactorial_table)
      if(allocated(this%Y_lm_prefactor)) deallocate(this%Y_lm_prefactor)

      this%initialised = .false.

   endsubroutine soap_express_finalise

!!!!!!!!
!!! Joint work of A. P. Bartok and Miguel Caro ENDS here
!!!!!!!!


  
