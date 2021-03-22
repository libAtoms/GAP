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

!
! Author: Alan Nichol
!

         case(DT_AN_MONOMER)
            allocate(distance_matrix(this%descriptor_AN_monomer%N,this%descriptor_AN_monomer%N), atom_permutations(this%descriptor_AN_monomer%N,np))

            if(this%descriptor_AN_monomer%do_atomic) then
               atom_permutations(1,:) = 0
               call generate_AN_permutations(atom_permutations(2:this%descriptor_AN_monomer%N,:))
               atom_permutations = atom_permutations + 1
            else
               call generate_AN_permutations(atom_permutations(:,:))
            endif

            i = 0
            distance_matrix = 0
            do n = 2, this%descriptor_AN_monomer%N
               i = i + 1
               distance_matrix(1,n) = i
               distance_matrix(n,1) = i
               do m = n+1, this%descriptor_AN_monomer%N
                  i = i + 1
                  distance_matrix(m,n) = i
                  distance_matrix(n,m) = i
               enddo
            enddo

            do ip = 1, np
               i = 0
               do n = 2, this%descriptor_AN_monomer%N
                  i = i + 1
                  permutations(i,ip) = distance_matrix(atom_permutations(1,ip),atom_permutations(n,ip))
                  do m = n+1, this%descriptor_AN_monomer%N
                     i = i + 1
                     permutations(i,ip) = distance_matrix(atom_permutations(m,ip),atom_permutations(n,ip))
                  enddo
               enddo
            enddo
            deallocate(distance_matrix,atom_permutations)


         case(DT_GENERAL_MONOMER)
            if (.not. this%descriptor_general_monomer%permutation_data%initialised) then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_monomer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_monomer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_GENERAL_DIMER)
            if (.not. this%descriptor_general_dimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_dimer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_dimer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_GENERAL_TRIMER)
            if (.not. this%descriptor_general_trimer%permutation_data%initialised)then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_general_trimer%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_general_trimer%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            permutations=my_permutation_data%dist_vec_permutations

         case(DT_MOLECULE_LO_D)
            if (.not. this%descriptor_molecule_lo_d%permutation_data%initialised) then
              RAISE_ERROR("descriptor_permutations: permutation_data not initialised "//this%descriptor_type,error)
            else if (this%descriptor_molecule_lo_d%permutation_data%perm_number /= 1) then
              RAISE_ERROR("descriptor_permutations: permutation_data%perm_number must be initialised to one"//this%descriptor_type,error)
            end if

            call permutation_data_copy(my_permutation_data, this%descriptor_molecule_lo_d%permutation_data)

            if (my_permutation_data%n_perms > 1) then
              call next(my_permutation_data, 1)
            end if

            allocate(sliced_permutations(size(this%descriptor_molecule_lo_d%included_components),my_permutation_data%n_perms))
            allocate(this_perm(size(this%descriptor_molecule_lo_d%included_components)))
            sliced_permutations =my_permutation_data%dist_vec_permutations(this%descriptor_molecule_lo_d%included_components,:)

            do j=1,my_permutation_data%n_perms
              this_perm=sliced_permutations(:,j)
              do i=1,size(this%descriptor_molecule_lo_d%included_components)
                unit_vec=maxloc(this%descriptor_molecule_lo_d%included_components, mask=this%descriptor_molecule_lo_d%included_components .eq. this_perm(i))
                if (unit_vec(1) == 0) then
                  RAISE_ERROR("descriptor_permutations: you have specified symmetries between atoms with different connectivity",error)
                end if
                permutations(i,j) =unit_vec(1)
              end do
            end do
! begin brau
            if(size(this%descriptor_molecule_lo_d%included_components) > maxval(this%descriptor_molecule_lo_d%included_components)) then
              permutations=my_permutation_data%dist_vec_permutations
            end if
! end brau