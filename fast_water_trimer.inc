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

! This file contains and adaptation of the general_trimer routine to rapidly calculate water trimers.
! The main issues of the general trimer (and dimer) routines are that large arrays (monomer_triplets, triplets_diffs, ...) are
! - reallocated at each found trimer (or dimer)
! - double-checked to avoid double inclusions
! This routines has (hopefully) equivalent results as long as 
! - the lattice is orthogonal
! - the shortest lattice vector is > 2*trimer_cutoff
! - there is only water in the simulation
!
! Work done by Jonatan Öström, jonatan.ostrom@gmail.com, @sujona
!

subroutine water_trimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error) !//////////////////////////////////////////////////////////////////////////////////
    type(water_trimer), intent(in) :: this
    type(atoms), intent(in) :: at
    type(descriptor_data), intent(out) :: descriptor_out
    logical, intent(in), optional :: do_descriptor, do_grad_descriptor
    character(len=*), intent(in), optional :: args_str
    integer, optional, intent(out) :: error
    logical :: my_do_descriptor, my_do_grad_descriptor
    integer :: n_descriptors, i, i_desc, n_mon_trip
    integer :: n_monomer1, n_monomer2, n_monomer3
    real(dp), dimension(:,:), allocatable :: triplets_diffs
    integer, dimension(:), allocatable :: triplets_diffs_map 
    integer, dimension(:,:), allocatable :: imonomers1 ,imonomers2, imonomers3, itriplets
    logical :: associated_to_monomer(at%N)
    integer, parameter :: trimer_size = 9, nd = 36
    INIT_ERROR(error)
    
   !  print*, "Entered water_trimer_calc"
    call system_timer('water_trimer_calc')

    if(.not. this%initialised) then
       RAISE_ERROR("water_trimer_calc: descriptor object not initialised", error)
    endif

    my_do_descriptor = optional_default(.false., do_descriptor)
    my_do_grad_descriptor = .true. !optional_default(.false., do_grad_descriptor) ! fix!!!


    if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

    call finalise(descriptor_out)

    associated_to_monomer=.false. 

    call find_general_monomer(at,imonomers1,&
         this%signature,associated_to_monomer,&
         this%monomer_cutoff,.false.,error) !this%atom_ordercheck=false
   ! Allocation on assignment
   block 
       integer s1,s2
       s1 = size(imonomers1,1)
       s2 = size(imonomers1,2)
       allocate(imonomers2(s1,s2))
       allocate(imonomers3(s1,s2))
       imonomers2   = imonomers1
       imonomers3 = imonomers1
   endblock
       
    n_monomer1 = size(imonomers1,2) !, imonomers2, imonomers3
    n_monomer2 = size(imonomers2,2)
    n_monomer3 = size(imonomers3,2)

   call system_timer('water_trimer_calc: find_monomer_triplets')
   
   print'(a,*(f10.5))', "at%lattice = ",at%lattice
   
   call find_monomer_triplets(at,itriplets,triplets_diffs,triplets_diffs_map,&
           imonomers1,imonomers2,imonomers3,&
           ! this%one_two_identical,this%one_three_identical,this%two_three_identical,&
           .true.,.true.,.true.,this%cutoff,error)
   
    n_descriptors = size(triplets_diffs_map)
    n_mon_trip = size(itriplets,2)
    
    if (n_descriptors /= n_mon_trip) then 
       RAISE_ERROR("water_triplets_calc: #diffs="//n_descriptors//" and #tiplets="//n_mon_trip//", make sure box is orthogonal and > 2*cutoff ",error)
    endif
    
    print*, "  wt calc: found #mono="//n_monomer1//" #trip="//n_mon_trip
    
    call system_timer('water_trimer_calc: find_monomer_triplets')
    call print("itriplets("//size(itriplets,1)//", "//size(itriplets,2)//")", PRINT_NERD)
    call print("triplets_diffs("//size(triplets_diffs,1)//", "//size(triplets_diffs,2)//")", PRINT_NERD)
    call print("triplets_diffs_map("//size(triplets_diffs_map)//")", PRINT_NERD)

    !call print("itriplets",PRINT_NERD)
    !call print(itriplets,PRINT_NERD)

    call print("ready to construct "//n_descriptors //" descriptors",PRINT_NERD)
    allocate(descriptor_out%x(n_descriptors))
    do i = 1, n_descriptors
       if(my_do_descriptor) then
          allocate(descriptor_out%x(i)%data(nd))
          allocate(descriptor_out%x(i)%ci(trimer_size))
          descriptor_out%x(i)%has_data = .false.
          descriptor_out%x(i)%data = 0.0_dp
          descriptor_out%x(i)%ci = 0
          descriptor_out%x(i)%covariance_cutoff = 1.0_dp
       endif
       if(my_do_grad_descriptor) then
          allocate(descriptor_out%x(i)%grad_data(nd,3,trimer_size))
          allocate(descriptor_out%x(i)%ii(trimer_size))
          allocate(descriptor_out%x(i)%pos(3,trimer_size))
          allocate(descriptor_out%x(i)%has_grad_data(trimer_size))
          allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,trimer_size))
          descriptor_out%x(i)%has_grad_data = .false.
          descriptor_out%x(i)%grad_data = 0.0_dp
          descriptor_out%x(i)%ii = 0
          descriptor_out%x(i)%pos = 0.0_dp
          descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
       endif
    enddo
    
   the_loop : do i_desc = 1, n_descriptors
            
        block 
            integer, parameter :: monomer_size = 3, trimer_size = 9, nd = 36
            real(dp) cut
            real(dp), dimension(3,9) :: pos, gcuts
            real(dp), dimension(nd) :: desc0, rlen
            ! real(dp), dimension(3,36) :: rvec
            ! real(dp), dimension(3,9,36) :: gdesc0
            real(dp), dimension(nd,3,9) :: gdesc0
            integer ii, jj
            integer i1, i2,i3
            real(dp) :: diffs1(3,monomer_size)
            real(dp) :: diffs2(3,monomer_size)
            real(dp) :: diffs3(3,monomer_size)
            integer :: iatoms1(monomer_size)
            integer :: iatmos2(monomer_size)
            integer :: iatoms3(monomer_size)
            integer :: iatoms123(trimer_size)
            real(dp), dimension(3) :: diff12, diff13
            real(dp), dimension(3) :: mean_pos1,mean_pos2,mean_pos3
            
            
            i1 = itriplets(1,i_desc)
            i2 = itriplets(2,i_desc)
            i3 = itriplets(3,i_desc)
            
            iatoms1 = imonomers1(:,i1) !store the indices of atoms in this monomer
            iatmos2 = imonomers2(:,i2) 
            iatoms3 = imonomers3(:,i3)
            
            iatoms123=(/iatoms1,iatmos2,iatoms3/)
            
            mean_pos1 = calc_mean_pos(at,iatoms1)
            mean_pos2 = calc_mean_pos(at,iatmos2)
            mean_pos3 = calc_mean_pos(at,iatoms3)
            
            
            diff12 = triplets_diffs(1:3,i_desc)
            diff13 = triplets_diffs(4:6,i_desc)
            
            
            do ii=1,3
                diffs1(:,ii) = diff_min_image(at,at%pos(:,iatoms1(ii)),mean_pos1)
                diffs2(:,ii) = diff_min_image(at,at%pos(:,iatmos2(ii)),mean_pos2)
                diffs3(:,ii) = diff_min_image(at,at%pos(:,iatoms3(ii)),mean_pos3)
            end do
            
            
            do ii=1,3
                pos(:,  ii) = mean_pos1          - diffs1(:,ii)
                pos(:,3+ii) = mean_pos1 + diff12 - diffs2(:,ii)
                pos(:,6+ii) = mean_pos1 + diff13 - diffs3(:,ii)
            end do
               
            ! Calcualate descriptor
            block 
                real(dp) a0, shft, pow
                real(dp), dimension(3) :: a3, b3
                integer ii, jj, ij

                shft  = this%dist_shift
                pow   = this%power
                gdesc0 = 0
                do ij = 1, nd
                    jj = this%component_atoms(ij,1)
                    ii = this%component_atoms(ij,2)
                    a3 = pos(:,ii)-pos(:,jj)
                    a0 = norm(a3)
                    ! rvec(:,ij) = a3
                    rlen(ij) = a0
                    desc0(ij) = (a0 + shft) ** pow
                    
                    b3 = - pow * (a0+shft)**(pow-1.0_dp)  * a3 / a0
                    gdesc0(ij,:,jj) = b3
                    gdesc0(ij,:,ii) = -b3
                    
                enddo
            endblock
            
            ! Calculate cutoff
            block 
                real(dp), dimension(3) :: oo12, oo13, oo23
                real(dp), dimension(3) :: gcut1, gcut2, gcut3
                real(dp) :: o12, o13, o23
                real(dp) :: cut12, cut13, cut23
                real(dp) :: dcut12, dcut13, dcut23
                real(dp) :: cutoff, width
                cutoff= this%cutoff
                width = this%cutoff_transition_width
                    
                oo12 = pos(:,4)-pos(:,1)
                oo13 = pos(:,7)-pos(:,1)
                oo23 = pos(:,7)-pos(:,4)
                o12 = norm(oo12)
                o13 = norm(oo13)
                o23 = norm(oo23)
                cut12 = coordination_function(o12,cutoff,width)
                cut13 = coordination_function(o13,cutoff,width)
                cut23 = coordination_function(o23,cutoff,width)
                cut = (cut12*cut13 + cut12*cut23 + cut13*cut23) / 3
                
                
                dcut12 = dcoordination_function(o12, cutoff, width)
                dcut13 = dcoordination_function(o13, cutoff, width)
                dcut23 = dcoordination_function(o23, cutoff, width)
                gcut1 = (- oo12 / o12 * dcut12 * (cut13+cut23) - oo13 / o13 * dcut13 * (cut12+cut23)) / 3
                gcut2 = (  oo12 / o12 * dcut12 * (cut13+cut23) - oo23 / o23 * dcut23 * (cut12+cut13)) / 3
                gcut3 = (  oo13 / o13 * dcut13 * (cut12+cut23) + oo23 / o23 * dcut23 * (cut12+cut13)) / 3
                
                gcuts = 0 
                gcuts(:,1) = gcut1
                gcuts(:,4) = gcut2
                gcuts(:,7) = gcut3
            endblock
            
            if(my_do_descriptor) then
                descriptor_out%x(i_desc)%has_data = .true.
                descriptor_out%x(i_desc)%data(:) = desc0
                descriptor_out%x(i_desc)%ci(:) = iatoms123
                descriptor_out%x(i_desc)%covariance_cutoff = cut
            end if
        
            if(my_do_grad_descriptor) then
                descriptor_out%x(i_desc)%has_grad_data = .true.
                descriptor_out%x(i_desc)%grad_data(:,:,:) = gdesc0
                descriptor_out%x(i_desc)%ii(:) = iatoms123
                descriptor_out%x(i_desc)%pos(:,:) = pos
                descriptor_out%x(i_desc)%grad_covariance_cutoff(:,:) = gcuts
            end if 
            
        endblock
        write(12,'(a,*(f10.5))') "  desc = ", descriptor_out%x(i_desc)%data
        write(12,'(a,*(f10.5))') "g desc = ", descriptor_out%x(i_desc)%grad_data
        write(12,'(a,*(f10.5))') "   cut = ", descriptor_out%x(i_desc)%covariance_cutoff
        write(12,'(a,*(f10.5))') "g  cut = ", descriptor_out%x(i_desc)%grad_covariance_cutoff
        stop "testing"
    
   enddo the_loop

    deallocate(imonomers1)
    deallocate(imonomers2)
    deallocate(imonomers3)
    if (allocated(itriplets)) deallocate(itriplets)
    if (allocated(triplets_diffs)) deallocate(triplets_diffs)
    if (allocated(triplets_diffs_map)) deallocate(triplets_diffs_map)



    call system_timer('water_trimer_calc')

endsubroutine water_trimer_calc !////////////////////////////////////////////////////////////////////////////////////////////////////////////

