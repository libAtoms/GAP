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

!!!!!!!!
!!! This file was written by Jonatan Öström (@sujona, jonatan.ostrom@gmail.com) and Lars G.M. Pettersson, Stockholm University
!!! Here are implementations of water-trimer search routines, that work only for orthogonal unit cells with the shortest dimension longer than 2 x (3-body cutoff)
!!!!!!!!

module find_water_triplets
    ! use backend
    implicit none
    integer, parameter :: dp = kind(0d0)
    
    interface insert 
        module procedure insert_i2, insert_r2, insert_i1
    end interface
    
    
contains

! ////////////////////////////////////////////////
! Dynamic insert in allocatable array
! insert(a,i,x) does a(i) = x but reallocates a to length 2*i if len(a) < i

subroutine insert_i1(array,ii,val)
    integer, intent(inout), allocatable :: array(:)
    integer, intent(in) :: ii, val
    integer, allocatable :: tmp(:)
    if (ii>size(array))then
        tmp = array
        deallocate(array)
        allocate(array(2*ii))
        array(:size(tmp)) = tmp
    endif
    array(ii) = val
end
subroutine insert_i2(array,ii,val)
    integer, intent(inout), allocatable :: array(:,:)
    integer, intent(in) :: ii, val(:)
    integer, allocatable :: tmp(:,:)
    if (ii>size(array,2))then
        tmp = array
        deallocate(array)
        allocate(array(size(tmp,1),2*ii))
        array(:,:size(tmp,2)) = tmp
    endif
    array(:,ii) = val
end
subroutine insert_r2(array,ii,val)
    real(dp), intent(inout), allocatable :: array(:,:)
    real(dp), intent(in) :: val(:)
    integer, intent(in) :: ii
    real(dp), allocatable :: tmp(:,:)
    if (ii>size(array,2))then
        tmp = array
        deallocate(array)
        allocate(array(size(tmp,1),2*ii))
        array(:,:size(tmp,2)) = tmp
    endif
    array(:,ii) = val
end

! ////////////////////////////////////////////////
! Utilitity routines

subroutine min_img(XO,NO,i1,i2,box,lsq,x12,s12)
    real(dp), intent(in)    :: XO(3,NO), box(3)
    integer, intent(in)     :: NO, i1, i2
    real(dp), intent(out)   :: lsq
    real(dp), intent(out)   :: x12(3)
    integer, intent(out)    :: s12(3)
    x12 = XO(:,i1)-XO(:,i2)
    s12 = nint(x12/box)
    x12 = x12 - s12*box
    lsq = sum(x12**2)
end

subroutine min_img_c(XO,XC,NO,i1,i2,box,lsq,x12,c12,s12)
    real(dp), intent(in)    :: XO(3,NO), XC(3,NO), box(3)
    integer, intent(in)     :: NO, i1, i2
    real(dp), intent(out)   :: lsq
    real(dp), intent(out)   :: x12(3), c12(3)
    integer, intent(out)    :: s12(3)
    c12 = XC(:,i2)-XC(:,i1)
    s12 = -nint(c12/box)
    x12 = XO(:,i2)-XO(:,i1)
    x12 = x12 + s12*box
    c12 = c12 + s12*box
    lsq = sum(x12**2)
end

subroutine average_position(NW,XW,XC,XO)
    integer, intent(in) :: NW
    real(dp), intent(in) :: XW(3,NW*3)
    real(dp), intent(out) :: XC(3,NW),XO(3,NW)
    integer ii,jj,kk
    do ii = 1,NW
        jj = 3*(ii-1) !+1,2,3 for O,H,H
        XO(:,ii) = XW(:,jj+1)
        XC(:,ii) = 0
        do kk = 1,3
            XC(:,ii) = XC(:,ii) + XW(:,jj+kk)/3
        enddo
    enddo
end

function wall_time() result(tt)
    integer count,count_rate
    real(dp) tt
    call system_clock(count,count_rate)
    tt  = dble(count)/count_rate
end
subroutine take(tt,text) 
    real(dp), intent(inout) :: tt
    real(dp) :: old
    character(*), intent(in) :: text
    old = tt
    tt = wall_time()
    print'(f6.3,a)',tt-old,"s "//text
end


! ////////////////////////////////////////////////
! Finding triplets

subroutine find_triplets_brute_force(XW,NW,box,rcut,n_trip)
    ! find water triplets with brute force
    integer, intent(in) :: NW
    real(dp), intent(in) :: XW(3,NW*3), box(3), rcut
    integer, intent(out) :: n_trip
    real(dp) :: rcut2, d2ij, d2ik,d2jk, XC(3,NW),XO(3,NW)
    integer ii,jj,kk, n_short !n_2, n_3, 
    real(dp),dimension(3) :: xij,xik,xjk ! Oxygen diff
    real(dp),dimension(3) :: cij,cik,cjk ! Center diff
    integer, dimension(3) :: sij,sik,sjk ! Shifts
    rcut2 = rcut**2
    
    call average_position(NW,XW,XC,XO)
    
    !!! $omp parallel do private(d2ij, d2ik,d2jk,n_short) reduction(+:n_2,n_3)
    write(13,'(a)') ' '
    write(13,'(a)') ' TRIPLETS'
    write(13,'(a)') ' '
    do ii = 1,NW
        do jj = ii+1,NW
            call min_img_c(XO,XC,NW,ii,jj,box,d2ij,xij,cij,sij)
            do kk = jj+1,NW
                call min_img_c(XO,XC,NW,ii,kk,box,d2ik,xik,cik,sik)
                call min_img_c(XO,XC,NW,jj,kk,box,d2jk,xjk,cjk,sjk)
                n_short = count([d2ij,d2ik,d2jk]<rcut2)
                ! Counting "linear" and "triangular" separately
                if (n_short>1)then
                    n_trip = n_trip+1
                    write(13,'(3i5,2(3i3,2x),2(3f10.5,2x))') ii,jj,kk, sij, sik, cij,cik
                endif
            enddo
        enddo
    enddo
end

subroutine find_pairs_jona(XW,NO,box,rcut,n_pair,id2, sh2, dx2)
    integer, intent(in) :: NO
    real(dp), intent(in) :: XW(3,NO*3), box(3), rcut
    integer, intent(out) :: n_pair
    real(dp) :: rcut2, d2ij, xij(3), cij(3), XO(3,NO*3), XC(3,NO*3)
    integer ii,jj, sij(3)
    integer,  allocatable, dimension(:,:), intent(out)  :: id2, sh2
    real(dp), allocatable, dimension(:,:), intent(out)  :: dx2
    allocate(id2(2,0), sh2(3,0), dx2(3,0))
    call average_position(NO,XW,XC,XO)
    rcut2 = rcut**2
    n_pair = 0
    do ii = 1,NO
        do jj = ii+1,NO
            call min_img_c(XO,XC,NO,ii,jj,box,d2ij,xij,cij,sij)
            if (d2ij<rcut2) then
                n_pair = n_pair + 1
                call insert(id2, n_pair, [ii,jj])
                call insert(sh2, n_pair, sij)
                call insert(dx2, n_pair, xij)
            endif
        enddo
    enddo
    
    ! Automatic reallocation to the correct size upon assignment
    id2 = id2(:,:n_pair)
    sh2 = sh2(:,:n_pair)
    dx2 = dx2(:,:n_pair)
    
end

subroutine find_triplets_jona(XW,NO,box,rcut,n_trip,id3,sh3,dx3)
    ! Finds water triplets using a full pair list
    integer , intent(in)    :: NO
    real(dp), intent(in)    :: XW(3,NO*3), box(3), rcut
    integer , intent(out)   :: n_trip
    
    ! Allocatable pair-lists for use in trimers (with trimer cutoff)
    integer , allocatable       :: nindex(:), sh2(:,:)
    real(dp), allocatable       :: dx2(:,:)
    
    real(dp), dimension(3,NO)   :: XC,XO
    real(dp), dimension(3)      :: xij , cij, cik
    integer , dimension(3)      :: sij, sik
    
    integer     :: ii, jj, kk, jl, kl
    integer     :: i0, j0, i1, n_pair
    real(dp)    :: rcut2, d2ij
    integer     :: ncount(NO), offset(NO+1)
    
    ! Allocatable output arrays (id=atom index, sh=shift, dx=diff)
    integer,  allocatable, dimension(:,:), intent(out) :: id3, sh3
    real(dp), allocatable, dimension(:,:), intent(out) :: dx3
    
    allocate(id3(3,0), sh3(6,0), dx3(6,0))
    allocate(nindex(0),sh2(3,0),dx2(3,0))
    
    ! Get molecular center (XC) and oxygen positions (XO) from water coordinates (XW)
    call average_position(NO,XW,XC,XO)
    
    ! Square cutoffs for speed
    rcut2 = rcut**2
    
    n_pair = 0
    n_trip = 0
    ncount = 0
    offset = 0
    
    ! PAIRS (with triplet cutoff)
    do ii = 1,NO
        do jj = 1,NO ! Double counting for full neighbor-lists
            if (ii==jj) cycle
            call min_img_c(XO,XC,NO,ii,jj,box,d2ij,xij,cij,sij)
            if (d2ij<rcut2) then
                n_pair      = n_pair + 1
                ncount(ii)  = ncount(ii) + 1
                call insert(nindex, n_pair, jj)
                call insert(sh2, n_pair, sij)
                call insert(dx2, n_pair, cij)
            endif
        enddo
        offset(ii+1) = n_pair
    enddo
    
    ! TRIPLETS
    ! a, b, c represent actual atomic indices
    ! All unique triplets with a<b<c are: 
    ! 1. c-a-b(-c)
    ! 2. a-b-c
    ! 3. a-c-b
    ! where "-" is a distance shorter than cutoff. 
    ! Loop indices ii, jj, kk run over molecules.
    ! G(ii) is the set of nighbors of ii, etc.
    do ii = 1,NO
        i0 = offset(ii)
        i1 = offset(ii+1)
        ! Take jj>ii in G(ii)
        do jl = 1,ncount(ii)
            jj  = nindex(i0+jl)
            j0  = offset(jj)
            sij = sh2(:,i0+jl)
            cij = dx2(:,i0+jl)
            if (ii<jj)then
                ! 1. c-a-b(-c): take kk>jj in G(ii) to get ii=a, jj=b, kk=c
                do kl = jl+1,ncount(ii)
                    kk = nindex(i0+kl)
                    n_trip  = n_trip + 1
                    sik     = sh2(:,i0+kl)
                    cik     = dx2(:,i0+kl)
                    call insert(id3, n_trip, [ii,jj,kk])
                    call insert(sh3, n_trip, [sij,sik])
                    call insert(dx3, n_trip, [cij,cik])
                enddo
                ! 2. a-b-c & 3. a-c-b: take kk>ii in G(jj)\G(ii) so that ii=a, jj=b/c, kk=c/b
                do kl = 1,ncount(jj)
                    kk = nindex(j0+kl)
                    if (ii<kk.and..not.any(nindex(i0+1:i1)==kk)) then
                        n_trip  = n_trip + 1
                        sik     = sh2(:,i0+jl) + sh2(:,j0+kl)
                        cik     = dx2(:,i0+jl) + dx2(:,j0+kl)
                        call insert(id3, n_trip, [ii,jj,kk])
                        call insert(sh3, n_trip, [sij,sik])
                        call insert(dx3, n_trip, [cij,cik])
                    endif
                enddo
            endif
        enddo
    enddo
    
    ! Automatic reallocation to the correct size upon assignment
    id3 = id3(:,:n_trip)
    sh3 = sh3(:,:n_trip)
    dx3 = dx3(:,:n_trip)
    
end

subroutine find_triplets_lars(XW,NO,box,rcut,id3,dx3,n_trip)!sh3,
    ! Finds water triplets using a pair list of only unique pairs
    integer, intent(in) :: NO
    real(dp),  intent(in) :: rcut, box(3), XW(3,3*NO)
    real(dp),  intent(out), allocatable :: dx3(:,:)
    integer, intent(out), allocatable :: id3(:,:)!, ma p3(:)!sh3(:,:)
    integer, intent(out) :: n_trip
    ! internal
    integer i0, i1, ii, jj, jl, kk, kl, jn, j0, k0
    integer n_pair, num(NO), ioff(NO)
    integer, allocatable :: neighbor(:) !sh2(:,:), 
    real(dp) , allocatable :: dx2(:,:)
    real(dp)  rcut2, dd, oij(3)
    real(dp)  XC(3,NO), XO(3,NO)
    real(dp),  dimension(3) :: cij,cik
    integer, dimension(3) :: sij!,sik
    
    allocate(id3(3,0),dx3(6,0))!,sh3(6,0)
    allocate(neighbor(0),dx2(3,0)) !sh2(3,0),
    
    call average_position(NO,XW,XC,XO)
    
    ! PAIRS
    rcut2 = rcut**2
    n_pair = 0
    ioff(1) = 0
    do ii = 1, NO - 1
        num(ii) = 0
        do jj = ii + 1, NO
            oij = XO(:,jj) - XO(:,ii)
            cij = XC(:,jj) - XC(:,ii)
            sij = -nint(oij/box)
            oij = oij + sij*box
            cij = cij + sij*box
            dd = sum(oij**2)
            if (dd .gt. rcut2) cycle
            num(ii) = num(ii) + 1
            n_pair = n_pair + 1
            call insert(neighbor, n_pair, jj)
            ! call insert(sh2, n_pair, sij)
            call insert(dx2, n_pair, cij)
        enddo
        ioff(ii+1) = n_pair
    enddo
    ! ioff(NO)  = npairs ! important
    num(NO)   = 0 ! important
    
    
    ! TRIPLETS
    ! With a < b < c being actual atomic indices of unique oxygen/molecule triplets
    ! all forms we have are three cases
    ! 1. c-a-b-(c-)
    ! 2. a-b-c
    ! 3. a-c-b
    ! j -> N(i) means j runs over the (N)eighbors of i
    ! case 1. c-a-b-(c-) : let i<j<k and j -> N(i) and k -> N(i)
    n_trip = 0
    do ii = 1, NO - 1
        i0 = ioff(ii)
        i1 = ioff(ii+1)
        do jl = 1, num(ii) - 1
            jj = neighbor(i0 + jl)
            ! sij = sh2(:,i0 + jl)
            cij = dx2(:,i0 + jl)
            do kl = jl + 1, num(ii)
                kk = neighbor(i0 + kl)
                n_trip = n_trip + 1
                ! sik = sh2(:,i0 + kl)
                cik = dx2(:,i0 + kl)
                
                call insert(id3, n_trip, [ii,jj,kk])
                ! call insert(sh3, n_trip, [sij,sik])
                call insert(dx3, n_trip, [cij,cik])
                
            enddo
        enddo
        ! case 2. a-b-c : let i<j<k and find j->N(i) and k -> N(j)\N(i)
        do jl = 1, num(ii)
            jj = neighbor(i0 + jl)
            j0 = ioff(jj)
            ! sij = sh2(:,i0 + jl)
            cij = dx2(:,i0 + jl)
            do kl = 1, num(jj) !(1)
                kk = neighbor(j0 + kl)   !(2)
                if (any(neighbor(i0+1:i1)==kk)) cycle
                n_trip = n_trip + 1
                
                ! sik = sh2(:,i0 + jl) + sh2(:,j0 + kl)
                cik = dx2(:,i0 + jl) + dx2(:,j0 + kl)
                
                call insert(id3, n_trip, [ii,jj,kk])
                ! call insert(sh3, n_trip, [sij,sik])
                call insert(dx3, n_trip, [cij,cik])
            enddo
            ! case 3. a-c-b : let i=a < k=b < j=c, so for j->N(i) find k->N(j)\N(i)
            ! since k<j is not listed in N(j) we take k->[i+1,j-1], exclude k->N(i) and include k when j->N(k)
            do kk = ii + 1, jj - 1
                k0 = ioff(kk)
                do jn = 1, num(kk)
                    if (jj /= neighbor(k0+jn)) cycle ! include j->N(k)
                    if (any(kk==neighbor(i0+1:i1))) cycle ! exclude k->N(i)
                    n_trip = n_trip + 1
                    
                    ! sik = sh2(:,i0 + jl) - sh2(:,k0 + jn)
                    cik = dx2(:,i0 + jl) - dx2(:,k0 + jn)
                    
                    ! insert in i<k<j order
                    call insert(id3, n_trip, [ii,kk,jj])
                    ! call insert(sh3, n_trip, [sik,sij])
                    call insert(dx3, n_trip, [cik,cij])
                enddo
            enddo
        enddo
    enddo
    ! Automatic reallocation to the correct size upon assignment
    id3 = id3(:,:n_trip)
    ! sh3 = sh3(:,:n_trip)
    dx3 = dx3(:,:n_trip)
    ! allocate(map3(n_trip))
    ! do ii = 1,n_trip
    !     map3(ii) = ii
    ! enddo
end

! ////////////////////////////////////////////////
! Printing

subroutine print_pairs_or_triplets(idx,shift,diff,NO)
    real(dp), allocatable,dimension(:,:), intent(inout) :: diff
    integer,  allocatable,dimension(:,:), intent(inout) :: idx,shift
    integer, intent(in) :: NO
    integer ii, nn
    character(:), allocatable :: form, phrase, name
    character(100) ctemp
    nn = size(idx,2)
    
    if (size(idx,1)==2)then
        form = '(2i5,3i3,3f10.5)'
        phrase = ' PAIRS Fast #'
        write(ctemp,'(a,i0,a)') "water",NO,"-2.dat"
    else
        form = '(3i5,2(3i3,2x),2(3f10.5,2x))'
        phrase = "Triplets Fast #"
        write(ctemp,'(a,i0,a)') "water",NO,"-3.dat"
    endif
    name = trim(ctemp)
    
    open(44,file=name,action="write")
    write(44,'(a,i0)') phrase,nn
    do ii = 1,nn
        write(44,form) idx(:,ii),shift(:,ii),diff(:,ii)
    enddo
    close(44)
end 


end module
