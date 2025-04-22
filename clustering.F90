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

#include "error.inc"

module clustering_module

  ! use libatoms_module
  use error_module
  use system_module ! , only : dp, optional_default, ran_uniform, reallocate
  use linearalgebra_module

  implicit none
  private

  public :: pivot, bisect_kmedoids, cluster_kmeans, select_uniform, cluster_fuzzy_cmeans, cur_decomposition

  integer, parameter  :: n_trial = 10
  integer, parameter  :: n_trial_k_med = 100
  real(dp), parameter :: cluster_jitter = 1.0e-7_dp
  real(dp), parameter :: KMEANS_THRESHOLD = 1.0e-6_dp

  type lst
     integer, dimension(:), allocatable :: object
     integer :: medoid
     real(dp) :: sse
     integer :: N
  endtype lst

  type clstr
     type(lst), dimension(:), allocatable :: cluster
     real(dp), dimension(:,:), pointer :: dm
     integer :: N
  endtype clstr

  contains

  subroutine distance_matrix(x,dm,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     real(dp), dimension(:,:), intent(out) :: dm
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), target, optional :: theta

     real(dp), dimension(:), pointer :: my_theta => null()
     real(dp) :: my_theta_fac
     integer :: i, j, d, n

     my_theta_fac = optional_default(1.0_dp, theta_fac)
     d = size(x,1)
     n = size(x,2)

     if( present(theta) ) then
        if( size(theta) == d) then
           my_theta => theta
        else
           allocate(my_theta(d))
           my_theta = theta(1)
        endif
     else
        allocate(my_theta(d))

        do i = 1, d
           my_theta(i) = ( maxval(x(i,:)) - minval(x(i,:)) )
   !        theta(i) = sqrt( & !take square root
   !                         & sum( x(i,:)**2 ) / size(x(i,:)) - &
   !                         & (sum( x(i,:) ) / size(x(i,:)))**2 )
           if( my_theta(i) .feq. 0.0_dp ) my_theta(i) = 1.0_dp
        enddo
        my_theta = my_theta * my_theta_fac
     endif

     do i = 1, n
        do j = i + 1, n
           dm(j,i) = cluster_jitter*ran_uniform()
        enddo
        dm(i,i) = 0.0_dp
     enddo

!$omp parallel do default(none) shared(dm,n,x,my_theta) private(i,j) schedule(dynamic)
     do i = 1, n
        do j = i + 1, n
           dm(j,i) = dm(j,i) + sqrt( sum( ( (x(:,j) - x(:,i)) / my_theta )**2 ) )
           dm(i,j) = dm(j,i)
        enddo
     enddo
!$omp end parallel do

     do i = 1, n
        do j = i + 1, n
           dm(i,j) = dm(j,i)
        enddo
     enddo

     if( present(theta) ) then
        my_theta => null()
     else
        deallocate(my_theta)
     endif

  endsubroutine distance_matrix

  subroutine pca(x,x_mean,v)

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:), intent(out) :: x_mean
    real(dp), dimension(:,:), intent(out) :: v

    real(dp), dimension(:), allocatable :: diag_c
    real(dp), dimension(:,:), allocatable :: cov
    integer :: i, j, d, n

    d = size(x,1)
    n = size(x,2)
    allocate(cov(d,d),diag_c(d))

    x_mean = sum(x,dim=2) / n ! empirical mean

    do i = 1, d
       do j = 1, d
          cov(j,i) = dot_product(x(i,:),x(j,:)) / n - x_mean(i)*x_mean(j)
       enddo
    enddo

    call diagonalise(cov,diag_c, evects=v)

    deallocate(cov, diag_c)

  endsubroutine pca

  subroutine pivot(x,pivout,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: pivout
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), optional :: theta

     real(dp), dimension(:,:), allocatable :: knn
     real(dp), dimension(:), allocatable :: ktmp
     integer, dimension(:), allocatable :: pivin

     integer :: stat, i, j, k, d, m, n, jtmp, jmax
     real(dp) :: dmax

     d = size(x,1)
     n = size(x,2)

     m = size(pivout)

     if( m > n ) call system_abort('pivot: required number of changes ('//m//') greater than possible number of changes ('//n//')')

     allocate(knn(n,n),stat=stat)
     if(stat /=0 ) call system_abort('pivot: could not allocate knn matrix.')

     allocate(pivin(n),ktmp(n))

     call distance_matrix(x,knn,theta_fac=theta_fac,theta=theta)
     do i = 1, n
        do j = 1, n
           knn(j,i) = exp(-0.5_dp*knn(j,i))
        enddo
     enddo

     pivin = (/ (i, i=1,n) /)

     do k = 1, m
        dmax = 0.0_dp
        do j = k, n
           if( dmax < knn(j,j) ) then
              jmax = j
              dmax = knn(j,j)
           endif
        enddo
        if( jmax /= k ) then
            jtmp = pivin(jmax)
            pivin(jmax) = pivin(k)
            pivin(k) = jtmp

            ktmp = knn(k,:)
            knn(k,:) = knn(jmax,:)
            knn(jmax,:) = ktmp

            ktmp = knn(:,k)
            knn(:,k) = knn(:,jmax)
            knn(:,jmax) = ktmp
         endif

         knn(k,k) = sqrt(knn(k,k))

         knn(k+1:n,k) = knn(k+1:n,k)/knn(k,k)
         do j = k+1, n
            knn(j:n,j) = knn(j:n,j) - knn(j:n,k)*knn(j,k)
         enddo

         do j = 1, n
            do i = j+1,n
               knn(j,i) = knn(i,j)
            enddo
         enddo
      enddo

      pivout = pivin(1:m)

      deallocate(knn,pivin,ktmp)

  endsubroutine pivot

  subroutine bisect_kmedoids(dat,n_clusters_in, c,med, theta_fac,theta, is_distance_matrix)
     real(dp), dimension(:,:), intent(in), target :: dat
     integer, intent(in) :: n_clusters_in
     integer, dimension(:), intent(out),optional :: c, med
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), optional :: theta
     logical, intent(in), optional :: is_distance_matrix

     type(clstr) :: my_cluster, tmp

     logical :: must_calculate_distance
     real(dp), dimension(:,:), allocatable, target :: dm

     real(dp), dimension(:), allocatable :: dv
     real(dp) :: max_sse, min_sse, sse

     integer, dimension(:), allocatable :: sub_cluster1, sub_cluster2, sub_cluster1_min, sub_cluster2_min
     integer, dimension(1) :: ml
     integer  :: stat, i, j, k, km, m, n, nc, &
     lo_med, hi_med, lo_med_new, hi_med_new, lo_med_min, hi_med_min, n1, n2, n1_min, n2_min, iter

     must_calculate_distance = .not. optional_default(.true., is_distance_matrix)

     n = size(dat,2)
     if (.not. must_calculate_distance) then
       if (size(dat,1) /= n) call system_abort('is_distance_matrix but not square')
     endif

     if( n_clusters_in > n ) call system_abort('bisect_kmedoids: required number of cluster greater than total number of data points')

     if(present(c) ) c = 0

     if (must_calculate_distance) then
       allocate(dm(n,n), stat=stat)
       if(stat /=0 ) call system_abort('bisect_kmedoids: could not allocate dm matrix.')

       call print('Started distance matrix calculation', verbosity=PRINT_NERD)
       call distance_matrix(dat, dm, theta_fac=theta_fac,theta=theta)
       call print('Finished distance matrix calculation', verbosity=PRINT_NERD)
       my_cluster%dm => dm
     else
       my_cluster%dm => dat
     endif

     ! start clustering
     my_cluster%N = 1                               ! start with one big cluster
     allocate( my_cluster%cluster(1) )
     my_cluster%cluster(1)%N = n                    ! put every object in the initial cluster
     allocate( my_cluster%cluster(1)%object(n) )
     my_cluster%cluster(1)%object = (/(i,i=1,n)/)

     allocate(dv(n)) ! distance vector, the sum of square of distances of points from central object
     dv = sum(my_cluster%dm,dim=1)
     my_cluster%cluster(1)%sse = minval( dv )  ! determine initial medoid, the object that is the
     ml = minloc( dv )                         ! closest to any other object in cluster
     my_cluster%cluster(1)%medoid = ml(1)
     deallocate(dv)

     ! main loop starts here, bisects initial clusters until desired number of
     ! clusters are found

     iter = 0
     do
        iter = iter + 1
        call print("Starting iteration "//iter,verbosity=PRINT_NERD)

        if( my_cluster%N == n_clusters_in )  exit
        max_sse = -1.0_dp                                 ! select cluster with greatest sse
        do j = 1, my_cluster%N
           if( max_sse < my_cluster%cluster(j)%sse ) then
              i = j
              max_sse = my_cluster%cluster(j)%sse
           endif
        enddo
        nc = my_cluster%cluster(i)%N
        if( nc==1 ) cycle
        allocate( sub_cluster1(nc), sub_cluster2(nc), sub_cluster1_min(nc),sub_cluster2_min(nc) )

        min_sse = huge(1.0_dp)
        do j = 1, n_trial
           m = ceiling( ran_uniform()*(nc-1) ) ! choose a bisecting point randomly
           ml = minloc( sum( my_cluster%dm( my_cluster%cluster(i)%object(:m), my_cluster%cluster(i)%object(:m) ), dim=1) )
           lo_med_new = my_cluster%cluster(i)%object(ml(1))

           ml = minloc( sum( my_cluster%dm( my_cluster%cluster(i)%object(m+1:), my_cluster%cluster(i)%object(m+1:) ), dim=1) )

           hi_med_new = my_cluster%cluster(i)%object(ml(1) + m)

           ! the median of the 2 subclusters determined
           lo_med = 0
           hi_med = 0

           ! perform k-medoid clustering on the two subclusters
           do km = 1, n_trial_k_med
              if( (lo_med_new == lo_med) .and. (hi_med_new == hi_med) ) exit
              lo_med = lo_med_new
              hi_med = hi_med_new
              n1 = 0
              n2 = 0
              !n1 = 1
              !n2 = 1
              !sub_cluster1(n1) = lo_med
              !sub_cluster1(n2) = hi_med

              do k = 1, my_cluster%cluster(i)%N
                 if( my_cluster%dm(lo_med,my_cluster%cluster(i)%object(k)) < &
                 & my_cluster%dm(hi_med,my_cluster%cluster(i)%object(k)) ) then
                    n1 = n1 + 1
                    sub_cluster1(n1) = my_cluster%cluster(i)%object(k)
                 else
                    n2 = n2 + 1
                    sub_cluster2(n2) = my_cluster%cluster(i)%object(k)
                 endif
              enddo

              ml = minloc( sum( my_cluster%dm( sub_cluster1(:n1), sub_cluster1(:n1) ), dim=1) )
              lo_med_new = sub_cluster1(ml(1))
              ml = minloc( sum( my_cluster%dm( sub_cluster2(:n2), sub_cluster2(:n2) ), dim=1) )
              hi_med_new = sub_cluster2(ml(1))
           enddo
           sse = sum( my_cluster%dm(lo_med_new,sub_cluster1(:n1)) ) + sum( my_cluster%dm(hi_med_new,sub_cluster2(:n2)) )

           ! choose the clustering that resulted the smallest sse
           if( sse < min_sse ) then
              min_sse = sse
              sub_cluster1_min = sub_cluster1
              sub_cluster2_min = sub_cluster2
              n1_min = n1
              n2_min = n2
              lo_med_min = lo_med_new
              hi_med_min = hi_med_new
           endif
        enddo

        ! now update the the clusters with the two new subclusters
        tmp = my_cluster

        do j = 1, my_cluster%N
           deallocate( my_cluster%cluster(j)%object )
        enddo
        deallocate( my_cluster%cluster )
        my_cluster%N = my_cluster%N + 1
        allocate( my_cluster%cluster( my_cluster%N ) )

        do j = 1, my_cluster%N - 1
           if( i == j ) then
              allocate( my_cluster%cluster(j)%object(n1_min) )
              my_cluster%cluster(j)%N = n1_min
              my_cluster%cluster(j)%object = sub_cluster1_min(:n1_min)
              my_cluster%cluster(j)%sse = sum( my_cluster%dm(lo_med_min,sub_cluster1_min(:n1_min)) )
              my_cluster%cluster(j)%medoid = lo_med_min
           else
              my_cluster%cluster(j) = tmp%cluster(j)
           endif
        enddo
        allocate( my_cluster%cluster(my_cluster%N)%object(n2_min) )
        my_cluster%cluster(my_cluster%N)%N = n2_min
        my_cluster%cluster(my_cluster%N)%object = sub_cluster2_min(:n2_min)
        my_cluster%cluster(my_cluster%N)%sse = sum( my_cluster%dm(hi_med_min,sub_cluster2_min(:n2_min)) )
        my_cluster%cluster(my_cluster%N)%medoid = hi_med_min

        do j = 1, tmp%N
           deallocate( tmp%cluster(j)%object )
        enddo
        deallocate( tmp%cluster, sub_cluster1, sub_cluster2, sub_cluster1_min, sub_cluster2_min )

        call kmedoid(my_cluster)
     enddo

     if( present(c) ) then
        do j = 1, my_cluster%N
           do k = 1, my_cluster%cluster(j)%N
              i = my_cluster%cluster(j)%object(k)
              c(i) = j
           enddo
        enddo
     endif

     if( present(med) ) then
        do j = 1, my_cluster%N
           med(j) =  my_cluster%cluster(j)%medoid
        enddo
     endif

     do j = 1, my_cluster%N
        deallocate( my_cluster%cluster(j)%object )
     enddo
     deallocate(my_cluster%cluster)
     if (allocated(dm)) deallocate(dm)

  endsubroutine bisect_kmedoids

  subroutine kmedoid(this)
     type(clstr), intent(inout) :: this

     type(clstr) :: tmp
     integer, dimension(:), allocatable :: medoids
     integer, dimension(1) :: ml
     integer :: n, j, k
     logical :: refined

     ! k-medoid-refinement
     n = size(this%dm,1)
     ! n: total number of objects

     tmp%N = this%N
     allocate( tmp%cluster(tmp%N), medoids(tmp%N) )
     do j = 1, tmp%N
        allocate( tmp%cluster(j)%object(n) )
        medoids(j) = this%cluster(j)%medoid
     enddo

     ! main loop starts here, perfom k-medoid clustering until medoids don't
     ! change anymore
     do
        do j = 1, tmp%N
           tmp%cluster(j)%N = 0
        enddo
        do j = 1, n
           ml = minloc( this%dm(j,medoids) ) ! determine to which medoid each object belongs
           k = ml(1)
           tmp%cluster(k)%N = tmp%cluster(k)%N + 1
           tmp%cluster(k)%object(tmp%cluster(k)%N) = j
        enddo

        ! re-determine the medoid in each cluster
        do j = 1, tmp%N
           ml = minloc( sum( this%dm( tmp%cluster(j)%object(:tmp%cluster(j)%N), &
           & tmp%cluster(j)%object(:tmp%cluster(j)%N) ), dim=1) )
           tmp%cluster(j)%medoid = tmp%cluster(j)%object(ml(1))
        enddo

        refined = .true.

        ! check whether medoids have changed
        do j = 1, tmp%N
           refined = refined .and. (tmp%cluster(j)%medoid == medoids(j))
           medoids(j) = tmp%cluster(j)%medoid
        enddo
        if(refined) exit
     enddo

     ! write results
     do j = 1, tmp%N
        deallocate( this%cluster(j)%object )
        allocate( this%cluster(j)%object( tmp%cluster(j)%N ) )
        this%cluster(j)%object = tmp%cluster(j)%object(:tmp%cluster(j)%N)
        this%cluster(j)%N = tmp%cluster(j)%N
        this%cluster(j)%medoid = tmp%cluster(j)%medoid
        this%cluster(j)%sse = sum( this%dm(this%cluster(j)%medoid,&
        & this%cluster(j)%object ) )

        deallocate( tmp%cluster(j)%object )
     enddo
     deallocate( tmp%cluster, medoids )

  endsubroutine kmedoid

  subroutine cluster_kmeans(x,cluster_index,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: cluster_index
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), target, optional :: theta

     real(dp), dimension(:), pointer :: my_theta => null()
     real(dp) :: my_theta_fac, d_min, d_ij, d_total, d_total_prev

     real(dp), dimension(:,:), allocatable :: cluster_centre
     integer, dimension(:), allocatable :: cluster_info
     integer :: d, n, m, i, j, k, cluster_info_old, iter, n_points_cluster_j
     logical :: cluster_same

     d = size(x,1)
     n = size(x,2)
     m = size(cluster_index)
     if( m > n ) call system_abort('cluster_kmeans: required number of clusters ('//m//') greater than total number of points ('//n//')')

     my_theta_fac = optional_default(1.0_dp, theta_fac)
     if( present(theta) ) then
        if( size(theta) == d) then
           my_theta => theta
        else
           allocate(my_theta(d))
           my_theta = theta(1)
        endif
     else
        allocate(my_theta(d))
        do i = 1, d
           my_theta(i) = ( maxval(x(i,:)) - minval(x(i,:)) )
           if( my_theta(i) .feq. 0.0_dp ) my_theta(i) = 1.0_dp
        enddo
        my_theta = my_theta * my_theta_fac
     endif

     allocate(cluster_centre(d,m),cluster_info(n))

     call fill_random_integer(cluster_index, n) !choose random points as cluster centres.

     cluster_centre = x(:,cluster_index)
     cluster_info = 0

     iter = 0
     d_total = huge(1.0_dp)
     do
        iter = iter + 1
        call print("iteration: "//iter,verbosity=PRINT_NERD)
        cluster_same = .true.

        d_total_prev = d_total
        d_total = 0.0_dp
!$omp parallel do default(none) shared(n,m,x,cluster_info,cluster_centre,my_theta) &
!$omp reduction(.and.:cluster_same) &
!$omp private(i,j,d_min,d_ij,cluster_info_old) reduction(+:d_total)
        do i = 1, n
           d_min = huge(0.0_dp)
           cluster_info_old = cluster_info(i)
           do j = 1, m
              d_ij = sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2)
              if( d_ij < d_min ) then
                 d_min = d_ij
                 cluster_info(i) = j
              endif
           enddo
           if( cluster_info_old /= cluster_info(i) ) cluster_same = cluster_same .and. .false.
           d_total = d_total + d_min
        enddo
!$omp end parallel do
        call print("cluster_kmeans iteration="//iter//" d_total="//d_total)

!$omp parallel do default(none) shared(x,cluster_centre,cluster_info,m,d,n) private(j,k,n_points_cluster_j)
        do j = 1, m
           n_points_cluster_j = count(cluster_info==j)
           if( n_points_cluster_j == 0 ) then
              cluster_centre(:,j) = x(:,ceiling(ran_uniform()*n))
           else
              do k = 1, d
                 cluster_centre(k,j) = sum(x(k,:),mask=(cluster_info==j)) / n_points_cluster_j
              enddo
           endif
        enddo
!$omp end parallel do
        if( cluster_same ) exit
        if( abs(d_total - d_total_prev) < KMEANS_THRESHOLD * d_total ) exit
     enddo

     do j = 1, m
        d_min = huge(0.0_dp)
        do i = 1, n
           d_ij = sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2)
           if( d_ij < d_min ) then
              d_min = d_ij
              cluster_index(j) = i
           endif
        enddo
     enddo

     deallocate(cluster_centre, cluster_info)

     if(present(theta)) then
        my_theta => null()
     else
        deallocate(my_theta)
     endif

  endsubroutine cluster_kmeans

  ! https://sites.google.com/site/dataclusteringalgorithms/fuzzy-c-means-clustering-algorithm
  subroutine cluster_fuzzy_cmeans(x,cluster_index,theta_fac,theta,fuzziness)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: cluster_index
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), target, optional :: theta
     real(dp), intent(in), optional :: fuzziness

     real(dp), dimension(:), pointer :: my_theta => null()
     real(dp) :: my_theta_fac, d_min, d_ij, d_total, d_total_prev

     real(dp), dimension(:,:), allocatable :: cluster_centre
     real(dp), dimension(:,:), allocatable :: w
     real(dp), dimension(:), allocatable, save :: wx_j, d_i
     real(dp) :: w_j, w_old, my_fuzziness, alpha
     integer :: d, n, m, i, j, iter
     logical :: cluster_same
!$omp threadprivate(d_i, wx_j)

     d = size(x,1)
     n = size(x,2)
     m = size(cluster_index)
     if( m > n ) call system_abort('cluster_fuzzy_cmeans: required number of clusters ('//m//') greater than total number of points ('//n//')')

     my_theta_fac = optional_default(1.0_dp, theta_fac)
     my_fuzziness = optional_default(4.0_dp, fuzziness)
     if( present(theta) ) then
        if( size(theta) == d) then
           my_theta => theta
        else
           allocate(my_theta(d))
           my_theta = theta(1)
        endif
     else
        allocate(my_theta(d))
        do i = 1, d
           my_theta(i) = ( maxval(x(i,:)) - minval(x(i,:)) )
           if( my_theta(i) .feq. 0.0_dp ) my_theta(i) = 1.0_dp
        enddo
        my_theta = my_theta * my_theta_fac
     endif

     allocate(cluster_centre(d,m), w(n,m))
!$omp parallel
     allocate(d_i(m), wx_j(d))
!$omp end parallel

     call fill_random_integer(cluster_index, n) !choose random points as cluster centres.

     cluster_centre = x(:,cluster_index)
     do i = 1, m
        do j = 1, d
           cluster_centre(j,i) = cluster_centre(j,i) + ( ran_uniform() - 0.5_dp ) * cluster_jitter
        enddo
     enddo

     w = 0.0_dp

     iter = 0
     d_total = huge(1.0_dp)
     do
        iter = iter + 1
        call print("iteration: "//iter,verbosity=PRINT_NERD)
        cluster_same = .true.

        d_total_prev = d_total
        d_total = 0.0_dp
        ! Calculate fuzzy membership
!$omp parallel do default(none) shared(n,m,my_theta,my_fuzziness,w,x,cluster_centre) &
!$omp private(i,j,alpha,w_old) reduction(.and.:cluster_same) reduction(+:d_total)
        do i = 1, n
           alpha = 0.0_dp
           do j = 1, m
              d_i(j) = sqrt(sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2))
              alpha = alpha + 1.0_dp / d_i(j)**(2.0_dp / (my_fuzziness - 1.0_dp))
           enddo

           do j = 1, m
              w_old = w(i,j)
              w(i,j) = 0.0_dp

              w(i,j) = 1.0_dp / d_i(j)**(2.0_dp / (my_fuzziness - 1.0_dp)) / alpha
              if( w_old .fne. w(i,j) ) cluster_same = cluster_same .and. .false.

              d_total = d_total + d_i(j)**2 * w(i,j)**my_fuzziness
           enddo
        enddo
!$omp end parallel do
        call print("cluster_fuzzy_cmeans iteration="//iter//" d_total="//d_total)

        ! Calculate fuzzy centres
!$omp parallel do default(none) shared(m,n,w,x,my_fuzziness,cluster_centre) &
!$omp private(i,j,w_j)
        do j = 1, m
           w_j = 0.0_dp
           wx_j = 0.0_dp

           do i = 1, n
              w_j = w_j + w(i,j)**my_fuzziness
              wx_j = wx_j + x(:,i) * w(i,j)**my_fuzziness
           enddo

           cluster_centre(:,j) = wx_j / w_j
        enddo
!$omp end parallel do

        call print("cluster_same: "//cluster_same,verbosity=PRINT_NERD)
        call print("d_total: "//d_total,verbosity=PRINT_NERD)
        call print("d_total_prev: "//d_total_prev,verbosity=PRINT_NERD)
        call print("d_total-d_total_prev: "//(d_total-d_total_prev),verbosity=PRINT_NERD)

        if( cluster_same ) exit
        if( abs(d_total - d_total_prev) < KMEANS_THRESHOLD * d_total ) exit
     enddo

     ! Allocate cluster centres to nearest points
     do j = 1, m
        d_min = huge(0.0_dp)
        do i = 1, n
           d_ij = sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2)
           if( d_ij < d_min ) then
              d_min = d_ij
              cluster_index(j) = i
           endif
        enddo
     enddo

     deallocate(cluster_centre, w)
!$omp parallel
     if(allocated(d_i)) deallocate(d_i)
     if(allocated(wx_j)) deallocate(wx_j)
!$omp end parallel

     if(present(theta)) then
        my_theta => null()
     else
        deallocate(my_theta)
     endif

  endsubroutine cluster_fuzzy_cmeans

  subroutine select_uniform(x,index_out)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: index_out

     integer :: i, d, n, m, n_grid, i_global, i_index_out, d_max
     integer, dimension(:), allocatable :: p_grid, i_hist, histogram, x_histogram, index_out_histogram
     real(dp), dimension(:), allocatable :: lower_bound, upper_bound, x_range

     d = size(x,1)
     n = size(x,2)
     m = size(index_out)

     if( n < m ) call system_abort('select_uniform: n = '//n//' < m = '//m)

     allocate(lower_bound(d), upper_bound(d), x_range(d), p_grid(d), i_hist(d))

     lower_bound = minval(x,dim=2)
     upper_bound = maxval(x,dim=2)
     x_range = upper_bound - lower_bound

     n_grid = ceiling( real(m, dp)**(1.0_dp / real(d, dp)) )
     d_max = floor( log(real(huge(1), dp)) / log(real(n_grid, dp)) )
     if (d > d_max) then
        call system_abort('select_uniform: Descriptor is too large ('//d//' > '//d_max//'). &
                          &Use another sparse method or descriptor.')
     end if
     p_grid = (/ ( n_grid**(i-1), i = 1, d ) /)

     allocate(histogram(n_grid**d))
     allocate(x_histogram(n),index_out_histogram(m))

     histogram = 0

     do i = 1, n
        ! for each datapoint x(:,i) compute the bin index in each d direction
        i_hist = nint( ( x(:,i) - lower_bound ) / x_range * (n_grid-1) ) + 1

        ! map the bin index to a flat histogram bin index
        i_global = sum((i_hist-1)*p_grid)+1
        histogram(i_global) = histogram(i_global) + 1

        ! the i-th datapoint belongs to the i_global-th index in the histogram
        x_histogram(i) = i_global
     enddo

     index_out = 0
     i_index_out = 0

     ! To monitor which bins the sparse points belong to.
     index_out_histogram = 0

     do i = 1, n
        ! That's the exit condition if all sparse points are assigned before we
        ! finish with the data points
        if( all(index_out /= 0) ) exit

        if( all(x_histogram(i) /= index_out_histogram) ) then
        ! We have just found a point which belongs to a bin that we haven't
        ! selected yet in the sparse points
           i_index_out = i_index_out + 1
           index_out(i_index_out) = i
           index_out_histogram(i_index_out) = x_histogram(i)
        endif
     enddo

     do while ( any(index_out == 0) )
        ! We haven't yet assigned all sparse points.

        ! Select a bin randomly
        i_global = ceiling( ran_uniform() * size(histogram) )

        ! cycle if the bin is empty
        if( histogram(i_global) == 0 ) cycle

        ! check if there are points belonging to this bin which we haven't
        ! selected yet
        if( count(x_histogram == i_global) == count(index_out_histogram == i_global) ) cycle

        do while (.true.)
           ! select a point from x which belongs to that bin and add it to
           ! the output.
           i = ceiling( ran_uniform() * n )
           if( x_histogram(i) /= i_global .or. any(index_out == i) ) then
              cycle
           else
              i_index_out = i_index_out + 1
              index_out(i_index_out) = i
              index_out_histogram(i_index_out) = x_histogram(i)
              exit
           endif
        enddo
     enddo

     deallocate(lower_bound, upper_bound, x_range, p_grid, i_hist, histogram, x_histogram,index_out_histogram)

     if (.not. all(index_out /= 0)) call system_abort('select_uniform: could not assign all sparse points')

  endsubroutine select_uniform

  subroutine cur_decomposition(this, index_out, rank, n_iter)
    ! based on 10.1073/pnas.0803205106

    real(dp), intent(in), dimension(:,:) :: this
    integer, dimension(:), intent(out) :: index_out
    integer, intent(in), optional :: rank, n_iter

    integer :: n
    integer :: expected_columns
    integer :: my_n_iter, my_rank
    type(LA_Matrix) :: LA_this
    real(dp), allocatable, dimension(:) :: p, s, p_minus_ran_uniform
    real(dp), allocatable, dimension(:,:) :: v
    integer :: j, l
    integer, allocatable, dimension(:), target :: p_index
    integer, pointer, dimension(:) :: tmp_index_out => null()
    real(dp), allocatable, dimension(:,:) :: C, Cp
    real(dp) :: err, min_err
    integer :: error

    expected_columns = size(index_out)

    if( expected_columns <= 0 ) then
       call print_warning("cur_decomposition: called with expected_columns "//expected_columns//", can't be zero or less")
       return
    endif

    call initialise(LA_this,this)

    my_n_iter = optional_default(1, n_iter)

    if (present(rank)) then
       call LA_Matrix_SVD_Allocate(LA_this,v=v,error=error)
       HANDLE_ERROR(error)
       call LA_Matrix_SVD(LA_this,v=v,error=error)
       HANDLE_ERROR(error)
       my_rank = rank
    else
       call LA_Matrix_SVD_Allocate(LA_this,s=s,v=v,error=error)
       HANDLE_ERROR(error)
       call LA_Matrix_SVD(LA_this,s=s,v=v,error=error)
       HANDLE_ERROR(error)
       my_rank = count(s > TOL_SVD) / 2
    endif

    n = size(v,1)
    allocate(p(n), p_minus_ran_uniform(n), p_index(n))
    allocate( C(size(this,1),expected_columns), Cp(expected_columns,size(this,1)) )

    p = sum(v(:,1:my_rank)**2, dim=2)
    p = p * expected_columns
    p = p / my_rank
    p = min(p,1.0_dp)

    if(my_n_iter <= 0) then ! do not do probabilistic selection of columns
       p_index = (/(j, j=1,n )/)
       p_minus_ran_uniform = -p
       call heap_sort(p_minus_ran_uniform,i_data=p_index)
       index_out = p_index(1:expected_columns)
    else
       min_err = huge(1.0_dp)
       do l = 1, my_n_iter

          ! randomly select columns according to the probabilities
          do j = 1, n
             p_minus_ran_uniform(j) = ran_uniform() - p(j)
             p_index(j) = j ! initialise index array
          end do

          call heap_sort(p_minus_ran_uniform,i_data=p_index)
          tmp_index_out => p_index(1:expected_columns)

          C = this(:,tmp_index_out)
          ! pinv: Moore-Penrose pseudo-inverse
          call pseudo_inverse(C,Cp)
          err = sum( (this - ( C .mult. Cp .mult. this))**2 )

          call print("cur_decomposition: iteration: "//l//", error: "//err)
          if(err < min_err) then        ! this happens at least once
             index_out = tmp_index_out
             min_err = err
          endif

       end do
    endif

    call finalise(LA_this)

    tmp_index_out => null()
    if(allocated(s)) deallocate(s)
    if(allocated(v)) deallocate(v)
    if(allocated(p)) deallocate(p)
    if(allocated(p_minus_ran_uniform)) deallocate(p_minus_ran_uniform)
    if(allocated(p_index)) deallocate(p_index)
    if(allocated(C)) deallocate(C)
    if(allocated(Cp)) deallocate(Cp)

  end subroutine cur_decomposition

endmodule clustering_module
