module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls(calc_source)
    implicit none

    logical(lgt), intent(in)                      :: calc_source
    integer(i4b)                                  :: i, j, l, l_num, x_num, n_spline, l_num_hires, n_k_full_size, n_x_full_size
    real(dp)                                      :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand, ls_hires
    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     allocatable, dimension(:)       :: x_hires, k_hires, cls_hires, eta0_min_eta

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    n_k_full_size = 5000
    n_x_full_size = 5000

    

    ! Task: Get source function from evolution_mod
    allocate(S(n_x_full_size,n_k_full_size))
    allocate(x_hires(n_x_full_size))
    allocate(k_hires(n_k_full_size))

    filename = "output/milestone4/source.bin"
    inquire(file=filename,exist=exist)
    if(exist .and. (.not. calc_source)) then
      write(*,*) "Reading Source Function"
      open(1,file=filename,form="unformatted")
      read(1) k_hires, x_hires, S
      close(1)
    else
      write(*,*) "Calculating Source Function"
      write(*,*) "Integrating pertubation"
      call initialize_perturbation_eqns
      call integrate_perturbation_eqns


      call get_hires_source_function(k_hires, x_hires, S)

      open(1,file=filename,form="unformatted")
      write(1)  k_hires, x_hires, S
      close(1)
    end if

    


    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    z_spline = [((3500.d0 -1.d0)/(n_spline-1.d0)*i,i=1,n_spline)]

    filename = "output/milestone4/jl.bin"
    inquire(file=filename,exist=exist)
    if(exist) then
      open(1,file=filename,form="unformatted")
      read(1) j_l, j_l2
      close(1)

    else
      write(*,*) "Could not find file for the Bessel Functions, will calculate and write"
      ! Calculating j_l and j_l2

    
      do l=1, l_num
        j_l(l,1) = 0.d0
        do i=2, n_spline
          call sphbes(ls(l), z_spline(i),j_l(i,l))
        end do
        !call spline(z_spline,j_l(:,l),1.d30,1.d30,j_l2(:,l))
      end do

      do l=1, l_num
        call spline(z_spline,j_l(:,l),1.d30,1.d30,j_l2(:,l))
      end do

      open(1,file=filename,form="unformatted")
      write(1) j_l, j_l2
      close(1)
    end if


    ! Overall task: Compute the C_l's for each given l

    allocate(Theta(l_num,n_k_full_size))
    allocate(integrand(n_x_full_size))
    allocate(eta0_min_eta(n_x_full_size))
    allocate(cls(l_num))

    eta0_min_eta = [(get_eta(0.d0)-get_eta(x_hires(i)), i=1,n_x_full_size)]


    ! ########################################################
    ! Before integrating, we start with a sanity check of S*Jl
    ! ########################################################


    j = locate_dp(k_hires, 340.d0*H_0/c)
    l = 17

    do i = 1, n_x_full_size
      integrand(i) = S(i,j)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(j)*(eta0_min_eta(i)))
    end do

    write(*,*) "Writing sanity check to file"
    open(5, file="output/milestone4/sjl.dat")
    write(5,*) integrand
    close(5)
    open(4, file="output/milestone4/x.dat")
    write(4,*) x_hires
    close(4)


    ! #########################################
    ! Integrating to get C_l
    ! #########################################

    write(*,*) "Finding C_l"
    do l = 1, l_num

      write(*,*) "Integrating for l=", l

       ! Task: Compute the transfer function, Theta_l(k)

      !Optimize this be precalculating k(eta0 - eta)!
      do j=1,n_k_full_size
        do i=1,n_x_full_size
          integrand(i) = S(i,j)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(j)*(eta0_min_eta(i)))
          
        end do
        
        Theta(l,j) = (sum(integrand) - 0.5d0*(integrand(1) - integrand(n_x_full_size)))*(x_hires(n_x_full_size) - x_hires(1))/real(n_x_full_size - 1.d0)
        
      end do

       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       ! Task: Store C_l in an array. Optionally output to file

      cls(l) = 0.d0
      do j=1,n_k_full_size
        cls(l) = cls(l) + (c*k_hires(j)/H_0)**(n_s - 1.d0)*(Theta(l,j)**2.d0)/k_hires(j)
      end do
      cls(l) = cls(l)*(k_hires(n_k_full_size) - k_hires(1))/(n_k_full_size-1.d0)*ls(l)*(ls(l) + 1.d0)
      write(*,*) cls(l)

    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l

    write(*,*) "Done finding C_l, now spling"
    l_num_hires = 1200
    allocate(cls2(l_num))
    allocate(cls_hires(l_num_hires))
    allocate(ls_hires(l_num_hires))
    allocate(ls_dp(l_num))
    write(*,*) "Making ls_hires"
    ls_hires = [(l,l=1,l_num_hires)]![(ls(1) + (ls(l_num)-ls(1))/(l_num_hires-1.d0)*float(i-1),i=1,l_num_hires)]
    write(*,*) "Making ls_dp"
    ls_dp = [(ls(i),i=1,l_num)]
    
    write(*,*) "Doing the Spline"
    call spline(ls_dp,cls,1.d30,1.d30,cls2)

    write(*,*) "Splinting"

    do l=1,l_num_hires
      cls_hires(l) = splint(ls_dp,cls,cls2,ls_hires(l))
    end do

    write(*,*) "Saving to file"
    open(1,file="output/milestone4/cls.dat")
    open(2,file="output/milestone4/ls.dat")

    write(1,*) cls_hires
    write(2,*) ls_hires

    close(1)
    close(2)





  end subroutine compute_cls
  
end module cl_mod
