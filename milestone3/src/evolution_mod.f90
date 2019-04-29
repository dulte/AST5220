module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 6
  integer(i4b), parameter          :: n_before = 1000             !Number of points before recombination
  integer(i4b), parameter, private :: lmax_int = 6


  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Time array
  real(dp), allocatable, dimension(:)   :: x

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    ks = [(k_min + ((k_max-k_min)/(n_k-1))*(i-1),i=1,n_k)]
    ks(1) = 0.1 * H_0 / c
    ks(2) = 8.36 * H_0 / c
    ks(3) = 85.9 * H_0 / c
    ks(4) = 245.1 * H_0 / c
    ks(5) = 636.8 * H_0 / c
    ks(6) = 1.d3 * H_0 / c


    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t+n_before, 0:lmax_int, n_k))
    allocate(delta(0:n_t+n_before, n_k))
    allocate(delta_b(0:n_t+n_before, n_k))
    allocate(v(0:n_t+n_before, n_k))
    allocate(v_b(0:n_t+n_before, n_k))
    allocate(Phi(0:n_t+n_before, n_k))
    allocate(Psi(0:n_t+n_before, n_k))
    allocate(dPhi(0:n_t+n_before, n_k))
    allocate(dPsi(0:n_t+n_before, n_k))
    allocate(dv_b(0:n_t+n_before, n_k))
    allocate(dTheta(0:n_t+n_before, 0:lmax_int, n_k))

    ! Allocates array for x
    allocate(x(0:n_t+n_before))


    ! Open files for writing
    open(1,file='output/milestone3/x_modes.dat')
    open(2,file='output/milestone3/v.dat')
    open(3,file='output/milestone3/v_b.dat')
    open(4,file='output/milestone3/delta.dat')
    open(5,file='output/milestone3/delta_b.dat')
    open(6,file='output/milestone3/Phi.dat')
    open(7,file='output/milestone3/Theta.dat')
    open(8,file='output/milestone3/Psi.dat')
    open(9,file='output/milestone3/k.dat')


    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    Psi(0,:)     = -1.d0
    delta(0,:)   = 1.5*Phi(0,1)
    delta_b(0,:) = 1.5*Phi(0,1)
       
    do i = 1, n_k
       v(0,i)       = (c*ks(i))/(get_H_p(log(a_init))*2.d0)*Phi(0,1)
       v_b(0,i)     = (c*ks(i))/(get_H_p(log(a_init))*2.d0)*Phi(0,1)
       Theta(0,0,i) = 0.5*Phi(0,1)
       Theta(0,1,i) = -(c*ks(i))/(get_H_p(log(a_init))*6.d0)*Phi(0,1)
       Theta(0,2,i) = -(20.d0*c*ks(i))/(get_H_p(log(a_init))*45.d0*get_dtau(log(a_init)))*Theta(0,1,i)
       do l = 3, lmax_int
          Theta(0,l,i) = -l/(2.d0*l + 1.d0)*(c*ks(i))/(get_H_p(log(a_init))*get_dtau(log(a_init)))*Theta(0,l-1,i)
       end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2

    real(dp)     :: x_start_rec, x_end_rec, dx1, dx2, dx
    integer(i4b) :: n1, n2, index_tc

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx
    real(dp), allocatable, dimension(:) :: x_before_tc,x_after
  


    ! #########################################
    ! This part differs from the skeleton,
    ! since I've defined the x array to go
    ! from the initial time until now, but
    ! with different step sizes for the
    ! different eras!
    ! #########################################
    x_init = log(a_init)
    eps    = 1.d-5
    hmin   = 0.d0

    n1          = 200                                   ! Number of grid points during recombination
    n2          = 300                                   ! Number of grid points after recombination

    x_start_rec = -log(1.d0 + 1630.4d0)                 ! x of start of recombination
    x_end_rec   = -log(1.d0 + 614.2d0)                  ! x of end of recombination
    
    dx1         = (x_end_rec - x_start_rec)/(n1+1)      ! Increment of x during rec
    dx2         = (- x_end_rec)/(n2+1)                  ! Increment of x after rec 
    dx          = (x_start_rec - x_init)/(n_before) 
    


    ! Fills the x array with values before recombination
    do i = 0, n_before
      x(i) = x_init + dx*i
    end do

    ! Fills the x array with values during recombination
    do i = 1, n1
      x(n_before + i) = x_start_rec + dx1*(i+1)
    end do

    ! Fills the x array with values after recombination
    do i = 1, n2
      x(n1+n_before+i) = x_end_rec + (i+1)*dx2
    end do
    
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))
    allocate(x_before_tc(1000))
    allocate(x_after(0:n_t))


    ! ################################
    ! Start of integration
    ! ################################

    ! Propagate each k-mode independently
    do k = 1, n_k

       write(*,*) "Starting to Integrate mode ", k, ". With k = ", ks(k)/(H_0/c)

       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-2


       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       
       x_tc =  get_tight_coupling_time(k_current,index_tc)
       write(*,*) "x_tc =", x_tc


       y(1) = delta(0,k)
       y(2) = delta_b(0,k)
       y(3) = v(0,k)
       y(4) = v_b(0,k)
       y(5) = Phi(0,k)
       y(6) = Theta(0,0,k)
       y(7) = Theta(0,1,k)


       write(*,*) "Integrating Tight Decoupling"
       
       ! Integration during tight coupling
       do i = 1, index_tc

       
          
          call odeint(y(1:7), x(i-1), x(i), eps, h1, hmin, derivs_tight_coupling, bsstep, output) 

 
          y(8)   = -(20.d0*c*k_current)/(45.d0*get_H_p(x(i))*get_dtau(x(i)))*y(7)
          do l = 3, lmax_int
              y(6+l) = -l/(2.d0*l + 1.d0)*c*k_current/(get_H_p(x(i))*get_dtau(x(i)))*y(6+l-1)
          end do

          call save_to_arrays(y,x(i),k,i)
  
        end do
        
        write(*,*) "Integrating Rest"
        ! Integration after tight coupling
        do i = index_tc +1, n_t + n_before
        
          call odeint(y, x(i-1), x(i),  eps, h1, hmin, derivs_after_decoupling, bsstep, output) 
          call save_to_arrays(y,x(i),k,i)

        end do

      

    end do


    ! Writes to array

    write(1,*) x
    write(2,*) v
    write(3,*) v_b
    write(4,*) delta
    write(5,*) delta_b
    write(6,*) Phi
    write(7,*) Theta
    write(8,*) Psi
    write(9,*) ks/(H_0/c)


    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dtau) > 0.1 or x > x(start of recombination)


  ! ######################################
  ! Function for finding where tight 
  ! coupling ends. This returns both
  ! the x value and index where it happens
  ! ######################################
  function get_tight_coupling_time(k,index)
    implicit none

    real(dp), intent(in)      :: k
    integer(i4b), intent(out) :: index
    real(dp)                  ::  x_init, x_end, dx, get_tight_coupling_time
    integer(i4b)              :: n, i
  
    x_end   = -log(1.d0 + 1630.4d0)
  
    get_tight_coupling_time = x_end  !Default values
    index = n_before

    do i =0, n_before
      
      if (abs(get_dtau(x(i)))<10) then 
        get_tight_coupling_time = x(i)
        index = n_before
        exit
      
      else if (abs((c*k_current)/(get_H_p(x(i))*get_dtau(x(i)))) > 0.1) then
        get_tight_coupling_time = x(i)
        index = n_before
        exit

      else if (x(i)>x_end) then
        get_tight_coupling_time = x_end
        index = n_before
        exit
      end if

    end do  

  end function get_tight_coupling_time


  ! ################################################
  ! Subroutines for saving to the main arrays. 
  ! Have its own subroutine to tidy up the main loop
  ! ################################################

  subroutine save_to_arrays(y,x,k,i)
    implicit none
    integer, intent(in)                 :: i,k
    real(dp), intent(in)                :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), allocatable, dimension(:) :: dydx
    integer l

    allocate(dydx(npar))

    delta(i,k)   = y(1)
    delta_b(i,k) = y(2)
    v(i,k)       = y(3)
    v_b(i,k)     = y(4)
    Phi(i,k)     = y(5)
    do l = 0, lmax_int
      Theta(i,l,k) = y(6+l)
    end do
    Psi(i,k)     = -y(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*x))*omega_r*y(8) !Temp: May have to use evolving version of Omega_r

    ! Task: Store derivatives that are required for C_l estimation

    call derivs_after_decoupling(x,y,dydx)

    dPhi(i,k)     = dydx(5)
    dv_b(i,k)     = dydx(4)
    dTheta(i,:,k) = dydx(6:lmax_int)
    dPsi(i,k)     = -dydx(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*x))*(dydx(8)+1.d0/((x)*exp(x))*y(8)) !Temp: Find out how to do this


  end subroutine save_to_arrays


  ! ##############################################
  ! Subroutines needed for the integration
  ! ##############################################

  subroutine derivs_tight_coupling(x, y, dydx)
    !use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp)                            :: Psi, R, q, Theta_2, H_p, dtau, ck, ck_H!ddelta, dv, ddelta_b, dv_b, dPhi, Psi, R, dTheta_0, dTheta_1, q, Theta_2, H_p, dtau
    

    ! Calculates the different derivatives
    H_p = get_H_p(x)
    dtau = get_dtau(x)
    ck = c*k_current
    ck_H = ck/H_p


    R           = (4.d0*omega_r)/(3.d0*omega_b*exp(x)) 

    Theta_2     = -(20.d0*ck_H)/(45.d0*dtau)*y(7)

    Psi         = -y(5) - (12.d0*H_0**2.d0)/(ck*ck*exp(2.d0*x))*omega_r*Theta_2 
    dydx(5)        = Psi - (ck*ck)/(3.d0*H_p**2.d0)*y(5) + H_0**2.d0/(2.d0*H_p**2.d0)*(omega_m*exp(-x)*y(1)+omega_b*exp(-x)*y(2)+4.d0*omega_r*exp(-2.d0*x)*y(6))

    dydx(1)      = ck_H*y(3) - 3.d0*dydx(5)
    dydx(2)    = ck_H*y(4) - 3.d0*dydx(5)

    
    
    dydx(6)    = -ck_H*y(7) - dydx(5)
    

    q           = (-(dtau*(1.d0-2.d0*R)+(1.d0 + R)*get_ddtau(x))*(3.d0*y(7)+y(4)) - ck_H*Psi + ck_H*(1-get_dH_p(x)/H_p)*(2.d0*Theta_2-y(6)) - ck_H*dydx(6))/((1.d0+R)*dtau - 1.d0 + get_dH_p(x)/H_p)

    dydx(3)          = -y(3) - ck_H*Psi
    dydx(4)        = (1.d0/(1.d0 + R))*(-y(4) - ck_H*Psi + R*(q + ck_H*(2.d0*Theta_2-y(6)) - ck_H*Psi))

    dydx(7)    = 1.d0/3.d0*(q - dydx(4))


 

    
  end subroutine derivs_tight_coupling

  subroutine derivs_after_decoupling(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    
    real(dp)                            :: Psi, R, H_p, dtau,ck,ck_H !ddelta, dv, ddelta_b, dv_b, dPhi, Psi, R, H_p, dtau
    integer(i4b)                        :: l
    
    
    

    ! Calculates the different derivatives

    H_p = get_H_p(x)
    dtau = get_dtau(x)
    ck = c*k_current
    ck_H = ck/H_p

    R                 = (4.d0*omega_r)/(3.d0*omega_b*exp(x)) !Temp: Will probabily need to use the evolving versions of the omegas
    Psi               = -y(5) - (12.d0*H_0**2.d0)/(ck*ck*exp(2.d0*x))*omega_r*y(8) !Temp: May have to use evolving version of Omega_r
    dydx(5)           = Psi - (ck*ck)/(3.d0*H_p**2.d0)*y(5) + H_0**2.d0/(2.d0*H_p**2.d0)*(omega_m*exp(-x)*y(1)+omega_b*exp(-x)*y(2)+4.d0*omega_r*exp(-2.d0*x)*y(6))

    dydx(2)           = ck_H*y(4) - 3.d0*dydx(5) 
    dydx(1)           = ck_H*y(3) - 3.d0*dydx(5) 

    dydx(3)           = -y(3) - ck_H*Psi
    dydx(4)           = -y(4) - ck_H*Psi + R*dtau*(3.d0*y(7)+y(4))

    dydx(6)           = -ck_H*y(7) - dydx(5) 
    dydx(7)           = (ck_H)/(3.d0)*y(6) - ck_H*(2.d0)/(3.d0)*y(8) +(ck_H)/(3.d0)*Psi + dtau*(y(7) + 1.d0/3.d0 * y(4))

    dydx(8)           = ck_H*(2.d0)/((2.d0*2.d0 + 1.d0))*y(7) - ck_H*((2.d0+1.d0))/((2.d0*2.d0 + 1.d0))*y(9)  + dtau*(y(8) - 1.d0/10.d0*y(8))

    do l = 3, lmax_int-1
      dydx(6+l)       = ck_H*(l)/((2.d0*l + 1.d0))*y(6+l-1) - ck_H*((l+1.d0))/((2.d0*l + 1.d0))*y(6+l+1) + dtau*(y(6+l))
    end do
     
    dydx(6+lmax_int)  = ck_H*y(6+lmax_int-1) - c*(lmax_int+1.d0)/(H_p*get_eta(x))*y(6+lmax_int) + dtau*y(6+lmax_int)



  end subroutine derivs_after_decoupling

end module evolution_mod
