module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  use rk_mod

  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains


  ! #################################################
  ! The main subroutine, used mainly to integrate eta
  ! #################################################

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, dx1, dx2, dx_eta

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    dx1         = (x_end_rec - x_start_rec)/(n1-1)     ! Increment of x during rec
    dx2         = (x_0 - x_end_rec)/(n2-1)  ! Increment of x after rec    

 
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation
    dx_eta      = (x_eta2-x_eta1)/(n_eta-1) ! Increments of x

    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    
    ! Fills the x_t array with values during recombination
    do i = 1, n1
        x_t(i) = x_start_rec + dx1*(i-1)
    end do

    ! Fills the x_t array with values after recombination
    do i = 1, n2
        x_t(n1+i-1) = x_end_rec + (i-1)*dx2
    end do

    ! Gets a from the definition of x (x = ln a -> a = exp(x))
    a_t = exp(x_t)

    
    
   

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))

    ! Fills the x array
    x_eta = [(x_eta1 + dx_eta*(i-1),i=1,n_eta)]

    ! Sets the initial value of eta given in the report
    eta(1) = c*a_init/(H_0*sqrt(Omega_r))


    ! Solving the differential equation to find eta. 
    ! This used the odesover from Numerical Recipes, with a Runge Kutta (qs) integrator.
    do i = 2, n_eta
      eta(i) = eta(i-1)
      call odeint(eta(i:i), x_eta(i-1), x_eta(i), 1.d-7, 1.d-2, 1.d-10, derivs, rkqs, output) 
    end do

    

    ! Gets the second derivatives used later to spline eta
    call spline(x_eta,eta,1.d30,1.d30,eta2)

  end subroutine initialize_time_mod





  ! ################################################################
  ! Subroutine that defines new arrays with z and x, which are used
  ! to write parameters to file, which are asked for in the exercise
  ! ################################################################

  subroutine write_arrays_to_file
    implicit none

    integer(i4b) :: i, n ! Defines the interation variable, and number of steps
    real(dp)     :: z_0, z_start, x, z, dz, dx, rho_crit
    real(dp),    allocatable, dimension(:) :: z_array, x_array


    ! Gives values to variables
    z_start       = 1.d10                                 ! Red-shift at start of time
    z_0           = 0.d0                                  ! Red-shift at present age
    n             = 5000                                  ! Number of steps
    dz            = (z_0 - z_start)/(n-1)                 ! Increment needed to fill z array
    dx            = (-log(1+z_0) + log(1+z_start))/(n-1)  ! Increment needed to fill x array
    

    ! Allocates and filles array with red-shift values and x
    ! Since there is a log difference in z and x, I have chosen to make
    ! two array. If one dont, one of the grids will be lacking points for certain scales
    allocate(z_array(n))
    allocate(x_array(n))
    z_array = [(z_start + dz*(i-1),i=1,n)]
    x_array = [(-log(1+z_start) + dx*(i-1),i=1,n)]


    ! Opens all the necessary files
    open(1,file='output/eta.dat')
    open(2,file='output/H.dat')
    open(3,file='output/H_z.dat')
    open(4,file='output/Omega_b.dat')
    open(5,file='output/Omega_r.dat')
    open(6,file='output/Omega_m.dat')
    open(7,file='output/Omega_L.dat')
    

    
    
    ! Loops over z and x, and writes the different fuctions to file
    do i=1,n
      z = z_array(i)              
      x = x_array(i)

      ! Writes the conformal time and H to file
      Write(1,*) x, get_eta(x)
      Write(2,*) x, get_H(x)
      Write(3,*) z, get_H(-log(1+z))

      ! Calculates the critical density
      rho_crit = get_critical_density(x)

      ! Calculates and writes density fractions to file
      Write(4,*) x, Omega_b*rho_c*exp(-3*x)/rho_crit
      Write(5,*) x, Omega_r*rho_c*exp(-4*x)/rho_crit
      Write(6,*) x, Omega_m*rho_c*exp(-3*x)/rho_crit
      Write(7,*) x, Omega_lambda*rho_c/rho_crit
    end do
  end subroutine write_arrays_to_file




  ! ###################################################################
  ! Subroutines and functions used to calculate various parameters
  ! ###################################################################

  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x 
    real(dp)             :: get_H

    get_H = H_0*sqrt((Omega_b+Omega_m)*exp(-3*x)+(Omega_r+Omega_nu)*exp(-4*x) + Omega_lambda)

  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p

    get_H_p = exp(x)*get_H(x)

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p

    
    get_dH_p = exp(x)*(get_H(x) + H_0**2/(get_H(x)*2.)*(-3*(Omega_b+Omega_m)*exp(-3*x)-4*(Omega_r+Omega_nu)*exp(-4*x)))

  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    get_eta = splint(x_eta,eta,eta2,x_in)

  end function get_eta


  ! Calculates the critical density for a given x, used in the calculations for the Omegas
  function get_critical_density(x)
    real(dp), intent(in) :: x
    real(16), parameter :: PI_16 = 4 * atan (1.0_16)
    real(dp)             :: get_critical_density

    get_critical_density = (3.*get_H(x)**(2.))/(8.*PI_16*G_grav)


  end function get_critical_density

  ! Subroutine for the derivative used in the odesolver
  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    dydx = c/get_H_p(x)
  end subroutine derivs

  ! Empty subroutine used for the odesolver
  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output


  

end module time_mod
