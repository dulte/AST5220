module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

  real(dp)                                     :: PI_16
  

contains


  ! ###################################################
  ! Init subrouting. Calculates most of the quantities
  ! ###################################################

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop, dx_rec, rhs, n_H
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstop
    dx_rec     = (xstop-xstart)/(n-1) ! Increment of x

    PI_16      = 4 * atan (1.0_16) ! Value of pi with double precition
    
    
    ! Allocation of the Arrays
    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! Task: Fill in x (rec) grid
    do i=1,n
      x_rec(i) = xstart + (i-1)*dx_rec
    end do


    


    ! Task: Compute X_e and n_e at all grid times

    ! Variable holding whether to use saha or peebles.
    use_saha = .true.
    
    do i = 1, n
      ! Used in both Saha and to calculate n_e
      n_H = (Omega_b*rho_c/(m_H*exp(3*x_rec(i))))
       
       if (use_saha) then

          ! Computes the right and side of the equation
          T_b = T_0/exp(x_rec(i))
          rhs =  (1.d0/(n_H)*((m_e*T_b*k_b)/(hbar*hbar*2.d0*PI_16))**(3.d0/2.d0)*exp(-epsilon_0/(T_b*k_b))) !1.d0/(1.d0/(Omega_b*rho_c/(m_H*exp(3.d0*x_rec(i))))*(m_e*(T_0*k_b/x_rec(i))/(2.d0*PI_16*hbar*hbar))**(3.d0/2.d0)*exp(-epsilon_0/(T_0*k_b/x_rec(i))))

          
          ! Citardauq formula
          X_e(i) = -2.d0*rhs/(-rhs-sqrt(rhs*rhs + 4*rhs)) 
          
          !Checks if one should switch to peebles.
          if (X_e(i) < saha_limit) then 
            use_saha = .false.
            Write(*,*) "Now starting to use Peebles' at", x_rec(i)
          end if
       else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), 1.d-5, dx_rec, 1.d-8, derivs_rec, rkqs, output_rec) 
       end if


       ! Gets n_e from X_e
       n_e(i) = X_e(i)*n_H
    end do

    



    ! Task: Compute splined (log of) electron density function
    call spline(x_rec,log(n_e),1.d30,1.d30,n_e2)

    ! Task: Compute optical depth at all grid points

    !Since we only know the end value, we have to integrate in reverse, from present to past
    tau(n) = 0.d0
    do i = n-1, 1, -1
      tau(i) = tau(i+1)
      call odeint(tau(i:i), x_rec(i+1), x_rec(i), 1.d-8, 1.d-5, 1.d-8, derivs_tau, rkqs, output_tau) 
    end do



    ! Task: Compute splined (log of) optical depth
    tau(n) = tau(n-1) ! Since tau(n) = 0, we will have problems with log
    call spline(x_rec,log(tau),1.d30,1.d30,tau2)


    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec,tau2,1.d30,1.d30,tau22)


    ! Task: Compute splined visibility function

    ! Integrates tau to get g
    do i = 1,n
      g(i) = -get_dtau(x_rec(i))*exp(-get_tau(x_rec(i)))
    end do

    call spline(x_rec,g,1.d30,1.d30,g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec,g2,1.d30,1.d30,g22)


  end subroutine initialize_rec_mod


  ! ##########################################################
  ! Function for saving the quantities with higher resolution
  ! ##########################################################

  subroutine write_rec_to_file
    implicit none

    integer(i4b) :: i, n                                ! Defines the interation variable, and number of steps
    real(dp)     :: x_0, x_start, x, dx, n_H            ! Boundary values and quatities needed to calculate things
    real(dp),    allocatable, dimension(:) ::  x_array  ! Array holding x


    ! Gives values to variables
    
    x_start       = -18.d0                                ! x at start of time
    x_0           = 0.d0                                  ! x at present age
    n             = 5000                                  ! Number of steps
    dx            = (x_0 + x_start)/(n-1)                 ! Increment needed to fill x array
    

    ! Allocates and filles array with x
    allocate(x_array(n))
    x_array = [(x_0 + dx*(i-1),i=1,n)]


    ! Opens all the necessary files
    open(1,file='output/X_e.dat')
    open(2,file='output/tau.dat')
    open(3,file='output/dtau.dat')
    open(4,file='output/ddtau.dat')
    open(5,file='output/g.dat')
    open(6,file='output/dg.dat')
    open(7,file='output/ddg.dat')
    

    
    
    ! Loops over x, and writes the different fuctions to file
    do i = 1, n

      x = x_array(i)
      n_H = (Omega_b*rho_c/(m_H*exp(3*x)))

      Write(1,*) x, get_n_e(x)/n_H
      Write(2,*) x, get_tau(x)
      Write(3,*) x, get_dtau(x)
      Write(4,*) x, get_ddtau(x)
      Write(5,*) x, get_g(x)
      Write(6,*) x, get_dg(x)
      Write(7,*) x, get_ddg(x)
    end do
  end subroutine write_rec_to_file

  ! ###################################################
  ! Spline/splint functions to get continous quantities
  ! ###################################################

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: log_n_e
    real(dp)             :: get_n_e

    log_n_e = splint(x_rec,log(n_e),n_e2,x)
    get_n_e = exp(log_n_e)

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: log_tau
    real(dp)             :: get_tau

    log_tau = splint(x_rec,log(tau),tau2,x)
    get_tau = exp(log_tau)

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: log_dtau
    real(dp)             :: get_dtau

    log_dtau = splint_deriv(x_rec,log(tau),tau2,x)
    get_dtau = get_tau(x)*log_dtau

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: log_ddtau
    real(dp)             :: get_ddtau

    log_ddtau = splint(x_rec,tau2,tau22,x)

    ! Formula from getting from (log tau)'' to tau''
    get_ddtau = get_dtau(x)**2.d0/get_tau(x) + get_tau(x)*log_ddtau

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(x_rec,g,g2,x)

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(x_rec,g,g2,x)

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

    get_ddg = splint(x_rec,g2,g22,x)

  end function get_ddg


  ! ################################################
  ! Function for dX_e/dx
  ! ################################################

  ! Functions for getting the terms in the Peebles equation

  function get_dX(x,X_e)
    implicit none

    real(dp), intent(in) :: x
    real(dp), intent(in) :: X_e
    real(dp)             :: get_dX
    real(dp)             :: lambda_2s_1s = 8.227

    real(dp)             :: Cr,lambda_a, n_1s, beta_2, beta, alpha_2, phi_2, T_b, n_H

    ! May be unnecessary since it is defined in init, (and in healpix)
    PI_16      = 4 * atan (1.0_16) ! Value of pi with double precition. 

    T_b = T_0/exp(x)
    n_H = (Omega_b*rho_c/(m_H*exp(3*x)))

    ! Various Functions used for dX_e/dx
    phi_2     = 0.448*log(epsilon_0/(T_b*k_b))
    alpha_2   = 64.d0*PI_16/sqrt(27.d0*PI_16)*(alpha*alpha)/(m_e*m_e)*sqrt(epsilon_0/(T_b*k_b))*phi_2   *hbar*hbar/c
    beta      = alpha_2*(m_e*T_b*k_b/(2.d0*PI_16*hbar*hbar))**(3.d0/2.d0)*exp(-epsilon_0/(T_b*k_b)) 
    beta_2    = alpha_2*(m_e*T_b*k_b/(2.d0*PI_16*hbar*hbar))**(3.d0/2.d0)*exp(-epsilon_0/(4.0*T_b*k_b))
    n_1s      = (1-X_e)*n_H
    lambda_a  = get_H(x)*(3.0*epsilon_0)**(3.d0)/((8.0*PI_16)**(2.d0)*n_1s) *1.d0/(hbar**(3.d0)*c**(3.d0)) 
    Cr        = (lambda_2s_1s + lambda_a)/(lambda_2s_1s + lambda_a + beta_2)


    ! Final Derivative
    get_dX = Cr/get_H(x)*(beta*(1-X_e) - n_H*alpha_2*X_e**(2.d0))



  end function get_dX



    

  ! ################################################
  ! Subroutines for the integration subroutines
  ! ################################################



  ! Subroutine for the derivative used in the odesolver
  subroutine derivs_rec(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    dydx = get_dX(x,y(1))
  end subroutine derivs_rec

  ! Empty subroutine used for the odesolver
  subroutine output_rec(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output_rec
  

   ! Subroutine for the derivative used in the odesolver
  subroutine derivs_tau(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    dydx = -get_n_e(x)*sigma_T*exp(x)*c/(get_H_p(x))
  end subroutine derivs_tau

  ! Empty subroutine used for the odesolver
  subroutine output_tau(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output_tau
  



end module rec_mod
