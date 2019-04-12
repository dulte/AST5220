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
  integer(i4b), parameter          :: n_k      = 100
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

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))

    ! Allocates array for x
    allocate(x(0:n_t))


    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1
    delta(0,:)   = 1.5*Phi(0,0)
    delta_b(0,:) = 1.5*Phi(0,0)
       
    do i = 1, n_k
       v(0,i)       = (c*ks(i))/(get_H_p(log(a_init))*2.d0)*Phi(0,0)
       v_b(0,i)     = (c*ks(i))/(get_H_p(log(a_init))*2.d0)*Phi(0,0)
       Theta(0,0,i) = 0.5*Phi(0,0)
       Theta(0,1,i) = -(c*ks(i))/(get_H_p(log(a_init))*6.d0)*Phi(0,0)
       Theta(0,2,i) = -(20.d0*c*ks(i))/(get_H_p(log(a_init))*45.d0*get_dtau(log(a_init)))*Phi(0,0)
       do l = 3, lmax_int
          Theta(0,l,i) = -l/(2.d0*l + 1)*(c*ks(i))/(get_H_p(log(a_init))*get_dtau(log(a_init)))*Theta(0,l-1,i)
       end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 0.d0

    

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k

       write(*,*) "Starting to Integrate mode ", k

       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5

       ! Initialize equation set for tight coupling7
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
       dt   = (-x_tc)/(n_t-1)

       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       write(*,*) "Integrating over tight coupling"
       call odeint(y_tight_coupling, x_init, x_tc, eps, h1, hmin, derivs_tight_coupling, bsstep, output) 

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling
       y(8)   = -(20.d0*c*k_current)/(45.d0*get_H_p(x_tc)*get_dtau(x_tc))*y_tight_coupling(7)
       do l = 3, lmax_int
          y(6+l) = -l/(2.d0*l + 1.d0)*c*k_current/(get_H_p(x_tc)*get_dtau(x_tc))*y(6+l-1)
       end do

       
       do i = 1, n_t
          ! Task: Integrate equations from tight coupling to today
          write(*,*) "Integrating Rest"
          call odeint(y, x_tc + (i-1)*dt, x_tc + (i)*dt,  eps, h1, hmin, derivs_after_decoupling, bsstep, output) 
          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do
          Psi(i,k)     = -y(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*(x_tc + (i)*dt)))*omega_r*y(8) !Temp: May have to use evolving version of Omega_r

          ! Task: Store derivatives that are required for C_l estimation

          call derivs_after_decoupling(x_tc + (i)*dt,y,dydx)

          dPhi(i,k)     = dydx(5)
          dv_b(i,k)     = dydx(4)
          dTheta(i,:,k) = dydx(6:lmax_int)
          dPsi(i,k)     = -dydx(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*(x_tc + (i)*dt)))*(dydx(8)+1.d0/((x_tc + (i)*dt)*exp(x_tc + (i)*dt))*y(8)) !Temp: Find out how to do this
       end do

    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns


  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time

    get_tight_coupling_time = -7.0 ! Temp: Recomb peaks at around x = -7, where dtau < 10

  end function get_tight_coupling_time


  ! ##############################################
  ! Subroutines needed for the integration
  ! ##############################################

  subroutine derivs_tight_coupling(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp)                            :: ddelta, dv, ddelta_b, dv_b, dPhi, Psi, R, dTheta_0, dTheta_1, q, Theta_2
    integer(i4b)                        :: l
    

    ! Calculates the different derivatives
    R           = (4.d0*omega_r)/(3.d0*omega_b*exp(x)) !Temp: Will probabily need to use the evolving versions of the omegas

    Psi         = -y(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*x))*omega_r*y(8) !Temp: May have to use evolving version of Omega_r
    dPhi        = Psi - (c*c*k_current*k_current)/(3.d0*get_H_p(x)**2.d0)*y(5) + H_0**2.d0/(2.d0*get_H_p(x)**2.d0)*(omega_m*exp(-x)*y(1)+omega_b*exp(-x)*y(2)+4.d0*omega_r*exp(-2.d0*x)*y(6))

    ddelta      = c*k_current/(get_H_p(x))*y(3)
    ddelta_b    = c*k_current/(get_H_p(x))*y(4)

    Theta_2     = -(20.d0*c*k_current)/(45.d0*get_H_p(x)*get_dtau(x))*y(7)
    
    dTheta_0    = -(c*k_current)/(get_H_p(x))*y(7) - dPhi
    

    q           = -((get_dtau(x)*(1.d0*2.d0*R)+(1.d0 + R)*get_ddtau(x))*(3.d0*y(7)+y(4)) - c*k_current/get_H_p(x)*Psi + c*k_current/get_H_p(x)*(1-get_dH_p(x)/get_H_p(x))*(2.d0*Theta_2-y(6)) - c*k_current/get_H_p(x)*dTheta_0)/((1.d0+R)*get_dtau(x) - 1 + get_dH_p(x)/get_H_p(x))
    
    dTheta_1    = 1.d0/3.d0*(q - y(4))

    dv          = -y(3) - c*k_current/(get_H_p(x))*Psi
    dv_b        = 1.d0/(1.d0 + R)*(-y(4) - c*k_current/(get_H_p(x))*Psi + R*(q + c*k_current/get_H_p(x)*(2*Theta_2-y(6)) - c*k_current/(get_H_p(x))*Psi))

    


    ! Places the derivatives into dydx for output
    dydx(1) = ddelta
    dydx(2) = ddelta_b
    dydx(3) = dv
    dydx(4) = dv_b
    dydx(5) = dPhi
    dydx(6) = dTheta_0
    dydx(7) = dTheta_1
  end subroutine derivs_tight_coupling

  subroutine derivs_after_decoupling(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp),allocatable, dimension(:)  :: dTheta
    real(dp)                            :: ddelta, dv, ddelta_b, dv_b, dPhi, Psi, R
    integer(i4b)                        :: l
    
    allocate(dTheta(0:lmax_int))


    ! Trenger ikke aa lage egen funksjon for hver av de deriverte. 
    ! y inneholder alle de forrige verdiene. Eneste tidligere dydx man trenger er Phi', hvor man nok kan bruke den man regnet ut
    ! saa
    

    ! Calculates the different derivatives
    R                 = (4.d0*omega_r)/(3.d0*omega_b*exp(x)) !Temp: Will probabily need to use the evolving versions of the omegas
    Psi               = -y(5) - (12.d0*H_0**2.d0)/(c*c*k_current*k_current*exp(2.d0*x))*omega_r*y(8) !Temp: May have to use evolving version of Omega_r
    dPhi              = Psi - (c*c*k_current*k_current)/(3.d0*get_H_p(x)**2.d0)*y(5) + H_0**2.d0/(2.d0*get_H_p(x)**2.d0)*(omega_m*exp(-x)*y(1)+omega_b*exp(-x)*y(2)+4.d0*omega_r*exp(-2.d0*x)*y(6))

    ddelta            = c*k_current/(get_H_p(x))*y(3)
    ddelta_b          = c*k_current/(get_H_p(x))*y(4)

    dv                = -y(3) - c*k_current/(get_H_p(x))*Psi
    dv_b              = -y(4) - c*k_current/(get_H_p(x))*Psi + R*get_dtau(x)*(3.d0*y(7)+y(4))

    dTheta(0)         = -(c*k_current)/(get_H_p(x))*y(7) - dPhi
    dTheta(1)         = (c*k_current)/(3.d0*get_H_p(x))*y(6) - (2.d0*c*k_current)/(3.d0*get_H_p(x))*y(8) +(c*k_current)/(3.d0*get_H_p(x))*Psi + get_dtau(x)*(y(7) + 1.d0/3.d0 * y(4))

    dTheta(2)         = (2.d0*k_current*c)/((2.d0*2 + 1)*get_H_p(x))*y(7) - ((2+1)*c*k_current)/((2.d0*2 + 1)*get_H_p(x))*y(7)  + get_dtau(x)*(y(8) + 1.d0/10.d0*y(8))

    do l = 3, lmax_int-1
      dTheta(l)       = (l*k_current*c)/((2.d0*l + 1)*get_H_p(x))*y(6+l-1) - ((l+1)*c*k_current)/((2.d0*l + 1)*get_H_p(x))*y(6+l-1) + get_dtau(x)*(y(6+l))
    end do
     
    dTheta(lmax_int)  = (c*k_current)/(get_H_p(x))*y(6+lmax_int-1) - c*(lmax_int+1)/(get_H_p(x)*get_eta(x))*y(6+lmax_int) + get_dtau(x)*y(6+lmax_int)


    ! Places the derivatives into dydx for output
    dydx(1) = ddelta
    dydx(2) = ddelta_b
    dydx(3) = dv
    dydx(4) = dv_b
    dydx(5) = dPhi
    do l = 0, lmax_int
      dydx(6+l) = dTheta(l)
    end do

  end subroutine derivs_after_decoupling

end module evolution_mod
