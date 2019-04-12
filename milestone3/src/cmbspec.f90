program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  
  implicit none

  ! Initialize time grids
  write(*,*) "Calculating eta"
  call initialize_time_mod

  
  ! Initialize rec_mod
  write(*,*) "Calculating tau"
  call initialize_rec_mod

  ! Initialize evolution mod and integrates
  write(*,*) "Integrating pertubation"
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns

  ! Output to file desired quantities here
  ! call write_arrays_to_file
  ! call write_rec_to_file

  
end program cmbspec
