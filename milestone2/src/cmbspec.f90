program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  
  implicit none

  ! Initialize time grids
  call initialize_time_mod

  ! Output to file desired quantities here
  call write_arrays_to_file

  call initialize_rec_mod

  
end program cmbspec
