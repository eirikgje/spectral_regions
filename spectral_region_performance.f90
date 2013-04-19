program spectral_region_performance
   use healpix_types
   use region_mod
   implicit none

   character(len=512)   :: paramfile
   integer(i4b) :: unit
   integer(i4b) :: test_iterations
   integer(i4b) :: i
   real(dp)     :: t1, t2

   unit = getlun()
   call getarg(1, paramfile)
   call get_parameter(unit, paramfile, 'NUM_TEST_ITERATIONS', par_int=test_iterations)
   call initialize_region_mod(unit, paramfile)


end program spectral_region_performance
