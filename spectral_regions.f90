program spectral_regions
   use healpix_types
   use region_mod
   use utils
   implicit none
   !Program that samples regions of common spectral indices on the sky given
   !a skymap and noise properties

   integer(i4b), allocatable, dimension(:, :) :: pixel_regions
   real(dp), allocatable, dimension(:, :) :: res_maps
   character(len=512)   :: paramfile
   integer(i4b) :: unit
   integer(i4b) :: mcmc_steps
   integer(i4b) :: i
   real(dp)     :: t1, t2, t1_int, t2_int

   unit = getlun()
   call getarg(1, paramfile)
   call get_parameter(unit, paramfile, 'NUM_MCMC_STEPS', par_int=mcmc_steps)
   call initialize_region_mod(unit, paramfile)

   call cpu_time(t1)
   do i = 1, mcmc_steps
      call cpu_time(t1_int)
      print *, 'mcmc_step : ', i
      call sample_marginalised_regions()
      call output_maps(i)
      call cpu_time(t2_int)
      print *, "Time for this step: ", t2_int - t1_int
   end do
   call cpu_time(t2)
   print *, 'Total time spent:'
   print *, t2-t1

end program spectral_regions
