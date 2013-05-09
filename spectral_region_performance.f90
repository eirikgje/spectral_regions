program spectral_region_performance
   use healpix_types
   use region_mod
   use utils
   implicit none

   character(len=512)   :: paramfile, testfunction
   integer(i4b) :: unit
   integer(i4b) :: test_iterations, num_init_reg
   integer(i4b) :: i, num_components, nside, npix
   real(dp)     :: t1, t2, dum
   integer(i4b), dimension(4)   :: state
   integer(i4b), dimension(:, :), allocatable   :: pixlist_scramble
   real(dp)     :: prior_low, prior_high, par
   real(dp), dimension(10) :: par_array, dum_array

   unit = getlun()
   call getarg(1, paramfile)
   call initialize_region_mod(unit, paramfile)

   call get_parameter(unit, paramfile, 'NUM_TEST_ITERATIONS', par_int=test_iterations)
   call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
   npix = 12 * nside ** 2
   call get_parameter(unit, paramfile, 'TEST_FUNCTION', par_string=testfunction)
   call get_parameter(unit, paramfile, 'NUM_SPECTRAL_BEHAVIORS', & 
      & par_int=num_components)
   allocate(pixlist_scramble(num_components, 0:npix-1))
   call get_scrambled_pixlist(pixlist_scramble)
   call get_parameter(unit, paramfile, 'PRIOR_LOW_01', par_dp=prior_low)
   call get_parameter(unit, paramfile, 'PRIOR_HIGH_01', par_dp=prior_high)
   call get_parameter(unit, paramfile, 'NUM_INIT_REGIONS_01', par_int=num_init_reg)

   state(2) = 1
   state(3:4) = -1
   par = (prior_high + prior_low) * 0.5
   call fill_par_array(prior_low, prior_high, par_array)
   select case (trim(testfunction))
      case ("singlepix")
         call cpu_time(t1)
         do i = 0, test_iterations - 1
            state(1) = pixlist_scramble(1, mod(i, npix))
            dum = get_single_pixel_chisq_singlepar(par, state)
         end do
         call cpu_time(t2)
      case ("region_singlepar")
         call cpu_time(t1)
         do i = 1, test_iterations
            state(1) = mod(i, num_init_reg)
            dum = get_region_like_singlepar(par, state, 0.d0)
         end do
         call cpu_time(t2)
      case ("region_multipar")
         call cpu_time(t1)
         do i = 1, test_iterations
            state(1) = mod(i, num_init_reg)
            dum_array = get_region_like(par_array, state, 0.d0)
         end do
         call cpu_time(t2)
      case default
         stop "Unknown test function."
   end select

   print *, "Total time: ", t2-t1

contains
   subroutine fill_par_array(low, high, array)
      implicit none
      real(dp) :: low, high, leng, dl
      real(dp), dimension(:)    :: array
      integer(i4b)      :: i
      leng = size(array)
      dl = (high - low) / (leng + 1)
      do i = 0, int(leng)
         array(i) = low + i * dl
      end do
   end subroutine
end program spectral_region_performance
