program spectral_region_performance
   use healpix_types
   use region_mod
   implicit none

   character(len=512)   :: paramfile, testfuntion
   integer(i4b) :: unit
   integer(i4b) :: test_iterations
   integer(i4b) :: i, num_components, nside, npix
   real(dp)     :: t1, t2, dum
   integer(i4b), dimension(4)   :: state
   integer(i4b), dimension(:, :), allocatable   :: pixlist_scramble
   real(dp)     :: prior_low, prior_high, par
   real(dp), dimension(10) :: par_array

   unit = getlun()
   call getarg(1, paramfile)
   call get_parameter(unit, paramfile, 'NUM_TEST_ITERATIONS', par_int=test_iterations)
   call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
   npix = 12 * nside ** 2
   call get_parameter(unit, paramfile, 'TEST_FUNCTION', par_string=testfunction)
   call initialize_region_mod(unit, paramfile)
   call get_parameter(unit, paramfile, 'NUM_SPECTRAL_BEHAVIORS', & 
      & par_int=num_components)
   allocate(pixlist_scramble(num_components, 0:npix-1))
   call get_scrambled_pixlist(pixlist_scramble)
   call get_parameter(unit, paramfile, 'PRIOR_LOW_01' // itext, par_dp=prior_low)
   call get_parameter(unit, paramfile, 'PRIOR_HIGH_01' // itext, par_dp=prior_high)
   call get_parameter(unit, paramfile, 'NUM_INIT_REGIONS_01' // itext, par_int=num_init_reg)

   state(2) = 1
   state(3:4) = -1
   par = (prior_high + prior_low) * 0.5
   call fill_par_array(prior_low, prior_high, par_array)
   select case (trim(testfunction))
      case ("singlepix")
         do i = 0, test_iterations - 1
            state(1) = pixlist_scramble(1, i % npix)
            dum = get_single_pixel_chisq_singlepar(par, state)
         end do
      case ("region_singlepar")
         do i = 1, test_iterations
            state(1) = i % num_init_reg
            dum = get_region_like_singlepar(par, state, 0)
         end do
      case ("region_multipar")
         do i = 1, test_iterations
            state(1) = i % num_init_reg
            dum = get_region_like(par_array, state, 0)
         end do
      case default
         stop "Unknown test function."
   end select

contains
   subroutine fill_par_array(low, high, array)
      implicit none
      real(dp) :: low, high, leng, dl
      real(dp), dimension(:)    :: array
      integer(i4b)      :: i
      leng = len(array)
      dl = (high - low) / (leng + 1)
      do i = 0, int(leng)
         array(i) = low + i * dl
      end do
   end subroutine
end program spectral_region_performance
