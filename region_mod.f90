module region_mod
   use healpix_types
   use rngmod
   use stackobject
   use quiet_mapfile_mod
   use optimization_mod
   use utils
   use udgrade_nr
   implicit none

   private

   public initialize_region_mod, sample_marginalised_regions, output_maps, get_single_pixel_chisq_singlepar, get_region_like, get_region_like_singlepar, get_region_like_singlepar_noomp, get_scrambled_pixlist

   interface output_likelihood
      module procedure output_likelihood_nodp, output_likelihood_dp
   end interface output_likelihood

   type region_data
      integer(i4b)      :: num_pix
      real(dp)  :: param
      real(dp)  :: currlnL
      logical(lgt)      :: first_eval
      real(dp), dimension(2)    :: param_bounds
   end type region_data

   integer(i4b), allocatable, dimension(:, :)      :: pixel_curr_region
   integer(i4b), allocatable, dimension(:)              :: num_init_reg
   integer(i4b), allocatable, dimension(:, :)           :: neighbours
   integer(i4b), allocatable, dimension(:)              :: numneighbours
   real(dp), allocatable, dimension(:)                  :: a2t
   real(dp), allocatable, dimension(:)                  :: prior_low, prior_high
   real(dp), allocatable, dimension(:)                  :: freq
   real(dp), allocatable, dimension(:, :)               :: invN, dat_invn
   real(dp), allocatable, dimension(:, :)       :: mixmat_base
   real(dp), allocatable, dimension(:)       :: ref_freq, initdisp, initpar
   integer(i4b), allocatable, dimension(:)   :: behavior
   character(len=512)   :: output_dir

   real(dp)     :: lambda
   type(planck_rng)     :: rng_handle
   integer(i4b)         :: npix, num_components, numband
   integer(i4b)              :: nside
   type(StackT),        allocatable, dimension(:) :: available_region_numbers
   type(region_data), allocatable, dimension(:, :) :: region

   real(dp), parameter          :: K_BOLTZMANN = 1.3806503d-23 !In SI-units
   real(dp), parameter          :: H_PLANCK = 6.626068d-34 !In SI-units
   real(dp), parameter          :: T_CMB = 2.72548 !In SI-units

contains
   subroutine initialize_region_mod(unit, paramfile)
      implicit none

      integer(i4b), intent(in)              :: unit
      character(len=512), intent(in)        :: paramfile

      integer(i4b)              :: i, k, j
      integer(i4b)              :: seed, ordering, curr_region_num, currpix
      integer(i4b)              :: nlist
      character(len=2)          :: itext
      real(dp)                  :: x, radius
      integer(i4b), allocatable, dimension(:)   :: listpix
      character(len=512), allocatable, dimension(:)     :: mapfname, invnfname
      real(dp), allocatable, dimension(:, :)    :: map_temp
      real(dp), dimension(3)                    :: pixvec
      real(dp), allocatable, dimension(:, :)    :: dat 

      call get_parameter(unit, paramfile, 'SEED', par_int=seed)
      call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
      call get_parameter(unit, paramfile, 'NUM_SPECTRAL_BEHAVIORS', & 
         & par_int=num_components)
      call get_parameter(unit, paramfile, 'NUMBAND', par_int=numband)
      call get_parameter(unit, paramfile, 'LAMBDA_PRIOR', par_dp=lambda)
      call get_parameter(unit, paramfile, 'OUTPUT_DIR', par_string=output_dir)
      allocate(freq(numband))
      allocate(mapfname(numband))
      allocate(invnfname(numband))
      do i = 1, numband
         call int2string(i, itext)
         call get_parameter(unit, paramfile, 'BAND_FREQ_' // itext, par_dp=freq(i))
         call get_parameter(unit, paramfile, 'MAP_' // itext,  par_string=mapfname(i))
         call get_parameter(unit, paramfile, 'INVN_' // itext,  par_string=invnfname(i))
      end do
      !Convert to Hz
      freq = freq * 1d9
      allocate(available_region_numbers(num_components))
      allocate(behavior(num_components))
      allocate(num_init_reg(num_components))
      allocate(ref_freq(num_components))
      allocate(prior_low(num_components))
      allocate(prior_high(num_components))
      npix = 12 * nside ** 2
      do i = 1, num_components
         call int2string(i, itext)
         !1: Power-law
         !2: template curve
         !etc.
         call get_parameter(unit, paramfile, 'BEHAVIOR_' // itext, par_int=behavior(i))
         call get_parameter(unit, paramfile, 'NUM_INIT_REGIONS_' // itext, par_int=num_init_reg(i))
         call get_parameter(unit, paramfile, 'REF_FREQ_' // itext, par_dp=ref_freq(i))
         call get_parameter(unit, paramfile, 'PRIOR_LOW_' // itext, par_dp=prior_low(i))
         call get_parameter(unit, paramfile, 'PRIOR_HIGH_' // itext, par_dp=prior_high(i))
         !Initialize stack. Its size will be npix, since that is the maximal
         !number of regions possible. We add all these numbers to the stack
         !after initializing
         available_region_numbers(i) = newstack(npix)
         do j = npix, 1, -1
            call push(j, available_region_numbers(i))
         end do
      end do
      ref_freq = ref_freq * 1d9
      !We can change the size of this later if needed
      allocate(region(maxval(num_init_reg) + 100, num_components))
      region(:, :)%num_pix = 0

      call rand_init(rng_handle, seed)

      allocate(a2t(numband))

      do i = 1, numband
         x = H_PLANCK * freq(i) / (K_BOLTZMANN * T_CMB)
         a2t(i) = (exp(x) - 1) ** 2 / (x ** 2 * exp(x))
      end do

      allocate(dat(numband, 0:npix-1))
      allocate(invN(numband, 0:npix-1))
      allocate(dat_invn(numband, 0:npix-1))
      do i = 1, numband
         call read_map(map_temp, ordering, trim(mapfname(i)))
         if (ordering == 1) then
            call convert_ring2nest(nside, map_temp(:, 1))
         end if
         dat(i, :) = map_temp(:, 1)
         deallocate(map_temp)
         call read_map(map_temp, ordering, trim(invnfname(i)))
         if (ordering == 1) then
            call convert_ring2nest(nside, map_temp(:, 1))
         end if
         invN(i, :) = map_temp(:, 1)
         deallocate(map_temp)
      end do
      dat_invn = dat * invN
      deallocate(dat)

      !Initialize regions
      allocate(pixel_curr_region(0:npix-1, num_components))
      allocate(listpix(0:npix-1))
      !First region covers the whole map. Then we randomly distribute the rest
      !of the regions
      do i = 1, num_components
         !Fetch the next available region number from the stack
         call pop(available_region_numbers(i), curr_region_num)
         pixel_curr_region(:, i) = curr_region_num
         region(curr_region_num, i)%num_pix = npix
         do j = 2, num_init_reg(i)
            call pop(available_region_numbers(i), curr_region_num)
            !Find a random pixel that will be the start of the new region
            !Not sure if this includes all pixels - but rand_uni gives numbers
            !from 0 to 1, 0 and 1 not inclusive, so I think this is right.
            currpix = int(rand_uni(rng_handle) * npix)
            !We make a disc around this pixel.
            call pix2vec_nest(nside, currpix, pixvec)
            !If somehow we don't get any pixels in this region, do it again
            do while (region(curr_region_num, i)%num_pix == 0)
               nlist = 0
               do while (nlist == 0)
                  !pi/2 was chosen arbitrarily, can be bigger
                  radius = rand_uni(rng_handle) * pi / 2
                  call query_disc(nside, pixvec, radius, listpix, nlist, nest=1)
               end do
               do k = 0, nlist - 1
                  if (region(pixel_curr_region(listpix(k), i), i)%num_pix == 1) then
                     !Last pixel - we cannot make this pixel a member of the new
                     !set
                     continue
                  else
                     !Remove a pixel from the old region, add it to the new
                     region(pixel_curr_region(listpix(k), i), i)%num_pix = & 
                        & region(pixel_curr_region(listpix(k), i), i)%num_pix - 1
                     region(curr_region_num, i)%num_pix = & 
                        & region(curr_region_num, i)%num_pix + 1
                     pixel_curr_region(listpix(k), i) = curr_region_num
                  end if
               end do
            end do
         end do
      end do

      !Initialize spectral indices
      allocate(initdisp(num_components))
      allocate(initpar(num_components))
      do i = 1, num_components
         !Again, if we have variable numbers of regions, this must be changed
         call int2string(i, itext)
         call get_parameter(unit, paramfile, 'INIT_SPECTRAL_INDEX_' // itext, par_dp=initpar(i))
         call get_parameter(unit, paramfile, 'INIT_SPECTRAL_INDEX_DISPERSION_' // itext, par_dp=initdisp(i))
         if (behavior(i) == 1) then
            do j = 1, num_init_reg(i)
               region(j, i)%param = initpar(i) + rand_gauss(rng_handle) * 0.01 * initdisp(i)
               region(j, i)%param_bounds(1) = initpar(i) - initdisp(i)
               region(j, i)%param_bounds(2) = initpar(i) + initdisp(i)
            end do
            !Add more behaviors here as needed!
         end if
      end do

      !Find neighbours (for use later)
      allocate(neighbours(8, 0:npix-1))
      allocate(numneighbours(0:npix-1))
      do i = 0, npix-1
         call neighbours_nest(nside, i, neighbours(:, i), numneighbours(i))
      end do

      allocate(mixmat_base(numband, num_components))
      !Initialize the 'base' of the mixing matrix - i.e. the part that does not
      !change with parameters
      do i = 1, num_components
         do j = 1, numband
            !Fill in with specific behaviors
            if (behavior(i) == 1) then
               !Power-law
               mixmat_base(j, i) = freq(j) / ref_freq(i)
            end if
         end do
      end do

   end subroutine initialize_region_mod

   subroutine sample_marginalised_regions
      implicit none

      integer(i4b), dimension(num_components, 0:npix-1)   :: pixlist_scramble
      integer(i4b), dimension(8)        :: neighbour_region
      integer(i4b)      :: i, j, k, l, pix
      integer(i4b)      :: nneighbour_region, reg_prop

      do i = 1, num_components
         region(:, i)%first_eval = .True.
      end do

      !During burn-in, the earlier pixels have a better chance of forming new
      !regions than the latter ones due to the region number prior. We therefore
      !scramble the list of pixels, and we do the inner loop over components
      call get_scrambled_pixlist(pixlist_scramble)
      do j = 0, npix-1
!         print *, 'count', j
         do i = 1, num_components
            pix = pixlist_scramble(i, j)
!            print *, 'curr_pixel', pix
            neighbour_region = -1
            nneighbour_region = 0
            !Find all neighbouring regions
            do k = 1, 8
               if (pixel_curr_region(neighbours(k, pix), i) /= & 
                  & pixel_curr_region(pix, i)) then
                  if (nneighbour_region == 0) then
                     nneighbour_region = nneighbour_region + 1
                     neighbour_region(nneighbour_region) = & 
                        & pixel_curr_region(neighbours(k, pix), i)
                  else
                     if (.not. any(pixel_curr_region(neighbours(k, pix), i) & 
                        & == neighbour_region(:))) then
                        !Means that the pixel is bordering a neighbouring
                        !region, and that this region is distinct from other
                        !neighbouring regions already added to the list
                        nneighbour_region = nneighbour_region + 1
                        neighbour_region(nneighbour_region) = &
                           & pixel_curr_region(neighbours(k, pix), i)
                     end if
                  end if
               end if
            end do
            !Now we choose whether to make current pixel a new region or propose
            !a switch to a neighbouring region. Currently, we allow every pixel
            !that does not border another region to form its own region,
            !allowing the prior on the number of regions to regulate this
            !entirely.
            if (nneighbour_region == 0) then
               call form_new_region(pix, i)
            else
               if (rand_uni(rng_handle) < 0.5) then
                  call form_new_region(pix, i)
               else
                  if (nneighbour_region > 1) then
                     reg_prop = int(rand_uni(rng_handle) * nneighbour_region) + 1
                  else
                     reg_prop = 1
                  end if
                  reg_prop = neighbour_region(reg_prop)
                  call switch_region(pix, i, reg_prop)
               end if
            end if
         end do
      end do
   end subroutine sample_marginalised_regions

   function gaussian(x, aux, var_dp)
      implicit none
      real(dp), dimension(:), intent(in) :: x
      integer(i4b), dimension(:), intent(in)        :: aux
      real(dp), intent(in)  :: var_dp
      real(dp), dimension(size(x))  :: gaussian

      gaussian = 1/ sqrt(2.d0 * pi) * exp(-0.5d0 * x ** 2)
   end function gaussian

   function gaussian_singlepar(x, aux, var_dp)
      implicit none
      real(dp),               intent(in) :: x
      integer(i4b), dimension(:), intent(in)        :: aux
      real(dp), intent(in)  :: var_dp
      real(dp)                      :: gaussian_singlepar

      gaussian_singlepar = 1/ sqrt(2.d0 * pi) * exp(-0.5d0 * x ** 2)
   end function gaussian_singlepar

   subroutine form_new_region(pixel, component)
      implicit none

      integer(i4b), intent(in)  :: pixel, component

      real(dp)          :: lnLcurr, lnLprop, lnLprop_new, lnLprop_curr
      real(dp)          :: currparam_maxlike, propparam_maxlike
      real(dp)          :: param_maxlike_new
      real(dp), dimension(2)    :: currparam_limits, propparam_limits
      real(dp), dimension(2)    :: limits_new
      logical(lgt)      :: accept
      integer(i4b)      :: reg_curr, reg_prop, i, j
      type(region_data), allocatable, dimension(:, :)      :: region_temp

!      print *, 'formnew'
      reg_curr = pixel_curr_region(pixel, component)
!      print *, 'reg_curr', reg_curr
      !If last pixel in region, do nothing
      if (region(reg_curr, component)%num_pix == 1) then 
         if (count(pixel_curr_region(:, component) == reg_curr) > 1) then
            print *, 'count - something is wrong', count(pixel_curr_region(:, component) == reg_curr)
         end if
!         print *, 'lastpix'
         return
      end if

      if (region(reg_curr, component)%first_eval) then
         lnLcurr = get_region_marginalised_likelihood(reg_curr, component, & 
            & maxlikepoint=currparam_maxlike, limits=currparam_limits)
      else
         lnLcurr = region(reg_curr, component)%currlnL
         currparam_maxlike = region(reg_curr, component)%param
         currparam_limits = region(reg_curr, component)%param_bounds
      end if
      lnLprop_new = get_pixel_marginalised_likelihood(pixel, component, &
         & maxlikepoint=param_maxlike_new, limits=limits_new)
      lnLprop_curr = & 
         & get_region_marginalised_likelihood(reg_curr, component, & 
         & subtract_pixel=pixel, maxlikepoint=propparam_maxlike, &
         & limits=propparam_limits)
      lnLprop = lnLprop_new + lnLprop_curr
      !The lambda term is due to the fact that we now have one more region and
      !we have a prior on the number of regions
      accept = exp(lnLprop - lnLcurr - lambda) > rand_uni(rng_handle)
      if (accept) then
!         print *, 'accept'
         region(reg_curr, component)%num_pix = & 
            & region(reg_curr, component)%num_pix - 1
         !Update bestfit-parameter and loglike
         region(reg_curr, component)%param = propparam_maxlike
         region(reg_curr, component)%currlnL = lnLprop_curr
         region(reg_curr, component)%param_bounds = propparam_limits
         !All the regions affected by this pixel should be reevaluated now
         do i = 1, num_components
            if (i /= component) then
               region(pixel_curr_region(pixel, i), i)%first_eval = .True.
            end if
         end do
         !Draw the next available region number from the stack
         call pop(available_region_numbers(component), & 
            & pixel_curr_region(pixel, component))
         reg_prop = pixel_curr_region(pixel, component)
         if (size(region, 1) < reg_prop) then
            allocate(region_temp(size(region, 1), num_components))
            do j = 1, num_components
               do i = 1, size(region, 1)
                  region_temp(i, j)%param = region(i, j)%param
                  region_temp(i, j)%param_bounds(:) = region(i, j)%param_bounds
                  region_temp(i, j)%num_pix = region(i, j)%num_pix
                  region_temp(i, j)%currlnL = region(i, j)%currlnL
                  region_temp(i, j)%first_eval = region(i, j)%first_eval
               end do
            end do
            deallocate(region)
            allocate(region(reg_prop + 100, num_components))
            do j = 1, num_components
               do i = 1, size(region_temp, 1)
                  region(i, j)%param = region_temp(i, j)%param
                  region(i, j)%param_bounds(:) = region_temp(i, j)%param_bounds
                  region(i, j)%num_pix = region_temp(i, j)%num_pix
                  region(i, j)%currlnL = region_temp(i, j)%currlnL
                  region(i, j)%first_eval = region_temp(i, j)%first_eval
               end do
            end do
            deallocate(region_temp)
         end if
         region(reg_prop, component)%num_pix = 1
         region(reg_prop, component)%param = param_maxlike_new
         region(reg_prop, component)%currlnL = lnLprop_new
         region(reg_prop, component)%param_bounds = limits_new
      else
!         print *, 'reject'
         !Just update the best-fit parameter and log-likelihood, even though
         !they shouldn't have changed much
         region(reg_curr, component)%param = currparam_maxlike
         region(reg_curr, component)%currlnL = lnLcurr
         region(reg_curr, component)%param_bounds = currparam_limits
      end if

   end subroutine form_new_region

   subroutine switch_region(pixel, component, reg_prop)
      implicit none

      integer(i4b), intent(in)  :: pixel, component, reg_prop

      real(dp)  :: lnLcurr, lnLprop, lnLcurr_curregion, lnLcurr_propregion
      real(dp)  :: lnLprop_curregion, lnLprop_propregion
      real(dp)  :: currparam_curregion, currparam_propregion
      real(dp)  :: propparam_curregion, propparam_propregion
      integer(i4b)      :: reg_curr, i
      real(dp), dimension(2)    :: currlimits_curregion, currlimits_propregion
      real(dp), dimension(2)    :: proplimits_curregion, proplimits_propregion
      logical(lgt)      :: accept
      real(dp)  :: parlow, parhigh

      !This looks a bit excessive, but in order to store the log-likelihood for
      !each region we need to save each part of the likelihood
!      print *, 'switch'
      reg_curr = pixel_curr_region(pixel, component)
      if (region(reg_curr, component)%first_eval) then
         lnLcurr_curregion = & 
            & get_region_marginalised_likelihood(reg_curr, component, & 
            & maxlikepoint=currparam_curregion, limits=currlimits_curregion)
      else
         lnLcurr_curregion = region(reg_curr, component)%currlnL
         currparam_curregion = region(reg_curr, component)%param
         currlimits_curregion = region(reg_curr, component)%param_bounds
      end if
      if (region(reg_prop, component)%first_eval) then
         lnLcurr_propregion = & 
            & get_region_marginalised_likelihood(reg_prop, component, & 
            & maxlikepoint=currparam_propregion, limits=currlimits_propregion)
      else
         lnLcurr_propregion = region(reg_prop, component)%currlnL
         currparam_propregion = region(reg_prop, component)%param
         currlimits_propregion = region(reg_prop, component)%param_bounds
      end if

      lnLcurr = lnLcurr_curregion + lnLcurr_propregion

      lnLprop_curregion = & 
         & get_region_marginalised_likelihood(reg_curr, component, & 
         & subtract_pixel=pixel, maxlikepoint=propparam_curregion, &
         & limits=proplimits_curregion)
      lnLprop_propregion = &
         & get_region_marginalised_likelihood(reg_prop, component, & 
         & add_pixel=pixel, maxlikepoint=propparam_propregion, &
         & limits=proplimits_propregion)
      lnLprop = lnLprop_curregion + lnLprop_propregion

      if (region(reg_curr, component)%num_pix == 1) then
         !We're potentially losing a region in this case, so we add the lambda
         !term
         accept = exp(lnLprop - lnLcurr + lambda) > rand_uni(rng_handle)
         if (accept) then
!            print *, 'accept'
            !Add the region number of the region that is about to disappear to
            !the stack of available region numbers
            call push(reg_curr, available_region_numbers(component))
            pixel_curr_region(pixel, component) = reg_prop
            region(reg_prop, component)%param = propparam_propregion
            region(reg_prop, component)%num_pix = & 
               & region(reg_prop, component)%num_pix + 1
            region(reg_prop, component)%currlnL = lnLprop_propregion
            !All the regions affected by this pixel should be reevaluated now
            do i = 1, num_components
               if (i /= component) then
                  region(pixel_curr_region(pixel, i), i)%first_eval = .True.
               end if
            end do
         else
!            print *, 'reject'
            region(reg_curr, component)%param = currparam_curregion
            region(reg_curr, component)%currlnL = lnLcurr_curregion
            region(reg_curr, component)%param_bounds = currlimits_curregion
            region(reg_prop, component)%param = currparam_propregion
            region(reg_prop, component)%currlnL = lnLcurr_propregion
            region(reg_prop, component)%param_bounds = currlimits_propregion
         end if
      else
         !No change in the number of regions
         accept = exp(lnLprop - lnLcurr)
         if (accept) then
!            print *, 'accept'
            region(reg_curr, component)%num_pix = & 
               & region(reg_curr, component)%num_pix - 1
            region(reg_prop, component)%num_pix = & 
               & region(reg_prop, component)%num_pix + 1
            pixel_curr_region(pixel, component) = reg_prop
            region(reg_curr, component)%param = propparam_curregion
            region(reg_curr, component)%currlnL = lnLprop_curregion
            region(reg_curr, component)%param_bounds = proplimits_curregion
            region(reg_prop, component)%param = propparam_propregion
            region(reg_prop, component)%currlnL = lnLprop_propregion
            region(reg_prop, component)%param_bounds = proplimits_propregion
         else
!            print *, 'reject'
            region(reg_curr, component)%param = currparam_curregion
            region(reg_curr, component)%currlnL = lnLcurr_curregion
            region(reg_curr, component)%param_bounds = currlimits_curregion
            region(reg_prop, component)%param = currparam_propregion
            region(reg_prop, component)%currlnL = lnLcurr_propregion
            region(reg_prop, component)%param_bounds = currlimits_propregion
         end if
      end if
   end subroutine switch_region

   function get_pixel_marginalised_likelihood(pixel, component, maxlikepoint, &
         & limits)
      implicit none
      real(dp)  :: get_pixel_marginalised_likelihood
      integer(i4b), intent(in)  :: pixel, component
      real(dp), intent(out)     :: maxlikepoint
      real(dp), dimension(:), intent(out)       :: limits

      real(dp)  :: parmin, chisqmin, dlnL, lnL0, delta, lnL, int_low, int_high
      real(dp)  :: res, parlow, parhigh
      integer(i4b), dimension(2)    :: pixel_lnL_state
      integer(i4b)      :: i
      logical(lgt)      :: err

      !We adjust the pixel_lnL_state array so that the
      !pixel-likelihood-per-parameter-routine needs only the current value of
      !the parameter as input
      pixel_lnL_state(1) = pixel
      pixel_lnL_state(2) = component

      !When calling this routine, it means we are starting a new region. The
      !parameter boundaries must thus be the ones set initially by user.

      call get_maxpoint_and_qrombint(1, pixel_lnL_state, maxlikepoint, res)

      get_pixel_marginalised_likelihood = res

   end function get_pixel_marginalised_likelihood

   function get_region_marginalised_likelihood(region_num, component, & 
         & subtract_pixel, add_pixel, maxlikepoint, limits)
      implicit none
      integer(i4b), intent(in)  :: region_num, component
      integer(i4b), intent(in), optional        :: subtract_pixel, add_pixel
      real(dp), intent(out), optional     :: maxlikepoint
      real(dp), intent(out), dimension(:), optional       :: limits
      real(dp)  :: get_region_marginalised_likelihood

      integer(i4b)      :: pix_sub, pix_add, npix_dum
      real(dp)  :: parlow, parhigh, parmin, chisqmin
      real(dp)  :: dlnL, lnL0, delta, lnL, int_low, int_high, res
      integer(i4b), dimension(4)      :: region_lnL_state
      logical(lgt)      :: err

      if (present(subtract_pixel)) then
         pix_sub = subtract_pixel
      else
         pix_sub = -1
      end if
      if (present(add_pixel)) then
         pix_add = add_pixel
      else
         pix_add = -1
      end if

      !Contains necessary information for the likelihood routine
      region_lnL_state(1) = region_num
      region_lnL_state(2) = component
      region_lnL_state(3) = pix_sub
      region_lnL_state(4) = pix_add

      !Do preliminary checks to see if we're evaluating an empty region
      if (region(region_num, component)%num_pix == 1 .and. pix_sub /= -1) then
         !This should mean that we want to evaluate the region with the only
         !remaining pixel subtracted from it. 

         !Do a test just to check whether this is actually true.
         npix_dum = count(pixel_curr_region(:, component) == region_num)
         if (npix_dum > 1) then
            !Something is wrong
            stop('Error in get_region_marginalised_likelihood:npix_dum>1')
         end if
         !Also check whether the last pixel actually is the same as the
         !subtracted one
         if (pixel_curr_region(pix_sub, component) /= region_num) then
            stop("Error in get_region_marginalised_likelihood:region numbers are not corresponding")
         end if
         !Everything should be fine, but we don't need to evaluate the
         !likelihood because the region is empty
         get_region_marginalised_likelihood = 0
         return
      end if

      call get_maxpoint_and_qrombint(2, region_lnL_state, maxlikepoint, res)

      get_region_marginalised_likelihood = res

   end function get_region_marginalised_likelihood

   function get_region_like(par, region_state, offset)
      implicit none
      real(dp), dimension(:), intent(in)      :: par
      real(dp), dimension(size(par))         :: get_region_like
      integer(i4b), dimension(:), intent(in) :: region_state
      real(dp), intent(in)      :: offset

      real(dp), dimension(size(par))    :: get_region_like_private
!      real(dp), dimension(size(par))  :: temp
      integer(i4b), dimension(2)        :: pix_state
      integer(i4b)      :: region_num, component, pix_sub, pix_add
      integer(i4b)      :: i, j

      region_num = region_state(1)
      component = region_state(2)
      pix_sub = region_state(3)
      pix_add = region_state(4)


      get_region_like = 0


      !$OMP PARALLEL PRIVATE(i, pix_state, get_region_like_private)
      get_region_like_private = 0
      pix_state(2) = component

      !$OMP DO SCHEDULE(GUIDED)
      do i = 0, npix - 1
         if ((pixel_curr_region(i, component) == region_num .and. .not. & 
            & i == pix_sub) .or. i == pix_add) then
            pix_state(1) = i
            get_region_like_private = get_region_like_private + get_single_pixel_chisq(par, pix_state)
         end if
      end do
      !$OMP END DO

      !$OMP CRITICAL
      get_region_like = get_region_like + get_region_like_private

      !$OMP END CRITICAL

      !$OMP END PARALLEL

      get_region_like = exp(-0.5d0 * get_region_like - offset)
!      print *, 'offset', offset
!      print *, 'reglike', get_region_like

   end function get_region_like

   function get_region_like_singlepar(par, region_state, offset)
      implicit none
      real(dp)          :: get_region_like_singlepar
      real(dp), intent(in)      :: par
      integer(i4b), dimension(:), intent(in) :: region_state
      real(dp), intent(in)      :: offset

      real(dp)  :: temp
      real(dp)  :: get_region_like_singlepar_private
      integer(i4b), dimension(2)        :: pix_state
      integer(i4b)      :: region_num, component, pix_sub, pix_add
      integer(i4b)      :: i, j

      region_num = region_state(1)
      component = region_state(2)
      pix_sub = region_state(3)
      pix_add = region_state(4)

      get_region_like_singlepar = 0

!      print *, 'singlepar'

      !$OMP PARALLEL PRIVATE(i, pix_state, get_region_like_singlepar_private)
      pix_state(2) = component
      get_region_like_singlepar_private = 0

      !$OMP DO SCHEDULE(GUIDED)
      do i = 0, npix - 1
         if ((pixel_curr_region(i, component) == region_num .and. .not. & 
            & i == pix_sub) .or. i == pix_add) then
            pix_state(1) = i
            get_region_like_singlepar_private = get_region_like_singlepar_private + get_single_pixel_chisq_singlepar(par, pix_state)
         end if
      end do
      !$OMP END DO

      !$OMP CRITICAL
      get_region_like_singlepar = get_region_like_singlepar + get_region_like_singlepar_private
      !$OMP END CRITICAL


      !$OMP END PARALLEL

      get_region_like_singlepar = exp(-0.5d0 * get_region_like_singlepar - offset)

   end function get_region_like_singlepar

   function get_region_like_singlepar_noomp(par, region_state, offset)
      implicit none
      real(dp)          :: get_region_like_singlepar_noomp
      real(dp), intent(in)      :: par
      integer(i4b), dimension(:), intent(in) :: region_state
      real(dp), intent(in)      :: offset

      real(dp)  :: temp
      integer(i4b), dimension(2)        :: pix_state
      integer(i4b)      :: region_num, component, pix_sub, pix_add
      integer(i4b)      :: i, j

      region_num = region_state(1)
      component = region_state(2)
      pix_sub = region_state(3)
      pix_add = region_state(4)

      get_region_like_singlepar_noomp = 0

      pix_state(2) = component

      do i = 0, npix - 1
         if ((pixel_curr_region(i, component) == region_num .and. .not. & 
            & i == pix_sub) .or. i == pix_add) then
            pix_state(1) = i
            get_region_like_singlepar_noomp = get_region_like_singlepar_noomp + get_single_pixel_chisq_singlepar(par, pix_state)
         end if
      end do

      get_region_like_singlepar_noomp = exp(-0.5d0 * get_region_like_singlepar_noomp - offset)

   end function get_region_like_singlepar_noomp


   function get_region_chisq(par, region_state)
      implicit none
      real(dp), dimension(:), intent(in)      :: par
      real(dp), dimension(size(par))         :: get_region_chisq
      integer(i4b), dimension(:), intent(in) :: region_state

      real(dp), dimension(size(par))    :: get_region_chisq_private
      integer(i4b), dimension(2)        :: pix_state
      integer(i4b)      :: region_num, component, pix_sub, pix_add
      integer(i4b)      :: i, j

      region_num = region_state(1)
      component = region_state(2)
      pix_sub = region_state(3)
      pix_add = region_state(4)

      get_region_chisq = 0

      !$OMP PARALLEL PRIVATE(i, pix_state, get_region_chisq_private)
      pix_state(2) = component
      get_region_chisq_private = 0

      !$OMP DO SCHEDULE(GUIDED)
      do i = 0, npix - 1
         if (pixel_curr_region(i, component) == region_num .and. .not. & 
            & i == pix_sub) then
            pix_state(1) = i
            get_region_chisq_private = get_region_chisq_private + get_single_pixel_chisq(par, pix_state)
         else if (i == pix_add) then
            pix_state(1) = i
            get_region_chisq_private = get_region_chisq_private + get_single_pixel_chisq(par, pix_state)
         end if
      end do
      !$OMP END DO

      !$OMP CRITICAL

      get_region_chisq = get_region_chisq + get_region_chisq_private
      !$OMP END CRITICAL

      !$OMP END PARALLEL
   end function get_region_chisq

   function get_region_chisq_singlepar(par, region_state)
      implicit none
      real(dp)          :: get_region_chisq_singlepar
      real(dp), intent(in)      :: par
      integer(i4b), dimension(:), intent(in) :: region_state

      real(dp)          :: get_region_chisq_singlepar_private
      integer(i4b), dimension(2)        :: pix_state
      real(dp)          :: chisq_tot
      integer(i4b)      :: region_num, component, pix_sub, pix_add
      integer(i4b)      :: i, j

      region_num = region_state(1)
      component = region_state(2)
      pix_sub = region_state(3)
      pix_add = region_state(4)

      get_region_chisq_singlepar = 0

      !$OMP PARALLEL PRIVATE(i, pix_state, get_region_chisq_singlepar_private)
      pix_state(2) = component
      get_region_chisq_singlepar_private = 0

      !$OMP DO SCHEDULE(GUIDED)
      do i = 0, npix - 1
         if (pixel_curr_region(i, component) == region_num .and. .not. & 
            & i == pix_sub) then
            pix_state(1) = i
            get_region_chisq_singlepar_private = get_region_chisq_singlepar_private + get_single_pixel_chisq_singlepar(par, pix_state)
         else if (i == pix_add) then
            pix_state(1) = i
            get_region_chisq_singlepar_private = get_region_chisq_singlepar_private + get_single_pixel_chisq_singlepar(par, pix_state)
         end if
      end do
      !$OMP END DO
      
      !$OMP CRITICAL
      get_region_chisq_singlepar = get_region_chisq_singlepar + get_region_chisq_singlepar_private 
      !$OMP END CRITICAL

      !$OMP END PARALLEL

   end function get_region_chisq_singlepar

   function get_single_pixel_like(par, pix_state, offset)
      implicit none
      real(dp), intent(in), dimension(:)        :: par
      real(dp), dimension(size(par))            :: get_single_pixel_like
      integer(i4b), dimension(:), intent(in)    :: pix_state
      real(dp), intent(in)      :: offset

      get_single_pixel_like = exp(-0.5 * get_single_pixel_chisq(par, pix_state) - offset)

   end function get_single_pixel_like

   function get_single_pixel_like_singlepar(par, pix_state, offset)
      implicit none
      real(dp), intent(in)                      :: par
      real(dp)  :: get_single_pixel_like_singlepar
      integer(i4b), dimension(:), intent(in)    :: pix_state
      real(dp), intent(in)      :: offset

      get_single_pixel_like_singlepar = exp(-0.5 * get_single_pixel_chisq_singlepar(par, pix_state) - offset)

   end function get_single_pixel_like_singlepar

   function get_single_pixel_chisq(par, pix_state)
      !It's not really chi-squared since this is not really something that is
      !gaussian. But it is the corresponding quantity.
      implicit none
      real(dp), intent(in), dimension(:)      :: par
      real(dp), dimension(size(par))          :: get_single_pixel_chisq
      integer(i4b), dimension(:), intent(in)      :: pix_state

      integer(i4b)      :: i, j, k
      integer(i4b)      :: pixnum, component
      real(dp), dimension(numband, num_components)      :: mixmat, mixmat_invn
      real(dp), dimension(num_components, num_components)       :: M, Minv, MinvLU
      real(dp)      :: lndet
      integer(i4b), dimension(num_components)   :: indx
      real(dp), dimension(num_components)       :: x, b, diag

      
      pixnum = pix_state(1)
      component = pix_state(2)
      where(par < prior_low(component) .or. par > prior_high(component))
         get_single_pixel_chisq = 1.d30
      end where
      do k = 1, size(par)
         call calc_mixing_matrix(par(k), pixnum, component, mixmat)
         do i = 1, numband
            mixmat_invn(i, :) = mixmat(i, :) * invN(i, pixnum)
         end do
         Minv = matmul(transpose(mixmat_invn), mixmat)

         MinvLU = Minv
         call choldc(MinvLU, diag)
         if (all(diag == 0)) then
             get_single_pixel_chisq(k) = 1.d30
         end if
         lndet = 0.d0
         do j = 1, num_components
            lndet = lndet + log(diag(j))
         end do
         lndet = 2.d0*lndet
         x = matmul(transpose(mixmat), dat_invn(:, pixnum))
         b = x
      !Solve Minvb = x for b
         call cholsl(MinvLU, diag, x, b)
         get_single_pixel_chisq(k) = lndet - dot_product(x, b)
      end do

   end function get_single_pixel_chisq

   function get_single_pixel_chisq_singlepar(par, pix_state)
      !It's not really chi-squared since this is not really something that is
      !gaussian. But it is the corresponding quantity.
      implicit none
      real(dp)          :: get_single_pixel_chisq_singlepar
      real(dp), intent(in)      :: par
      integer(i4b), dimension(:), intent(in)      :: pix_state

      integer(i4b)      :: i, j
      integer(i4b)      :: pixnum, component
      real(dp), dimension(numband, num_components)      :: mixmat, mixmat_invn
      real(dp), dimension(num_components, num_components)       :: M, Minv, MinvLU
      real(dp)      :: det
      integer(i4b), dimension(num_components)   :: indx
      real(dp), dimension(num_components)       :: x, b
      real(dp), dimension(num_components)       :: diag


      pixnum = pix_state(1)
      component = pix_state(2)

      if (par < prior_low(component) .or. par > prior_high(component)) then
         get_single_pixel_chisq_singlepar = 1.d30
         return
      end if
      call calc_mixing_matrix(par, pixnum, component, mixmat)
      do i = 1, numband
         mixmat_invn(i, :) = mixmat(i, :) * invN(i, pixnum)
      end do
      Minv = matmul(transpose(mixmat_invn), mixmat)
      MinvLU = Minv
      call choldc(MinvLU, diag)
      if (all(diag == 0)) then
         get_single_pixel_chisq_singlepar = 1.d30
         return
      end if
      det = 1
      do j = 1, num_components
         det = det * diag(j) ** 2
      end do
      x = matmul(transpose(mixmat), dat_invn(:, pixnum))
      b = x
      !Solve Minvb = x for b
      call cholsl(MinvLU, diag, x, b)
      get_single_pixel_chisq_singlepar = log(det) - dot_product(x, b)

   end function get_single_pixel_chisq_singlepar

   subroutine calc_mixing_matrix(par, pixel, component, mixmat)
      implicit none
      real(dp), intent(in)  :: par
      integer(i4b), intent(in)  :: pixel, component
      real(dp), intent(out), dimension(:, :)    :: mixmat

      integer(i4b)      :: i, j
      type(region_data)      :: curregion

      do i = 1, num_components
         if (i /= component) then
            curregion = region(pixel_curr_region(pixel, i), i)
         end if
         do j = 1, numband
            !Fill in with specific behaviors
            if (behavior(i) == 1) then
               !Power-law
               if (i == component) then
                  !Use the parameter provided by the caller
                  mixmat(j, i) = a2t(j) * mixmat_base(j, i) ** par
               else
                  !Use the best-fit parameter corresponding to that component
                  mixmat(j, i) = a2t(j) * mixmat_base(j, i) ** curregion%param
               end if
            end if
         end do
      end do

   end subroutine calc_mixing_matrix

   subroutine get_scrambled_pixlist(pixlist)
      implicit none
      integer(i4b), dimension(:, 0:), intent(inout)   :: pixlist

      integer(i4b)      :: i, j, temp, currind, currnumpix

      do j = 0, npix-1
         pixlist(:, j) = j
      end do

      currnumpix = npix - 1
      do while (currnumpix >= 0)
         do i = 1, num_components
            currind = int(rand_uni(rng_handle) * (currnumpix + 1))
            temp = pixlist(i, currind)
            pixlist(i, currind) = pixlist(i, currnumpix)
            pixlist(i, currnumpix) = temp
         end do
         currnumpix = currnumpix - 1
      end do

   end subroutine get_scrambled_pixlist

   subroutine make_result_map(res_map, comp_num)
      implicit none

      real(dp), allocatable, dimension(:, :), intent(inout) :: res_map
      integer(i4b), intent(in)  :: comp_num
      integer(i4b)              :: i, j, k

      if (behavior(comp_num) == 1) then
         !Only one parameter to output
         allocate(res_map(0:npix - 1, 2))
         do i = 0, npix-1
            res_map(i, 1) = region(pixel_curr_region(i, comp_num), comp_num)%param
            res_map(i, 2) = pixel_curr_region(i, comp_num)
         end do
      end if
   end subroutine make_result_map

   subroutine output_maps(iteration)
      implicit none

      integer(i4b), intent(in)  :: iteration

      integer(i4b)              :: i
      character(len=2)          :: itext
      character(len=5)          :: ittext
      character(len=512)        :: fname
      real(dp), allocatable, dimension(:, :)    :: res_map


      do i = 1, num_components
         call int2string(i, itext)
         call int2string(iteration, ittext)
         call make_result_map(res_map, i)
         fname = trim(output_dir) // 'output_map_comp_' // itext // '_iter_' // ittext //  '.fits'
         call write_map(res_map, 2, trim(fname))
         deallocate(res_map)
      end do
   end subroutine output_maps

   subroutine output_likelihood_nodp(parmin, parmax, nsteps, func, f_state, fname_in)
      implicit none
      real(dp), intent(in)      :: parmin, parmax
      integer(i4b)      :: nsteps
      integer(i4b), dimension(:)      :: f_state
      character(len=*), optional  :: fname_in
      interface 
         function func(x, s)
            use healpix_types
            implicit none
            real(dp)    :: func
            real(dp), intent(in)        :: x
            integer(i4b), dimension(:), intent(in)  :: s
         end function func
      end interface 

      integer(i4b)      :: i, j, unit
      real(dp)  :: dx, currx, res
      character(len=512)        :: fname

      if (present(fname_in)) then
         fname = fname_in
      else
         fname = 'like_out.dat'
      end if

      unit = getlun()

      open(unit, file=trim(fname))

      dx = (parmax - parmin) / nsteps

      currx = parmin
      do i = 1, nsteps + 1
         res = func(currx, f_state)
         write(unit, *) currx, res
         currx = currx + dx
      end do
      close(unit)

   end subroutine output_likelihood_nodp

   subroutine output_likelihood_dp(parmin, parmax, nsteps, func, f_state, aux_dp, fname_in)
      implicit none
      real(dp), intent(in)      :: parmin, parmax
      integer(i4b), intent(in)      :: nsteps
      integer(i4b), dimension(:), intent(in)      :: f_state
      real(dp), intent(in) :: aux_dp
      character(len=*), optional  :: fname_in
      interface 
         function func(x, s, v)
            use healpix_types
            implicit none
            real(dp), intent(in), dimension(:)        :: x
            integer(i4b), dimension(:), intent(in)      :: s
            real(dp), intent(in)        :: v
            real(dp), dimension(size(x))    :: func
         end function func
      end interface 

      integer(i4b)      :: i, j, unit
      real(dp), dimension(nsteps)      :: res, currx
      real(dp)  :: dx
      character(len=512)        :: fname

      if (present(fname_in)) then
         fname = fname_in
      else
         fname = 'like_out.dat'
      end if

      unit = getlun()

      open(unit, file=trim(fname))

      dx = (parmax - parmin) / nsteps

      do i = 1, nsteps + 1
         currx(i) = parmin + (i - 1) * dx
      end do
      res = func(currx, f_state, aux_dp)
      do i = 1, nsteps + 1
         write(unit, *) currx(i), res(i)
      end do
      close(unit)

   end subroutine output_likelihood_dp

   subroutine get_maxpoint_and_qrombint(mode, state, maxpoint, intres)
      implicit none
      integer(i4b), intent(in)  :: mode
      integer(i4b), intent(in), dimension(:)    :: state
      real(dp), intent(out)     :: maxpoint
      real(dp), intent(out)     :: intres

      real(dp), allocatable, dimension(:)       :: parvals, lnLvals
      integer(i4b)      :: currsize, component, i
      real(dp)  :: parlow, parhigh, lnL0, chisqmin, delta, lnL
      real(dp)  :: int_low, int_high
      real(dp), dimension(2)    :: limits
      real(dp)  :: dumchisq, dumchisq2
      logical(lgt)      :: err

      !This holds whether it is a pixel or region likelihood
      component = state(2)
      parlow = prior_low(component)
      parhigh = prior_high(component)

      if (mode == 1) then
         !Single pixel-likelihood
         
         call minimize_brent(parlow, parhigh, maxpoint, chisqmin, & 
            & get_single_pixel_chisq_singlepar, state)

         lnL0 = -0.5d0 * chisqmin
         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            delta = 2.d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         int_low = maxpoint - delta
         if (int_low < prior_low(component)) then
            int_low = prior_low(component)
         end if
         !Upper boundary
         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            delta = 2.d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         int_high = maxpoint + delta
         if (int_high > prior_high(component)) then
            int_high = prior_high(component)
         end if
         limits(1) = int_low
         limits(2) = int_high
         intres =  qromb(get_single_pixel_like, int_low, int_high, &
            & state, lnL0, err)

      else if (mode == 2) then
         !Region-likelihood
         call minimize_brent(parlow, parhigh, maxpoint, chisqmin, & 
            & get_region_chisq_singlepar, state)
         lnL0 = -0.5d0 * chisqmin

         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_region_like_singlepar(maxpoint-delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_region_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            delta = 2.d0 * delta
            lnL = get_region_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         int_low = maxpoint - delta
         if (int_low < prior_low(component)) then
            int_low = prior_low(component)
         end if
         !Upper boundary

         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            delta = 2.d0 * delta
            lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         int_high = maxpoint + delta
         if (int_high > prior_high(component)) then
            int_high = prior_high(component)
            lnL = get_region_like_singlepar(int_high, state, lnL0)
         end if
         limits(1) = int_low
         limits(2) = int_high
         !Normalize relative to parameter space
         intres = qromb(get_region_like, int_low, int_high, state, lnL0, &
            & err)
      end if
      if (isnan(intres)) stop
      intres = log(intres) + lnL0 - log(int_high - int_low)
   end subroutine get_maxpoint_and_qrombint

   subroutine get_maxpoint_and_splinedint(mode, state, maxpoint, intres)
      implicit none
      integer(i4b), intent(in)  :: mode
      integer(i4b), intent(in), dimension(:)    :: state
      real(dp), intent(out)     :: maxpoint
      real(dp), intent(out)     :: intres

      real(dp), allocatable, dimension(:)       :: parvals, lnLvals
      integer(i4b)      :: currsize, component, i
      real(dp)  :: parlow, parhigh, lnL0, chisqmin, delta, lnL
      real(dp)  :: int_low, int_high
      real(dp), dimension(2)    :: limits
      real(dp)  :: dumchisq, dumchisq2
      logical(lgt)      :: err

      !This holds whether it is a pixel or region likelihood
      component = state(2)
      parlow = prior_low(component)
      parhigh = prior_high(component)

      allocate(parvals(500))
      allocate(lnLvals(500))
      parvals = 0
      lnLvals = 0
      currsize = 0

      if (mode == 1) then
         !Single pixel-likelihood
         call minimize_brent(parlow, parhigh, maxpoint, chisqmin, & 
            & get_single_pixel_chisq_singlepar, state)

         lnL0 = -0.5d0 * chisqmin
         do i = 1, currsize
            lnLvals(i) = exp(-0.5d0 * lnLvals(i) - lnL0)
         end do
         call add_to_vals(maxpoint, 1.d0, parvals, lnLvals, currsize)
         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            call add_to_vals(maxpoint-delta, lnL, parvals, lnLvals, currsize)
            delta = 2.d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         int_low = maxpoint - delta
         if (int_low < prior_low(component)) then
            int_low = prior_low(component)
            lnL = get_single_pixel_like_singlepar(int_low, state, lnL0)
         end if
         call add_to_vals(int_low, lnL, parvals, lnLvals, currsize)
         !Upper boundary
         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            call add_to_vals(maxpoint+delta, lnL, parvals, lnLvals, currsize)
            delta = 2.d0 * delta
            lnL = get_single_pixel_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         int_high = maxpoint + delta
         if (int_high > prior_high(component)) then
            int_high = prior_high(component)
            lnL = get_single_pixel_like_singlepar(int_high, state, lnL0)
         end if
         call add_to_vals(int_high, lnL, parvals, lnLvals, currsize)
         limits(1) = int_low
         limits(2) = int_high
         call spline_int_refine(get_single_pixel_like_singlepar, int_low, int_high, state, lnL0, err, parvals, lnLvals, currsize, intres)

      else if (mode == 2) then
         !Region-likelihood
         call minimize_brent(parlow, parhigh, maxpoint, chisqmin, & 
            & get_region_chisq_singlepar, state)
         lnL0 = -0.5d0 * chisqmin
         do i = 1, currsize
            lnLvals(i) = exp(-0.5d0 * lnLvals(i) - lnL0)
         end do
         call add_to_vals(maxpoint, 1.d0, parvals, lnLvals, currsize)

         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_region_like_singlepar(maxpoint-delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_region_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            call add_to_vals(maxpoint-delta, lnL, parvals, lnLvals, currsize)
            delta = 2.d0 * delta
            lnL = get_region_like_singlepar(maxpoint - delta, state, lnL0)
         end do
         int_low = maxpoint - delta
         if (int_low < prior_low(component)) then
            int_low = prior_low(component)
            lnL = get_region_like_singlepar(int_low, state, lnL0)
         end if
         call add_to_vals(int_low, lnL, parvals, lnLvals, currsize)
         !Upper boundary

         delta = min(0.0001, 0.01 * abs(maxpoint))
         lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         do while (lnL < 1.d-5)
            delta = 0.5d0 * delta
            lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         do while (lnL > 1.d-5)
            call add_to_vals(maxpoint+delta, lnL, parvals, lnLvals, currsize)
            delta = 2.d0 * delta
            lnL = get_region_like_singlepar(maxpoint + delta, state, lnL0)
         end do
         int_high = maxpoint + delta
         if (int_high > prior_high(component)) then
            int_high = prior_high(component)
            lnL = get_region_like_singlepar(int_high, state, lnL0)
         end if
         call add_to_vals(int_high, lnL, parvals, lnLvals, currsize)
         limits(1) = int_low
         limits(2) = int_high
         call spline_int_refine(get_region_like_singlepar, int_low, int_high, state, lnL0, err, parvals, lnLvals, currsize, intres)
      end if
      if (isnan(intres)) stop
      intres = log(intres) + lnL0 - log(int_high - int_low)
   end subroutine get_maxpoint_and_splinedint

end module region_mod
