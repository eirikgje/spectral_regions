program find_bestfit_si
   use healpix_types
   use utils
   use rngmod
   use quiet_mapfile_mod
   use optimization_mod
   implicit none

   character(len=512)   :: paramfile
   integer(i4b) :: unit
   integer(i4b) :: iterations, nbins, nside, npix, seed
   integer(i4b) :: numband, num_components
   integer(i4b) :: i, j, k, ordering, curreg
   real(dp)     :: x, chisq
   real(dp), allocatable, dimension(:, :)       :: dat, map_temp, invN
   real(dp), allocatable, dimension(:, :)       :: par_map, amp_map
   real(dp), allocatable, dimension(:, :)       :: chisq_map
   integer(i4b), allocatable, dimension(:, :)   :: reg_map, scrambled_regions
   integer(i4b), allocatable, dimension(:, :)   :: representative_pixels
   integer(i4b), allocatable, dimension(:, :)   :: region_mapping
   real(dp), allocatable, dimension(:)  :: freq, prior_low, prior_high
   real(dp), allocatable, dimension(:)  :: a2t, ref_freq
   integer(i4b), allocatable, dimension(:)      :: behavior, num_regions
   type(planck_rng)     :: rng_handle
   character(len=512)   :: output_dir
   character(len=2)   :: itext
   character(len=512), allocatable, dimension(:)   :: regionfname, mapfname
   character(len=512), allocatable, dimension(:)   :: invnfname

   real(dp), parameter          :: K_BOLTZMANN = 1.3806503d-23 !In SI-units
   real(dp), parameter          :: H_PLANCK = 6.626068d-34 !In SI-units
   real(dp), parameter          :: T_CMB = 2.72548 !In SI-units

   call getarg(1, paramfile)

   unit = getlun()
   !Get parameters
   call get_parameter(unit, paramfile, 'ITERATIONS', par_int=iterations)
   call get_parameter(unit, paramfile, 'SEED', par_int=seed)
   call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
   call get_parameter(unit, paramfile, 'NUM_SPECTRAL_BEHAVIORS', & 
      & par_int=num_components)
   call get_parameter(unit, paramfile, 'NUMBAND', par_int=numband)
   call get_parameter(unit, paramfile, 'OUTPUT_DIR', par_string=output_dir)
   call get_parameter(unit, paramfile, 'NBINS', par_int=nbins)
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
   allocate(behavior(num_components))
   allocate(ref_freq(num_components))
   allocate(prior_low(num_components))
   allocate(prior_high(num_components))
   allocate(regionfname(num_components))
   npix = 12 * nside ** 2
   do i = 1, num_components
      call int2string(i, itext)
      !1: Power-law
      !2: template curve
      !etc.
      call get_parameter(unit, paramfile, 'BEHAVIOR_' // itext, par_int=behavior(i))
      call get_parameter(unit, paramfile, 'REF_FREQ_' // itext,par_dp=ref_freq(i))
      call get_parameter(unit, paramfile, 'PRIOR_LOW_' // itext, par_dp=prior_low(i))
      call get_parameter(unit, paramfile, 'PRIOR_HIGH_' // itext, par_dp=prior_high(i))
      call get_parameter(unit, paramfile, 'REGIONFNAME_' // itext, par_string=regionfname(i))
   end do
   call rand_init(rng_handle, seed)

   allocate(a2t(numband))

   do i = 1, numband
      x = H_PLANCK * freq(i) / (K_BOLTZMANN * T_CMB)
      a2t(i) = (exp(x) - 1) ** 2 / (x ** 2 * exp(x))
   end do

   allocate(dat(0:npix-1, numband))
   allocate(invN(0:npix-1, numband))
   do i = 1, numband
      call read_map(map_temp, ordering, trim(mapfname(i)))
      if (ordering == 1) then
         call convert_ring2nest(nside, map_temp(:, 1))
      end if
      dat(:, i) = map_temp(:, 1)
      deallocate(map_temp)
      call read_map(map_temp, ordering, trim(invnfname(i)))
      if (ordering == 1) then
         call convert_ring2nest(nside, map_temp(:, 1))
      end if
      invN(:, i) = map_temp(:, 1)
      deallocate(map_temp)
   end do

   allocate(num_regions(num_components))
   allocate(par_map(0:npix-1, num_components))
   allocate(reg_map(0:npix-1, num_components))
   allocate(region_mapping(npix, num_components))
   do i = 1, num_components
      call read_map(map_temp, ordering, trim(regionfname(i)))
      if (ordering == 1) then
         do j = 1, 2
            call convert_ring2nest(nside, map_temp(:, j))
         end do
      end if
      par_map(:, i) = map_temp(:, 1)
      reg_map(:, i) = map_temp(:, 2)
      call get_num_regions(map_temp(:, 2), num_regions(i), region_mapping(:, i))
      deallocate(map_temp)
   end do
   allocate(representative_pixels(num_components, maxval(num_regions)))
   do i = 1, num_components
      do j = 1, num_regions(i)
         do k = 0, npix - 1
            if (int(reg_map(k, i)) == region_mapping(j, i)) then
               representative_pixels(i, j) = k
               exit
            end if
            if (k == npix - 1) then
               print *, "No pixels found for region ", j, " component ", i
            end if
         end do
      end do
   end do

   !Find initial bestfit amplitudes
   allocate(amp_map(0:npix-1, num_components))
   do i = 0, npix - 1
      call solve_amp(amp_map(i, :), par_map(i, :), i)
   end do

   !Find initial chisq
   call eval_chisq(chisq_map, par_map, amp_map)
   chisq = sum(chisq_map)

   call output_maps(0)

   !Loop through
   do i = 1, iterations
      print *, 'iter', i
      !Scramble region lists
      call scramble_region_list(scrambled_regions)
      do j = 1, num_components
         do k = 1, num_regions(j)
            curreg = region_mapping(scrambled_regions(k, j), j)
            print *, 'region', curreg
            call find_bestfit_si_and_amp(curreg, scrambled_regions(k, j), j, & 
              & chisq)
         end do
      end do
      print *, 'curr chisq:', chisq
      call output_maps(i)
   end do

contains
   subroutine scramble_region_list(list)
      implicit none
      integer(i4b), allocatable, dimension(:, :), intent(out)       :: list

      integer(i4b)      :: i, j, currnumreg, currind, temp

      if (allocated(list)) deallocate(list)
      allocate(list(maxval(num_regions), num_components))

      do i = 1, num_components
         do j = 1, num_regions(i)
            list(j, i) = j
         end do
      end do

      do i = 1, num_components
         currnumreg = num_regions(i)
         do while (currnumreg > 0)
            currind = int(rand_uni(rng_handle) * (currnumreg)) + 1
            temp = list(currind, i)
            list(currind, i) = list(currnumreg, i)
            list(currnumreg, i) = temp
            currnumreg = currnumreg - 1
         end do
      end do

   end subroutine scramble_region_list

   subroutine find_bestfit_si_and_amp(curreg, reg_num, component, chisq)
      implicit none

      integer(i4b), intent(in)      :: curreg, component, reg_num
      real(dp), intent(inout)   :: chisq

      real(dp), allocatable, dimension(:, :)    :: curr_amp_map, prop_amp_map
      real(dp), allocatable, dimension(:, :)    :: chisq_map
      real(dp)  :: dv, curr_bf_par, currchisq, currval, prop_chisq

      integer(i4b)      :: n, i, j

      allocate(curr_amp_map(0:npix-1, num_components))
      allocate(prop_amp_map(0:npix-1, num_components))

      dv = (prior_high(component) - prior_low(component)) / (nbins - 1)
      curr_bf_par = par_map(representative_pixels(component, reg_num), component)
      curr_amp_map = amp_map
      prop_amp_map = amp_map
      currchisq = chisq

      do n = 0, nbins - 1
         currval = prior_low(component) + n * dv
         prop_amp_map = curr_amp_map
         !$OMP PARALLEL PRIVATE(i)
         !$OMP DO SCHEDULE(STATIC)
         do i = 0, npix-1
            if (reg_map(i, component) == curreg) then
               par_map(i, component) = currval
               call solve_amp(prop_amp_map(i, :), par_map(i, :), i)
            end if
         end do
         !$OMP END DO
         !$OMP END PARALLEL
         call eval_chisq(chisq_map, par_map, prop_amp_map)
         prop_chisq = sum(chisq_map)
!         print *, 'prop_chisq', prop_chisq
!         print *, 'curr_chisq', currchisq      
         deallocate(chisq_map)
         if (prop_chisq < currchisq) then
            print *, 'accept'
            currchisq = prop_chisq
            curr_amp_map = prop_amp_map
            curr_bf_par = currval
         end if
      end do
      do i = 0, npix - 1
         if (reg_map(i, component) == curreg) then
            par_map(i, component) = curr_bf_par
         end if
      end do
      amp_map = curr_amp_map
      chisq = currchisq

   end subroutine find_bestfit_si_and_amp

   subroutine solve_amp(amp_map, par_map, pixnum)
      implicit none
      real(dp), dimension(:), intent(in)        :: par_map
      integer(i4b), intent(in)      :: pixnum
      real(dp), dimension(:), intent(out)       :: amp_map

      real(dp), allocatable, dimension(:, :)       :: mixmat
      real(dp), allocatable, dimension(:, :)    :: temp, minv
      real(dp), allocatable, dimension(:)       :: diag, x
      integer(i4b)      :: i, j

      allocate(mixmat(numband, num_components))
      allocate(temp(numband, num_components))
      allocate(minv(num_components, num_components))
      allocate(diag(num_components))
      allocate(x(num_components))
      mixmat = 0

      do i = 1, num_components
         do j = 1, numband
            if (behavior(i) == 1) then
               mixmat(j, i) = a2t(j) * (freq(j) / ref_freq(i)) ** par_map(i)
            end if
         end do
      end do

      do i = 1, numband
         temp(i, :) = mixmat(i, :) * invN(pixnum, i)
      end do

      minv = matmul(transpose(mixmat), temp)
      call choldc(minv, diag)

      x = matmul(transpose(temp), dat(pixnum, :))
      amp_map = x
      call cholsl(minv, diag, x, amp_map)

   end subroutine solve_amp

   subroutine eval_chisq(chisq_map, par_map, amp_map)
      implicit none

      real(dp), dimension(0:, :), intent(in)    :: par_map, amp_map
      real(dp), allocatable, dimension(:, :)  :: chisq_map

      real(dp), allocatable, dimension(:, :) :: signal_maps
      integer(i4b)      :: i, j

      allocate(signal_maps(0:npix-1, numband))
      signal_maps = 0

      do i = 1, num_components
         if (behavior(i) == 1) then
            do j = 1, numband
               signal_maps(:, j) = signal_maps(:, j) + amp_map(:, i) * a2t(j) &
                  & * (freq(j) / ref_freq(i)) ** par_map(:, i)
            end do
         end if
      end do

      signal_maps = dat - signal_maps
      allocate(chisq_map(0:npix-1, 1))
      chisq_map(:, 1) = sum(signal_maps ** 2 * invN, 2)
      
   end subroutine eval_chisq

   subroutine output_maps(iteration)
      implicit none

      integer(i4b), intent(in)  :: iteration

      integer(i4b)              :: i
      character(len=2)          :: itext
      character(len=5)          :: ittext
      character(len=512)        :: fname
      real(dp), allocatable, dimension(:, :)    :: res_map

      allocate(res_map(0:npix - 1, 1))
      do i = 1, num_components
         call int2string(i, itext)
         call int2string(iteration, ittext)
         res_map = amp_map(:, i:i)
         fname = trim(output_dir) // 'amp_map_comp_' // itext // '_iter_' // ittext //  '.fits'
         call write_map(res_map, 2, trim(fname))
         res_map = par_map(:, i:i)
         fname = trim(output_dir) // 'par_map_comp_' // itext // '_iter_' // ittext //  '.fits'
         call write_map(res_map, 2, trim(fname))
      end do
      deallocate(res_map)
      call eval_chisq(res_map, par_map, amp_map)
      fname = trim(output_dir) // 'chisq_map_iter_' // ittext // '.fits'
      call write_map(res_map, 2, trim(fname))
   end subroutine output_maps

   subroutine get_num_regions(regmap, num_regions, mapping)
      implicit none
      real(dp), dimension(:), intent(in)    :: regmap
      integer(i4b), dimension(:), intent(out)   :: mapping
      integer(i4b)      :: num_regions

      integer(i4b)      :: i, j, foundreg_ind

!      allocate(foundregs(size(regmap)))
      mapping = -1

      foundreg_ind = 0

      do i = 1, size(regmap)
         if (.not. any(mapping(1:foundreg_ind) == regmap(i))) then
            foundreg_ind = foundreg_ind + 1
            mapping(foundreg_ind) = regmap(i)
         end if
      end do

      num_regions = foundreg_ind

   end subroutine get_num_regions
end program find_bestfit_si
