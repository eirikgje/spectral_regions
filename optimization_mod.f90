module optimization_mod
   use healpix_types
   use nr_mod
   use utils
   implicit none

contains
   subroutine spline_int_refine(func, int_low, int_high, aux_var, var_dp, err, xarr, yarr, currsize, res)
      implicit none
      real(dp), intent(in)      :: int_low, int_high
      integer(i4b), intent(in), dimension(:)      :: aux_var
      real(dp), intent(in)      :: var_dp
      logical(lgt), intent(out) :: err
      real(dp), intent(inout), allocatable, dimension(:)     :: xarr, yarr
      integer(i4b), intent(inout)       :: currsize
      real(dp), intent(out)     :: res
      interface
         function func(x, aux, var_dp)
            use healpix_types
            real(dp), intent(in)        :: x
            real(dp)                     :: func
            integer(i4b), dimension(:), intent(in)   :: aux
            real(dp), intent(in) :: var_dp
         end function func
      end interface

      real(dp)  :: curres, newres
      real(dp)  :: TOL = 1d-6
      real(dp), allocatable, dimension(:)       :: y2

      allocate(y2(currsize))
      call quicksort_simarrays(xarr, yarr, 1, currsize)
      if (xarr(1) /= int_low) stop("xarr not equal to int_low")
      if (xarr(currsize) /= int_high) stop("xarr not equal to int_high")
      call spline(xarr(1:currsize), yarr(1:currsize), 1.d30, 1.d30, y2)
      call spline_int(xarr(1:currsize), yarr(1:currsize), y2, curres)
      call populate_arrays(xarr, yarr, currsize)
      deallocate(y2)
      allocate(y2(currsize))
      call spline(xarr(1:currsize), yarr(1:currsize), 1.d30, 1.d30, y2)
      call spline_int(xarr(1:currsize), yarr(1:currsize), y2, newres)
      do while (abs((curres - newres) / newres) > TOL)
         print *, (curres - newres) / newres
         print *, 'currsize', currsize
         print *, 'curres', curres
         print *, 'newres', newres
         curres = newres
         call populate_arrays(xarr, yarr, currsize)
         deallocate(y2)
         allocate(y2(currsize))
         call spline(xarr(1:currsize), yarr(1:currsize), 1.d30, 1.d30, y2)
         call spline_int(xarr(1:currsize), yarr(1:currsize), y2, newres)
      end do
      res = newres
      err = .False.

      contains
         subroutine populate_arrays(xarr, yarr, currsize)
            !Assumes xarr is sorted (and that yarr corresponds to xarr)
            implicit none
            real(dp), intent(inout), allocatable, dimension(:)   :: xarr, yarr
            integer(i4b), intent(inout) :: currsize

            real(dp), allocatable, dimension(:) :: temparr1, temparr2
            integer(i4b)        :: i

            allocate(temparr1(currsize + currsize - 1))
            allocate(temparr2(currsize + currsize - 1))
            do i = 1, currsize - 1
               temparr1(2*i-1) = xarr(i)
               temparr1(2*i) = 0.5d0 * (xarr(i + 1) + xarr(i))
               temparr2(2*i-1) = yarr(i)
               temparr2(2*i) =  func(temparr1(2*i), aux_var, var_dp)
            end do
            temparr1(2 * currsize - 1) = xarr(currsize)
            temparr2(2 * currsize - 1) = yarr(currsize)
            currsize = 2 * currsize - 1
            if (currsize > size(xarr)) then
               deallocate(xarr, yarr)
               allocate(xarr(2 * currsize), yarr(2 * currsize))
            end if
            xarr(1:currsize) = temparr1
            yarr(1:currsize) = temparr2
         end subroutine populate_arrays
   end subroutine spline_int_refine

   subroutine spline_int(x, y, y2, res)
      !Computes the analytic integral of a splined function
      implicit none
      real(dp), dimension(:), intent(in)        :: x, y, y2
      real(dp), intent(out)     :: res

      real(dp)  :: z, onebyz, d0, d1, d2, d3
      integer(i4b)      :: n, i, j, k
      real(dp)  :: onebytwentyfour
      onebytwentyfour = 1.d0 / 24.d0
      n = size(x)

      res = 0
      do i = 1, n - 1
         z = x(i + 1) - x(i)
         res = res + 0.5d0 * z * (y(i+1) + y(i)) - onebytwentyfour * z ** 3 & 
            & * (y2(i + 1) + y2(i)) 
      end do
   end subroutine spline_int

   recursive subroutine quicksort_simarrays(arr, dumarr, left, right)
      implicit none
      real(dp), intent(inout), dimension(:)     :: arr, dumarr
      integer(i4b), intent(in)      :: left, right

      integer(i4b)      :: piv_ind, piv_new_ind

      if (left < right) then
         piv_ind = (right + left) / 2
         piv_new_ind = partition(arr, dumarr, left, right, piv_ind)
         call quicksort_simarrays(arr, dumarr, left, piv_new_ind - 1)
         call quicksort_simarrays(arr, dumarr, piv_new_ind + 1, right)
      end if
      contains
         function partition(arr, dumarr, left, right, piv_ind) result(store_ind)
            implicit none
            real(dp), intent(inout), dimension(:)       :: arr, dumarr
            integer(i4b), intent(in)    :: left, right
            integer(i4b), intent(in)    :: piv_ind
            integer(i4b)        :: store_ind

            real(dp)    :: piv_val
            integer(i4b)        :: i
            
            piv_val = arr(piv_ind)

            call swap(arr(piv_ind), arr(right))
            call swap(dumarr(piv_ind), dumarr(right))
            store_ind = left
            do i = left, right - 1
               if (arr(i) < piv_val) then
                  call swap(arr(i), arr(store_ind))
                  call swap(dumarr(i), dumarr(store_ind))
                  store_ind = store_ind + 1
               end if
            end do
            call swap(arr(store_ind), arr(right))
            call swap(dumarr(store_ind), dumarr(right))
         end function partition
         
         subroutine swap(val1, val2)
            implicit none
            real(dp), intent(inout)       :: val1, val2

            real(dp)    :: swapval

            swapval = val1
            val1 = val2
            val2 = swapval
         end subroutine swap

   end subroutine quicksort_simarrays

   subroutine dump_splint(x, y, y2, xmin, xmax, num_points, fname)
      implicit none
      real(dp), intent(in), dimension(:)  :: x, y, y2
      real(dp), intent(in)      :: xmin, xmax
      integer(i4b), intent(in)      :: num_points
      character(len=*), intent(in)      :: fname

      integer(i4b)      :: n
      real(dp), dimension(num_points)   :: xarr, yarr
      real(dp)  :: dx
      integer(i4b)      :: i, unit

      n = size(x)

      unit = getlun()
      open(unit, file=trim(fname))
      dx = (xmax - xmin) / (num_points - 1)
      do i = 1, num_points
         xarr(i) = xmin + (i - 1) * dx
         yarr(i) = splint(x, y, y2, n, xarr(i))
         write(unit, *), xarr(i), yarr(i)
      end do

      close(unit)

   end subroutine dump_splint
end module optimization_mod
