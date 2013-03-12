module utils
   use healpix_types
   implicit none

contains
   subroutine int2string(integer, string)
      implicit none

      integer(i4b),     intent(in) :: integer
      character(len=*), intent(out)     :: string

      integer(i4b)      :: temp_int, i, k

      temp_int = integer
      do i = 1, len(string)
         k = temp_int / 10 ** (len(string) - i)
         write(string(i:i), '(I1)') k
         temp_int = temp_int - k * 10 ** (len(string) - i)
      end do

   end subroutine int2string

   function getlun()
      implicit none
      integer(i4b) :: getlun
      logical(lgt) :: exists, isopen
      getlun = 9
      do
         getlun = getlun+1
         inquire(unit=getlun,exist=exists)
         if(exists) then
            inquire(unit=getlun,opened=isopen)
            if(.not. isopen) return
         end if
      end do
   end function getlun

   subroutine get_parameter(unit, parfile, parname, par_int, par_char, &
      & par_string, par_sp, par_dp, par_lgt)
    implicit none

    integer(i4b),      intent(in)  :: unit
    character(len=*),  intent(in)  :: parfile, parname
    integer(i4b),      intent(out), optional :: par_int
    logical(lgt),      intent(out), optional :: par_lgt
    character(len=1),  intent(out), optional :: par_char
    character(len=*),  intent(out), optional :: par_string
    real(sp),          intent(out), optional :: par_sp
    real(dp),          intent(out), optional :: par_dp


    integer(i4b)        :: i
    character(len=128)  :: string, variable, value
    character(len=1)    :: equals


    open(unit, file=trim(parfile))

    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*) variable, equals, value

       if (trim(variable) == trim(parname)) then

          if (present(par_int)) then
             read(value,*) par_int
          else if (present(par_char)) then
             read(value,*) par_char
          else if (present(par_string)) then
             read(value,'(a)') par_string
          else if (present(par_sp)) then
             read(value,*) par_sp
          else if (present(par_dp)) then
             read(value,*) par_dp
          else if (present(par_lgt)) then
             read(value,*) par_lgt
          end if

          close(unit)
          return

       end if

    end do

1   write(*,*) 'GET_PARAMETER:    Critical error -- parameter not found'
    write(*,*) 'GET_PARAMETER:       Parameter file = ', trim(parfile)
    write(*,*) 'GET_PARAMETER:       Parameter name = ', trim(parname)

    close(unit)
    stop

  end subroutine get_parameter

  subroutine dump_arrays(x, y, fname)
     implicit none
     
     real(dp), dimension(:), intent(in) :: x, y
     character(len=*), intent(in)   :: fname

     integer(i4b)       :: n, i, unit

     n = size(x)
     unit = getlun()
     open(unit, file=trim(fname))
     do i = 1, n
        write(unit, *), x(i), y(i)
     end do

     close(unit)

  end subroutine dump_arrays

end module utils
