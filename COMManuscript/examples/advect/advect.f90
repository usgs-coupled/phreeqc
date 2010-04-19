module Subs
  integer    (kind=4), dimension(7) :: vt
  real       (kind=8), dimension(7) :: dv
  character (len=100), dimension(7) :: sv 
  contains
  subroutine ExtractWrite(cell)
    include "IPhreeqc.f90.inc"
    integer    (kind=4), intent(in) :: cell
    do j = 1, 7
      ! Headings are on row 0
      if(GetSelectedOutputValue(1,j,vt(j),dv(j),sv(j)).ne.VR_OK) call EHandler()
    enddo
    write(*,"(a,i2/2(5x,a,f7.2))") "Cell",cell,"pH:",dv(6),"SR(calcite):",dv(7) 
  end subroutine ExtractWrite
  subroutine EHandler()
    include "IPhreeqc.f90.inc"
!    integer i
!    character(len=100) errstr
!    do i = 1, GetErrorLineCount()
!       call GetErrorLine(i, errstr)
!       write (*,*) trim(errstr)
!    enddo
     call OutputLastError()
    stop
  end subroutine EHandler    
end module Subs
program Advect
  use Subs
  include "IPhreeqc.f90.inc"
  character(len=1024) Istring
  
!Load database, define initial conditions and selected output
  if (LoadDatabase("phreeqc.dat") .ne. 0) call EHandler()
  If (RunFile("ic") .ne. 0) call EHandler()

!Run cell 1, extract/write result
  if (RunString("RUN_CELLS; -cells; 1; END") .ne. 0) call EHandler()
  call ExtractWrite(1)

!Advect cell 1 solution to cell 2, run cell 2, extract/write results
  Ierr = AccumulateLine("SOLUTION_MODIFY 2")    
  Ierr = AccumulateLine("   -cb      " // sv(1))
  Ierr = AccumulateLine("   -total_h " // sv(2))
  Ierr = AccumulateLine("   -total_o " // sv(3))
  Ierr = AccumulateLine("   -xtotals  ")
  Ierr = AccumulateLine("      C     " // sv(4))
  Ierr = AccumulateLine("      Ca    " // sv(5))
  Ierr = AccumulateLine("RUN_CELLS; -cells; 2; END")
  if (RunAccumulated() .ne. 0) call EHandler()
  call ExtractWrite(2)
end program Advect