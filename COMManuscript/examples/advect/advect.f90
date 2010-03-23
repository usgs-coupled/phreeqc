program Advect
  include "IPhreeqc.f90.inc"
  interface
    subroutine extractwrite(i, vtype, dv, sv)
      integer (kind=4),               intent(in)    :: i
      integer (kind=4), dimension(:), intent(inout) :: vtype
      real    (kind=8), dimension(:), intent(inout) :: dv
      character(len=*), dimension(:), intent(inout) :: sv   
    end subroutine
  end interface
  character(len=1024) Istring
  integer,            dimension(7) :: vtype
  real      (kind=8), dimension(7) :: dvalue
  character(len=100), dimension(7) :: svalue

!Load database, define initial conditions and selected output
  if (LoadDatabase("phreeqc.dat") < 0) call ErrorHandler()
  If (RunFile("ic") < 0) call ErrorHandler()

!Run cell 1, extract/write result
  if (RunString("RUN_CELLS; -cells; 1; END") < 0) call ErrorHandler()
  call ExtractWrite(1, vtype, dvalue, svalue)

!Transfer cell 1 solution to cell 2,run cell 2, extract/write results
  Ierr = AccumulateLine("SOLUTION_MODIFY 2")    
  Ierr = AccumulateLine("   -cb      " // svalue(1))
  Ierr = AccumulateLine("   -total_h " // svalue(2))
  Ierr = AccumulateLine("   -total_o " // svalue(3))
  Ierr = AccumulateLine("   -totals  ")
  Ierr = AccumulateLine("      C     " // svalue(4))
  Ierr = AccumulateLine("      Ca    " // svalue(5))
  Ierr = AccumulateLine("RUN_CELLS; -cells; 2; END")
  if (RunAccumulated() < 0) call ErrorHandler()
  call ExtractWrite(2, vtype, dvalue, svalue)
end program Advect

subroutine ErrorHandler()
  call OutputLastError()
  stop
end subroutine ErrorHandler

subroutine ExtractWrite(i, vtype, dv, sv)
  include "IPhreeqc.f90.inc"
  integer (kind=4),               intent(in)    :: i
  integer (kind=4), dimension(:), intent(inout) :: vtype
  real    (kind=8), dimension(:), intent(inout) :: dv
  character(len=*), dimension(:), intent(inout) :: sv    
  do j = 1, 7
    ! Headings are on row 0
    if (GetSelectedOutputValue(1,j,vtype(j),dv(j),sv(j)) < 0) call ErrorHandler()
  enddo
  write(*,"(a,i2/2(5x,a,f7.2))") "Cell",i,"pH:",dv(6),"SR(calcite):",dv(7)
end subroutine ExtractWrite