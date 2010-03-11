program advect
  include "IPhreeqc.f90.inc"
  integer Errors, rows, columns, i  
  character(len=1024) Istring
  integer, dimension(7) :: vtype
  character(len=30), dimension(7) :: svalue
  real(kind=8), dimension(7) :: dvalue

! Load database, define initial conditions and selected output
  if (LoadDatabase("phreeqc.dat") .ne. 0) call ErrorHandler()
  If (RunFile("ic",.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()

! Run cell 1, extract results, and write result
  Errors = AccumulateLine("USE solution 1; USE equilibrium_phases 1; SAVE solution 1; END")
  If (Run(.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  !if (RunString("RUN_CELLS; -cells; 1; END") .ne. 0) call ErrorHandler()
  call ExtractSelectedOutput(vtype, dvalue, svalue)
  call WriteResults(1, dvalue)
  
! Transfer cell 1 solution to cell 2,run cell 2, extract results, and write result
!  Errors = AccumulateLine("SOLUTION_MODIFY 2")    
!  Errors = AccumulateLine("   -cb      " // svalue(1))
!  Errors = AccumulateLine("   -total_h " // svalue(2))
!  Errors = AccumulateLine("   -total_o " // svalue(3))
!  Errors = AccumulateLine("   -totals ")
!  Errors = AccumulateLine("      C     " // svalue(4))
!  Errors = AccumulateLine("      Ca    " // svalue(5))
!  Errors = AccumulateLine("RUN_CELLS; -cells; 2; END")
  Errors = AccumulateLine("COPY solution 1 2; END")
  Errors = AccumulateLine("USE solution 2; USE equilibrium_phases 2; END")
  if (Run(.true.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  call ExtractSelectedOutput(vtype, dvalue, svalue)
  call WriteResults(2, dvalue)
end program advect
subroutine ErrorHandler()
  call OutputLastError()
  stop
end subroutine ErrorHandler
subroutine WriteResults(cell, dvalue)
  integer cell
  real(kind=8), dimension(:) :: dvalue
  write(*,"(a,i2/a,f7.2,a,f7.2)") &
    "Cell", cell, "    pH: ", dvalue(6), "    SR(calcite): ", dvalue(7)
end subroutine WriteResults
subroutine ExtractSelectedOutput(vtype, dvalue, svalue)
  integer, dimension(:) :: vtype
  real(kind=8), dimension(:) :: dvalue  
  character(*), dimension(:) :: svalue
  do j = 1, 7
    Iresult = GetSelectedOutputValue(1, j, vtype(j), dvalue(j), svalue(j))
    if (vtype(j) .eq. 3) write(svalue(j),"(1pE21.14)") dvalue(j)
  enddo
end subroutine ExtractSelectedOutput