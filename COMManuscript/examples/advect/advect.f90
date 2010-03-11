program advect
  implicit none
  include "IPhreeqc.f90.inc"
  integer Errors, rows, columns, i  
  character(len=1024) Istring
  integer, dimension(7) :: vtype
  character(len=30), dimension(7) :: svalue
  real(kind=8), dimension(7) :: dvalue
! Load database
  Errors = LoadDatabase("phreeqc.dat")
  if (Errors .ne. 0) call ErrorHandler()
! Initial conditions and selected output

    !Errors = RunFile("ic",.TRUE.,.FALSE.,.FALSE.,.FALSE.)
    ! ----to be replaced by runfile
    Errors = AccumulateLine("Solution 1-2; END")
    Errors = AccumulateLine("EQUILIBRIUM_PHASES 1; CO2(g) -1.5 10")
    Errors = AccumulateLine("EQUILIBRIUM_PHASES 2; Calcite 0   10")
    Errors = AccumulateLine("END")
    Errors = AccumulateLine("SELECTED_OUTPUT; -reset false")
    Errors = AccumulateLine("USER_PUNCH")
    Errors = AccumulateLine("   -Heading  charge    H   O   C   Ca  pH  SR(calcite)")
    Errors = AccumulateLine("10 w = TOT(""water"")")
    Errors = AccumulateLine("20 PUNCH charge_balance, TOT(""H"")*w, TOT(""O"")*w, TOT(""C"")*w, TOT(""Ca"")*w")
    Errors = AccumulateLine("30 PUNCH -LA(""H+""), SR(""calcite"")")
    Errors = AccumulateLine("END")
    ! ----to be replaced by runfile
    !Errors = AccumulateLine("RUN_CELLS; -cells 1; END")
    Errors = AccumulateLine("USE solution 1; USE equilibrium_phases 1; END")
    !call OutputLines 
    If (Run(.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  !If (RunFile("ic",.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
! Run cell 1
  !write(Istring,"(a/)")   "RUN_CELLS; -cells; 1; END"
  !Errors = AccumulateLine(Istring)    
  !if (Run(.FALSE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  !if (RunString("RUN_CELLS; -cells; 1; END") .ne. 0) call ErrorHandler()
! Extract results and write    
  call ExtractSelectedOutput(vtype, dvalue, svalue)
  call WriteResults(1, dvalue)
! Transfer cell 1 solution to cell 2 and run cell 2
  Errors = AccumulateLine("SOLUTION_MODIFY 2")    
  Errors = AccumulateLine("   -cb      " // svalue(1))
  Errors = AccumulateLine("   -total_h " // svalue(2))
  Errors = AccumulateLine("   -total_o " // svalue(3))
  Errors = AccumulateLine("   -totals ")
  Errors = AccumulateLine("      C     " // svalue(4))
  Errors = AccumulateLine("      Ca    " // svalue(5))
  !Errors = AccumulateLine("RUN_CELLS; -cells; 2; END")
  Errors = AccumulateLine("USE solution 2; USE equilibrium_phases 2; END")
  call OutputLines()
  if (Run(.true.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  !if (RunString(Istring) .ne. 0) call ErrorHandler()
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
  write(*,"(a,i/a,f7.2,a,f7.2)") &
    "Cell: ", cell, "    pH: ", dvalue(6), "    SR(calcite): ", dvalue(7)
end subroutine WriteResults
subroutine ExtractSelectedOutput(vtype, dvalue, svalue)
  integer, dimension(:) :: vtype
  real(kind=8), dimension(:) :: dvalue  
  character(len=30), dimension(:) :: svalue
  integer :: j, Iresult
  do j = 1, 7
    Iresult = GetSelectedOutputValue(1, j, vtype(j), dvalue(j), svalue(j))
    !write (*,*) "Extracting: ", i, j, "iresult: ", Iresult, dvalue(j)
    !if (GetSelectedOutputValue(1, j, vtype(j), dvalue(j), svalue(j)) .ne. 0) call ErrorHandler()
    if (vtype(j) .eq. 3) write(svalue(j),"(1pE21.14)") dvalue(j)
  enddo
end subroutine ExtractSelectedOutput