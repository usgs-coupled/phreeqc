program advect
  include "IPhreeqc.f90.inc"
  character(len=1024) Istring
  integer,           dimension(7) :: vtype
  character(len=30), dimension(7) :: svalue
  real(kind=8),      dimension(7) :: dvalue

! Load database, define initial conditions and selected output
  if (LoadDatabase("phreeqc.dat") .ne. 0) call ErrorHandler()
  If (RunFile("ic",.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()

! Run cell 1, extract/write result
  Ierr = AccumulateLine("USE solution 1; USE equilibrium_phases 1; SAVE solution 1; END")
  If (Run(.TRUE.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  !if (RunString("RUN_CELLS; -cells; 1; END") .ne. 0) call ErrorHandler()
  call ExtractWrite(1, vtype, dvalue, svalue)
  
! Transfer cell 1 solution to cell 2,run cell 2, extract/write results
!  Ierr = AccumulateLine("SOLUTION_MODIFY 2")    
!  Ierr = AccumulateLine("   -cb      " // svalue(1))
!  Ierr = AccumulateLine("   -total_h " // svalue(2))
!  Ierr = AccumulateLine("   -total_o " // svalue(3))
!  Ierr = AccumulateLine("   -totals ")
!  Ierr = AccumulateLine("      C     " // svalue(4))
!  Ierr = AccumulateLine("      Ca    " // svalue(5))
!  Ierr = AccumulateLine("RUN_CELLS; -cells; 2; END")
  Ierr = AccumulateLine("COPY solution 1 2; END")
  Ierr = AccumulateLine("USE solution 2; USE equilibrium_phases 2; END")
  if (Run(.true.,.FALSE.,.FALSE.,.FALSE.) .ne. 0) call ErrorHandler()
  call ExtractWrite(2, vtype, dvalue, svalue)
end program advect

subroutine ErrorHandler()
  call OutputLastError()
  stop
end subroutine ErrorHandler

subroutine ExtractWrite(icell, vtype, dvalue, svalue)
  integer,      dimension(:) :: vtype
  real(kind=8), dimension(:) :: dvalue  
  character(*), dimension(:) :: svalue
  do j = 1, 7
    Iresult = GetSelectedOutputValue(1, j, vtype(j), dvalue(j), svalue(j))
    if (vtype(j) .eq. 3) write(svalue(j),"(1pE21.14)") dvalue(j)
  enddo
  write(*,"(a,i2/a,f7.2,a,f7.2)") &
    "Cell", icell, "    pH: ", dvalue(6), "  SR(calcite): ", dvalue(7)
end subroutine ExtractWrite