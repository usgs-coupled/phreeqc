Dim Phreeqc, Errors, Arr, Ostring

' Create module 
Set Phreeqc = CreateObject("PhreeqcCOM.Phreeqc")

' Initiate module
Phreeqc.LoadDatabase("C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database\phreeqc.dat")

' Initial conditions
Phreeqc.AccumulateLine("Solution 1-2")
Phreeqc.AccumulateLine("END")
Phreeqc.AccumulateLine("EQUILIBRIUM_PHASES 1")
Phreeqc.AccumulateLine("   CO2(g) -1.5 10")
Phreeqc.AccumulateLine("EQUILIBRIUM_PHASES 2")
Phreeqc.AccumulateLine("   Calcite 0   10")
Phreeqc.AccumulateLine("END")
Errors = Phreeqc.Run()
If Errors <> 0 Then Call ErrorMessage(Errors)

' Write selected output and run cell 1
Phreeqc.AccumulateLine("SELECTED_OUTPUT")
Phreeqc.AccumulateLine("   -reset false")
Phreeqc.AccumulateLine("END")
Phreeqc.AccumulateLine("USER_PUNCH")
Phreeqc.AccumulateLine("   -Heading  charge" & vbTab & "H" & vbTab & "O" & vbTab & "C" & vbTab & "Ca" & vbTab & "pH" & vbTab & "SR(calcite)")
Phreeqc.AccumulateLine("10 w = TOT(""water"")")
Phreeqc.AccumulateLine("20 PUNCH charge_balance, TOT(""H"")*w, TOT(""O"")*w, TOT(""C"")*w, TOT(""Ca"")*w")
Phreeqc.AccumulateLine("30 PUNCH -LA(""H+""), SR(""calcite"")")
Phreeqc.AccumulateLine("END")
Phreeqc.AccumulateLine("RUN_CELLS" & vbCRLF & "-cells" & vbCRLF &"   1" & vbCRLF & "END")
Phreeqc.OutputOn = True
Errors = Phreeqc.Run()
If Errors <> 0 Then Call ErrorMessage(Errors)
Call WriteResults(1, Ostring, Phreeqc.GetSelectedOutputArray())

' Transfer to cell 2 and run cell 2
Phreeqc.AccumulateLine("SOLUTION_MODIFY 2")
Phreeqc.AccumulateLine("   -cb      " & Phreeqc.GetSelectedOutputValue(1, 0))
Phreeqc.AccumulateLine("   -total_h " & Phreeqc.GetSelectedOutputValue(1, 1))
Phreeqc.AccumulateLine("   -total_o " & Phreeqc.GetSelectedOutputValue(1, 2))
Phreeqc.AccumulateLine("   -totals ")
Phreeqc.AccumulateLine("      C     " & Phreeqc.GetSelectedOutputValue(1, 3))
Phreeqc.AccumulateLine("      Ca    " & Phreeqc.GetSelectedOutputValue(1, 4))
Phreeqc.AccumulateLine("RUN_CELLS" & vbCRLF & "-cells" & vbCRLF & "   2" & vbCRLF & "END")
Errors = Phreeqc.Run()
Call WriteResults(2, Ostring, Phreeqc.GetSelectedOutputArray())

' Print results
MsgBox(Ostring)




'' Write SELECTED_OUTPUT and USER_PUNCH
'Dim clist, ccount, i
'ccount = Phreeqc.list_components(clist)
'var sInput = "SELECTED_OUTPUT\n";
'sInput +=    "   -reset false\n"
'sInput +=    "USER_PUNCH\n"
'sInput +=    "   -heading\t"
'For i = 0 To ccount - 1
'    sInput += clist(i) + "\t"
'Next
'sInput +=    "\n"
'sInput +=    "10 PUNCH "
'Next
'For i = 0 To ccount - 1
'    sInput += "TOT(""" + clist(i) + """)*m\t"
'Next
'sInput +=    "\n"
'sInput +=    "END"
'Errors = RunString(sInput)
''


Sub ErrorMessage(ByRef Errors) 
  MsgBox("Number of errors: " + STR(err) & vbCRLF & Phreeqc.GetLastErrorString()) 
  STOP
End Sub
Sub WriteResults(ByVal cell, ByRef Ostring, ByRef Arr) 
  Ostring = Ostring & "Cell: " & cell & vbCRLF
  For i = 0 To Phreeqc.RowCount - 1
      For j = Phreeqc.ColumnCount - 2 To Phreeqc.ColumnCount - 1
          If (IsNumeric(Arr(i,j))) Then Tstring = FormatNumber(Arr(i,j), 2) Else Tstring = Arr(i,j)
          Ostring = ostring & Tstring & vbTab
      Next
      Ostring = Ostring & vbCRLF
  Next

End Sub