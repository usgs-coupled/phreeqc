# Microsoft Developer Studio Project File - Name="phreeqc_console" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=phreeqc_console - Win32 Test
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "phreeqc_console.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "phreeqc_console.mak" CFG="phreeqc_console - Win32 Test"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "phreeqc_console - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "phreeqc_console - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE "phreeqc_console - Win32 Test" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "phreeqc_console - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Release/"
# ADD F90 /include:"Release/"
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "DOS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /version:2.5 /subsystem:console /debug /machine:I386 /out:"Release/phreeqc.exe" /RELEASE
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/"
# ADD F90 /browser /include:"Debug/"
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "DOS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"Debug/phreeqc.exe" /pdbtype:sept

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Test"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "phreeqc_console___Win32_Test"
# PROP BASE Intermediate_Dir "phreeqc_console___Win32_Test"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Test"
# PROP Intermediate_Dir "Test"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /browser /include:"Debug/"
# ADD F90 /browser /include:"Debug/"
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "DOS" /FR /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\..\src" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "DOS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"Debug/phreeqc.exe" /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"Test/phreeqc.exe" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "phreeqc_console - Win32 Release"
# Name "phreeqc_console - Win32 Debug"
# Name "phreeqc_console - Win32 Test"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\advection.c
# End Source File
# Begin Source File

SOURCE=..\..\src\basic.c
# End Source File
# Begin Source File

SOURCE=..\..\src\basicsubs.c
# End Source File
# Begin Source File

SOURCE=..\..\src\cl1.c
# End Source File
# Begin Source File

SOURCE=..\..\src\cvdense.c
# End Source File
# Begin Source File

SOURCE=..\..\src\cvode.c
# End Source File
# Begin Source File

SOURCE=..\..\src\dense.c
# End Source File
# Begin Source File

SOURCE=..\..\src\input.c
# End Source File
# Begin Source File

SOURCE=..\..\src\integrate.c
# End Source File
# Begin Source File

SOURCE=..\..\src\inverse.c
# End Source File
# Begin Source File

SOURCE=..\..\src\isotopes.c
# End Source File
# Begin Source File

SOURCE=..\..\src\kinetics.c
# End Source File
# Begin Source File

SOURCE=..\..\src\main.c

!IF  "$(CFG)" == "phreeqc_console - Win32 Release"

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Debug"

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Test"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\mainsubs.c
# End Source File
# Begin Source File

SOURCE=..\..\src\model.c
# End Source File
# Begin Source File

SOURCE=..\..\src\nvector.c
# End Source File
# Begin Source File

SOURCE=..\..\src\nvector_serial.c
# End Source File
# Begin Source File

SOURCE=..\..\src\output.c
# End Source File
# Begin Source File

SOURCE=..\..\src\p2clib.c
# End Source File
# Begin Source File

SOURCE=..\..\src\parse.c
# End Source File
# Begin Source File

SOURCE=..\..\src\phqalloc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\phreeqc_files.c
# End Source File
# Begin Source File

SOURCE=..\..\src\prep.c
# End Source File
# Begin Source File

SOURCE=..\..\src\print.c
# End Source File
# Begin Source File

SOURCE=..\..\src\read.c
# End Source File
# Begin Source File

SOURCE=..\..\src\readtr.c
# End Source File
# Begin Source File

SOURCE=..\..\src\smalldense.c
# End Source File
# Begin Source File

SOURCE=..\..\src\spread.c
# End Source File
# Begin Source File

SOURCE=..\..\src\step.c
# End Source File
# Begin Source File

SOURCE=..\..\src\structures.c
# End Source File
# Begin Source File

SOURCE=..\..\src\sundialsmath.c
# End Source File
# Begin Source File

SOURCE=..\..\src\tally.c
# End Source File
# Begin Source File

SOURCE=.\test.cxx

!IF  "$(CFG)" == "phreeqc_console - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "phreeqc_console - Win32 Test"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\tidy.c
# End Source File
# Begin Source File

SOURCE=..\..\src\transport.c
# End Source File
# Begin Source File

SOURCE=..\..\src\utilities.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\cvdense.h
# End Source File
# Begin Source File

SOURCE=..\..\src\cvode.h
# End Source File
# Begin Source File

SOURCE=..\..\src\dense.h
# End Source File
# Begin Source File

SOURCE=..\..\src\global.h
# End Source File
# Begin Source File

SOURCE=..\..\src\input.h
# End Source File
# Begin Source File

SOURCE=..\..\src\kinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\src\nvector.h
# End Source File
# Begin Source File

SOURCE=..\..\src\nvector_serial.h
# End Source File
# Begin Source File

SOURCE=..\..\src\output.h
# End Source File
# Begin Source File

SOURCE=..\..\src\p2c.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phqalloc.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phrqproto.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phrqtype.h
# End Source File
# Begin Source File

SOURCE=..\..\src\smalldense.h
# End Source File
# Begin Source File

SOURCE=..\..\src\sundialsmath.h
# End Source File
# Begin Source File

SOURCE=..\..\src\sundialstypes.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\resource.rc
# End Source File
# End Group
# End Target
# End Project
