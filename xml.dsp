# Microsoft Developer Studio Project File - Name="xml" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=xml - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "xml.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "xml.mak" CFG="xml - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "xml - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "xml - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "xml - Win32 Release"

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
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /browser /compile_only /nologo /optimize:4 /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "C:\xerces-c-windows_2000-msvc_60\include" /I "$(DEV_GMP_INC)" /D "INVERSE_CL1MP" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 xerces-c_2.lib gmp.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  /nologo /subsystem:console /machine:I386 /libpath:"$(DEV_GMP_LIB)" /libpath:"C:\xerces-c-windows_2000-msvc_60\lib" /libpath:"C:\xerces-c-windows_2000-msvc_60\bin"
# SUBTRACT LINK32 /incremental:yes /nodefaultlib

!ELSEIF  "$(CFG)" == "xml - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "xml___Win32_Debug"
# PROP BASE Intermediate_Dir "xml___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "C:\xerces-c-windows_2000-msvc_60\include" /I "$(DEV_GMP_INC)" /D "WIN32_MEMORY_DEBUG" /D "INVERSE_CL1MP" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /i "C:\xerces-c-windows_2000-msvc_60\bin" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 xerces-c_2.lib gmpDebug.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept /libpath:"$(DEV_GMP_LIB)" /libpath:"C:\xerces-c-windows_2000-msvc_60\lib" /libpath:"C:\xerces-c-windows_2000-msvc_60\bin"

!ENDIF 

# Begin Target

# Name "xml - Win32 Release"
# Name "xml - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\advection.c
# End Source File
# Begin Source File

SOURCE=.\basic.c
# End Source File
# Begin Source File

SOURCE=.\basicsubs.c
# End Source File
# Begin Source File

SOURCE=.\cl1.c
# End Source File
# Begin Source File

SOURCE=.\cl1mp.c
# End Source File
# Begin Source File

SOURCE=.\cvdense.c
# End Source File
# Begin Source File

SOURCE=.\cvode.c
# End Source File
# Begin Source File

SOURCE=.\dense.c
# End Source File
# Begin Source File

SOURCE=.\dw.c
# End Source File
# Begin Source File

SOURCE=.\input.c
# End Source File
# Begin Source File

SOURCE=.\integrate.c
# End Source File
# Begin Source File

SOURCE=.\inverse.c
# End Source File
# Begin Source File

SOURCE=.\isotopes.c
# End Source File
# Begin Source File

SOURCE=.\kinetics.c
# End Source File
# Begin Source File

SOURCE=.\main.c
# End Source File
# Begin Source File

SOURCE=.\mainsubs.c
# End Source File
# Begin Source File

SOURCE=.\model.c
# End Source File
# Begin Source File

SOURCE=.\nvector.c
# End Source File
# Begin Source File

SOURCE=.\nvector_serial.c
# End Source File
# Begin Source File

SOURCE=.\output.c
# End Source File
# Begin Source File

SOURCE=.\p2clib.c
# End Source File
# Begin Source File

SOURCE=.\parse.c
# End Source File
# Begin Source File

SOURCE=.\phqalloc.c
# End Source File
# Begin Source File

SOURCE=.\phreeqc_files.c
# End Source File
# Begin Source File

SOURCE=.\pitzer.c
# End Source File
# Begin Source File

SOURCE=.\pitzer_structures.c
# End Source File
# Begin Source File

SOURCE=.\prep.c
# End Source File
# Begin Source File

SOURCE=.\print.c
# End Source File
# Begin Source File

SOURCE=.\read.c
# End Source File
# Begin Source File

SOURCE=.\readtr.c
# End Source File
# Begin Source File

SOURCE=.\SAXPhreeqc.cpp
# End Source File
# Begin Source File

SOURCE=.\smalldense.c
# End Source File
# Begin Source File

SOURCE=.\spread.c
# End Source File
# Begin Source File

SOURCE=.\step.c
# End Source File
# Begin Source File

SOURCE=.\structures.c
# End Source File
# Begin Source File

SOURCE=.\sundialsmath.c
# End Source File
# Begin Source File

SOURCE=.\tally.c
# End Source File
# Begin Source File

SOURCE=.\tidy.c
# End Source File
# Begin Source File

SOURCE=.\transport.c
# End Source File
# Begin Source File

SOURCE=.\utilities.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\cvdense.h
# End Source File
# Begin Source File

SOURCE=.\cvode.h
# End Source File
# Begin Source File

SOURCE=.\dense.h
# End Source File
# Begin Source File

SOURCE=.\global.h
# End Source File
# Begin Source File

SOURCE=.\input.h
# End Source File
# Begin Source File

SOURCE=.\kinetics.h
# End Source File
# Begin Source File

SOURCE=.\nvector.h
# End Source File
# Begin Source File

SOURCE=.\nvector_serial.h
# End Source File
# Begin Source File

SOURCE=.\output.h
# End Source File
# Begin Source File

SOURCE=.\p2c.h
# End Source File
# Begin Source File

SOURCE=.\phqalloc.h
# End Source File
# Begin Source File

SOURCE=.\phrqproto.h
# End Source File
# Begin Source File

SOURCE=.\phrqtype.h
# End Source File
# Begin Source File

SOURCE=.\pitzer.h
# End Source File
# Begin Source File

SOURCE=.\SAXPhreeqc.h
# End Source File
# Begin Source File

SOURCE=.\SaxPhreeqcHandlers.h
# End Source File
# Begin Source File

SOURCE=.\smalldense.h
# End Source File
# Begin Source File

SOURCE=.\sundialsmath.h
# End Source File
# Begin Source File

SOURCE=.\sundialstypes.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project