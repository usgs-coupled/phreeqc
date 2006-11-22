[Info]
Name=
Type=CompDef
Version=2.10.000

[Source Code]
OBJECT=No
DESCRIPTION=This component contains the program source code.
STATUS=
VISIBLE=Yes
DISK=ANYDISK
FILENEED=STANDARD
INCLUDEINBUILD=Yes
PASSWORD=
ENCRYPT=No
COMPRESSIFSEPARATE=No
UNINSTALLABLE=Yes
COMMENT=
DEFSELECTION=Yes
SELECTED=Yes
IMAGE=
TARGETDIRCDROM=
DISPLAYTEXT=
HTTPLOCATION=
FTPLOCATION=
MISC=
GUID=2758654F-F591-4E6E-8408-EE41AA83DB7D
_SPLIT_BEFORE=
_SPLIT_AFTER=
_DATAASFILES=
_NO_SPLIT=
_NO_SPLIT_BEFORE=
VOLATILE=
filegroup0=src
filegroup1=top-level
HANDLERONInstalling=
HANDLERONInstalled=
HANDLERONUnInstalling=
HANDLERONUnInstalled=

[Components]
component0=Source Code
component1=Program Files
component2=Verification Tests
component3=Documentation

[Program Files]
OBJECT=No
DESCRIPTION=This component contains the program executable, files required during execution, and minimal informal documentation.  It cannot be deselected unless the Verification Tests component is deselected.
STATUS=
VISIBLE=Yes
DISK=ANYDISK
FILENEED=STANDARD
INCLUDEINBUILD=Yes
PASSWORD=
ENCRYPT=No
COMPRESSIFSEPARATE=No
UNINSTALLABLE=Yes
COMMENT=
DEFSELECTION=Yes
SELECTED=Yes
IMAGE=
TARGETDIRCDROM=
DISPLAYTEXT=
HTTPLOCATION=
FTPLOCATION=
MISC=
GUID=63F0BED5-148A-475C-AD4C-7D93D57B1438
_SPLIT_BEFORE=
_SPLIT_AFTER=
_DATAASFILES=
_NO_SPLIT=
_NO_SPLIT_BEFORE=
VOLATILE=
requiredby0=Verification Tests
filegroup0=batch-file
filegroup1=bin
filegroup2=top-level
HANDLERONInstalling=ProgramFiles1_Installing
HANDLERONInstalled=ProgramFiles_Installed
HANDLERONUnInstalling=ProgramFiles_UnInstalling
HANDLERONUnInstalled=ProgramFiles_UnInstalled

[Verification Tests]
OBJECT=No
DESCRIPTION=This component contains the standard data sets and batch files used in running the verification tests.
STATUS=
VISIBLE=Yes
DISK=ANYDISK
FILENEED=STANDARD
INCLUDEINBUILD=Yes
PASSWORD=
ENCRYPT=No
COMPRESSIFSEPARATE=No
UNINSTALLABLE=Yes
COMMENT=
DEFSELECTION=Yes
SELECTED=Yes
IMAGE=
TARGETDIRCDROM=
DISPLAYTEXT=
HTTPLOCATION=
FTPLOCATION=
MISC=
GUID=23AC7AAD-D3C3-44AB-9FBE-039FFBB6E464
_SPLIT_BEFORE=
_SPLIT_AFTER=
_DATAASFILES=
_NO_SPLIT=
_NO_SPLIT_BEFORE=
VOLATILE=
required0=Program Files
filegroup0=examples
filegroup1=test
HANDLERONInstalling=
HANDLERONInstalled=VerificationTests_Installed
HANDLERONUnInstalling=
HANDLERONUnInstalled=

[Documentation]
OBJECT=No
DESCRIPTION=This component contains brief, summary documentation as well as the full, formal program documentation.
STATUS=
VISIBLE=Yes
DISK=ANYDISK
FILENEED=STANDARD
INCLUDEINBUILD=Yes
PASSWORD=
ENCRYPT=No
COMPRESSIFSEPARATE=No
UNINSTALLABLE=Yes
COMMENT=
DEFSELECTION=Yes
SELECTED=Yes
IMAGE=
TARGETDIRCDROM=
DISPLAYTEXT=
HTTPLOCATION=
FTPLOCATION=
MISC=
GUID=34045E58-6070-4FEB-8F2A-3C42E97C3DD9
_SPLIT_BEFORE=
_SPLIT_AFTER=
_DATAASFILES=
_NO_SPLIT=
_NO_SPLIT_BEFORE=
VOLATILE=
filegroup0=doc
filegroup1=top-level
HANDLERONInstalling=Documentation_Installing
HANDLERONInstalled=
HANDLERONUnInstalling=
HANDLERONUnInstalled=

[TopComponents]
component0=Program Files
component1=Documentation
component2=Verification Tests
component3=Source Code

[SetupType]
setuptype0=Compact
setuptype1=Typical
setuptype2=Custom

[SetupTypeItem-Compact]
Comment=
Descrip=
DisplayText=
item0=Program Files
item1=Tutorial
item2=Main App
item3=Examples

[SetupTypeItem-Typical]
Comment=
Descrip=
DisplayText=
item0=Source Code
item1=Program Files
item2=Verification Tests
item3=Tutorial
item4=Main App
item5=Examples
item6=Documentation

[SetupTypeItem-Custom]
Comment=
Descrip=
DisplayText=
item0=Source Code
item1=Program Files
item2=Verification Tests
item3=Tutorial
item4=Main App
item5=Examples
item6=Documentation

