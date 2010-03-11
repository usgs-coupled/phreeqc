
#import win32com.client as w32c 

from win32com.client import Dispatch
Wateq4f = Dispatch('PhreeqcCOM.Phreeqc');

from pylab import *
#
Istring = 'SOLUTION 1; END; REACTION; NaCl 1.0; 0 1 2 3 4 5 6 moles; EQUILIBRIUM_PHASES; gypsum; USE solution 1;'
Istring += 'SELECTED_OUTPUT; -reset false; -rxn; -total S(6); END'

print '\nWateq4f'
Wateq4f.LoadDatabase('C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database\wateq4f.dat');
Wateq4f.OutputOn = True
Wateq4f.RunString(Istring);
arr = Wateq4f.GetSelectedOutputArray()
for i in range(1, Wateq4f.RowCount):
	print arr[i][0], arr[i][1]

print '\nSit'	
Sit = Dispatch('PhreeqcCOM.Phreeqc');
Sit.LoadDatabase('C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database\minteq.dat');
Sit.RunString(Istring);
arrsit = Sit.GetSelectedOutputArray()
for i in range(1, Sit.RowCount):
	print arrsit[i][0], arrsit[i][1]


print '\nPitzer'
Pitzer = Dispatch('PhreeqcCOM.Phreeqc');
Pitzer.LoadDatabase('C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database\pitzer.dat');
Pitzer.RunString(Istring);
arrpitz = Pitzer.GetSelectedOutputArray()
for i in range(1, Pitzer.RowCount):
	print arrpitz[i][0], arrpitz[i][1]