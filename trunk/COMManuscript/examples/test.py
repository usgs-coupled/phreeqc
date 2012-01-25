from win32com.client import Dispatch
from pylab import *

Istring =  'SOLUTION 1; END; INCREMENTAL_REACTIONS;'
Istring += 'REACTION; NaCl 1.0; 0 60*0.1 moles;'
Istring += 'EQUILIBRIUM_PHASES; gypsum; USE solution 1;'
Istring += 'SELECTED_OUTPUT; -reset false; -total Na S(6); END'

dbdir='C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database'
Wateq4f = Dispatch('PhreeqcCOM.Phreeqc');
Wateq4f.LoadDatabase(dbdir + '\wateq4f.dat');
Wateq4f.OutputOn = True
Wateq4f.RunString(Istring);
arrw = Wateq4f.GetSelectedOutputArray()

Pitzer = Dispatch('PhreeqcCOM.Phreeqc');
Pitzer.LoadDatabase(dbdir + '\pitzer.dat');
Pitzer.OutputOn = True
Pitzer.RunString(Istring);
arrp = Pitzer.GetSelectedOutputArray()
	
x = []
yw = []
yp = []
for i in range(1, Wateq4f.RowCount):
	x.append(arrw[i][0])
	yw.append(arrw[i][1])
	yp.append(arrp[i][1])
	
import matplotlib.pyplot as plt
plt.plot(x, yp, 'r', x, yw,'b--')
plt.axis([0, 6, 0, .06])
plt.legend(('PITZER','WATEQ4F'), loc = (0.4, 0.4))
plt.ylabel('GYPSUM SOLUBILITY, MOLES PER KILOGRAM WATER')
plt.xlabel('NaCl, MOLES PER KILOGRAM WATER')
plt.show()

	
