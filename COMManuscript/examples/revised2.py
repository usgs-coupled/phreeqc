"""Compares gypsum solubility for WATEQ4F and Pitzer databases.
"""
# Import standard library modules first.
import os
# Then get third party modules.
from win32com.client import Dispatch
import matplotlib.pyplot as plt

def selected_array(db_path, input_string):
    """Load database via COM and run input string.
    """
    dbase = Dispatch('IPhreeqcCOM.Object')
    dbase.LoadDatabase(db_path)
    dbase.RunString(input_string)
    return dbase.GetSelectedOutputArray()

def show_results(dbdir, input_string):
    """Get results for different databases
    """
    # Use os.path.join to build platform-independent path names.
    # Make a short function name for convenience.
    join = os.path.join
    wateq4f_result = selected_array(join(dbdir, 'wateq4f.dat'), input_string)
    pitzer_result  = selected_array(join(dbdir, 'pitzer.dat'), input_string)
    # Get data from the arrays.
    nacl_conc      = [entry[0] for entry in wateq4f_result][1:]
    wateq4f_values = [entry[1] for entry in wateq4f_result][1:]
    pitzer_values  = [entry[1] for entry in pitzer_result][1:]
    # Plot
    plt.plot(nacl_conc, pitzer_values, 'r', nacl_conc, wateq4f_values,'b--')
    plt.axis([0, 6, 0, .06])
    plt.legend(('PITZER','WATEQ4F'), loc = (0.4, 0.4))
    plt.ylabel('GYPSUM SOLUBILITY, MOLES PER KILOGRAM WATER')
    plt.xlabel('NaCl, MOLES PER KILOGRAM WATER')
    #plt.backend      : PS
    plt.show()
    
if __name__ == '__main__':
    # This will only run when called as script from the command line
    # and not when imported from another script.
    DBDIR = r'C:\Program Files (x86)\USGS\Phreeqc Interactive 2.17.4137\database'
    INPUT_STRING = """
    SOLUTION 1
    END
    INCREMENTAL_REACTIONS
    REACTION
    	NaCl 1.0
    	0 60*0.1 moles
    EQUILIBRIUM_PHASES
    	Gypsum
    USE solution 1
    SELECTED_OUTPUT
    	-reset false
    	-total Na S(6)
    END"""
    show_results(DBDIR, INPUT_STRING)