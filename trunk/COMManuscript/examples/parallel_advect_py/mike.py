"""Advection with COM server. Object oruente solution.

Using MODFIY we update the concentration on every
time step. We shift by one cell per step.
"""

import multiprocessing
import os
import time

import matplotlib.pyplot as plt
from win32com.client import Dispatch

USE_MULTIRPOCESSING = True

class StreamTube(object):
    """Stream tube for 1D-advection.

    This is just a row of cells that can do
    reactions. The actual advection  is done
    "from outside" by shifting the cell contents.
    A tube has intial state a in flow and an ouflow.
    Several tubes can be chained to run the calcultions
    in parallel. For this purpose the outflow from
    the last time step of the connected tube would
    serve as inflow.
    """

    def __init__(self, cells, initial_conditions, inflow):
        """
        cells - first and last cell inclusive, e.g. (10, 25) would be
                translate into a cell range 10-25 (16 cells)
        intial_conditions - string containing PHREEQC input for
                            solution and exchange, see example below
        inflow - dictionary with inflow cocentrations, e.g.
                 {'Ca': 0.002, 'cb': 3e-12, 'Na': 0.001}
        """
        self.cells = cells
        self.initial_conditions = initial_conditions
        self.inflow = inflow
        self.phreeqc = Dispatch('IPhreeqcCOM.Object')
        self.phreeqc.LoadDatabase(r"phreeqc.dat")
        self.tools = PhreeqcTools(self.phreeqc)
        self.components = []
        self.create_initial_state()

    def create_initial_state(self):
        """Copy solution to all cells.
        """
        self.phreeqc.RunString(self.initial_conditions)
        self.components = self.phreeqc.GetComponentList()
        start = self.cells[0]
        end = self.cells[1]
        code = ''
        code += "COPY solution 1 %d-%d\n" % (start, end)
        code += "COPY exchange 1 %d-%d\n" % (start, end)
        code += "END\n"
        code += "RUN_CELLS; -cells %d-%d\n" % (start, end)
        code += self.tools.make_selected_output(self.components)
        self.phreeqc.RunString(code)
        self.conc = self.tools.get_selected_output()

    def advect_step(self):
        """Advect by shifting concentrations from previous time step.
        """
        conc = self.conc
        all_names = conc.keys()
        names = [name for name in all_names if name not in ('cb', 'H', 'O')]
        for name in conc:
            # shift one cell
            conc[name][1:] = conc[name][:-1]
            conc[name][0] = self.inflow[name]
        modify = []
        for index, cell in enumerate(xrange(self.cells[0], self.cells[1] + 1)):
            modify.append("SOLUTION_MODIFY %d" % cell)
            modify.append("\t-cb      %e" % conc['cb'][index])
            modify.append("\t-total_h %f" % conc['H'][index])
            modify.append("\t-total_o %f" % conc['O'][index])
            modify.append("\t-totals")
            for name in names:
                modify.append("\t\t%s\t%f" % (name, conc[name][index]))
        modify.append("RUN_CELLS; -cells %d-%d\n" % (self.cells[0],
                                                     self.cells[1]))
        code = '\n'.join(modify)
        self.phreeqc.RunString(code)
        outflow = {}
        conc = self.tools.get_selected_output()
        for name in all_names:
            outflow[name] = conc[name][-1]
        self.conc = conc
        return outflow


class Advection(object):
    """Use plug with shifts for advection.

    We use one or more instance of StreamTube to
    calculate advection.
    Each strem tube might run in its own process
    using multiprocessing
    """

    def __init__(self, ncells, shifts, initial_conditions, processes):
        if processes > ncells:
            raise ValueError('Number of proceeses needs to be less or equal '
                             'than number of cells. %d processes %d cells.'
                             % (processes, ncells))
        if processes < 1:
            raise ValueError('Need at least one process got %d' % processes)
        self.ncells = ncells
        self.shifts = shifts
        self.initial_conditions = initial_conditions
        self.processes = processes
        self.phreeqc = Dispatch('IPhreeqcCOM.Object')
        self.phreeqc.LoadDatabase(r"phreeqc.dat")
        self.tools = PhreeqcTools(self.phreeqc)
        self.inflow = {}
        self.initial = {}
        self.components = []
        self.make_initial_state()

    def make_initial_state(self):
        """Run PHREEQC to calculate initial concentations.
        """
        self.phreeqc.RunString(self.initial_conditions)
        self.components = self.phreeqc.GetComponentList()
        code = ''
        code += self.tools.make_selected_output(self.components)
        code += "RUN_CELLS; -cells 0-1\n"
        self.phreeqc.RunString(code)
        conc = self.tools.get_selected_output()
        for name in conc:
            self.inflow[name] = conc[name][0]
            self.initial[name] = conc[name][1]

    def run_single(self):
        """Do one run in one process.
        """
        cells = (1, self.ncells)
        tube = StreamTube(cells, self.initial_conditions, self.inflow)
        outflows = {}
        for name in self.components:
            outflows[name] = []
        for _ in xrange(self.shifts):
            # advect
            outflow = tube.advect_step()
            for name in self.components:
                outflows[name].append(outflow[name])
        return outflows

    def run_multiprocessing(self):
        """Do parallel runs with multiprocessing.
        """
        if self.processes == 1:
            return self.run_single()
        # Domain decomposition.
        slave_ncells, reminder = divmod(self.ncells, self.processes)
        root_ncells = slave_ncells + reminder
        root_cells = (1, root_ncells)
        root_tube = StreamTube(root_cells, self.initial_conditions,
                               self.inflow)
        current_cell = root_ncells
        slaves = []
        for process in xrange(self.processes - 1):
            slave_cells = (current_cell + 1, current_cell + slave_ncells)
            slaves.append(StreamTubeProxy(slave_cells, self.initial_conditions,
                                          self.initial))
            current_cell += slave_ncells
        assert current_cell == self.ncells
        # We want them starting from the last.
        slaves.reverse()
        outflows = {}
        for name in self.components:
            outflows[name] = []
        for _ in xrange(self.shifts):
            # advect
            inflows = []
            outflow = slaves[0].advect_step()
            for name in self.components:
                outflows[name].append(outflow[name])
            for slave in slaves[1:]:
                inflows.append(slave.advect_step())
            inflows.append(root_tube.advect_step())
            for index, slave in enumerate(slaves):
                slave.inflow = inflows[index]
        for slave in slaves:
            slave.finish()
        return outflows

class PhreeqcTools(object):
    """Collection of useful functions.

    We put them in a class to avoid moving arround
    the PHREEQC-COM-server. It cannot be pickle and
    hence not be moved between processes.
    """

    def __init__(self, phreeqc):
        """Set the COM-server.
        """
        self.phreeqc = phreeqc

    @ staticmethod # this is just a function but belongs here
    def make_selected_output(components):
        """
        Build SELECTED_OUTPUT data block
        """
        headings = "-headings    cb    H    O    "
        headings += '\t'.join(components)
        selected_output = """
        SELECTED_OUTPUT
            -reset false
        USER_PUNCH
        """
        selected_output += headings + "\n"
        # charge balance, H, and O
        code = '10 w = TOT("water")\n'
        code += '20 PUNCH CHARGE_BALANCE, TOTMOLE("H"), TOTMOLE("O")\n'
        # All other elements
        lino = 30
        for component in components:
            code += '%d PUNCH w*TOT(\"%s\")\n' % (lino, component)
            lino += 10
        selected_output += code
        return selected_output

    def get_selected_output(self):
        """Return calculation result as dict.

        Header entries are the keys and the columns
        are the values as lists of numbers.
        """
        output = self.phreeqc.GetSelectedOutputArray()
        header = output[0]
        conc = {}
        for head in header:
            conc[head] = []
        for row in output[1:]:
            for col, head in enumerate(header):
                conc[head].append(row[col])
        return conc


def plot(ncells, outflow, specie_names):
    """Plot the results.
    """
    colors = {'Ca': 'r', 'Cl': 'b', 'K': 'g', 'N': 'y', 'Na': 'm'}
    x = [i / float(ncells) for i in xrange(1,
                                           len(outflow[specie_names[0]]) + 1)]
    args = []
    for name in specie_names:
        args.extend([x, outflow[name], colors[name]])
    # pylint: disable-msg=W0142
    # we do want *
    plt.plot(*args)
    plt.legend(specie_names, loc=(0.8, 0.5))
    plt.ylabel('MILLIMOLES PER KILOGRAM WATER')
    plt.xlabel('PORE VOLUME')
    plt.show()

if not USE_MULTIRPOCESSING:
    class StreamTubeProxy(StreamTube):
        """Proxy for non-multiprocessing test.

        We can use this for testing. There is no multiprocessing
        but we can use instances of StreamTubeProxy that
        should produce the same result as the multiprocessing
        version.
        """
        # pylint: disable-msg=R0903
        # one method is enought
        def finish(self):
            """Mock to make it work just like multiprocessing version.
            """
            pass
else:
    class StreamTubeProxy(object):
        """Proxy that communicates with other process.
        """
        def __init__(self, cells, initial_conditions, inflow):
            """Go parallel.
            """
            self.inflow = inflow
            self.inflow_queue = multiprocessing.JoinableQueue()
            self.outflow_queue = multiprocessing.JoinableQueue()
            name = 'tube_%d_%d' % cells
            self.process = multiprocessing.Process(
                target=process_worker, name=name,
                args=(cells, initial_conditions, inflow, self.inflow_queue,
                      self.outflow_queue))
            self.process.start()
            self.inflow_queue.put(self.inflow)
            
        def advect_step(self):
            """Do advection in other process.
            """
            self.inflow_queue.put(self.inflow)
            outflow = self.outflow_queue.get()
            return outflow

        def finish(self):
            """Terminate the process.
            """
            self.inflow_queue.put(None)
            self.process.join()

def process_worker(cells, initial_conditions, inflow, inflow_queue,
                   outflow_queue):
    """This runs in another process.
    """
    print 'started', os.getpid()
    tube = StreamTube(cells, initial_conditions, inflow)
    while True:
        inflow = inflow_queue.get()
        # None is the sentinel. We are done
        if inflow is None:
            break
        tube.inflow = inflow
        outflow_queue.put(tube.advect_step())

def measure_time(func, *args, **kwargs):
    """Convinience function to measure run times.
    """
    import sys
    if sys.platform == 'win32':
        # time.clock is more accurate on Windows
        timer_func = time.clock
    else:
        # but behaves differently on other platfroms
        timer_func = time.time
    start = timer_func()
    result = func(*args, **kwargs)
    return result, time.clock() - start

if __name__ == '__main__':

    def main(ncells, shifts, processes=2):
        """
        Specify initial conditions data blocks.

        Uniform initial conditions are assumed.
        """
        initial_conditions = """
        TITLE Example 11.--Transport and ion exchange.
        SOLUTION 0  CaCl2
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Ca               0.6
            Cl               1.2
        SOLUTION 1  Initial solution for column
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Na               1.0
            K                0.2
            N(5)             1.2
            END
        EXCHANGE 1
            equilibrate 1
            X                0.0011
        END
            """
        def run():
            """Do the work.
            """
            advect = Advection(ncells, shifts, initial_conditions, processes)
            outflow = advect.run_multiprocessing()
            return advect, outflow
        (advect, outflow), run_time = measure_time(run)
        print 'Statistics'
        print '=========='
        print 'number of cells:    ', ncells
        print 'number of shifts:   ', shifts
        print 'number of processes:', processes
        print 'run_time:           ', run_time 
        plot(ncells, outflow, advect.components)
    main(ncells=40, shifts=120, processes=2)
