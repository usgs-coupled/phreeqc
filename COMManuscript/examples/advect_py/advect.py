from win32com.client import Dispatch
import time 
from  multiprocessing import *
from numpy import *
import matplotlib.pyplot as plt

def set_initial_conditions(CELLS):
	#
	# Specify initial conditions data blocks
	# Uniform initial conditions are assumed
	#
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
        return initial_conditions

def set_selected_output(components):  
	#
	# Build SELECTED_OUTPUT data block
	#
	headings = "-headings	H	O	cb	"
	for i in range(len(components)):
		headings += components[i] + "\t"
	selected_output = """
	SELECTED_OUTPUT
		-reset false
	USER_PUNCH
	"""
	selected_output += headings + "\n"
	#
	# charge balance, H, and O
	#
	code  = '10 w = TOT("water")\n'
	code += '20 PUNCH CHARGE_BALANCE, TOTMOLE("H"), TOTMOLE("O")\n'
	#
	# All other elements
	#
	no = 30
	for i in xrange(len(components)):
		code += str(no) + ' PUNCH w*TOT(\"' + components[i] + '\")\n'
		no += 10
	selected_output += code 
	return selected_output

def initialize(CELLS, PROCESSES, task_queue, done_queue):
	#
	# Initialize IPhreeqc module
	#
	Phreeqc = Dispatch('IPhreeqcCOM.Object')
	Phreeqc.LoadDatabase(r"phreeqc.dat")
	initial_conditions = set_initial_conditions(Phreeqc)
	Phreeqc.RunString(initial_conditions)
	components = Phreeqc.GetComponentList()
	#
	# Send components to processes
	# No data returned in done queue
	#
	for j in range(PROCESSES):
		args = []
		args.append("components")
		args.append(components)
		task_queue[j].put(args)
	#
	# Create, SELECTED_OUTPUT data block
	# and run in root process
	#
	selected_output = set_selected_output(components)
	Phreeqc.RunString(selected_output)
	#
	# Calculate and distribute cell ranges for processes
	# No data returned in done queue
	#
	ranges = []
	for j in range(PROCESSES):
		#
		# Calculate cell range for process
		# Send to process
		#	
		min = j*CELLS/PROCESSES + 1
		max = (j+1)*CELLS/PROCESSES
		if j == PROCESSES - 1:
			max = CELLS
		range_string = str(min) + "-" + str(max)
		ranges.append([min, max])
		r = [min, max]
		args = []
		args.append("range")
		args.append(r)
		task_queue[j].put(args)
		#
		# Distribute initial conditions to process
		#
        	task = initial_conditions + "\n"
	        task += "COPY solution 1 " + range_string + "\n"
		task += "COPY exchange 1 " + range_string + "\n"
		task += "END\n"
		task += "RUN_CELLS; -cells " + range_string + "\n"
		task += selected_output
		args = []
		args.append("string")
		args.append(task)
		task_queue[j].put(args)
	#
	# Clear done queue from initial conditions
	#
	for j in range(PROCESSES):
		done_queue[j].get()
	
	#
	# Generate initial and infilling solution compositions
	#	
	Phreeqc.RunString("RUN_CELLS; -cells 0-1\n")
	infilling = zeros([Phreeqc.ColumnCount], double)	
	for j in range(Phreeqc.ColumnCount):
		infilling[j] = Phreeqc.GetSelectedOutputValue(1, j)
	initial = zeros([Phreeqc.ColumnCount], double)	
	for j in range(Phreeqc.ColumnCount):
		initial[j] = Phreeqc.GetSelectedOutputValue(2, j)
	#
	# Generate initial concentration array
	#
	conc = zeros([CELLS, Phreeqc.ColumnCount], double)
	for i in range(CELLS):
		conc[i] = initial
	results = []
	results.append(ranges)
	results.append(infilling)
	results.append(conc)
	return results

def IP_run(task_queue, done):
	initial_run = True
	IP_worker = Dispatch('IPhreeqcCOM.Object')
	IP_worker.LoadDatabase(r"phreeqc.dat")
	cell_range = []
	components = []
	while True:
		next_task = task_queue.get()
		if next_task is None:
			break
		elif next_task[0] == "components":
			components = next_task[1]
		elif next_task[0] == "range":
			cell_range = next_task[1]
		elif next_task[0] == "string":
			IP_worker.RunString(next_task[1])
			done.put(IP_worker.GetSelectedOutputArray())
		elif next_task[0] == "array":
			task = ""
			eol ="\n"
			conc = next_task[1]
			k = 0
			for i in xrange(cell_range[0], cell_range[1] + 1):
				soln = conc[k]
				k += 1
				modify = "SOLUTION_MODIFY " + str(i) + eol
				modify += "\t-cb      " + str(soln[0]) + eol
				modify += "\t-total_h " + str(soln[1]) + eol
				modify += "\t-total_o " + str(soln[2]) + eol
				modify += "\t-totals  " + eol
				for j in range(len(components)):
					modify += "\t\t" + components[j] + "\t" + str(soln[3 + j]) + eol
				task += modify
			task += "RUN_CELLS; -cells " + str(cell_range[0]) + "-" + str(cell_range[1]) + eol
				
			IP_worker.RunString(task)
			conc = []
			for i in xrange(1, IP_worker.RowCount):
				soln = []
				for j in range(IP_worker.ColumnCount):
					soln.append(IP_worker.GetSelectedOutputValue(i, j))
				conc.append(soln)
			done.put(conc)
		else:
			print "Error in IP_run"
	return
	
if __name__ == '__main__':

	PROCESSES = 8
	CELLS = 400
	SHIFTS = 1200
	
	start = time.clock()
	#
	# Create queues and processes
	#
	task_queue = []
	done_queue = []
	for i in range(PROCESSES):
		task_queue.append(Queue())
		done_queue.append(Queue())

	for i in range(PROCESSES):
		Process(target=IP_run, args=(task_queue[i], done_queue[i])).start()
	print "Done starting processes."

	#
	# Defined and process initial conditions for all cells
	#
	initialize_results = initialize(CELLS, PROCESSES, task_queue, done_queue)
	ranges = initialize_results[0]
	infilling = initialize_results[1]
	conc = initialize_results[2]
	print "Done initialize."
	#
	# Loop for advection + reaction
	#
	outflow = zeros([SHIFTS, len(conc[0])], double)
	for t in xrange(SHIFTS):
	 	#
 		# Advect
 		#
 		for i in range(CELLS - 1):
 			conc[CELLS - i - 1] = conc[CELLS - i - 2]
 		conc[0] = infilling
 			
 		#
		# Fill queues for new solution compositions
		#
		for i in range(PROCESSES):
			args = []
			args.append("array")
			args.append(conc[(ranges[i][0] - 1):(ranges[i][1]) , :])
			task_queue[i].put(args)
		#
		# Retrieve reacted solution compositions
		#
		for i in range(PROCESSES):
			segment = done_queue[i].get()
			k = 0
			for j in xrange(ranges[i][0] - 1, ranges[i][1]):
				conc[j] = segment[k]
				k += 1
		if (t % 10 == 0):
			print "Finished step %d." % t
		outflow[t] = conc[CELLS - 1]

	#
	# Finalize processes
	#
	for i in range(PROCESSES):
		task_queue[i].put(None)
	print 'Time: ', time.clock() - start

	#
	# Plot results
	#
	
	Ca = outflow[:,3]
	Cl = outflow[:,4]
	K  = outflow[:,5]
	N  = outflow[:,6]
	Na = outflow[:,7]
	x = []
	for i in xrange(1, len(outflow) + 1):
		a = double(i)
		x.append(double(i)/double(CELLS))

	plt.plot(x, Ca, 'r', x, Cl, 'b', x, K, 'g', x, N, 'y', x, Na, 'm')
	plt.legend(('Ca','Cl','K','N','Na'), loc = (0.8, 0.5))
	plt.ylabel('MILLIMOLES PER KILOGRAM WATER')
	plt.xlabel('PORE VOLUME')
	plt.show()
    

		
	print "Finished simulation\n"
	
	#
	# Process done_queue
	#
	
	print 'Time: ', time.clock() - start	
