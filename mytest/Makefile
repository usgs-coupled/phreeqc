FILES=alkalinity \
	aj1 \
	all_llnl \
	as.ex \
	ball \
	brown \
	co2 \
	comments \
	dissolve_only \
	exch_related\
	exch_kin_related\
	fixed_volume_gas \
	fluorite \
	gas \
	hard.tests\
	jgb2col\
	jgbcol\
	joe\
	krest \
	lam7 \
	lam10 \
	llnl \
	mahoney \
	mu \
	np.test\
	org_gas\
	oxygen \
	phantom_surface \
	phreeqc.test\
	prob8\
	saver\
	saunders \
	solid_solution \
	solution_s \
	spread_solution \
	iso_example \
	surf3\
	surf_related\
	surf_kin_related\
	ss_r \
	t922 \
	transport\
	transport.dl\
	user_print \
	using \
	zn.test \
	ex1 \
	ex2 \
	ex3 \
	ex4 \
	ex5 \
	ex6 \
	ex7 \
	ex8 \
	ex9 \
	ex10 \
	ex11 \
	ex12 \
	colt1 \
	colt2 \
	colt3 \
	colt3_mi \
	colt4 \
	colt5 \
#	denitrify


DATAFILE=../database/phreeqc.dat
DATAFILE1=../database/wateq4f.dat
DATAFILE2=tony.dat
DATAFILE3=../database/iso.dat
PHREEQC=../bin/phreeqc

all:
	for FILE in $(FILES); do \
		echo ; \
		echo "====>" $$FILE; \
		echo ; \
		rm -f $$FILE.out $$FILE.log ; \
		if [ $$FILE = zn.test -o $$FILE = hard.tests -o $$FILE = ex10 -o $$FILE = brown \
			-o $$FILE = lam7 -o $$FILE = lam10 -o $$FILE = saunders -o $$FILE = as.ex ]; \
			then time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE1); \
			elif [ $$FILE = mahoney ]; \
			then time $(PHREEQC) $$FILE $$FILE.out mahoney.dat; \
			elif [ $$FILE = iso_example -o $$FILE = co2 ]; \
			then time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE3); \
			elif [ $$FILE = colt1 -o $$FILE = colt2 -o $$FILE = colt3 -o \
				$$FILE = colt4 -o $$FILE = colt5 ]; \
			then time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE2); \
			else time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE); fi; \
		if [ -f logfile ]; then mv logfile $$FILE.logfile; fi; \
		mv phreeqc.log $$FILE.log; \
		echo "====>" $$FILE.out; \
			rcsdiff -bw $$FILE.out; \
		if [ $$FILE = using ]; then rcsdiff -bw $$FILE.sel; fi; \
		echo "====>" $$FILE.log; \
		done

#			rcsdiff -bw $$FILE.log; \

#$(FILE).out: $(FILE)
#		echo ; \
#		echo "====>" $$FILE; \
#		echo ; \
#		if [ $$FILE = zn.test -o $$FILE = hard.tests ]; \
#			then DATAFILE=../data/phreeqc.dat; \
#			else DATAFILE=../phreeqc.dat; fi; \
#		
#		rm -f $$FILE.out $$FILE.log ; \
#		/home/dlpark/bin/phreeqc $$FILE $$DATAFILE; \
#		if [ -f logfile ]; then mv logfile $$FILE.logfile; fi; \
#		mv phreeqc.log $$FILE.log; \
#		echo "====>" $$FILE.out; \
#			rcsdiff $$FILE.out; \
#		echo "====>" $$FILE.log; \
#			rcsdiff $$FILE.log; \




clean:
	rm *.out *.log *.pun *.prn *.flat *.sel
