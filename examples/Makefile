FILES=ex1 \
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
	ex12a \
	ex13a \
	ex13b \
	ex13c \
	ex14 \
	ex15 \
	ex16 \
	ex17 \
	ex18

DATAFILE=../database/phreeqc.dat
DATAFILE1=../database/wateq4f.dat
PHREEQC=../bin/phreeqc

all:
	for FILE in $(FILES); do \
		echo ; \
		echo "====>" $$FILE; \
		echo ; \
		rm -f $$FILE.out $$FILE.log ; \
		if [ $$FILE = ex14 ]; \
			then time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE1); \
			elif [ $$FILE = ex15 ]; \
			then time $(PHREEQC) $$FILE $$FILE.out ex15.dat; \
			else time $(PHREEQC) $$FILE $$FILE.out $(DATAFILE); fi; \
		if [ -f logfile ]; then mv logfile $$FILE.logfile; fi; \
		mv phreeqc.log $$FILE.log; \
		done
	echo "Now diffing files"
#	for FILE in ex2.sel ex5.sel ex6A-B.sel ex6C.sel ex7.sel ex8.sel ex9.sel ex10.sel ex11adv.sel ex11trn.sel ex12.sel ex12a.sel ex13a.sel ex13b.sel ex13c.sel ex14.sel ex15.sel; do \
#		echo "====>" $$FILE; \
#			rcsdiff -ibw $$FILE; \
#		done
	svn diff .

clean:
	rm *.out *.log *.sel
