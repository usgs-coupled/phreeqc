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

diff:
	svn diff .

clean:
	rm -f *.out *.log *.sel
