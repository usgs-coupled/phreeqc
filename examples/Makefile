PHREEQCDAT=../database/phreeqc.dat
WATEQ4FDAT=../database/wateq4f.dat
PHREEQC=../bin/phreeqc

all: ex1.out ex2.out ex3.out ex4.out ex5.out ex6.out ex7.out ex8.out ex9.out \
	ex10.out ex11.out ex12.out ex12a.out ex13a.out ex13b.out ex13c.out ex14.out \
	ex15.out ex16.out ex17.out ex18.out 

ex1.out: ex1 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex1 ex1.out $(PHREEQCDAT)

ex2.out: ex2 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex2 ex2.out $(PHREEQCDAT)

ex3.out: ex3 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex3 ex3.out $(PHREEQCDAT)

ex4.out: ex4 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex4 ex4.out $(PHREEQCDAT)

ex5.out: ex5 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex5 ex5.out $(PHREEQCDAT)

ex6.out: ex6 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex6 ex6.out $(PHREEQCDAT)

ex7.out: ex7 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex7 ex7.out $(PHREEQCDAT)

ex8.out: ex8 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex8 ex8.out $(PHREEQCDAT)

ex9.out: ex9 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex9 ex9.out $(PHREEQCDAT)

ex10.out: ex10 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex10 ex10.out $(PHREEQCDAT)

ex11.out: ex11 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex11 ex11.out $(PHREEQCDAT)

ex12.out: ex12 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex12 ex12.out $(PHREEQCDAT)

ex12a.out: ex12a $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex12a ex12a.out $(PHREEQCDAT)

ex13a.out: ex13a $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex13a ex13a.out $(PHREEQCDAT)

ex13b.out: ex13b $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex13b ex13b.out $(PHREEQCDAT)

ex13c.out: ex13c $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex13c ex13c.out $(PHREEQCDAT)

ex14.out: ex14 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex14 ex14.out $(WATEQ4FDAT)

ex15.out: ex15 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex15 ex15.out ex15.dat

ex16.out: ex16 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex16 ex16.out $(PHREEQCDAT)

ex17.out: ex17 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex17 ex17.out $(PHREEQCDAT)

ex18.out: ex18 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) ex18 ex18.out $(PHREEQCDAT)

diff:
	svn diff .

clean:
	rm -f *.out *.log *.sel
