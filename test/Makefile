DATA=../examples
PHREEQCDAT=../database/phreeqc.dat
WATEQ4FDAT=../database/wateq4f.dat
PHREEQC=../bin/phreeqc

all: ex1.out ex2.out ex3.out ex4.out ex5.out ex6.out ex7.out ex8.out ex9.out \
	ex10.out ex11.out ex12.out ex12a.out ex13a.out ex13b.out ex13c.out ex14.out \
	ex15.out ex16.out ex17.out ex18.out 

ex1.out: $(DATA)/ex1 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex1 ex1.out $(PHREEQCDAT)

ex2.out: $(DATA)/ex2 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex2 ex2.out $(PHREEQCDAT)

ex3.out: $(DATA)/ex3 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex3 ex3.out $(PHREEQCDAT)

ex4.out: $(DATA)/ex4 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex4 ex4.out $(PHREEQCDAT)

ex5.out: $(DATA)/ex5 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex5 ex5.out $(PHREEQCDAT)

ex6.out: $(DATA)/ex6 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex6 ex6.out $(PHREEQCDAT)

ex7.out: $(DATA)/ex7 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex7 ex7.out $(PHREEQCDAT)

ex8.out: $(DATA)/ex8 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex8 ex8.out $(PHREEQCDAT)

ex9.out: $(DATA)/ex9 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex9 ex9.out $(PHREEQCDAT)

ex10.out: $(DATA)/ex10 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex10 ex10.out $(PHREEQCDAT)

ex11.out: $(DATA)/ex11 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex11 ex11.out $(PHREEQCDAT)

ex12.out: $(DATA)/ex12 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex12 ex12.out $(PHREEQCDAT)

ex12a.out: $(DATA)/ex12a $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex12a ex12a.out $(PHREEQCDAT)

ex13a.out: $(DATA)/ex13a $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex13a ex13a.out $(PHREEQCDAT)

ex13b.out: $(DATA)/ex13b $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex13b ex13b.out $(PHREEQCDAT)

ex13c.out: $(DATA)/ex13c $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex13c ex13c.out $(PHREEQCDAT)

ex14.out: $(DATA)/ex14 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex14 ex14.out $(WATEQ4FDAT)

ex15.out: $(DATA)/ex15 $(PHREEQC) $(DATA)/ex15.dat
	$(PHREEQC) $(DATA)/ex15 ex15.out $(DATA)/ex15.dat

ex16.out: $(DATA)/ex16 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex16 ex16.out $(PHREEQCDAT)

ex17.out: $(DATA)/ex17 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex17 ex17.out $(PHREEQCDAT)

ex18.out: $(DATA)/ex18 $(PHREEQC) $(PHREEQCDAT)
	$(PHREEQC) $(DATA)/ex18 ex18.out $(PHREEQCDAT)

diff:
	svn diff .

clean:
	rm -f *.out *.log *.sel
