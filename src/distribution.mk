
# Locations to save compressed tar file for distribution
DIST_DIR=~dlpark/temp
VERSION=2.11

# list of files for distribution
FILES=  \
	src/Makefile \
	src/advection.c \
	src/basic.c \
	src/basicsubs.c \
	src/cl1.c \
	src/cvdense.c \
	src/cvode.c \
	src/dense.c \
	src/input.c \
	src/integrate.c \
	src/inverse.c \
	src/isotopes.c \
	src/kinetics.c \
	src/main.c \
	src/mainsubs.c \
	src/model.c \
	src/nvector.c \
	src/nvector_serial.c \
	src/output.c \
	src/p2clib.c \
	src/parse.c \
	src/phqalloc.c \
	src/phreeqc_files.c \
	src/prep.c \
	src/print.c \
	src/read.c \
	src/readtr.c \
	src/smalldense.c \
	src/spread.c \
	src/step.c \
	src/structures.c \
	src/sundialsmath.c \
	src/tally.c \
	src/tidy.c \
	src/transport.c \
	src/utilities.c \
	src/cvdense.h \
	src/cvode.h \
	src/dense.h \
	src/global.h \
	src/input.h \
	src/kinetics.h \
	src/nvector.h \
	src/nvector_serial.h \
	src/output.h \
	src/p2c.h \
	src/phqalloc.h \
	src/phrqproto.h \
	src/phrqtype.h \
	src/smalldense.h \
	src/sundialsmath.h \
	src/sundialstypes.h \
	database/llnl.dat \
	database/minteq.dat \
	database/phreeqc.dat \
	database/wateq4f.dat \
	database/iso.dat \
	examples/ex1 examples/ex1.out \
	examples/ex2 examples/ex2.out examples/ex2.sel \
	examples/ex3 examples/ex3.out \
	examples/ex4 examples/ex4.out \
	examples/ex5 examples/ex5.out examples/ex5.sel \
	examples/ex6 examples/ex6.out examples/ex6A-B.sel examples/ex6C.sel \
	examples/ex7 examples/ex7.out examples/ex7.sel \
	examples/ex8 examples/ex8.out examples/ex8.sel \
	examples/ex9 examples/ex9.out examples/ex9.sel \
	examples/ex10 examples/ex10.out examples/ex10.sel \
	examples/ex11 examples/ex11.out examples/ex11adv.sel examples/ex11trn.sel \
	examples/ex12 examples/ex12.out examples/ex12.sel \
	examples/ex13a examples/ex13a.out examples/ex13a.sel \
	examples/ex13b examples/ex13b.out examples/ex13b.sel \
	examples/ex13c examples/ex13c.out examples/ex13c.sel \
	examples/ex14 examples/ex14.out examples/ex14.sel \
	examples/ex15 examples/ex15.dat examples/ex15.out examples/ex15.sel \
	examples/ex16 examples/ex16.out \
	examples/ex17 examples/ex17.out \
	examples/ex18 examples/ex18.out \
	doc/NOTICE.TXT \
	doc/README.TXT \
	doc/RELEASE.TXT \
	doc/manual.pdf \
	doc/wrir02-4172.pdf \
	doc/phreeqc.txt \
	bin/phreeqc.orig \
	test/test.sh \
	test/clean.sh \
	test/check.sh 

OUTPUT_FILES = \
	ex1.out \
	ex2.out ex2.sel \
	ex3.out \
	ex4.out \
	ex5.out ex5.sel \
	ex6.out ex6A-B.sel ex6C.sel \
	ex7.out ex7.sel \
	ex8.out ex8.sel \
	ex9.out ex9.sel \
	ex10.out ex10.sel \
	ex11.out ex11adv.sel ex11trn.sel \
	ex12.out ex12.sel \
	ex13a.out ex13a.sel \
	ex13b.out ex13b.sel \
	ex13c.out ex13c.sel \
	ex14.out ex14.sel \
	ex15.out ex15.sel \
	ex16.out \
	ex17.out \
	ex18.out 

all_dist: linux sun dist win

sun: clean_sun 
	cp Makefile *.c *.h $(TOPDIR)/sun
	ssh u450rcolkr "cd ~dlpark/programs/phreeqc/sun; make -j 4 EXE=$(TOPDIR)/bin/$(PROGRAM).sun"
	cp ../bin/phreeqc ../bin/phreeqc.linux
	cp ../bin/phreeqc.sun ../bin/phreeqc
	ssh u450rcolkr "cd ~dlpark/programs/phreeqc/test.sun; ./clean.sh; ./test.sh"
	mv ../bin/phreeqc.linux ../bin/phreeqc

clean_sun:
	rm -f $(TOPDIR)/sun/*.c $(TOPDIR)/sun/*.h $(TOPDIR)/sun/*.o $(TOPDIR)/bin/phreeqc.sun

linux: clean_linux 
	ssh srv2rcolkr 'cd ~dlpark/programs/phreeqc/src; make -f Makefile EXE=$(TOPDIR)/bin/$(PROGRAM) CCFLAGS="-O3 -Wall -ansi -pedantic"'
	cd ../examples; make clean; make >& make.out
	cd ../test; ./clean.sh; ./test.sh

clean_linux:
	rm -f $(TOPDIR)/bin/$(PROGRAM) $(TOPDIR)/src/*.o 

dist: 
	error=0; for FILE in $(OUTPUT_FILES); do \
	   if [ ! -f ../test/$$FILE ]; then echo test/$$FILE is missing.; error=1; fi; \
	   if [ ! -f ../test.sun/$$FILE ]; then echo test.sun/$$FILE is missing.; error=1; fi; \
	   done; if [ $$error -eq 1 ]; then echo Stopping because of missing files; exit 4; fi;
	rm -f Makefile.internaldist; cp Makefile Makefile.internaldist
	make -C $(TOPDIR) -f src/Makefile.internaldist internaldist 
	rm -f Makefile.internaldist
	echo Done with distribution.
	echo

internaldist: clean_dist dist_linux dist_sun dist_unix
	cp $(DIST_DIR)/$(PROGRAM)$(VERSION).Linux.tar.gz ~/programs/unix/phreeqc/ftp
	cp $(DIST_DIR)/$(PROGRAM)$(VERSION).SunOS.tar.gz ~/programs/unix/phreeqc/ftp
	cp $(DIST_DIR)/$(PROGRAM)$(VERSION).source.tar.gz ~/programs/unix/phreeqc/ftp

dist_linux: $(FILES)
	rm -f $(PROGRAM).tar
	for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	tar -rf $(PROGRAM).tar bin/$(PROGRAM)
	rm -rf $(PROGRAM)$(VERSION)
	mkdir $(PROGRAM)$(VERSION)
	mv $(PROGRAM).tar $(PROGRAM)$(VERSION)
	cd $(PROGRAM)$(VERSION); tar -xf $(PROGRAM).tar; rm -f $(PROGRAM).tar
	for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)$(VERSION)/examples; done;
	tar -czf $(PROGRAM).Linux.tar.gz $(PROGRAM)$(VERSION)
	mv $(PROGRAM).Linux.tar.gz $(DIST_DIR)/$(PROGRAM)$(VERSION).Linux.tar.gz
	echo $(PROGRAM)$(VERSION).Linux.tar.gz saved in $(DIST_DIR).
	rm -rf $(PROGRAM)$(VERSION)

dist_sun: $(FILES)
	rm -f $(PROGRAM).tar
	for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	rm -rf $(PROGRAM)$(VERSION)
	mkdir $(PROGRAM)$(VERSION)
	mv $(PROGRAM).tar $(PROGRAM)$(VERSION)
	cd $(PROGRAM)$(VERSION); tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
	cp bin/$(PROGRAM).sun $(PROGRAM)$(VERSION)/bin/$(PROGRAM)
	rm -f $(PROGRAM)$(VERSION)/src/Makefile
	cp src/Makefile $(PROGRAM)$(VERSION)/src/Makefile
	for FILE in $(OUTPUT_FILES); do cp test.sun/$$FILE $(PROGRAM)$(VERSION)/examples; done;
	tar -czf $(PROGRAM).SunOS.tar.gz $(PROGRAM)$(VERSION)
	mv $(PROGRAM).SunOS.tar.gz $(DIST_DIR)/$(PROGRAM)$(VERSION).SunOS.tar.gz
	echo $(PROGRAM)$(VERSION).SunOS.tar.gz saved in $(DIST_DIR).
	rm -rf $(PROGRAM)$(VERSION)

dist_unix: $(FILES)
	rm -f $(PROGRAM).tar
	for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	rm -rf $(PROGRAM)$(VERSION)
	mkdir $(PROGRAM)$(VERSION)
	mv $(PROGRAM).tar $(PROGRAM)$(VERSION)
	cd $(PROGRAM)$(VERSION); tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
	for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)$(VERSION)/examples; done;
	tar -czf $(PROGRAM).source.tar.gz $(PROGRAM)$(VERSION)
	mv $(PROGRAM).source.tar.gz $(DIST_DIR)/$(PROGRAM)$(VERSION).source.tar.gz
	echo $(PROGRAM)$(VERSION).source.tar.gz saved in $(DIST_DIR).
	rm -rf $(PROGRAM)$(VERSION)

clean_dist:
	rm -f $(DIST_DIR)/README.TXT
	rm -f $(DIST_DIR)/$(PROGRAM)*gz
	echo Removed README.TXT and tar.gz files for $(PROGRAM) from $(DIST_DIR).

win: 
	rm -f $(PROGRAM).tar
	cd ..; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd ..; rm -rf $(PROGRAM).win
	cd ..; mkdir $(PROGRAM).win
	cd ..; mv $(PROGRAM).tar $(PROGRAM).win
	cd ../$(PROGRAM).win; tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
# remove example output
	cd ..; rm -f $(PROGRAM).win/examples/*.out $(PROGRAM).win/examples/*.sel
# remove database, bin directories
	cd ..; rm -rf $(PROGRAM).win/database $(PROGRAM).win/bin
# copy Windows files
	cd ..; cp database/phreeqc.dat database/llnl.dat database/wateq4f.dat database/minteq.dat database/iso.dat $(PROGRAM).win/
	cd ..; cp win/phreeqc.bat $(PROGRAM).win/
	cd ..; rm -f $(PROGRAM).win/test/*
	cd ..; cp win/clean.bat win/check.bat win/test.bat $(PROGRAM).win/test
	cd ..; rm -f $(PROGRAM).win/README.TXT
	cd ..; cp win/README.TXT $(PROGRAM).win/
	cd ..; rm -f $(PROGRAM).Windows.tar.gz
	cd ..; tar -czf $(PROGRAM).Windows.tar.gz $(PROGRAM).win
	cd ..; echo $(PROGRAM).Windows.tar.gz created.
	cd ..; rm -rf $(PROGRAM).win
