# Locations to save compressed tar file for distribution
DIST_DIR=~dlpark/temp
DEBUG_DIR=phreeqc_debug
DEBUG_EXE=$(TOPDIR)/src/phreeqc
VERSION=2.11
REVISION=`svnversion .`
ROOTNAME=$(PROGRAM)-$(VERSION)-$(REVISION)

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
	database/minteq.v4.dat \
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

all_dist: ../doc/RELEASE.TXT linux sun mytest dist

../doc/RELEASE.TXT: revisions
	cp revisions ../doc/RELEASE.TXT

sun: clean_sun 
	cp Makefile *.c *.h $(TOPDIR)/sun
	ssh u450rcolkr "cd ~dlpark/programs/phreeqc/sun; make -j 4 EXE=/z/parkplace/home/dlpark/programs/phreeqc/bin/$(PROGRAM).sun"
	cp ../bin/phreeqc ../bin/phreeqc.linux
	cp ../bin/phreeqc.sun ../bin/phreeqc
	ssh u450rcolkr "cd ~dlpark/programs/phreeqc/test.sun; ./clean.sh; ./test.sh"
	mv ../bin/phreeqc.linux ../bin/phreeqc

clean_sun:
	rm -f $(TOPDIR)/sun/*.c $(TOPDIR)/sun/*.h $(TOPDIR)/sun/*.o $(TOPDIR)/bin/phreeqc.sun

linux: clean_linux 
	make -k 
	cd ../examples; make clean; make >& make.out
	cd ../test; ./clean.sh; ./test.sh

clean_linux:
	rm -f $(TOPDIR)/bin/$(PROGRAM) $(TOPDIR)/src/*.o 

dist: svn_update clean_dist dist_linux dist_sun dist_source win
	echo Done with distribution.

dist_linux: 
	cd ..; rm -f $(PROGRAM).tar
	cd ..; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd ..; tar -rf $(PROGRAM).tar bin/$(PROGRAM)
	cd ..; rm -rf $(PROGRAM)-$(VERSION)
	cd ..; mkdir $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd ..; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm -f $(PROGRAM).tar
	cd ..; for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)-$(VERSION)/examples; done;
	cd ..; tar -czf $(PROGRAM).Linux.tar.gz $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).Linux.tar.gz $(DIST_DIR)/$(ROOTNAME).Linux.tar.gz
	cd ..; echo $(ROOTNAME).Linux.tar.gz saved in $(DIST_DIR).
	cd ..; rm -rf $(PROGRAM)-$(VERSION)

dist_sun: 
	cd ..; rm -f $(PROGRAM).tar
	cd ..; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd ..; rm -rf $(PROGRAM)-$(VERSION)
	cd ..; mkdir $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd ..; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
	cd ..; cp bin/$(PROGRAM).sun $(PROGRAM)-$(VERSION)/bin/$(PROGRAM)
	cd ..; for FILE in $(OUTPUT_FILES); do cp test.sun/$$FILE $(PROGRAM)-$(VERSION)/examples; done;
	cd ..; tar -czf $(PROGRAM).SunOS.tar.gz $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).SunOS.tar.gz $(DIST_DIR)/$(ROOTNAME).SunOS.tar.gz
	cd ..; echo $(ROOTNAME).SunOS.tar.gz saved in $(DIST_DIR).
	cd ..; rm -rf $(PROGRAM)$(VERSION)

dist_source: 
	cd ..; rm -f $(PROGRAM).tar
	cd ..; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd ..; rm -rf $(PROGRAM)-$(VERSION)
	cd ..; mkdir $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd ..; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
	cd ..; for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)-$(VERSION)/examples; done;
	cd ..; tar -czf $(PROGRAM).source.tar.gz $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).source.tar.gz $(DIST_DIR)/$(ROOTNAME).source.tar.gz
	cd ..; echo $(ROOTNAME).source.tar.gz saved in $(DIST_DIR).
	cd ..; rm -rf $(PROGRAM)-$(VERSION)

clean_dist:
	rm -f $(DIST_DIR)/README.TXT
	rm -f $(DIST_DIR)/phreeqc-$(VERSION)-*gz
	echo Removed README.TXT and tar.gz files for $(PROGRAM) from $(DIST_DIR).

win: 
	cd ..; rm -f $(PROGRAM).tar
	cd ..; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd ..; rm -rf $(PROGRAM)-$(VERSION)
	cd ..; mkdir $(PROGRAM)-$(VERSION)
	cd ..; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd ../$(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm $(PROGRAM).tar
# remove example output
	cd ..; rm -f $(PROGRAM)-$(VERSION)/examples/*.out $(PROGRAM)-$(VERSION)/examples/*.sel
# copy Windows files
o	cd ..; cp database/*.dat $(PROGRAM)-$(VERSION)/
# remove database, bin directories
	cd ..; rm -rf $(PROGRAM)-$(VERSION)/database $(PROGRAM)-$(VERSION)/bin
# copy bat file
	cd ..; cp win/phreeqc.bat $(PROGRAM)-$(VERSION)/
	cd ..; rm -f $(PROGRAM)-$(VERSION)/test/*
	cd ..; cp win/clean.bat win/check.bat win/test.bat $(PROGRAM)-$(VERSION)/test
	cd ..; rm -f $(PROGRAM)-$(VERSION)/README.TXT
	cd ..; cp win/README.TXT $(PROGRAM)-$(VERSION)/
	cd ..; rm -f $(PROGRAM).Windows.tar.gz
	cd ..; tar -czf $(PROGRAM).$(VERSION).$(REVISION).Windows.tar.gz $(PROGRAM)-$(VERSION)
	cd ..; echo $(PROGRAM).$(VERSION).$(REVISION).Windows.tar.gz created.
#	cd ..; rm -rf $(PROGRAM)-$(VERSION)

mytest:
	cd ../mytest; make clean; make >& make.out 

debug: 
	mkdir -p $(DEBUG_DIR)
	cd $(DEBUG_DIR); make -f $(TOPDIR)/src/debug.mk; make -f $(TOPDIR)/src/Makefile CCFLAGS="-Wall -ansi -g" EXE=$(DEBUG_EXE)

test_dist: test_linux test_sunos test_source

test_linux:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).Linux
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.Linux.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).Linux
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/test; ./test.sh
	rm -f $(DIST_DIR)/phreeqc-$(VERSION).Linux/bin/phreeqc
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/src; make -k
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/test; ./clean.sh; ./test.sh

test_sunos:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).SunOS
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.SunOS.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).SunOS
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).SunOS/test; ./test.sh"
	ssh u450rcolkr "rm -f $(DIST_DIR)/phreeqc-$(VERSION).SunOS/bin/phreeqc"
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).SunOS/src; make -k"
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).SunOS/test; ./clean.sh; ./test.sh"

test_source:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).source
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.source.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).source
	cd $(DIST_DIR)/phreeqc-$(VERSION).source/src; make -k
	cd $(DIST_DIR)/phreeqc-$(VERSION).source/test; ./test.sh

web:
	cp $(DIST_DIR)/phreeqc-$(VERSION)*.tar.gz /var/anonymous/ftp/dlpark/geochem/unix/phreeqc
	cp $(TOPDIR)/doc/README.TXT /var/anonymous/ftp/dlpark/geochem/unix/phreeqc/README.TXT 
	cp $(TOPDIR)/doc/README.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/README.Unix.TXT
	cp $(TOPDIR)/win/README.TXT /var/anonymous/ftp/dlpark/geochem/pc/phreeqc/README.TXT 
	cp $(TOPDIR)/win/README.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/README.Win.TXT
	cp $(TOPDIR)/doc/phreeqc.txt /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/phreeqc.txt
	cp $(TOPDIR)/doc/RELEASE.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/RELEASE.TXT

svn_update:
	cd ..; svn update .

