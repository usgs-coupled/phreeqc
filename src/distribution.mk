# Locations to save compressed tar file for distribution
EXPORT=$(TOPDIR)/src/phreeqc_export
EXPORT_DIR=$(HOME)/../..$(EXPORT)
WIN_DIR=$(TOPDIR)/win
DIST_DIR=$(EXPORT_DIR)
DEBUG_DIR=phreeqc_debug
DEBUG_EXE=$(TOPDIR)/src/phreeqc
VERSION=2.11
VER_DATE=February 7, 2005
REVISION=$(shell svnversion .)
ROOTNAME=$(PROGRAM)-$(VERSION)-$(REVISION)
TEXTCP=textcp DOS

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

# make sure program is compiles, run examples and mytest
output_files: all examples mytest

examples:
	cd ../examples; make clean; make >& make.out 

mytest:
	cd ../mytest; make clean; make >& make.out 

all_dist:  clean_dist linux sun source win

test_dist: linux_test source_test sun_test

#
#Linux
#
linux: linux_export linux_clean linux_sed_files linux_compile linux_output linux_dist source_dist

linux_export:
	mkdir -p $(EXPORT_DIR)
	rm -rf $(EXPORT_DIR)/Linux
	svn export .. $(EXPORT_DIR)/Linux

linux_clean:
	rm -f $(EXPORT_DIR/Linux/bin/$(PROGRAM) $(EXPORT_DIR/Linux/src/*.o 

linux_sed_files:
	sed -e "s/VERSION/$(VERSION)/" \
	    -e "s/VER_DATE/$(VER_DATE)/" \
	    -e "s/VERSION_DATE/$(VERSION)/" \
	    -e "s/REVISION/$(REVISION)/" < $(EXPORT_DIR)/Linux/src/revisions > $(EXPORT_DIR)/Linux/doc/RELEASE.TXT
	for FILE in $(EXPORT_DIR)/Linux/doc/README.TXT; do \
		sed  -e "s/VER_DATE/$(VER_DATE)/" -e "s/VERSION/$(VERSION)/" -e "s/REVISION/$(REVISION)/" < $$FILE > t; mv t $$FILE; done

linux_compile:
	make -C $(EXPORT_DIR)/Linux/src

linux_output:
	cd $(EXPORT_DIR)/Linux/examples; make clean; make >& make.out
	cd $(EXPORT_DIR)/Linux/test; ./clean.sh; ./test.sh

linux_dist: 
	cd $(EXPORT_DIR)/Linux; rm -f $(PROGRAM).tar
	cd $(EXPORT_DIR)/Linux; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd $(EXPORT_DIR)/Linux; tar -rf $(PROGRAM).tar bin/$(PROGRAM)
	cd $(EXPORT_DIR)/Linux; rm -rf $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Linux; mkdir $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Linux; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Linux; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm -f $(PROGRAM).tar
	cd $(EXPORT_DIR)/Linux; for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)-$(VERSION)/examples; done;
	cd $(EXPORT_DIR)/Linux; tar -czf $(PROGRAM).Linux.tar.gz $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Linux; mv $(PROGRAM).Linux.tar.gz $(DIST_DIR)/$(ROOTNAME).Linux.tar.gz
	cd $(EXPORT_DIR)/Linux; echo $(ROOTNAME).Linux.tar.gz saved in $(DIST_DIR).
	cd $(EXPORT_DIR)/Linux; rm -rf $(PROGRAM)-$(VERSION)

source_dist:
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

linux_test:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).Linux
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.Linux.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).Linux
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/test; ./test.sh
	rm -f $(DIST_DIR)/phreeqc-$(VERSION).Linux/bin/phreeqc
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/src; make -k
	cd $(DIST_DIR)/phreeqc-$(VERSION).Linux/test; ./clean.sh; ./test.sh

source_test:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).source
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.source.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).source
	cd $(DIST_DIR)/phreeqc-$(VERSION).source/src; make -k
	cd $(DIST_DIR)/phreeqc-$(VERSION).source/test; ./test.sh

#
#Sun
#
sun: sun_export sun_clean sun_sed_files sun_compile sun_output sun_dist 

sun_export:
	mkdir -p $(EXPORT_DIR)
	rm -rf $(EXPORT_DIR)/Sun
	svn export .. $(EXPORT_DIR)/Sun

sun_clean:
	rm -f $(EXPORT_DIR/Sun/bin/$(PROGRAM) $(EXPORT_DIR/Sun/src/*.o 

sun_sed_files:
	sed -e "s/VERSION/$(VERSION)/" \
	    -e "s/VER_DATE/$(VER_DATE)/" \
	    -e "s/VERSION_DATE/$(VERSION)/" \
	    -e "s/REVISION/$(REVISION)/" < $(EXPORT_DIR)/Sun/src/revisions > $(EXPORT_DIR)/Sun/doc/RELEASE.TXT
	for FILE in $(EXPORT_DIR)/Sun/doc/README.TXT; do \
		sed -e "s/VER_DATE/$(VER_DATE)/" -e "s/VERSION/$(VERSION)/" -e "s/REVISION/$(REVISION)/" < $$FILE > t; mv t $$FILE; done

sun_compile:
	ssh u450rcolkr "make -C $(EXPORT_DIR)/Sun/src"

sun_output:
	ssh u450rcolkr "cd $(EXPORT_DIR)/Sun/examples; make clean; make >& make.out"
	ssh u450rcolkr "cd $(EXPORT_DIR)/Sun/test; ./clean.sh; ./test.sh"

sun_dist: 
	cd $(EXPORT_DIR)/Sun; rm -f $(PROGRAM).tar
	cd $(EXPORT_DIR)/Sun; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd $(EXPORT_DIR)/Sun; tar -rf $(PROGRAM).tar bin/$(PROGRAM)
	cd $(EXPORT_DIR)/Sun; rm -rf $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Sun; mkdir $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Sun; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Sun; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm -f $(PROGRAM).tar
	cd $(EXPORT_DIR)/Sun; for FILE in $(OUTPUT_FILES); do cp test/$$FILE $(PROGRAM)-$(VERSION)/examples; done;
	cd $(EXPORT_DIR)/Sun; tar -czf $(PROGRAM).Sun.tar.gz $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Sun; mv $(PROGRAM).Sun.tar.gz $(DIST_DIR)/$(ROOTNAME).Sun.tar.gz
	cd $(EXPORT_DIR)/Sun; echo $(ROOTNAME).Sun.tar.gz saved in $(DIST_DIR).
	cd $(EXPORT_DIR)/Sun; rm -rf $(PROGRAM)-$(VERSION)

sun_test:
	rm -rf $(DIST_DIR)/phreeqc-$(VERSION).Sun
	cd $(DIST_DIR); tar -xzf phreeqc-$(VERSION)-*.Sun.tar.gz; mv phreeqc-$(VERSION) phreeqc-$(VERSION).Sun
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).Sun/test; ./test.sh"
	rm -f $(DIST_DIR)/phreeqc-$(VERSION).Sun/bin/phreeqc
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).Sun/src; make -k"
	ssh u450rcolkr "cd $(DIST_DIR)/phreeqc-$(VERSION).Sun/test; ./clean.sh; ./test.sh"

clean_dist:
	rm -rf $(EXPORT_DIR)

#
#Win
#
win: win_export win_sed_files win_dist 

win_export:
	mkdir -p $(EXPORT_DIR)
	rm -rf $(EXPORT_DIR)/Win
	svn export .. $(EXPORT_DIR)/Win

win_sed_files:
	sed -e "s/VERSION/$(VERSION)/" \
	    -e "s/VER_DATE/$(VER_DATE)/" \
	    -e "s/REVISION/$(REVISION)/" < $(EXPORT_DIR)/Win/src/revisions > $(EXPORT_DIR)/Win/doc/RELEASE.TXT
	sed -e "s/VERSION/$(VERSION)/" \
	    -e "s/VER_DATE/$(VER_DATE)/" \
	    -e "s/REVISION/$(REVISION)/" < $(WIN_DIR)/README.TXT > $(EXPORT_DIR)/Win/doc/README.TXT

win_dist: 
	cd $(EXPORT_DIR)/Win; rm -f $(PROGRAM).tar
# Translate cr/lf
	cd $(EXPORT_DIR)/Win; for FILE in $(FILES); do \
		if [ $$FILE = doc/manual.pdf -o $$FILE = doc/wrir02-4172.pdf ]; then cp $$FILE t; mv t $$FILE; \
		else $(TEXTCP) $$FILE t; mv t $$FILE; fi; done
	cd $(EXPORT_DIR)/Win; for FILE in $(FILES); do tar -rf $(PROGRAM).tar $$FILE; done
	cd $(EXPORT_DIR)/Win; rm -rf $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Win; mkdir $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Win; mv $(PROGRAM).tar $(PROGRAM)-$(VERSION)
	cd $(EXPORT_DIR)/Win; cd $(PROGRAM)-$(VERSION); tar -xf $(PROGRAM).tar; rm -f $(PROGRAM).tar
# remove example output
	cd $(EXPORT_DIR)/Win; rm -f $(PROGRAM)-$(VERSION)/examples/*.out $(PROGRAM)-$(VERSION)/examples/*.sel
# remove bin directory
	cd $(EXPORT_DIR)/Win; rm -rf $(PROGRAM)-$(VERSION)/bin
# remove test directory files
	cd $(EXPORT_DIR)/Win; rm -f $(PROGRAM)-$(VERSION)/test/*
	cd $(EXPORT_DIR)/Win; $(TEXTCP) $(WIN_DIR)/clean.bat $(PROGRAM)-$(VERSION)/test/clean.bat
	cd $(EXPORT_DIR)/Win; $(TEXTCP) $(WIN_DIR)/check.bat $(PROGRAM)-$(VERSION)/test/check.bat
	cd $(EXPORT_DIR)/Win; $(TEXTCP) $(WIN_DIR)/test.bat $(PROGRAM)-$(VERSION)/test/test.bat
# copy bat file
	cd $(EXPORT_DIR)/Win; $(TEXTCP) $(WIN_DIR)/phreeqc.bat $(PROGRAM)-$(VERSION)/phreeqc.bat
	cd $(EXPORT_DIR); rm -f $(PROGRAM).Windows.tar.gz
	cd $(EXPORT_DIR)/Win/$(PROGRAM)-$(VERSION); tar -czf $(PROGRAM).Windows.tar.gz .
	cd $(EXPORT_DIR)/Win/$(PROGRAM)-$(VERSION); mv $(PROGRAM).Windows.tar.gz $(DIST_DIR)/$(ROOTNAME).Windows.tar.gz
	echo $(ROOTNAME).Windows.tar.gz saved in $(DIST_DIR).
#	cd $(EXPORT_DIR)/Win; rm -rf $(PROGRAM)-$(VERSION)

debug: 
	mkdir -p $(DEBUG_DIR)
	cd $(DEBUG_DIR); make -f $(TOPDIR)/src/debug.mk; make -f $(TOPDIR)/src/Makefile CCFLAGS="-Wall -ansi -g" EXE=$(DEBUG_EXE)

test_dist: linux_test sunos_test source_test


web:
	cp $(DIST_DIR)/phreeqc-$(VERSION)*.tar.gz /var/anonymous/ftp/dlpark/geochem/unix/phreeqc
	cp $(TOPDIR)/doc/README.TXT /var/anonymous/ftp/dlpark/geochem/unix/phreeqc/README.TXT 
	cp $(TOPDIR)/doc/README.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/README.Unix.TXT
	cp $(TOPDIR)/win/README.TXT /var/anonymous/ftp/dlpark/geochem/pc/phreeqc/README.TXT 
	cp $(TOPDIR)/win/README.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/README.Win.TXT
	cp $(TOPDIR)/doc/phreeqc.txt /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/phreeqc.txt
	cp $(TOPDIR)/doc/RELEASE.TXT /z/linarcolkr/home/www/projects/GWC_coupled/phreeqc/RELEASE.TXT


