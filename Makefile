# 
# Make file for PHREEQC
#
TOPDIR=~dlpark/programs/phreeqc
PROGRAM=phreeqc
EXE=$(TOPDIR)/bin/$(PROGRAM)

# Do not print commands before executing
#.SILENT:

# Provides compatibility with GNU make
#.SUFFIXES:

# Change to pawd if using automounter
PWD=pwd

# Change to C compiler on your system
CC=gcc

# Change to C compiler options on your system
CCFLAGS=-O3 -Wall -ansi -pedantic # -pg

.c.o :
	${CC} ${CCFLAGS} -c -o $@ $<

LOADFLAGS= -lm # -pg

# Location to copy scripts on installation
BINDIR=$(HOME)/bin

OBJECTS=	main.o \
		advection.o \
		basic.o \
		basicsubs.o \
		cl1.o \
		input.o \
		integrate.o \
		inverse.o \
		isotopes.o \
		kinetics.o \
		mainsubs.o \
		output.o \
		model.o \
		p2clib.o \
		parse.o \
		phreeqc_files.o \
		phqalloc.o \
		prep.o \
		print.o \
		read.o \
		readtr.o \
		spread.o \
		step.o \
		structures.o \
		tally.o \
		tidy.o \
		transport.o \
		utilities.o \
		cvdense.o \
		cvode.o \
		dense.o \
		nvector.o \
		nvector_serial.o \
		smalldense.o \
		sundialsmath.o 

all: $(EXE)

install: 
#
#    Create directory for binary and scripts if necessary
#
	if [ ! -d $(BINDIR) ]; \
           then \
              mkdir $(BINDIR); \
              echo Created directory $(BINDIR); \
        fi
#
#    Put path name of current directory into script for
#    locating data files, put script in top directory,
#    put symbolic link in BINDIR
#
	cd $(TOPDIR); dir1=`$(PWD)`/bin; cd $(BINDIR); if [ `$(PWD)` = $$dir1 ]; \
			then \
				echo "Can not install to $(BINDIR). Choose another directory."; \
				exit 4 ; \
			fi
	cd $(TOPDIR); \
		rm -f $(BINDIR)/$(PROGRAM); \
		rm -f $(PROGRAM); \
		sed "s?TOPDIR=.\{0,80\}?TOPDIR=`$(PWD)`?" bin/$(PROGRAM).orig > $(PROGRAM); \
		chmod 755 $(PROGRAM)
	cd $(TOPDIR); dir1=`$(PWD)`; cd $(BINDIR); if [ `$(PWD)` != $$dir1 ]; then \
			ln -s $$dir1/$(PROGRAM) $(BINDIR); \
			echo Symbolic link for $(PROGRAM) has been placed in $(BINDIR).	; \
		fi
#
#     Check that all necessary files are in place.
#
	if [ -f $(BINDIR)/$(PROGRAM) -a \
             -f $(TOPDIR)/bin/$(PROGRAM) -a \
             -f $(TOPDIR)/$(PROGRAM) ]; \
           then echo "Installation complete."; \
           else echo "Installation incomplete."; \
              for FILE in $(BINDIR)/$(PROGRAM) \
                 $(TOPDIR)/bin/$(PROGRAM) $(TOPDIR)/$(PROGRAM) ; \
              do \
                 if [ ! -f $$FILE ]; then echo $$FILE is missing.; fi; \
              done; \
           fi
	echo "Add directory $(BINDIR) to PATH."

clean:
	rm -f $(BINDIR)/$(PROGRAM)
	rm -f $(TOPDIR)/$(PROGRAM)
	rm -f $(TOPDIR)/bin/$(PROGRAM) *.o 
	echo Removed files generated by make.

$(EXE): $(OBJECTS) 
	$(CC) -o $(EXE) $(OBJECTS) $(LOADFLAGS) # -L/z/estespark/home/dlpark/packages/efence -lefence
	echo Compilation complete, $(EXE).

advection.o: advection.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
basic.o: basic.c global.h phrqtype.h phqalloc.h output.h phrqproto.h \
  p2c.h
basicsubs.o: basicsubs.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
cl1.o: cl1.c phqalloc.h output.h phrqtype.h
cvdense.o: cvdense.c cvdense.h cvode.h sundialstypes.h phrqtype.h \
  nvector.h dense.h smalldense.h sundialsmath.h output.h phqalloc.h
cvode.o: cvode.c cvode.h sundialstypes.h phrqtype.h nvector.h \
  sundialsmath.h output.h kinetics.h phqalloc.h
dense.o: dense.c sundialstypes.h phrqtype.h sundialsmath.h dense.h \
  smalldense.h output.h phqalloc.h
input.o: input.c global.h phrqtype.h input.h output.h phrqproto.h \
  phqalloc.h
integrate.o: integrate.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
inverse.o: inverse.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
isotopes.o: isotopes.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
kinetics.o: kinetics.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h sundialstypes.h cvode.h nvector.h cvdense.h dense.h \
  smalldense.h nvector_serial.h kinetics.h
main.o: main.c global.h phrqtype.h output.h phrqproto.h input.h
mainsubs.o: mainsubs.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h input.h
model.o: model.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
nvector.o: nvector.c nvector.h sundialstypes.h phrqtype.h output.h
nvector_serial.o: nvector_serial.c nvector_serial.h nvector.h \
  sundialstypes.h phrqtype.h sundialsmath.h output.h phqalloc.h
output.o: output.c global.h phrqtype.h output.h phrqproto.h phqalloc.h
p2clib.o: p2clib.c p2c.h output.h
parse.o: parse.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
phqalloc.o: phqalloc.c global.h phrqtype.h output.h
phreeqc_files.o: phreeqc_files.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h input.h
prep.o: prep.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
print.o: print.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
read.o: read.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
readtr.o: readtr.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
smalldense.o: smalldense.c smalldense.h sundialstypes.h phrqtype.h \
  sundialsmath.h output.h phqalloc.h
spread.o: spread.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
step.o: step.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
structures.o: structures.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
sundialsmath.o: sundialsmath.c sundialsmath.h sundialstypes.h phrqtype.h \
  output.h
tally.o: tally.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
tidy.o: tidy.c global.h phrqtype.h phqalloc.h output.h phrqproto.h
transport.o: transport.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h
utilities.o: utilities.c global.h phrqtype.h phqalloc.h output.h \
  phrqproto.h

-include $(TOPDIR)/src/distribution.mk


