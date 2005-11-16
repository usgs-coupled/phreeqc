# 
# Make file for PHREEQC
#
# $(CURDIR) is current directory
TOPDIR:=$(CURDIR)/..
PROGRAM=phreeqc
EXE=$(TOPDIR)/bin/$(PROGRAM)
SRC=.

# Do not print commands before executing
#.SILENT:

# Provides compatibility with GNU make
#.SUFFIXES:

# Change to pawd if using automounter
PWD=pwd

# Change to C compiler on your system
CC=gcc

# Change to C compiler options on your system
CCFLAGS=-O3 -Wall -ansi -pedantic # -frounding-math  # -pg
CCFLAGS_MODEL=-O2 -Wall -ansi -pedantic  # -pg

# Remove the following definition if you do not have 
# gmp (Gnu Multiple Precision) package on your system
INVERSE_CL1MP=TRUE

LOADFLAGS= -lm  # -pg

#.c.o : 
#	${CC} ${CCFLAGS} -c -o $@ $<
%.o : $(SRC)/%.c
	${CC} ${CCFLAGS} -c -o $@ $<

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
		sundialsmath.o \
		dw.o \
		pitzer.o \
		pitzer_structures.o \

ifdef INVERSE_CL1MP
	LOADFLAGS += /z/parkplace/usr/lib/libgmp.a 
	CCFLAGS += -DINVERSE_CL1MP
	OBJECTS += cl1mp.o
endif

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
	rm -f $(TOPDIR)/bin/$(PROGRAM) 
	rm -f $(TOPDIR)/src/*.o 
	rm -f $(SUN_DIR)/bin/$(PROGRAM)
	rm -f $(SUN_DIR)/src/*.o
	rm -f $(TOPDIR)/src/$(PROGRAM) 
	rm -f $(DEBUG_DIR)/src/*.o
	echo Removed object and executable files generated by make.

$(EXE): $(OBJECTS) 
	echo $(TOPDIR)
	$(CC) -o $(EXE) $(OBJECTS) $(LOADFLAGS) # -L/z/estespark/home/dlpark/packages/efence -lefence
	echo Compilation complete, $(EXE).

advection.o: $(SRC)/advection.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

basic.o: $(SRC)/basic.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/p2c.h

basicsubs.o: $(SRC)/basicsubs.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

cl1.o: $(SRC)/cl1.c $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqtype.h

cl1mp.o: $(SRC)/cl1mp.c $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqtype.h

cvdense.o: $(SRC)/cvdense.c $(SRC)/cvdense.h $(SRC)/cvode.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/nvector.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

cvode.o: $(SRC)/cvode.c $(SRC)/cvode.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/nvector.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/kinetics.h $(SRC)/phqalloc.h

dense.o: $(SRC)/dense.c $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/output.h $(SRC)/phqalloc.h

input.o: $(SRC)/input.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/input.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/phqalloc.h

integrate.o: $(SRC)/integrate.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

inverse.o: $(SRC)/inverse.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

isotopes.o: $(SRC)/isotopes.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

kinetics.o: $(SRC)/kinetics.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/sundialstypes.h $(SRC)/cvode.h $(SRC)/nvector.h $(SRC)/cvdense.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/nvector_serial.h $(SRC)/kinetics.h

main.o: $(SRC)/main.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

mainsubs.o: $(SRC)/mainsubs.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

model.o: $(SRC)/model.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h
	${CC} model.c ${CCFLAGS_MODEL} -c -o model.o #-ffloat-store 

nvector.o: $(SRC)/nvector.c $(SRC)/nvector.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/output.h

nvector_serial.o: $(SRC)/nvector_serial.c $(SRC)/nvector_serial.h $(SRC)/nvector.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

output.o: $(SRC)/output.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/phqalloc.h

p2clib.o: $(SRC)/p2clib.c $(SRC)/p2c.h $(SRC)/output.h

parse.o: $(SRC)/parse.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

phqalloc.o: $(SRC)/phqalloc.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h

phreeqc_files.o: $(SRC)/phreeqc_files.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

pitzer.o: $(SRC)/pitzer.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/pitzer.h

dw.o: $(SRC)/dw.c $(SRC)/pitzer.h

pitzer_structures.o: $(SRC)/pitzer_structures.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/pitzer.h

prep.o: $(SRC)/prep.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

prep.o: $(SRC)/prep.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

print.o: $(SRC)/print.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

read.o: $(SRC)/read.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

readtr.o: $(SRC)/readtr.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

smalldense.o: $(SRC)/smalldense.c $(SRC)/smalldense.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

spread.o: $(SRC)/spread.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

step.o: $(SRC)/step.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

structures.o: $(SRC)/structures.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

sundialsmath.o: $(SRC)/sundialsmath.c $(SRC)/sundialsmath.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/output.h

tally.o: $(SRC)/tally.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

tidy.o: $(SRC)/tidy.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

transport.o: $(SRC)/transport.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

utilities.o: $(SRC)/utilities.c $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

-include $(TOPDIR)/src/distribution.mk


