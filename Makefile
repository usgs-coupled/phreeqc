# 
# Make file for PHREEQC
#
# $(CURDIR) is current directory
TOPDIR:=$(CURDIR)/..
PROGRAM=phreeqcsax
EXE=$(TOPDIR)/bin/$(PROGRAM)
EXE=$(PROGRAM)
SRC:=$(CURDIR)

# Do not print commands before executing
#.SILENT:

# Provides compatibility with GNU make
#.SUFFIXES:

# Change to pawd if using automounter
PWD=pwd

# Change to C compiler on your system
CC=gcc

USE_XML=TRUE
XERCESCROOT=/z/parkplace/home/dlpark/packages/xerces-c-src_2_7_0

# Change to C compiler options on your system
ifdef OPTIMIZE
  CCFLAGS=-O3 -Wall -ansi -pedantic -I${XERCESCROOT}/include # -frounding-math  # -pg
  CCFLAGS_MODEL=-O2 -Wall -ansi -pedantic  # -pg
else
  CCFLAGS=-g -Wall -ansi -pedantic -I${XERCESCROOT}/include # -frounding-math  # -pg
  CCFLAGS_MODEL=-g -Wall -ansi -pedantic  # -pg
endif
# Remove the following definition if you do not have 
# gmp (Gnu Multiple Precision) package on your system
INVERSE_CL1MP=TRUE

LOADFLAGS= -lm -lxerces-c # -pg

PLATFORM= LINUX
CXX= g++ -c -D${PLATFORM} -D_REENTRANT -fpic 
ifdef OPTIMIZE
  CXXFLAGS= -O3
else
  CXXFLAGS= -Wall -g
endif
LINK= g++ -D${PLATFORM} -fpic
PLATFORM_LIB_LINK_OPTIONS=-L/usr/lib -L/usr/local/lib
EXTRA_LINK_OPTIONS=-lc 
LIBRARY_NAMES= -lxerces-c
LIBRARY_SEARCH_PATHS= -L${XERCESCROOT}/lib # -L/home/dlpark/lib #
INCLUDES= -I${XERCESCROOT}/include 

#.c.o : 
#	${CC} ${CCFLAGS} -c -o $@ $<
#%.o : $(SRC)/%.c
#	${CC} ${CCFLAGS} -c -o $@ $<
%.o : $(SRC)/%.cpp
	${CXX} ${CXXFLAGS} $(INCLUDES) -c -o $@ $<
%.o : $(SRC)/%.cxx
	${CXX} ${CXXFLAGS} $(INCLUDES) -c -o $@ $<

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

CLASS_OBJECTS=  Conc.o \
		Exchange.o \
		ExchComp.o \
		GasPhase.o \
		ISolution.o \
		Isotope.o \
		KineticsCxx.o \
		KineticsComp.o \
		NameDouble.o \
		NumKeyword.o \
		Parser.o \
		PPassemblage.o \
		PPassemblageComp.o \
		Reaction.o \
		ReadClass.o \
		Solution.o \
		SSassemblage.o \
		SSassemblageSS.o \
		Surface.o \
		SurfComp.o \
		SurfCharge.o \
		Utils.o

OBJECTS += $(CLASS_OBJECTS)

ifdef USE_XML
	OBJECTS += SAXPhreeqc.o
endif

ifdef INVERSE_CL1MP
	LOADFLAGS += /z/parkplace/usr/lib/libgmp.a 
	CCFLAGS += -DINVERSE_CL1MP
	CXXFLAGS += -DINVERSE_CL1MP
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
ifdef USE_XML
	${LINK} ${PLATFORM_LIB_LINK_OPTIONS} ${OBJECTS} -o $(EXE) ${LIBRARY_SEARCH_PATHS} ${LIBRARY_NAMES} ${EXTRA_LINK_OPTIONS} ${LOADFLAGS}
else	
	$(CC) -o $(EXE) $(OBJECTS) $(LOADFLAGS) # -L/z/estespark/home/dlpark/packages/efence -lefence
endif
	echo Compilation complete, $(EXE).

advection.o: $(SRC)/advection.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

basic.o: $(SRC)/basic.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/p2c.h

basicsubs.o: $(SRC)/basicsubs.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

cl1.o: $(SRC)/cl1.cpp $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqtype.h

cl1mp.o: $(SRC)/cl1mp.cpp $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqtype.h

cvdense.o: $(SRC)/cvdense.cpp $(SRC)/cvdense.h $(SRC)/cvode.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/nvector.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

cvode.o: $(SRC)/cvode.cpp $(SRC)/cvode.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/nvector.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/kinetics.h $(SRC)/phqalloc.h

dense.o: $(SRC)/dense.cpp $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/output.h $(SRC)/phqalloc.h

input.o: $(SRC)/input.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/input.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/phqalloc.h

integrate.o: $(SRC)/integrate.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

inverse.o: $(SRC)/inverse.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

isotopes.o: $(SRC)/isotopes.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

kinetics.o: $(SRC)/kinetics.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/sundialstypes.h $(SRC)/cvode.h $(SRC)/nvector.h $(SRC)/cvdense.h $(SRC)/dense.h $(SRC)/smalldense.h $(SRC)/nvector_serial.h $(SRC)/kinetics.h

main.o: $(SRC)/main.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

mainsubs.o: $(SRC)/mainsubs.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

model.o: $(SRC)/model.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h
	${CC} $(SRC)/model.cpp ${CCFLAGS_MODEL} -c -o model.o #-ffloat-store 

nvector.o: $(SRC)/nvector.cpp $(SRC)/nvector.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/output.h

nvector_serial.o: $(SRC)/nvector_serial.cpp $(SRC)/nvector_serial.h $(SRC)/nvector.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

output.o: $(SRC)/output.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/phqalloc.h

p2clib.o: $(SRC)/p2clib.cpp $(SRC)/p2c.h $(SRC)/output.h

parse.o: $(SRC)/parse.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

phqalloc.o: $(SRC)/phqalloc.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/output.h

phreeqc_files.o: $(SRC)/phreeqc_files.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/input.h

pitzer.o: $(SRC)/pitzer.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/pitzer.h

dw.o: $(SRC)/dw.cpp $(SRC)/pitzer.h

pitzer_structures.o: $(SRC)/pitzer_structures.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h $(SRC)/pitzer.h

prep.o: $(SRC)/prep.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

prep.o: $(SRC)/prep.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

print.o: $(SRC)/print.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

read.o: $(SRC)/read.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

readtr.o: $(SRC)/readtr.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

SAXPhreeqc.o: $(SRC)/SAXPhreeqc.cpp $(SRC)/SAXPhreeqc.h $(SRC)/SaxPhreeqcHandlers.h
	${CXX} ${CXXFLAGS} $(INCLUDES) -o SAXPhreeqc.o $(SRC)/SAXPhreeqc.cpp

smalldense.o: $(SRC)/smalldense.cpp $(SRC)/smalldense.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/sundialsmath.h $(SRC)/output.h $(SRC)/phqalloc.h

spread.o: $(SRC)/spread.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

step.o: $(SRC)/step.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

structures.o: $(SRC)/structures.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

sundialsmath.o: $(SRC)/sundialsmath.cpp $(SRC)/sundialsmath.h $(SRC)/sundialstypes.h $(SRC)/phrqtype.h $(SRC)/output.h

tally.o: $(SRC)/tally.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

tidy.o: $(SRC)/tidy.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

transport.o: $(SRC)/transport.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

utilities.o: $(SRC)/utilities.cpp $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/output.h $(SRC)/phrqproto.h

Conc.o: $(SRC)/Conc.cxx $(SRC)/Conc.h $(SRC)/Utils.h $(SRC)/char_star.h $(SRC)/ISolution.h $(SRC)/NumKeyword.h \
  $(SRC)/Parser.h $(SRC)/Solution.h $(SRC)/Isotope.h $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h \
  $(SRC)/phrqproto.h $(SRC)/phqalloc.h
Exchange.o: $(SRC)/Exchange.cxx $(SRC)/Utils.h $(SRC)/Exchange.h $(SRC)/NumKeyword.h Parser.h \
  $(SRC)/char_star.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/ExchComp.h $(SRC)/NameDouble.h $(SRC)/phqalloc.h \
  $(SRC)/phrqproto.h
ExchComp.o: $(SRC)/ExchComp.cxx $(SRC)/Utils.h $(SRC)/ExchComp.h $(SRC)/NameDouble.h $(SRC)/global.h \
  $(SRC)/phrqtype.h $(SRC)/char_star.h $(SRC)/Parser.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
GasPhase.o: $(SRC)/GasPhase.cxx $(SRC)/Utils.h $(SRC)/GasPhase.h $(SRC)/NumKeyword.h $(SRC)/Parser.h \
  $(SRC)/char_star.h $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
ISolution.o: $(SRC)/ISolution.cxx $(SRC)/ISolution.h $(SRC)/NumKeyword.h $(SRC)/Parser.h $(SRC)/char_star.h \
  $(SRC)/Solution.h $(SRC)/Isotope.h $(SRC)/Conc.h $(SRC)/Utils.h $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h \
  $(SRC)/phqalloc.h $(SRC)/phrqproto.h
Isotope.o: $(SRC)/Isotope.cxx $(SRC)/Isotope.h $(SRC)/Parser.h $(SRC)/char_star.h $(SRC)/Utils.h $(SRC)/global.h \
  $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
KineticsComp.o: $(SRC)/KineticsComp.cxx $(SRC)/Utils.h $(SRC)/KineticsComp.h $(SRC)/NameDouble.h \
  $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/char_star.h $(SRC)/Parser.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
KineticsCxx.o: $(SRC)/KineticsCxx.cxx $(SRC)/Utils.h $(SRC)/KineticsCxx.h $(SRC)/NumKeyword.h $(SRC)/Parser.h \
  $(SRC)/char_star.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/KineticsComp.h $(SRC)/NameDouble.h $(SRC)/phqalloc.h \
  $(SRC)/phrqproto.h
NameDouble.o: $(SRC)/NameDouble.cxx $(SRC)/Utils.h $(SRC)/Conc.h $(SRC)/char_star.h $(SRC)/NameDouble.h \
  $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/Parser.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
NumKeyword.o: $(SRC)/NumKeyword.cxx $(SRC)/NumKeyword.h $(SRC)/Parser.h $(SRC)/char_star.h
Reaction.o: $(SRC)/Reaction.cxx $(SRC)/Utils.h $(SRC)/Reaction.h $(SRC)/NumKeyword.h $(SRC)/Parser.h \
  $(SRC)/char_star.h $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
ReadClass.o: $(SRC)/ReadClass.cpp
Parser.o: $(SRC)/Parser.cxx $(SRC)/Parser.h $(SRC)/char_star.h $(SRC)/Utils.h
PPassemblageComp.o: $(SRC)/PPassemblageComp.cxx $(SRC)/Utils.h $(SRC)/NameDouble.h \
  $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/char_star.h $(SRC)/Parser.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
PPassemblage.o: $(SRC)/PPassemblage.cxx $(SRC)/Utils.h $(SRC)/PPassemblage.h $(SRC)/NumKeyword.h \
  $(SRC)/Parser.h $(SRC)/char_star.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/PPassemblageComp.h \
  $(SRC)/NameDouble.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
Solution.o: $(SRC)/Solution.cxx $(SRC)/Utils.h $(SRC)/Solution.h $(SRC)/NumKeyword.h $(SRC)/Parser.h \
  $(SRC)/char_star.h $(SRC)/Isotope.h $(SRC)/Conc.h $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h \
  $(SRC)/phqalloc.h $(SRC)/phrqproto.h $(SRC)/ISolution.h
SSassemblage.o: $(SRC)/SSassemblage.cxx $(SRC)/Utils.h $(SRC)/SSassemblage.h $(SRC)/NumKeyword.h \
  $(SRC)/Parser.h $(SRC)/char_star.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/SSassemblageSS.h \
  $(SRC)/NameDouble.h $(SRC)/SSassemblageSS.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h
SSassemblageSS.o: $(SRC)/SSassemblageSS.cxx $(SRC)/Utils.h $(SRC)/SSassemblageSS.h \
  $(SRC)/NameDouble.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/char_star.h $(SRC)/Parser.h $(SRC)/phqalloc.h \
  $(SRC)/phrqproto.h
Surface.o: $(SRC)/Surface.cxx $(SRC)/Utils.h $(SRC)/Surface.h $(SRC)/NumKeyword.h $(SRC)/Parser.h \
  $(SRC)/char_star.h $(SRC)/global.h $(SRC)/phrqtype.h $(SRC)/SurfComp.h $(SRC)/NameDouble.h $(SRC)/phqalloc.h \
  $(SRC)/phrqproto.h
SurfComp.o: $(SRC)/SurfComp.cxx $(SRC)/Utils.h $(SRC)/SurfComp.h $(SRC)/NameDouble.h $(SRC)/global.h \
  $(SRC)/phrqtype.h $(SRC)/char_star.h $(SRC)/Parser.h $(SRC)/phqalloc.h $(SRC)/phrqproto.h

Utils.o: $(SRC)/Utils.cxx $(SRC)/Utils.h $(SRC)/Parser.h $(SRC)/char_star.h

-include $(SRC)/distribution.mk


