TOPDIR=..

all: advection.c basic.c basicsubs.c cl1.c cvdense.c cvode.c dense.c input.c integrate.c inverse.c isotopes.c kinetics.c main.c mainsubs.c model.c nvector.c nvector_serial.c output.c p2clib.c parse.c phqalloc.c phreeqc_files.c prep.c print.c read.c readtr.c smalldense.c spread.c step.c structures.c sundialsmath.c tally.c tidy.c transport.c utilities.c cvdense.h cvode.h dense.h global.h input.h kinetics.h nvector.h nvector_serial.h output.h p2c.h phqalloc.h phrqproto.h phrqtype.h smalldense.h sundialsmath.h sundialstypes.h

advection.c: $(TOPDIR)/advection.c
	cp $(TOPDIR)/advection.c .
basic.c: $(TOPDIR)/basic.c
	cp $(TOPDIR)/basic.c .
basicsubs.c: $(TOPDIR)/basicsubs.c
	cp $(TOPDIR)/basicsubs.c .
cl1.c: $(TOPDIR)/cl1.c
	cp $(TOPDIR)/cl1.c .
cvdense.c: $(TOPDIR)/cvdense.c
	cp $(TOPDIR)/cvdense.c .
cvode.c: $(TOPDIR)/cvode.c
	cp $(TOPDIR)/cvode.c .
dense.c: $(TOPDIR)/dense.c
	cp $(TOPDIR)/dense.c .
input.c: $(TOPDIR)/input.c
	cp $(TOPDIR)/input.c .
integrate.c: $(TOPDIR)/integrate.c
	cp $(TOPDIR)/integrate.c .
inverse.c: $(TOPDIR)/inverse.c
	cp $(TOPDIR)/inverse.c .
isotopes.c: $(TOPDIR)/isotopes.c
	cp $(TOPDIR)/isotopes.c .
kinetics.c: $(TOPDIR)/kinetics.c
	cp $(TOPDIR)/kinetics.c .
main.c: $(TOPDIR)/main.c
	cp $(TOPDIR)/main.c .
mainsubs.c: $(TOPDIR)/mainsubs.c
	cp $(TOPDIR)/mainsubs.c .
model.c: $(TOPDIR)/model.c
	cp $(TOPDIR)/model.c .
nvector.c: $(TOPDIR)/nvector.c
	cp $(TOPDIR)/nvector.c .
nvector_serial.c: $(TOPDIR)/nvector_serial.c
	cp $(TOPDIR)/nvector_serial.c .
output.c: $(TOPDIR)/output.c
	cp $(TOPDIR)/output.c .
p2clib.c: $(TOPDIR)/p2clib.c
	cp $(TOPDIR)/p2clib.c .
parse.c: $(TOPDIR)/parse.c
	cp $(TOPDIR)/parse.c .
phqalloc.c: $(TOPDIR)/phqalloc.c
	cp $(TOPDIR)/phqalloc.c .
phreeqc_files.c: $(TOPDIR)/phreeqc_files.c
	cp $(TOPDIR)/phreeqc_files.c .
prep.c: $(TOPDIR)/prep.c
	cp $(TOPDIR)/prep.c .
print.c: $(TOPDIR)/print.c
	cp $(TOPDIR)/print.c .
read.c: $(TOPDIR)/read.c
	cp $(TOPDIR)/read.c .
readtr.c: $(TOPDIR)/readtr.c
	cp $(TOPDIR)/readtr.c .
smalldense.c: $(TOPDIR)/smalldense.c
	cp $(TOPDIR)/smalldense.c .
spread.c: $(TOPDIR)/spread.c
	cp $(TOPDIR)/spread.c .
step.c: $(TOPDIR)/step.c
	cp $(TOPDIR)/step.c .
structures.c: $(TOPDIR)/structures.c
	cp $(TOPDIR)/structures.c .
sundialsmath.c: $(TOPDIR)/sundialsmath.c
	cp $(TOPDIR)/sundialsmath.c .
tally.c: $(TOPDIR)/tally.c
	cp $(TOPDIR)/tally.c .
tidy.c: $(TOPDIR)/tidy.c
	cp $(TOPDIR)/tidy.c .
transport.c: $(TOPDIR)/transport.c
	cp $(TOPDIR)/transport.c .
utilities.c: $(TOPDIR)/utilities.c
	cp $(TOPDIR)/utilities.c .
cvdense.h: $(TOPDIR)/cvdense.h
	cp $(TOPDIR)/cvdense.h .
cvode.h: $(TOPDIR)/cvode.h
	cp $(TOPDIR)/cvode.h .
dense.h: $(TOPDIR)/dense.h
	cp $(TOPDIR)/dense.h .
global.h: $(TOPDIR)/global.h
	cp $(TOPDIR)/global.h .
input.h: $(TOPDIR)/input.h
	cp $(TOPDIR)/input.h .
kinetics.h: $(TOPDIR)/kinetics.h
	cp $(TOPDIR)/kinetics.h .
nvector.h: $(TOPDIR)/nvector.h
	cp $(TOPDIR)/nvector.h .
nvector_serial.h: $(TOPDIR)/nvector_serial.h
	cp $(TOPDIR)/nvector_serial.h .
output.h: $(TOPDIR)/output.h
	cp $(TOPDIR)/output.h .
p2c.h: $(TOPDIR)/p2c.h
	cp $(TOPDIR)/p2c.h .
phqalloc.h: $(TOPDIR)/phqalloc.h
	cp $(TOPDIR)/phqalloc.h .
phrqproto.h: $(TOPDIR)/phrqproto.h
	cp $(TOPDIR)/phrqproto.h .
phrqtype.h: $(TOPDIR)/phrqtype.h
	cp $(TOPDIR)/phrqtype.h .
smalldense.h: $(TOPDIR)/smalldense.h
	cp $(TOPDIR)/smalldense.h .
sundialsmath.h: $(TOPDIR)/sundialsmath.h
	cp $(TOPDIR)/sundialsmath.h .
sundialstypes.h: $(TOPDIR)/sundialstypes.h
	cp $(TOPDIR)/sundialstypes.h .

