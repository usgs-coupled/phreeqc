all: advection.c basic.c basicsubs.c cl1.c cvdense.c cvode.c dense.c input.c integrate.c inverse.c isotopes.c kinetics.c main.c mainsubs.c model.c nvector.c nvector_serial.c output.c p2clib.c parse.c phqalloc.c phreeqc_files.c prep.c print.c read.c readtr.c smalldense.c spread.c step.c structures.c sundialsmath.c tally.c tidy.c transport.c utilities.c cvdense.h cvode.h dense.h global.h input.h kinetics.h nvector.h nvector_serial.h output.h p2c.h phqalloc.h phrqproto.h phrqtype.h smalldense.h sundialsmath.h sundialstypes.h

advection.c: ../advection.c
	cp ../advection.c .
basic.c: ../basic.c
	cp ../basic.c .
basicsubs.c: ../basicsubs.c
	cp ../basicsubs.c .
cl1.c: ../cl1.c
	cp ../cl1.c .
cvdense.c: ../cvdense.c
	cp ../cvdense.c .
cvode.c: ../cvode.c
	cp ../cvode.c .
dense.c: ../dense.c
	cp ../dense.c .
input.c: ../input.c
	cp ../input.c .
integrate.c: ../integrate.c
	cp ../integrate.c .
inverse.c: ../inverse.c
	cp ../inverse.c .
isotopes.c: ../isotopes.c
	cp ../isotopes.c .
kinetics.c: ../kinetics.c
	cp ../kinetics.c .
main.c: ../main.c
	cp ../main.c .
mainsubs.c: ../mainsubs.c
	cp ../mainsubs.c .
model.c: ../model.c
	cp ../model.c .
nvector.c: ../nvector.c
	cp ../nvector.c .
nvector_serial.c: ../nvector_serial.c
	cp ../nvector_serial.c .
output.c: ../output.c
	cp ../output.c .
p2clib.c: ../p2clib.c
	cp ../p2clib.c .
parse.c: ../parse.c
	cp ../parse.c .
phqalloc.c: ../phqalloc.c
	cp ../phqalloc.c .
phreeqc_files.c: ../phreeqc_files.c
	cp ../phreeqc_files.c .
prep.c: ../prep.c
	cp ../prep.c .
print.c: ../print.c
	cp ../print.c .
read.c: ../read.c
	cp ../read.c .
readtr.c: ../readtr.c
	cp ../readtr.c .
smalldense.c: ../smalldense.c
	cp ../smalldense.c .
spread.c: ../spread.c
	cp ../spread.c .
step.c: ../step.c
	cp ../step.c .
structures.c: ../structures.c
	cp ../structures.c .
sundialsmath.c: ../sundialsmath.c
	cp ../sundialsmath.c .
tally.c: ../tally.c
	cp ../tally.c .
tidy.c: ../tidy.c
	cp ../tidy.c .
transport.c: ../transport.c
	cp ../transport.c .
utilities.c: ../utilities.c
	cp ../utilities.c .
cvdense.h: ../cvdense.h
	cp ../cvdense.h .
cvode.h: ../cvode.h
	cp ../cvode.h .
dense.h: ../dense.h
	cp ../dense.h .
global.h: ../global.h
	cp ../global.h .
input.h: ../input.h
	cp ../input.h .
kinetics.h: ../kinetics.h
	cp ../kinetics.h .
nvector.h: ../nvector.h
	cp ../nvector.h .
nvector_serial.h: ../nvector_serial.h
	cp ../nvector_serial.h .
output.h: ../output.h
	cp ../output.h .
p2c.h: ../p2c.h
	cp ../p2c.h .
phqalloc.h: ../phqalloc.h
	cp ../phqalloc.h .
phrqproto.h: ../phrqproto.h
	cp ../phrqproto.h .
phrqtype.h: ../phrqtype.h
	cp ../phrqtype.h .
smalldense.h: ../smalldense.h
	cp ../smalldense.h .
sundialsmath.h: ../sundialsmath.h
	cp ../sundialsmath.h .
sundialstypes.h: ../sundialstypes.h
	cp ../sundialstypes.h .

