AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_srcdir)/src/phreeqcpp

noinst_LTLIBRARIES=libphreeqc.la

# dummy C++ source to cause C++ linking
nodist_EXTRA_libphreeqc_la_SOURCES = dummy.cxx

# library sources
libphreeqc_la_SOURCES=\
	advection.c\
	basic.c\
	basic.h\
	basicsubs.c\
	cl1.c\
	cvdense.c\
	cvdense.h\
	cvode.c\
	cvode.h\
	dense.c\
	dense.h\
	dw.c\
	global.h\
	global_structures.h\
	input.c\
	input.h\
	integrate.c\
	inverse.c\
	isotopes.c\
	kinetics.c\
	kinetics.h\
	mainsubs.c\
	model.c\
	nvector.c\
	nvector.h\
	nvector_serial.c\
	nvector_serial.h\
	output.c\
	output.h\
	p2c.h\
	p2clib.c\
	parse.c\
	phqalloc.c\
	phqalloc.h\
	phreeqc_files.c\
	phrqproto.h\
	phrqtype.h\
	pitzer.c\
	pitzer.h\
	pitzer_structures.c\
	pitzer_structures.h\
	prep.c\
	print.c\
	read.c\
	readtr.c\
	sit.c\
	smalldense.c\
	smalldense.h\
	spread.c\
	step.c\
	structures.c\
	sundialsmath.c\
	sundialsmath.h\
	sundialstypes.h\
	tally.c\
	tidy.c\
	transport.c\
	utilities.c

COMPILE = $(CXXCOMPILE)
LTCOMPILE = $(LTCXXCOMPILE)
