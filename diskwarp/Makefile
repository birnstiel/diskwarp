all: Makefile fort.f90
	f2py -m fort -c fort.f90

test: fort.f90
	gfortran fort.f90 -o test

run: test
	./test

clean:
	-rm *.o *.mod *.cpython*.so test x.dat y.dat z.dat v.dat