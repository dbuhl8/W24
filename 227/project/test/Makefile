OBJECTS = LinAl.o shearSolve2.o
	
DOBJECTS = LinAl.o 

MODULES = LinAl.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas -ffree-line-length-0 #-g3 -fbacktrace -fbounds-check

.PHONY: clean

output.txt: shearSolve2.exe
	./shearSolve2.exe > output.txt

shearSolve2.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o shearSolve2.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) shearSolve2.exe *.dat
