
COMPILER = mpif90

OBJECTS = LinAl.o shearSolve2.o
	
DOBJECTS = LinAl.o 

MODULES = LinAl.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas -ffree-line-length-0 #-g3 -fbacktrace -fbounds-check

.PHONY: clean

output.txt: shearSolve2
	./shearSolve2

shearSolve2: $(MODULES) $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(FFLAGS) -o shearSolve2

%.o: %.f90
	$(COMPILER) -c $(FFLAGS) $<

%.mod: %.f90
	$(COMPILER) -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) shearSolve2 *.dat
