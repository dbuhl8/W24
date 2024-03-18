OBJECTS = LinAl.o Driver_LinAl.o
	
DOBJECTS = LinAl.o 

MODULES = LinAl.mod


FFLAGS = -fdefault-real-8 -fdefault-double-8 -llapack -lblas #-fbounds-check -g

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt
	python3 plot.py

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(FFLAGS) -o LinAl.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) LinAl.exe Image_appn_10* output.txt
