OBJECTS = LinAl.o Driver_LinAl.o
	
DOBJECTS = LinAl.o debug.o

MODULES = LinAl.mod

FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbounds-check -llapack -lblas

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

debug: debug.exe
	./debug.exe > debug.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

debug.exe: $(MODULES) $(DOBJECTS)
	gfortran $(DOBJECTS) -o debug.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(DOBJECTS) $(MODULES) LinAl.exe debug.exe
