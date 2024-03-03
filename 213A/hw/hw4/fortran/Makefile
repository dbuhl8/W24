OBJECTS = LinAl.o Driver_LinAl.o

MODULES = LinAl.mod

OUTPUTS = qr.png chol.png plot.dat qr.dat chol.dat

FFLAGS = #-fdefault-real-8 -fdefault-double-8 -fbounds-check 

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt && gnuplot -c gplot.plt && gnuplot -c qplot.plt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o LinAl.exe

%.o: %.f90
	gfortran -c $(FFLAGS) $<

%.mod: %.f90
	gfortran -c $(FFLAGS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(OUTPUTS) LinAl.exe
