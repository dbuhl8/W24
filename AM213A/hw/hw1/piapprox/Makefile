
#Directory Stucture for Makefile
#Plot Directory
#PFOLD = ~/Documents/Fluids/gaussian/plotData
#TFOLD = ~/Documents/Fluids/gaussian/plots
#Log File Directory
#LFOLD = ~/Documents/Fluids/gaussian/logs


FC = gfortran
FFLAGS = -fdefault-real-8 -fdefault-double-8
#LFLAGS = -L/usr/local/lib -L/usr/lib

#Module File
#OMOD = gaussian_mod.o
#FMOD = ../src/gaussian_mod.f90

.PHONY : clean

computepi: computepi.e
		./$< >> exerciselog.txt

%.e : %.o
		$(FC) -o $@ $(FFLAGS) $< 

%.o : %.f90
		$(FC) -c $(FFLAGS) $<

clean: 
		rm -rf *.e *.o


