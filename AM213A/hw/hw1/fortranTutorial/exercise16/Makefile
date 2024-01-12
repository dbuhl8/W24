
#Directory Stucture for Makefile
#Plot Directory
#PFOLD = ~/Documents/Fluids/gaussian/plotData
#TFOLD = ~/Documents/Fluids/gaussian/plots
#Log File Directory
#LFOLD = ~/Documents/Fluids/gaussian/logs


FC = gfortran
#FFLAGS = -llapack -lblas
#What do LFLAGS do??????
#LFLAGS = -L/usr/local/lib -L/usr/lib

#Module File
OMOD = bisect.o fcosx.o
FMOD = bisect.f90 fcosx.f90

.PHONY : clean

solbybisection: solbybisection.e 
		./$< 

%.e : %.o
		$(FC) -o $@ $(OMOD) $< 

%.o : %.f90
		$(FC) -c $(FMOD) $<

clean: 
		rm -rf *.e *.o


