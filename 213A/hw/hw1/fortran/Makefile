


FC = gfortran
FFLAGS = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0

#Module File
OMOD = bisect.o fcosx.o
FMOD = bisect.f90 fcosx.f90

.PHONY : clean

all: computepi.e addmats.e betteraddmats.e solbybisection.e
		$(shell echo "Successfully made all executables")

%.e : %.o
		$(FC) -o $@ $(FFLAGS) $(OMOD) $< 

%.o : %.f90
		$(FC) -c $(FFLAGS) $(FMOD) $<

clean: 
		rm -rf *.e *.o


