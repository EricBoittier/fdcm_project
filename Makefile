F90 = gfortran
#F90 = /usr/local/bin/mpif90

FFLAGS_PAR = -fopenmp -O3 -fbacktrace 
LFLAGS_PAR = -fopenmp -fbacktrace 
FFLAGS     = -fbacktrace  -O3
LFLAGS     = -fbacktrace 

##########################
# Object Files for build #
##########################

OBJS = differential_evolution.o m_strings.o main.o

OBJS_PAR = pdifferential_evolution.o m_strings.o main.o

all : fdcm.x pfdcm.x clean

serial : fdcm.x clean
parallel : pfdcm.x clean

fdcm.x : $(OBJS)
	 ${F90}  -o $@ $(LFLAGS) $(OBJS)
	 
pfdcm.x : $(OBJS_PAR)
	 ${F90}  -o $@ $(LFLAGS_PAR) $(OBJS_PAR)

#######################################
# Object dependencies and compilation #
#######################################
m_strings.o : src/m_strings.f90 \
    $(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/m_strings.f90

differential_evolution.o : src/differential_evolution.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/differential_evolution.f90
	
pdifferential_evolution.o : src/differential_evolution.f90
	$(F90) -c $(FFLAGS_PAR) $(INCLUDES) -o $@ src/differential_evolution.f90

main.o : src/main.f90 \
differential_evolution.o m_string.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/main.f90


.PHONY: clean veryclean
clean:
	rm *.o *.mod

veryclean:
	rm *.x *.o *.mod



