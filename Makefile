
#OBJS = StackObject.o nrutil.o optimization_mod.o utils.o quiet_hdf_mod.o quiet_mapfile_mod.o region_mod.o
OBJS = StackObject.o utils.o nrutil.o nr_mod.o optimization_mod.o quiet_hdf_mod.o quiet_mapfile_mod.o region_mod.o

all : spectral_regions spectral_region_performance

C_LIBDIR      = /mn/stornext/u2/eirikgje/.local/lib

FC  = ifort
F90FLAGS = -O3 -vec_report0 -assume byterecl -openmp
#F90FLAGS = -O0 -g -check bounds -check format -check pointers -check uninit -check output_conversion -traceback -openmp

FFLAGS := -I$(HOME)/.local/include
LDFLAGS =-L$(C_LIBDIR) -lhealpix -lcfitsio -lm -lhdf5_fortran -lhdf5 -openmp

spectral_regions  : $(OBJS) spectral_regions.o
	$(FC) -o spectral_regions spectral_regions.o $(OBJS) $(LDFLAGS)

spectral_region_performance : $(OBJS) spectral_region_performance.o
	$(FC) -o spectral_region_performance spectral_region_performance.o $(OBJS) $(LDFLAGS)

quiet_mapfile_mod.o: quiet_hdf_mod.o

spectral_regions.o: $(OBJS)

spectral_region_performance.o: $(OBJS)

%.o : %.f90
	$(FC) $(F90FLAGS) $(FFLAGS) -c $<

.PHONY: clean
clean: 
	rm -f *.o

