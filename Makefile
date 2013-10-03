
#OBJS = StackObject.o nrutil.o optimization_mod.o utils.o quiet_hdf_mod.o quiet_mapfile_mod.o region_mod.o
COMMONOBJS = StackObject.o utils.o nrutil.o nr_mod.o optimization_mod.o quiet_hdf_mod.o quiet_mapfile_mod.o 
TIMEDOBJS = $(COMMONOBJS) region_mod_timed.o
MAINOBJS = $(COMMONOBJS) region_mod.o
PERFORMANCEOBJS = $(COMMONOBJS) region_mod.o
BFOBJS = utils.o nrutil.o nr_mod.o quiet_hdf_mod.o quiet_mapfile_mod.o


all : spectral_regions spectral_region_performance spectral_regions_timed find_bestfit_si

C_LIBDIR      = /mn/stornext/u2/eirikgje/.local/lib

FC  = ifort
F90FLAGS = -O3 -vec_report0 -assume byterecl -openmp
#F90FLAGS = -O0 -g -check bounds -check format -check pointers -check uninit -check output_conversion -traceback -openmp

FFLAGS := -I$(HOME)/.local/include
LDFLAGS =-L$(C_LIBDIR) -lhealpix -lcfitsio -lm -lhdf5_fortran -lhdf5 -openmp

spectral_regions  : $(MAINOBJS) spectral_regions.o
	$(FC) -o spectral_regions spectral_regions.o $(MAINOBJS) $(LDFLAGS)

spectral_region_performance : $(PERFORMANCEOBJS) spectral_region_performance.o
	$(FC) -o spectral_region_performance spectral_region_performance.o $(PERFORMANCEOBJS) $(LDFLAGS)

spectral_regions_timed  : $(TIMEDOBJS) spectral_regions_timed.o
	$(FC) -o spectral_regions_timed spectral_regions_timed.o $(TIMEDOBJS) $(LDFLAGS)

find_bestfit_si	: $(BFOBJS) find_bestfit_si.o
	$(FC) -o find_bestfit find_bestfit_si.o $(BFOBJS) $(LDFLAGS)

quiet_mapfile_mod.o: quiet_hdf_mod.o

spectral_regions.o: $(MAINOBJS)

spectral_region_performance.o: $(PERFORMANCEOBJS)

spectral_regions_timed.o: $(TIMEDOBJS)

find_bestfit_si.o: $(BFOBJS)
	$(FC) $(F90FLAGS) $(FFLAGS) -c find_bestfit_si/find_bestfit_si.f90

%.o : %.f90
	$(FC) $(F90FLAGS) $(FFLAGS) -c $<

.PHONY: clean
clean: 
	rm -f *.o
