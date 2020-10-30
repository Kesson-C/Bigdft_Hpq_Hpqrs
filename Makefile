
MKL_PATH=/home/dell/intel/compilers_and_libraries_2017.8.262/linux/mkl/lib/intel64

BIGDFT_BUILDDIR=/home/chou/software/BigDFT-1.8.3/bigdft/
BIGDFT_BUILDDIR=/home/chou/software/BigDFT-1.8.3/bigdft-revised/

         libdir=/home/chou/software/BigDFT-1.8.3/build/install/lib/

     FC = mpiifort -scalar-rep -fopenmp -m64 -g -Wl,--no-as-needed -ldl # ifort (IFORT) 17.0.4 20170411
FCFLAGS = -Ofast -xSSE4.2 -axAVX,CORE-AVX2

     FC = mpiifort             -fopenmp 
FCFLAGS = 

INCLUDES = -I$(BIGDFT_BUILDDIR)/includes -I$(BIGDFT_BUILDDIR)/../build/install/include
LIBS_DEPENDENCIES = 

LIBS_DEPENDENCIES = -L$(BIGDFT_BUILDDIR)/../bigdft/src         -lbigdft-1 -labinit -lxcf90 -lxc 
LIBS_DEPENDENCIES = -L$(BIGDFT_BUILDDIR)/../bigdft-revised/src -lbigdft-1 -labinit -lxcf90 -lxc 

LIBS = -L${libdir} -lbigdft-1 -lbabel  -labinit -lxcf90 -lxc -lGaIn -larchive -lCheSS-1 -lPSolver-1 \
       -L/home/chou/software/BigDFT-1.8.3/build/install/lib \
       -L${MKL_PATH} -lfutile-1 -lmkl_rt -lpthread -lm -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lyaml -ldl -lrt

all: toy_model.x
	${MAKE} toy_model.x

toy_model.x: toy_model.f90
	$(FC) -o $@ $(FCFLAGS) $(INCLUDES) toy_model.f90 $(LIBS_DEPENDENCIES) $(LIBS)
