#----------------------MAKEFILE----------------------------------------#
#----------------------------------------------------------------------#
#
#  Nuno de Sousa
#
#----------------------------------------------------------------------#
#
#Number of cores grep -c processor /proc/cpuinfo
# cat /proc/cpuinfo
#ldconfig -p | grep libjpeg
#-lcomplex_bessel for the inclusion of the a
#Number of cores
#SYS_CORES=$(shell system_profiler SPHardwareDataType | grep 'Total Number of Cores' | cut -d: -f2 | tr -d ' ')
#----------------------------------------------------------------------#


#Compilation and system commands
CC=g++
CCi=icc
CCHADES=/opt/bin/g++
CCPC185=g++
INPUT=main.cpp
OUTPUT=main
UNAME:=$(shell uname)
HOSTNAME:=$(shell hostname)
FLAGS:=
SYSTEM_INFO=$(shell system_profiler SPHardwareDataType)


#——————————————————————Cluster————————————————————————
#—————————————————————————————————————————————————————
CLUSTER_ARMADILLO_LIBS=-L/home/nsantos/libs/armadillo/lib/
CLUSTER_ARMADILLO_INCLUDE=-I/home/nsantos/libs/armadillo/include/
CLUSTER_ARMADILLO_FLAGS=-larmadillo
CLUSTER_GENERAL_FLAGS=-fopenmp  -std=c++11 -Wl,--start-group
CLUSTER_INTEL_STATIC=/opt/intel/composer_xe_2015.1.133/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2015.1.133/mkl/lib/intel64/libmkl_core.a /opt/intel/composer_xe_2015.1.133/mkl/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -L/opt/intel/composer_xe_2015.1.133/compiler/lib/intel64
CLUSTER_INTEL_FLAGS=-liomp5 -ldl -lpthread -lm -msse2
#—————————————————————————————————————————————————————



#————————————————————————Pc185————————————————————————
#—————————————————————————————————————————————————————
PC185_LIBS=-L/home/nsantos/libs/armadillo/lib/
PC185_ARMADILLO_INCLUDE=-I/home/nsantos/libs/armadillo/include/
PC185_ARMADILLO_FLAGS=-larmadillo
PC185_GENERAL_FLAGS=-fopenmp  -std=c++11 -Wl,--start-group
PC185_INTEL_STATIC=/opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_core.a /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -L/opt/intel/composer_xe_2013.3.163/mkl/lib/intel64
PC185_INTEL_FLAGS=-liomp5 -ldl -lpthread -lm -msse2

#——————————————————————————————————————————————————————

ifeq ($(UNAME),Linux)
    CC=g++
    SYSTEM_FLAGS=-O3 -march=native -std=c++11
	LLAPACK_LIBS=-llapack -larmadillo -lcomplex_bessel -fopenmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
	SYS_CORES=$(shell grep -c processor /proc/cpuinfo)
	SYS_PROCESSOR=
endif

ifeq ($(UNAME),Darwin)
    CC=g++
    SYSTEM_FLAGS=-ffast-math -march=native -std=c++11 -funroll-loops
	LLAPACK_LIBS=-larmadillo -framework Accelerate -march=native
    LCBESSEL_LIBS=
	LARMADILLO=-I/usr/include/ -I/usr/local/include/
	SYS_CORES=$(shell sysctl -n hw.ncpu)
	SYS_PROCESSOR=$(shell sysctl -n machdep.cpu.brand_string)
endif

all: main

main:
	$(info ************* COMPILATION INITIATED ************)
	$(info System processor: $(SYS_PROCESSOR))
	$(info Number of logical cores: $(SYS_CORES))
	$(info Operating system: $(UNAME))
	$(CC) -o $(OUTPUT) $(INPUT) $(LLAPACK_LIBS) $(LARMADILLO) $(LCBESSEL_LIBS) $(SYSTEM_FLAGS)
	$(info ************* COMPILATION FINISHED *************)

linux:
	$(info ************  COMPILATION INITIATED ************)
	$(info Number of cores available: $(SYS_CORES))
	$(info Operating system: $(UNAME))
	$(CCi) -o $(OUTPUT) $(INPUT) $(LLAPACK_LIBS_LINUX) $(CPU_NUMBER)
	$(info ************* COMPILATION FINISHED *************)

hades:
	$(info ************  COMPILATION INITIATED ************)
	$(info ********************HADES***********************)
	$(CCHADES) -o $(OUTPUT) $(INPUT) $(CLUSTER_ARMADILLO_LIBS) $(CLUSTER_ARMADILLO_INCLUDE) $(CLUSTER_ARMADILLO_FLAGS) $(CLUSTER_GENERAL_FLAGS) $(CLUSTER_INTEL_STATIC) $(CLUSTER_INTEL_FLAGS)
	$(info ************* COMPILATION FINISHED *************)

pc185:
	$(info ************  COMPILATION INITIATED ************)
	$(info ********************Pc185***********************)
	$(CCPC185) -o $(OUTPUT) $(INPUT) $(PC185_ARMADILLO_LIBS) $(PC185_ARMADILLO_INCLUDE) $(PC185_ARMADILLO_FLAGS) $(PC185_GENERAL_FLAGS) $(PC185_INTEL_STATIC) $(PC185_INTEL_FLAGS)
	$(info ************* COMPILATION FINISHED *************)
	

info:
	$(info ************** System Information **************)
	$(info $(SYSTEM_INFO))
	$(info ************************************************)

static:

clean:
	rm main results/*.dat screen.out script.* results/*.png

start:
	rm results/*Epfield.dat results/*sections_*.dat results/*projection* script.* results/ results/*polarizations_* results/*.png results/*.dat
