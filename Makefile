# all user-specific paths (headers, libraries) should be set in env. variables, i.e.
# export CPATH=$CPATH:/path/to/non-standars/headers/include
# export LIBRARY_PATH=/path/to/non-standars/libraries/lib

PRECISION = 2

CXXFLAGS =-std=c++11 -pipe
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp -flto=jobserver -fPIC
#CXXFLAGS +=-D CORR
CXXFLAGS +=-D NOISE_HALF
CXXFLAGS +=-D OPENMP # for multigrid_solver

# used boost libraries
CXXLIB +=-lboost_program_options -lboost_filesystem -lboost_system

# fftw and omp libraries depend on define precision (single / double / long double)
ifeq ($(PRECISION),1)
	CXXLIBP =-lfftw3f -lfftw3f_omp
else ifeq ($(PRECISION),2)
	CXXLIBP =-lfftw3 -lfftw3_omp
else ifeq ($(PRECISION),3)
	CXXLIBP =-lfftw3l -lfftw3l_omp
endif

# gsl and ccl libraries
CXXLIB +=-lgsl -lgslcblas
CXXLIB +=-lccl

COMPILE.cc = $(CXX) $(CXXFLAGS) -c -I./include# -D_GLIBCXX_USE_CXX11_ABI=0
COMPILE.fin = $(CXX) $(CXXFLAGS) $(CXXLIB_PATH)

# get list of object files for each cpp file in src/, src/**/
OBJ_FILES := $(wildcard src/*.cpp) $(wildcard src/**/*.cpp)
OBJ_FILES := $(OBJ_FILES:.cpp=.o)

# get list of object files for each cpp file in tests/, tests/**/
TEST_OBJ_FILES := $(wildcard tests/*.cpp) $(wildcard tests/**/*.cpp)
TEST_OBJ_FILES := $(TEST_OBJ_FILES:.cpp=.o)

# add the rest of object files that are not already included
TEST_OBJ_FILES += $(filter-out $(subst test_,,$(TEST_OBJ_FILES)), $(patsubst src/%.o, tests/%.o, $(OBJ_FILES)))

LIB = lib/fastsim.a
PCH = include/stdafx.h
PCH_O = $(PCH).gch

adh_app: CXXFLAGS +=-Ofast -march=native -D PRECISION=$(PRECISION)
adh_app: CXXLIB += $(CXXLIBP)
adh_app: $(OBJ_FILES)
	+$(COMPILE.fin) -o $@ $^ $(CXXLIB)
	
debug: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion -D PRECISION=$(PRECISION)
debug: CXXLIB += $(CXXLIBP)
debug: $(OBJ_FILES)
	+$(COMPILE.fin) -o $@ $^ $(CXXLIB)

# SWIG wrapper build always in double precision, single precision not working now with scipy optimization
swig: CXXFLAGS +=-Ofast -march=native -D LESSINFO -D PRECISION=2
swig: CXXLIB += -lfftw3 -lfftw3_omp
swig: $(LIB) swig/*.i
	swig -python -c++ -DPRECISION=2 -I/usr/local/include/ -I./include -o swig/swig_wrap.cpp swig/all.i
	$(COMPILE.cc) -c -I/usr/include/python2.7 -o swig/swig_wrap.o swig/swig_wrap.cpp
	+$(COMPILE.fin) -fuse-linker-plugin -shared -o swig/_fastsim.so swig/swig_wrap.o $(LIB) $(CXXLIB)
	ln -sf ${CURDIR}/swig/_fastsim.so simpy/
	ln -sf ${CURDIR}/swig/_fastsim.so lib/
	ln -sf ${CURDIR}/swig/fastsim.py simpy/

$(LIB): $(OBJ_FILES)
	gcc-ar rcs $@ $^

doc:
	cd doc/ && doxygen Doxyfile && ln -sf html/index.html main.html

mpi: CXX = mpic++ -D MPI_ENABLED
mpi: CXXLIB += -lboost_mpi
mpi: adh_app

check: test
	./tests/test

test: tests/test
test: CXXFLAGS +=-Og -g -w -D PRECISION=$(PRECISION)
test: CXXLIB += $(CXXLIBP)
test: COMPILE.cc += -I./src


tests/test: $(TEST_OBJ_FILES)
	+$(COMPILE.fin) -o tests/test $^ $(CXXLIB)

# generic rule for building object files
%.o: %.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

# specific rule for building test object files which DO have its own source file
tests/%.o: tests/%.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

# specific rule for building test object files which do NOT have its own source file
tests/%.o: src/%.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

$(PCH_O): $(PCH)
	$(CXX) $(CXXFLAGS) -nostartfiles -o $@ $<

# keep libraries
clean:
	find . -type f \( -name *~ -o -name *.o -o -name *.d -o -name *.gch \) -exec rm -f {} \;
	rm -f swig/*.cpp adh_app debug tests/test

cleanall: clean
		rm -f swig/*.py swig/*.pyc swig/*.so lib/*
		rm -rf doc/html doc/latex doc/main.html

-include $(OBJ_FILES:.o=.d)
-include $(TEST_OBJ_FILES:.o=.d)
-include $(PCH_O:.gch=.d)

.PHONY: swig check clean test doc mpi
