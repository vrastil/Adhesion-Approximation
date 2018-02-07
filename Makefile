# all user-specific paths (headers, libraries) should be set in env. variables, i.e.
# export CPATH=$CPATH:/path/to/non-standars/headers/include
# export LIBRARY_PATH=/path/to/non-standars/libraries/lib

PRECISION = 1

CXXFLAGS =-std=c++11 -pipe
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp -flto -fPIC
#CXXFLAGS +=-D CORR
CXXFLAGS +=-D NOISE_HALF

CXXLIB +=-lboost_program_options -lboost_filesystem -lboost_system
ifeq ($(PRECISION),1)
	CXXLIBP =-lfftw3f -lfftw3f_omp
else ifeq ($(PRECISION),2)
	CXXLIBP =-lfftw3 -lfftw3_omp
else ifeq ($(PRECISION),3)
	CXXLIBP =-lfftw3l -lfftw3l_omp
endif
CXXLIB +=-lgsl -lgslcblas
CXXLIB +=-lccl

COMPILE.cc = $(CXX) $(CXXFLAGS) -c -I./include
COMPILE.fin = $(CXX) $(CXXFLAGS) $(CXXLIB_PATH)

OBJ_FILES = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp) $(wildcard src/**/*.cpp))
TEST_OBJ_FILES = $(patsubst src/%,tests/%, $(filter-out src/main.o,$(OBJ_FILES))) tests/test_main.o
LIB = lib/fastsim.a
PCH = include/stdafx.h
PCH_O = $(PCH).gch

adh_app: CXXFLAGS +=-Ofast -march=native -D PRECISION=$(PRECISION)
adh_app: CXXLIB += $(CXXLIBP)
adh_app: $(OBJ_FILES)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)
	
debug: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion -D PRECISION=$(PRECISION)
debug: CXXLIB += $(CXXLIBP)
debug: $(OBJ_FILES)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

# SWIG wrapper build always in double precision, single precision not working now with scipy optimization
swig: CXXFLAGS +=-Ofast -march=native -D LESSINFO -D PRECISION=2
swig: CXXLIB += -lfftw3 -lfftw3_omp
swig: $(LIB) swig/*.i
	swig -python -c++ -DPRECISION=2 -I/usr/local/include/ -I./include -o swig/swig_wrap.cpp swig/all.i
	$(COMPILE.cc) -c -I/usr/include/python2.7 -o swig/swig_wrap.o swig/swig_wrap.cpp
	$(COMPILE.fin) -shared -o swig/_fastsim.so swig/swig_wrap.o $(LIB) $(CXXLIB)
	ln -sf ${CURDIR}/swig/_fastsim.so simpy/
	ln -sf ${CURDIR}/swig/_fastsim.so lib/
	ln -sf ${CURDIR}/swig/fastsim.py simpy/

$(LIB): $(OBJ_FILES)
	gcc-ar rcs $@ $^

check: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion -D TEST
check: COMPILE.cc += -I./tests
check: tests/test
	./tests/test

tests/test: $(TEST_OBJ_FILES)
	$(COMPILE.fin) -o tests/test $^ $(CXXLIB)

%.o: %.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

tests/%.o: src/%.cpp
	$(COMPILE.cc) -o $@ $<

$(PCH_O): $(PCH)
	$(CXX) $(CXXFLAGS) -o $@ $<

# keep libraries
clean:
	find . -type f \( -name *~ -o -name *.o -o -name *.d -o -name *.gch \) -exec rm -f {} \;
	rm -f swig/*.cpp adh_app debug tests/test

cleanall: clean
		rm -f swig/*.py swig/*.pyc swig/*.so lib/*

-include $(OBJ_FILES:.o=.d)
-include $(TEST_OBJ_FILES:.o=.d)
-include $(PCH_O:.gch=.d)

.PHONY: swig check clean
