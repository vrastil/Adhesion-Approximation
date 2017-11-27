# all user-specific paths (headers, libraries) should be set in env. variables, i.e.
# export CPATH=$CPATH:/path/to/non-standars/headers/include
# export LIBRARY_PATH=/path/to/non-standars/libraries/lib

CXXFLAGS =-std=c++11 -pipe
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp -flto -fPIC
#CXXFLAGS +=-D CORR
CXXFLAGS +=-D NOISE_HALF

CXXLIB +=-lboost_program_options -lboost_filesystem -lboost_system
CXXLIB +=-lfftw3 -lfftw3_omp
CXXLIB +=-lgsl -lgslcblas
CXXLIB +=-lccl

COMPILE.cc = $(CXX) $(CXXFLAGS) -c -I./include
COMPILE.fin = $(CXX) $(CXXFLAGS) $(CXXLIB_PATH)

OBJ_FILES = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
TEST_OBJ_FILES = $(patsubst src/%,tests/%, $(filter-out src/main.o,$(OBJ_FILES))) tests/test_main.o
LIB = lib/fastsim.a
PCH = include/stdafx.h
PCH_O = $(PCH).gch

adh_app: CXXFLAGS +=-Ofast -march=native
adh_app: $(LIB)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

debug: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion
debug: $(LIB)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

swig: CXXFLAGS +=-Ofast -march=native
swig: $(LIB) swig/swig.i
	swig -python -c++ -I./include swig/swig.i
	$(COMPILE.cc) -c -I/usr/include/python2.7 -o swig/swig_wrap.o swig/swig_wrap.cxx
	$(COMPILE.fin) -shared -o swig/_fastsim.so swig/swig_wrap.o $(LIB) $(CXXLIB)
	ln -sf ${CURDIR}/swig/_fastsim.so simpy/
	ln -sf ${CURDIR}/swig/_fastsim.so lib/
	ln -sf ${CURDIR}/swig/fastsim.py simpy/

$(LIB): $(OBJ_FILES)
	gcc-ar rcs $@ $^

src/%.o: src/%.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

check: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion
check: $(TEST_OBJ_FILES)
	$(COMPILE.fin) -o tests/test $^ $(CXXLIB)
	./tests/test

tests/%.o: src/%.cpp
	$(COMPILE.cc) -I./tests -D TEST -o $@ $<

$(PCH_O): $(PCH)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f src/*~ src/*.o src/*.d \
		swig/*~ swig/*.o swig/*.d \
		tests/*~ tests/*.o tests/*.d \
		include/*~ include/*.o include/*.d include/*.gch \
		lib/* adh_app

-include $(OBJ_FILES:.o=.d)
-include $(TEST_OBJ_FILES:.o=.d)
-include $(PCH_O:.gch=.d)
