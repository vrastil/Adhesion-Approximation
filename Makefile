# all user-specific paths (headers, libraries) should be set in env. variables, i.e.
# export CPATH=$CPATH:/path/to/non-standars/headers/include
# export LIBRARY_PATH=/path/to/non-standars/libraries/lib

CXXFLAGS =-std=c++11 -pipe
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp
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
PCH = include/stdafx.h
PCH_O = $(PCH).gch

all: CXXFLAGS +=-Ofast -march=native
all: adh_app

debug: CXXFLAGS +=-Og -g -Wall -Wunused-parameter -Wfloat-conversion
debug: adh_app

adh_app: $(OBJ_FILES)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

src/%.o: src/%.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

check: CXXFLAGS +=-Og -g -Wall
check: $(TEST_OBJ_FILES)
	$(COMPILE.fin) -o tests/test $^ $(CXXLIB)
	./tests/test

tests/%.o: src/%.cpp
	$(COMPILE.cc) -I./tests -D TEST -o $@ $<

$(PCH_O): $(PCH)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f src/*.o src/*~ src/*.d tests/*.o tests/*~ tests/*.d include/*.gch include/*~ include/*.o include/*.d

-include $(OBJ_FILES:.o=.d)
-include $(TEST_OBJ_FILES:.o=.d)
-include $(PCH_O:.gch=.d)
