
#CXX = g++
CXX = g++-6.4
#CXX = g++-7.2

CXXFLAGS =-std=c++11 -pipe
#CXXFLAGS +=-Og -g -Wall
CXXFLAGS +=-Ofast -march=native
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp

CXXFLAGS +=-D CORR
#CXXFLAGS +=-D REAL_NOISE
#CXXFLAGS +=-D FFTW_SYM
#CXXFLAGS +=-D OLD_NORM

CXXLIB_PATH +=-L/usr/local/lib/
# -L/lib/

CXXLIB +=-lboost_program_options -lboost_filesystem -lboost_system
CXXLIB +=-lfftw3 -lfftw3_omp
CXXLIB +=-lgsl -lgslcblas
CXXLIB +=-lccl

OBJ_FILES = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
TEST_OBJ_FILES = $(patsubst src/%,tests/%, $(filter-out src/main.o,$(OBJ_FILES))) tests/test_main.o
PCH = include/stdafx.h
PCH_O = $(PCH).gch

COMPILE.cc = $(CXX) $(CXXFLAGS) -c -I./include
COMPILE.fin = $(CXX) $(CXXFLAGS) $(CXXLIB_PATH)

adh_app: $(OBJ_FILES)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

src/%.o: src/%.cpp $(PCH_O)
	$(COMPILE.cc) -o $@ $<

$(PCH_O): $(PCH)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f src/*.o src/*~ src/*.d tests/*.o tests/*~ tests/*.d include/*.gch include/*~ include/*.o

check : test

test: $(TEST_OBJ_FILES)
	$(COMPILE.fin) -o tests/test $^ $(CXXLIB)
	./tests/test

tests/%.o: src/%.cpp
	$(COMPILE.cc) -I./tests -D TEST -o $@ $<

-include $(OBJ_FILES:.o=.d)
-include $(TEST_OBJ_FILES:.o=.d)
-include $(PCH_O:.gch=.d)
