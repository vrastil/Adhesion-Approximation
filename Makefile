
CXX = g++

CXXFLAGS +=-g -std=c++11 -Wall
#CXXFLAGS +=-std=c++11 -O3
CXXFLAGS +=-MMD
CXXFLAGS +=-fopenmp
CXXLIB_PATH +=-L/usr/local/lib/
# -L/lib/

CXXLIB +=-lboost_program_options -lboost_filesystem -lboost_system
CXXLIB +=-lfftw3 -lfftw3_omp
CXXLIB +=-lgsl -lgslcblas

OBJ_FILES = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

COMPILE.cc = $(CXX) $(CXXFLAGS) -c -I./include
COMPILE.fin = $(CXX) $(CXXFLAGS) $(CXXLIB_PATH)

adh_app: $(OBJ_FILES)
	$(COMPILE.fin) -o $@ $^ $(CXXLIB)

src/%.o: src/%.cpp
	$(COMPILE.cc) -o $@ $<

precompiled: src/stdafx.h
	$(COMPILE.cc) src/stdafx.h -o src/stdafx.h.gch

clean:
	rm -f src/*.o src/*~ src/*.d

-include $(OBJ_FILES:.o=.d)
