
CXX = g++
CXXFLAGS = -Wall -std=c++17

OBJECTS = $(SOURCES:.cc=.o) 

TARGETS = boundary isohedral

BOUNDARY_SOURCES = boundary.cpp
BOUNDARY_OBJECTS = $(BOUNDARY_SOURCES:.cpp=.o)
BOUNDARY_DEPENDS = $(BOUNDARY_OBJECTS:.o=.d)

ISOHEDRAL_SOURCES = isohedral.cpp isohedral_test.cpp
ISOHEDRAL_OBJECTS = $(ISOHEDRAL_SOURCES:.cpp=.o)
ISOHEDRAL_DEPENDS = $(ISOHEDRAL_OBJECTS:.o=.d)

SOURCES = boundary.cpp isohedral.cpp isohedral_test.cpp
OBJECTS = $(SOURCES:.cpp=.o)
DEPENDS = $(OBJECTS:.o=.d)


all: $(TARGETS)

boundary: $(BOUNDARY_OBJECTS)
	$(CXX) $(CXXFLAGS) $(BOUNDARY_OBJECTS) -o boundary

isohedral: $(ISOHEDRAL_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ISOHEDRAL_OBJECTS) -o isohedral

.PHONY: clean

clean:
	rm ${OBJECTS} ${DEPENDS} ${TARGETS}
