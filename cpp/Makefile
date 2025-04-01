
CXX = g++
OPT = -O3 -DNDEBUG
CXXFLAGS = -Wall -std=c++17 $(OPT)

OBJECTS = $(SOURCES:.cc=.o) 

TARGETS = boundary isohedral isohedral_e2e isohedral_hex_e2e

BOUNDARY_SOURCES = boundary_test.cpp boundary.cpp
BOUNDARY_OBJECTS = $(BOUNDARY_SOURCES:.cpp=.o)
BOUNDARY_DEPENDS = $(BOUNDARY_OBJECTS:.o=.d)

ISOHEDRAL_SOURCES = isohedral.cpp isohedral_test.cpp boundary.cpp
ISOHEDRAL_OBJECTS = $(ISOHEDRAL_SOURCES:.cpp=.o)
ISOHEDRAL_DEPENDS = $(ISOHEDRAL_OBJECTS:.o=.d)

ISOHEDRAL_E2E_SOURCES = isohedral.cpp isohedral_e2e_test.cpp boundary.cpp
ISOHEDRAL_E2E_OBJECTS = $(ISOHEDRAL_E2E_SOURCES:.cpp=.o)
ISOHEDRAL_E2E_DEPENDS = $(ISOHEDRAL_E2E_OBJECTS:.o=.d)

ISOHEDRAL_HEX_E2E_SOURCES = isohedral.cpp isohedral_hex_e2e_test.cpp boundary.cpp
ISOHEDRAL_HEX_E2E_OBJECTS = $(ISOHEDRAL_HEX_E2E_SOURCES:.cpp=.o)
ISOHEDRAL_HEX_E2E_DEPENDS = $(ISOHEDRAL_HEX_E2E_OBJECTS:.o=.d)

ISOHEDRAL_IAMOND_E2E_SOURCES = isohedral.cpp isohedral_iamond_e2e_test.cpp boundary.cpp
ISOHEDRAL_IAMOND_E2E_OBJECTS = $(ISOHEDRAL_IAMOND_E2E_SOURCES:.cpp=.o)
ISOHEDRAL_IAMOND_E2E_DEPENDS = $(ISOHEDRAL_IAMOND_E2E_OBJECTS:.o=.d)

SOURCES = boundary.cpp boundary_test.cpp isohedral.cpp isohedral_test.cpp isohedral_e2e_test.cpp isohedral_hex_e2e_test.cpp isohedral_iamond_e2e_test.cpp
OBJECTS = $(SOURCES:.cpp=.o)
DEPENDS = $(OBJECTS:.o=.d)


all: $(TARGETS)

boundary: $(BOUNDARY_OBJECTS)
	$(CXX) $(CXXFLAGS) $(BOUNDARY_OBJECTS) -o boundary

isohedral: $(ISOHEDRAL_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ISOHEDRAL_OBJECTS) -o isohedral

isohedral_e2e: $(ISOHEDRAL_E2E_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ISOHEDRAL_E2E_OBJECTS) -o isohedral_e2e

isohedral_hex_e2e: $(ISOHEDRAL_HEX_E2E_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ISOHEDRAL_HEX_E2E_OBJECTS) -o isohedral_hex_e2e

isohedral_iamond_e2e: $(ISOHEDRAL_IAMOND_E2E_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ISOHEDRAL_IAMOND_E2E_OBJECTS) -o isohedral_iamond_e2e

.PHONY: clean

clean:
	rm ${OBJECTS} ${DEPENDS} ${TARGETS} *.txt
