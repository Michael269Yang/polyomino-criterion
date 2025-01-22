CXX = g++
CXXFLAGS = -Wall -std=c++17

all: isohedral

isohedral: isohedral.o isohedral.h
	$(CXX) $(CXXFLAGS) isohedral.o -o isohedral

isohedral_test: isohedral_test.o isohedral.h 
	$(CXX) $(CXXFLAGS) isohedral_test.o -o isohedral_test

main: main.o isohedral.h 
	$(CXX) $(CXXFLAGS) main.o -o main

clean:
	rm -f isohedral.o main.o isohedral_test.o

