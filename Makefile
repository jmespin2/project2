EXENAME = isingmodel
OBJS = main.o

CXX = clang++
CXXFLAGS = -g -std=c++11 -Wall -pedantic


isingmodel: main.o
	$(CXX) -o isingmodel main.o

main.o : main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean :
	-rm -f *.o $(EXENAME)