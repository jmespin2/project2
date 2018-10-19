EXENAME = isingmodel
OBJS = main.o
HEADERS = input.h

CXX = g++
CXXFLAGS = -g


$(EXENAME): $(OBJS)
	$(CXX) -o isingmodel main.o

main.o : main.cpp $(HEADERS)
	$(CXX) -g -c main.cpp

clean :
	-rm -f *.o $(EXENAME)