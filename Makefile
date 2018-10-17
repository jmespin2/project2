EXENAME = isingmodel
OBJS = main.o

CXX = g++ 
CXXFLAGS = -c


isingmodel: main.o input.h
	$(CXX) -o isingmodel main.o

main.o : main.cpp input.h
	$(CXX) $(CXXFLAGS) main.cpp

clean :
	-rm -f *.o $(EXENAME)