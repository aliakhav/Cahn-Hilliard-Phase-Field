CXX=g++

CXXFLAGS=-O3 -Wall -c
#CXXFLAGS=-O0 -Wall -c
#CXXFLAGS=-O0 -g -Wall -c

OBJS = PhaseField.o Derivatives.o RND.o


all: PhaseField


PhaseField: $(OBJS) 
	$(CXX) -o $@ $(OBJS)

PhaseField.o: PhaseField.C
	$(CXX) $(CXXFLAGS) -o $@ $?

Derivatives.o: Derivatives.C
	$(CXX) $(CXXFLAGS) -o $@ $?

RND.o: RND.C
	$(CXX) $(CXXFLAGS) -o $@ $?

clean:
	rm -rf *.o PhaseField

	
