CXX = g++
RM=rm -f
CPPFLAGS=-I/usr/include/ -O2
LDFLAGS=-O2
LDLIBS=-lblas

SRCS=main.cpp \
	Edge.cpp \
	Simulation.cpp \
	Triangle.cpp \
	Triangulation.cpp \
	Vertex.cpp \
	Diffusion.cpp \
	LinearAlgebra.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: dt-diffusion

dt-diffusion: $(OBJS)
	g++ $(LDFLAGS) -o dt-diffusion $(OBJS) $(LDLIBS) 

Edge.o: Edge.cpp Edge.h
Simulation.o: Simulation.cpp Simulation.h
Triangle.o: Triangle.cpp Triangle.h
Triangulation.o: Triangulation.cpp Triangulation.h
Vertex.o: Vertex.cpp Vertex.h
Diffusion.o: Diffusion.cpp Diffusion.h
LinearAlgebra.o: LinearAlgebra.cpp LinearAlgebra.h

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) dt