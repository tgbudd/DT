CXX = g++
RM=rm -f
CPPFLAGS=-I/usr/include/ -O2 -D BOOST_DISABLE_ASSERTS
LDFLAGS=-O2
LDLIBS=-lblas

SRCS=main.cpp \
	Edge.cpp \
	Simulation.cpp \
	Triangle.cpp \
	Triangulation.cpp \
	Vertex.cpp \
	Diffusion.cpp \
	CirclePacking.cpp \
	CohomologyBasis.cpp \
	DualCohomologyBasis.cpp \
	ConnectivityRestrictor.cpp \
	ConformalDistribution.cpp \
	LinearAlgebra.cpp \
	BabyUniverseDetector.cpp \
	TriangulationProperties.cpp \
	Embedding.cpp \
	LaplacianMatrix.cpp \
	ShortestLoop.cpp \
	CMinusTwoBuilder.cpp \
	BabyUniverseRemover.cpp \
	HarmonicEmbedding.cpp




OBJS=$(subst .cpp,.o,$(SRCS))

all: dt-conf

dt-conf: $(OBJS)
	g++ $(LDFLAGS) -o ./hetws9/dt-conf $(OBJS) $(LDLIBS) 

Edge.o: Edge.cpp Edge.h
Simulation.o: Simulation.cpp Simulation.h
Triangle.o: Triangle.cpp Triangle.h
Triangulation.o: Triangulation.cpp Triangulation.h
Vertex.o: Vertex.cpp Vertex.h
Diffusion.o: Diffusion.cpp Diffusion.h
CohomologyBasis.o: CohomologyBasis.cpp CohomologyBasis.h
DualCohomologyBasis.o: DualCohomologyBasis.cpp DualCohomologyBasis.h
ConnectivityRestrictor.o: ConnectivityRestrictor.cpp ConnectivityRestrictor.h
ConformalDistribution.o: ConformalDistribution.cpp ConformalDistribution.h
LinearAlgebra.o: LinearAlgebra.cpp LinearAlgebra.h
BabyUniverseDetector.o: BabyUniverseDetector.cpp BabyUniverseDetector.h
TriangulationProperties.o: TriangulationProperties.cpp TriangulationProperties.h
Embedding.o: Embedding.cpp Embedding.h
LaplacianMatrix.o: LaplacianMatrix.cpp LaplacianMatrix.h
ShortestLoop.o: ShortestLoop.cpp ShortestLoop.h
CMinusTwoBuilder.o: CMinusTwoBuilder.cpp CMinusTwoBuilder.h
BabyUniverseRemover.o: BabyUniverseRemover.cpp BabyUniverseRemover.h
HarmonicEmbedding.o: HarmonicEmbedding.cpp HarmonicEmbedding.h

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) ./hetws9/dt-conf