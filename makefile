CXX = g++
RM=rm -f
CPPFLAGS=-I/usr/include/ -O2
LDFLAGS=-O2
LDLIBS=-lblas

SRCS=main.cpp \
	AngleHistogram.cpp \
	BabyUniverseDistribution.cpp \
	CirclePattern.cpp \
	CohomologyBasis.cpp \
	ConjugateGradient.cpp \
	ConnectivityRestrictor.cpp \
	DualCohomologyBasis.cpp \
	Edge.cpp \
	Embedding.cpp \
	HarmonicEmbedding.cpp \
	LaplacianDeterminant.cpp \
	ModuliObservable.cpp \
	PottsModel.cpp \
	ShortestLoop.cpp \
	Simulation.cpp \
	ThetaHistogram.cpp \
	ThetaModel.cpp \
	Triangle.cpp \
	Triangulation.cpp \
	Vertex.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: dt

dt: $(OBJS)
	g++ $(LDFLAGS) -o dt $(OBJS) $(LDLIBS) 

AngleHistogram.o: AngleHistogram.cpp AngleHistogram.h 
BabyUniverseDistribution.o: BabyUniverseDistribution.cpp BabyUniverseDistribution.h
CirclePattern.o: CirclePattern.cpp CirclePattern.h
CohomologyBasis.o: CohomologyBasis.cpp CohomologyBasis.h
ConjugateGradient.o: ConjugateGradient.cpp ConjugateGradient.h
ConnectivityRestrictor.o: ConnectivityRestrictor.cpp ConnectivityRestrictor.h
DualCohomologyBasis.o: DualCohomologyBasis.cpp DualCohomologyBasis.h
Edge.o: Edge.cpp Edge.h
Embedding.o: Embedding.cpp Embedding.h
HarmonicEmbedding.o: HarmonicEmbedding.cpp HarmonicEmbedding.h
LaplacianDeterminant.o: LaplacianDeterminant.cpp LaplacianDeterminant.h
ModuliObservable.o: ModuliObservable.cpp ModuliObservable.h
PottsModel.o: PottsModel.cpp PottsModel.h
ShortestLoop.o: ShortestLoop.cpp ShortestLoop.h
Simulation.o: Simulation.cpp Simulation.h
ThetaHistogram.o: ThetaHistogram.cpp ThetaHistogram.h
ThetaModel.o: ThetaModel.cpp ThetaModel.h
Triangle.o: Triangle.cpp Triangle.h
Triangulation.o: Triangulation.cpp Triangulation.h
Vertex.o: Vertex.cpp Vertex.h

tool.o: tool.cc support.hh

support.o: support.hh support.cc

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) dt