CXX = g++
RM=rm -f
CPPFLAGS=-I/usr/include/ -O2 -D BOOST_DISABLE_ASSERTS
LDFLAGS=-O2
LDLIBS=-lblas

CORESCRS=AngleHistogram.cpp \
	BabyUniverseDetector.cpp \
	BabyUniverseDistribution.cpp \
	BabyUniverseRemover.cpp \
	BallSizeDistribution.cpp \
	BitmapDrawer.cpp \
	BoundaryMeasure.cpp \
	CMinusTwoBuilder.cpp \
	CirclePacking.cpp \
	CirclePattern.cpp \
	CohomologyBasis.cpp \
	ConformalDistribution.cpp \
	ConnectivityRestrictor.cpp \
	Diffusion.cpp \
	DiskCirclePacking.cpp \
	DistanceDistribution.cpp \
	DualCohomologyBasis.cpp \
	DualDiffusion.cpp \
	DualRandomWalk.cpp \
	DualScalarField.cpp \
	Edge.cpp \
	Embedding.cpp \
	HarmonicDiffusion.cpp \
	HarmonicEmbedding.cpp \
	Histogram.cpp \
	HyperbolicStructure.cpp \
	LaplacianDeterminant.cpp \
	LaplacianMatrix.cpp \
	LaplacianSpectrum.cpp \
	LinearAlgebra.cpp \
	LoopLength.cpp \
	MetricGraphObservable.cpp \
	ModuliObservable.cpp \
	PeelingProcedure.cpp \
	PottsModel.cpp \
	PottsModelOnVertices.cpp \
	ShortestLoop.cpp \
	Simulation.cpp \
	SpanningTree.cpp \
	ThetaHistogram.cpp \
	ThetaModel.cpp \
	Triangle.cpp \
	Triangulation.cpp \
	TriangulationProperties.cpp \
	Vertex.cpp

COREOBJS=$(subst .cpp,.o,$(CORESRCS))

all: core dt-ball

core: $(COREOBJS)

ball_sizes_torus: Examples/ball_sizes_torus.o core
	 $(CXX) $(LDFLAGS) -o ./bin/ball_sizes_torus ./Examples/ball_sizes_torus.o $(COREOBJS) $(LDLIBS)

Examples/ball_sizes_torus.o:Examples/ball_sizes_torus.cpp

boundary_measure_on_cminus2_disk: Examples/_disk.o core
	 $(CXX) $(LDFLAGS) -o ./bin/_disk ./Examples/_disk.o $(COREOBJS) $(LDLIBS)

Examples/_disk.o:Examples/_disk.cpp

circle_packing_conformal_distribution: Examples/circle_packing_conformal_distribution.o core
	 $(CXX) $(LDFLAGS) -o ./bin/circle_packing_conformal_distribution ./Examples/circle_packing_conformal_distribution.o $(COREOBJS) $(LDLIBS)

Examples/circle_packing_conformal_distribution.o:Examples/circle_packing_conformal_distribution.cpp

circle_packing_peeling_process_sphere: Examples/circle_packing_peeling_process_sphere.o core
	 $(CXX) $(LDFLAGS) -o ./bin/circle_packing_peeling_process_sphere ./Examples/circle_packing_peeling_process_sphere.o $(COREOBJS) $(LDLIBS)

Examples/circle_packing_peeling_process_sphere.o:Examples/circle_packing_peeling_process_sphere.cpp

circle_packing_peeling_process_torus: Examples/circle_packing_peeling_process_torus.o core
	 $(CXX) $(LDFLAGS) -o ./bin/circle_packing_peeling_process_torus ./Examples/circle_packing_peeling_process_torus.o $(COREOBJS) $(LDLIBS)

Examples/circle_packing_peeling_process_torus.o:Examples/circle_packing_peeling_process_torus.cpp

conformal_distribution_circle_packing_disk: Examples/conformal_distribution_circle_packing_disk.o core
	 $(CXX) $(LDFLAGS) -o ./bin/conformal_distribution_circle_packing_disk ./Examples/conformal_distribution_circle_packing_disk.o $(COREOBJS) $(LDLIBS)

Examples/conformal_distribution_circle_packing_disk.o:Examples/conformal_distribution_circle_packing_disk.cpp

conformal_distribution_torus: Examples/conformal_distribution_torus.o core
	 $(CXX) $(LDFLAGS) -o ./bin/conformal_distribution_torus ./Examples/conformal_distribution_torus.o $(COREOBJS) $(LDLIBS)

Examples/conformal_distribution_torus.o:Examples/conformal_distribution_torus.cpp

diffusion_on_cminus2_sphere: Examples/_sphere.o core
	 $(CXX) $(LDFLAGS) -o ./bin/_sphere ./Examples/_sphere.o $(COREOBJS) $(LDLIBS)

Examples/_sphere.o:Examples/_sphere.cpp

diffusion_on_dt_coupled_to_scalar: Examples/diffusion_on_dt_coupled_to_scalar.o core
	 $(CXX) $(LDFLAGS) -o ./bin/diffusion_on_dt_coupled_to_scalar ./Examples/diffusion_on_dt_coupled_to_scalar.o $(COREOBJS) $(LDLIBS)

Examples/diffusion_on_dt_coupled_to_scalar.o:Examples/diffusion_on_dt_coupled_to_scalar.cpp

diffusion_on_dt_coupled_to_various_matter: Examples/diffusion_on_dt_coupled_to_various_matter.o core
	 $(CXX) $(LDFLAGS) -o ./bin/diffusion_on_dt_coupled_to_various_matter ./Examples/diffusion_on_dt_coupled_to_various_matter.o $(COREOBJS) $(LDLIBS)

Examples/diffusion_on_dt_coupled_to_various_matter.o:Examples/diffusion_on_dt_coupled_to_various_matter.cpp

diffusion_simulation: Examples/diffusion_simulation.o core
	 $(CXX) $(LDFLAGS) -o ./bin/diffusion_simulation ./Examples/diffusion_simulation.o $(COREOBJS) $(LDLIBS)

Examples/diffusion_simulation.o:Examples/diffusion_simulation.cpp

harmonic_measure_of_peeling_process: Examples/harmonic_measure_of_peeling_process.o core
	 $(CXX) $(LDFLAGS) -o ./bin/harmonic_measure_of_peeling_process ./Examples/harmonic_measure_of_peeling_process.o $(COREOBJS) $(LDLIBS)

Examples/harmonic_measure_of_peeling_process.o:Examples/harmonic_measure_of_peeling_process.cpp

metric_graph_geodesics: Examples/metric_graph_geodesics.o core
	 $(CXX) $(LDFLAGS) -o ./bin/metric_graph_geodesics ./Examples/metric_graph_geodesics.o $(COREOBJS) $(LDLIBS)

Examples/metric_graph_geodesics.o:Examples/metric_graph_geodesics.cpp

modulus_constrained_shortest_cycle: Examples/modulus_constrained_shortest_cycle.o core
	 $(CXX) $(LDFLAGS) -o ./bin/modulus_constrained_shortest_cycle ./Examples/modulus_constrained_shortest_cycle.o $(COREOBJS) $(LDLIBS)

Examples/modulus_constrained_shortest_cycle.o:Examples/modulus_constrained_shortest_cycle.cpp

save_torus_triangulation: Examples/save_torus_triangulation.o core
	 $(CXX) $(LDFLAGS) -o ./bin/save_torus_triangulation ./Examples/save_torus_triangulation.o $(COREOBJS) $(LDLIBS)

Examples/save_torus_triangulation.o:Examples/save_torus_triangulation.cpp

test_attributes: Examples/test_attributes.o core
	 $(CXX) $(LDFLAGS) -o ./bin/test_attributes ./Examples/test_attributes.o $(COREOBJS) $(LDLIBS)

Examples/test_attributes.o:Examples/test_attributes.cpp

thetamodel_cosine_power: Examples/thetamodel_cosine_power.o core
	 $(CXX) $(LDFLAGS) -o ./bin/thetamodel_cosine_power ./Examples/thetamodel_cosine_power.o $(COREOBJS) $(LDLIBS)

Examples/thetamodel_cosine_power.o:Examples/thetamodel_cosine_power.cpp

thetamodel_cosine_snapshots: Examples/thetamodel_cosine_snapshots.o core
	 $(CXX) $(LDFLAGS) -o ./bin/thetamodel_cosine_snapshots ./Examples/thetamodel_cosine_snapshots.o $(COREOBJS) $(LDLIBS)

Examples/thetamodel_cosine_snapshots.o:Examples/thetamodel_cosine_snapshots.cpp

thetamodel_moduli: Examples/thetamodel_moduli.o core
	 $(CXX) $(LDFLAGS) -o ./bin/thetamodel_moduli ./Examples/thetamodel_moduli.o $(COREOBJS) $(LDLIBS)

Examples/thetamodel_moduli.o:Examples/thetamodel_moduli.cpp

thetamodel_simulation_babyuniverses: Examples/thetamodel_simulation_babyuniverses.o core
	 $(CXX) $(LDFLAGS) -o ./bin/thetamodel_simulation_babyuniverses ./Examples/thetamodel_simulation_babyuniverses.o $(COREOBJS) $(LDLIBS)

Examples/thetamodel_simulation_babyuniverses.o:Examples/thetamodel_simulation_babyuniverses.cpp

thetamodel_snapshots_geodesics: Examples/thetamodel_snapshots_geodesics.o core
	 $(CXX) $(LDFLAGS) -o ./bin/thetamodel_snapshots_geodesics ./Examples/thetamodel_snapshots_geodesics.o $(COREOBJS) $(LDLIBS)

Examples/thetamodel_snapshots_geodesics.o:Examples/thetamodel_snapshots_geodesics.cpp

video_builder: Examples/video_builder.o core
	 $(CXX) $(LDFLAGS) -o ./bin/video_builder ./Examples/video_builder.o $(COREOBJS) $(LDLIBS)

Examples/video_builder.o:Examples/video_builder.cpp

video_builder_semiclassical: Examples/video_builder_semiclassical.o core
	 $(CXX) $(LDFLAGS) -o ./bin/video_builder_semiclassical ./Examples/video_builder_semiclassical.o $(COREOBJS) $(LDLIBS)

Examples/video_builder_semiclassical.o:Examples/video_builder_semiclassical.cpp

video_dt_coupled_to_spanning_tree_and_scalar: Examples/video_dt_coupled_to_spanning_tree_and_scalar.o core
	 $(CXX) $(LDFLAGS) -o ./bin/video_dt_coupled_to_spanning_tree_and_scalar ./Examples/video_dt_coupled_to_spanning_tree_and_scalar.o $(COREOBJS) $(LDLIBS)

Examples/video_dt_coupled_to_spanning_tree_and_scalar.o:Examples/video_dt_coupled_to_spanning_tree_and_scalar.cpp

video_peeling_procedure: Examples/video_peeling_procedure.o core
	 $(CXX) $(LDFLAGS) -o ./bin/video_peeling_procedure ./Examples/video_peeling_procedure.o $(COREOBJS) $(LDLIBS)

Examples/video_peeling_procedure.o:Examples/video_peeling_procedure.cpp

zoom_torus_video: Examples/zoom_torus_video.o core
	 $(CXX) $(LDFLAGS) -o ./bin/zoom_torus_video ./Examples/zoom_torus_video.o $(COREOBJS) $(LDLIBS)

Examples/zoom_torus_video.o:Examples/zoom_torus_video.cpp



AngleHistogram.o: AngleHistogram.cpp AngleHistogram.h

BabyUniverseDetector.o: BabyUniverseDetector.cpp BabyUniverseDetector.h

BabyUniverseDistribution.o: BabyUniverseDistribution.cpp BabyUniverseDistribution.h

BabyUniverseRemover.o: BabyUniverseRemover.cpp BabyUniverseRemover.h

BallSizeDistribution.o: BallSizeDistribution.cpp BallSizeDistribution.h

BitmapDrawer.o: BitmapDrawer.cpp BitmapDrawer.h

BoundaryMeasure.o: BoundaryMeasure.cpp BoundaryMeasure.h

CMinusTwoBuilder.o: CMinusTwoBuilder.cpp CMinusTwoBuilder.h

CirclePacking.o: CirclePacking.cpp CirclePacking.h

CirclePattern.o: CirclePattern.cpp CirclePattern.h

CohomologyBasis.o: CohomologyBasis.cpp CohomologyBasis.h

ConformalDistribution.o: ConformalDistribution.cpp ConformalDistribution.h

ConnectivityRestrictor.o: ConnectivityRestrictor.cpp ConnectivityRestrictor.h

Diffusion.o: Diffusion.cpp Diffusion.h

DiskCirclePacking.o: DiskCirclePacking.cpp DiskCirclePacking.h

DistanceDistribution.o: DistanceDistribution.cpp DistanceDistribution.h

DualCohomologyBasis.o: DualCohomologyBasis.cpp DualCohomologyBasis.h

DualDiffusion.o: DualDiffusion.cpp DualDiffusion.h

DualRandomWalk.o: DualRandomWalk.cpp DualRandomWalk.h

DualScalarField.o: DualScalarField.cpp DualScalarField.h

Edge.o: Edge.cpp Edge.h

Embedding.o: Embedding.cpp Embedding.h

HarmonicDiffusion.o: HarmonicDiffusion.cpp HarmonicDiffusion.h

HarmonicEmbedding.o: HarmonicEmbedding.cpp HarmonicEmbedding.h

Histogram.o: Histogram.cpp Histogram.h

HyperbolicStructure.o: HyperbolicStructure.cpp HyperbolicStructure.h

LaplacianDeterminant.o: LaplacianDeterminant.cpp LaplacianDeterminant.h

LaplacianMatrix.o: LaplacianMatrix.cpp LaplacianMatrix.h

LaplacianSpectrum.o: LaplacianSpectrum.cpp LaplacianSpectrum.h

LinearAlgebra.o: LinearAlgebra.cpp LinearAlgebra.h

LoopLength.o: LoopLength.cpp LoopLength.h

MetricGraphObservable.o: MetricGraphObservable.cpp MetricGraphObservable.h

ModuliObservable.o: ModuliObservable.cpp ModuliObservable.h

PeelingProcedure.o: PeelingProcedure.cpp PeelingProcedure.h

PottsModel.o: PottsModel.cpp PottsModel.h

PottsModelOnVertices.o: PottsModelOnVertices.cpp PottsModelOnVertices.h

ShortestLoop.o: ShortestLoop.cpp ShortestLoop.h

Simulation.o: Simulation.cpp Simulation.h

SpanningTree.o: SpanningTree.cpp SpanningTree.h

ThetaHistogram.o: ThetaHistogram.cpp ThetaHistogram.h

ThetaModel.o: ThetaModel.cpp ThetaModel.h

Triangle.o: Triangle.cpp Triangle.h

Triangulation.o: Triangulation.cpp Triangulation.h

TriangulationProperties.o: TriangulationProperties.cpp TriangulationProperties.h

Vertex.o: Vertex.cpp Vertex.h

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) ./hetws9/dt-ball