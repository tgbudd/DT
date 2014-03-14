#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "ConnectivityRestrictor.h"
#include "PeelingProcedure.h"
#include "CMinusTwoBuilder.h"
#include "LinearAlgebra.h"

#include <Eigen/Sparse>

class PeelingMeasurement : public Observable {
public:
	PeelingMeasurement( Triangulation * triangulation );
	void Measure();
	std::string OutputData() const;
	void SetVolumeWindow(int vmin, int vmax);
	void SetLengthWindow(int lmin, int lmax);
private:
	void MeasureBoundary();
	Triangle * SelectTriangle();

	Triangulation * triangulation_;
	PeelingProcedure peelingprocedure_;
	BabyUniverseDetector babyuniversedetector_;

	bool use_window_;
	int vmin_, vmax_;
	int lmin_, lmax_;

	std::list<const Edge *> boundary_;
	std::vector<Triangle *> triangles_;

	Eigen::MatrixXd invLaplacian_;
	bool DetermineInverseLaplacian();
	void UpdateHistogram();

	int max_delta_;
	double minlogmeasure_, maxlogmeasure_;
	int measure_bins_;
	std::vector<std::vector<int> > measure_histogram_;
	std::vector<int> measurements_;

	int max_volume_, d_volume_;
	std::vector<int> volume_histogram_;
	int max_length_, d_length_;
	std::vector<int> length_histogram_;
};

PeelingMeasurement::PeelingMeasurement(Triangulation * triangulation ) 
	: triangulation_(triangulation), 
	  peelingprocedure_(triangulation),
	  babyuniversedetector_(triangulation),
	  max_delta_(100),
	  minlogmeasure_(-15.0),
	  maxlogmeasure_(0.0),
	  measure_bins_(300),
	  max_volume_(triangulation->NumberOfTriangles()),
	  d_volume_(20),
	  max_length_(2000),
	  d_length_(1),
	  use_window_(false)
{
	measure_histogram_.resize(max_delta_,std::vector<int>(measure_bins_,0));
	measurements_.resize(max_delta_,0);
	volume_histogram_.resize(max_volume_/d_volume_,0);
	length_histogram_.resize(max_length_/d_length_,0);
}

void PeelingMeasurement::SetVolumeWindow(int vmin, int vmax)
{
	vmin_ = vmin;
	vmax_ = vmax;
	use_window_=true;
}

void PeelingMeasurement::SetLengthWindow(int lmin, int lmax)
{
	lmin_ = lmin;
	lmax_ = lmax;
	use_window_=true;
}

void PeelingMeasurement::Measure()
{
	Triangle * triangle = SelectTriangle();
		
	peelingprocedure_.PreparePeelingOnSphere( triangle );
	Edge * centerEdge = peelingprocedure_.final_vertex_->getParent()->getPrevious();
	do {
		if( peelingprocedure_.volume_within_frontier_ > (use_window_?triangulation_->NumberOfTriangles()-vmax_:triangulation_->NumberOfTriangles()/2))
		{
			if( use_window_ && peelingprocedure_.volume_within_frontier_ > triangulation_->NumberOfTriangles()-vmin_ )
			{
				break;
			}
			boundary_ = peelingprocedure_.frontier_;
			if( use_window_ && (static_cast<int>(boundary_.size()) < lmin_ || static_cast<int>(boundary_.size()) > lmax_) )
			{
				break;
			}
			boundary_.reverse();
			for(std::list<const Edge *>::iterator it = boundary_.begin(); it != boundary_.end(); it++)
			{
				(*it) = (*it)->getAdjacent();
			}
			MeasureBoundary();

			break;
		}
	} while( !peelingprocedure_.DoPeelingStepOnSphere() );
}

void PeelingMeasurement::MeasureBoundary()
{
	triangles_.clear();
	babyuniversedetector_.EnclosedTriangles(boundary_,triangles_,true);

	if( DetermineInverseLaplacian() )
	{
		UpdateHistogram();
		std::cout << "Measurement succesful.\n";
	}
}

bool PeelingMeasurement::DetermineInverseLaplacian()
{
	int size = triangles_.size();

	std::vector<Eigen::Triplet<double> > rules;
	rules.reserve(4*size);
	std::vector<int> trianglePosition(triangulation_->NumberOfTriangles(),-1);
	for(int i=0;i<size;i++)
	{
		trianglePosition[triangles_[i]->getId()] = i;
	}
	double minusonethird = -1/3.0;
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<3;j++)
		{
			int pos = trianglePosition[triangles_[i]->getEdge(j)->getAdjacent()->getParent()->getId()];
			if( pos >= 0 )
			{
				rules.push_back(Eigen::Triplet<double>(i,pos,minusonethird));
			}
		}
		rules.push_back(Eigen::Triplet<double>(i,i,1.0));
	}

	Eigen::SparseMatrix<double> matrix(size,size);
	matrix.setFromTriplets(rules.begin(),rules.end());

	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(matrix);

	if( chol.info() != Eigen::Success )
	{
		return false;
	}

	Eigen::SparseMatrix<double> boundaryMatrix(size,boundary_.size());
	boundaryMatrix.reserve(Eigen::VectorXi::Constant(boundary_.size(),1));
	int index = 0;
	for(std::list<const Edge *>::const_iterator it=boundary_.begin();it!=boundary_.end();it++,index++)
	{
		boundaryMatrix.insert(trianglePosition[(*it)->getParent()->getId()],index) = 1.0;
	}

	Eigen::SparseMatrix<double> invLapl = chol.solve(boundaryMatrix);
	invLaplacian_ = invLapl.toDense();
	if( chol.info() != Eigen::Success )
	{
		return false;
	}

	return true;
}

void PeelingMeasurement::UpdateHistogram()
{
	int bSize = boundary_.size();
	int tSize = triangles_.size();
	for(int n=1;n <= max_delta_ && n <= bSize;n++)
	{
		for(int t=0;t<tSize;t++)
		{
			double measure = 0.0;
			for(int i=0;i<n;i++)
			{
				measure += invLaplacian_(t,i);
			}
			for(int i=0;i<bSize;i++)
			{
				if( measure >= 0.0 )
				{
					int bin = static_cast<int>(std::floor((std::log(measure/3.0) - minlogmeasure_)/(maxlogmeasure_ - minlogmeasure_)*measure_bins_));
					if( bin >= 0 && bin < measure_bins_ )
					{
						measure_histogram_[n-1][bin]++;
					}
				}
				measure += invLaplacian_(t,(i+n)%bSize) - invLaplacian_(t,i);
			}
		}
		measurements_[n-1] += tSize;
	}

	if( tSize/d_volume_ < static_cast<int>(volume_histogram_.size()) )
	{
		volume_histogram_[tSize/d_volume_]++;
	}
	if( bSize/d_length_ < static_cast<int>(length_histogram_.size()) )
	{
		length_histogram_[bSize/d_length_]++;
	}
}

Triangle * PeelingMeasurement::SelectTriangle() 
{
	return triangulation_->getRandomTriangle();
}

std::string PeelingMeasurement::OutputData() const 
{ 
	std::ostringstream stream;
	stream << std::fixed << "peelingmeasurement -> {minlogmeasure -> ";
	stream << minlogmeasure_ << ", maxlogmeasure -> " << maxlogmeasure_;
	stream << ",measurebins -> " << measure_bins_ << ", measurements ->";
	PrintToStream(stream, measurements_.begin(), measurements_.end());
	stream << ", histogram -> ";
	PrintToStream2D(stream, measure_histogram_.begin(), measure_histogram_.end());
	stream << ", volumehistogram -> ";
	PrintToStream(stream, volume_histogram_.begin(), volume_histogram_.end());
	stream << ", lengthhistogram -> ";
	PrintToStream(stream, length_histogram_.begin(), length_histogram_.end());
	stream << "}";
	return stream.str();
}


int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	int volumewindow = param.Read<int>("volume window size");
	int lengthmin = param.Read<int>("min length");
	int lengthmax = param.Read<int>("max length");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, 0, secondsperoutput, output );

	PeelingMeasurement peelingmeasurement( &triangulation );
	peelingmeasurement.SetVolumeWindow(static_cast<int>(0.27*n)-volumewindow/2,static_cast<int>(0.27*n)+volumewindow/2);
	peelingmeasurement.SetLengthWindow(lengthmin,lengthmax);

	simulation.AddObservable( &peelingmeasurement, 1 );

	simulation.Run();
	return 0;
}