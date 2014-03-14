#include "BoundaryMeasure.h"

BoundaryMeasure::BoundaryMeasure(Triangulation * triangulation ) 
	: triangulation_(triangulation), 
	  babyuniversedetector_(triangulation),
	  max_delta_(20),
	  minlogmeasure_(-15.0),
	  maxlogmeasure_(0.0),
	  measure_bins_(300),
	  samples_(100)
{
	measure_histogram_.resize(max_delta_,std::vector<int>(measure_bins_,0));
	measurements_.resize(max_delta_,0);
}

void BoundaryMeasure::Measure()
{
	DetermineBoundary();

	triangles_.clear();
	babyuniversedetector_.EnclosedTriangles(boundary_,triangles_,true);

	if( DetermineInverseLaplacian() )
	{
		UpdateHistogram();
		std::cout << "Measurement succesful.\n";
	}
}

bool BoundaryMeasure::DetermineInverseLaplacian()
{
	int size = triangles_.size();

	std::vector<Eigen::Triplet<double> > rules;
	rules.reserve(4*size);
	std::vector<int> trianglePosition(triangulation_->NumberOfTriangles(),-1);
	boost::array<bool,3> allFalse = {false,false,false};
	std::vector<boost::array<bool,3> > inBoundary(size,allFalse);
	for(int i=0;i<size;i++)
	{
		trianglePosition[triangles_[i]->getId()] = i;
	}
	for(std::list<const Edge *>::const_iterator it=boundary_.begin();it!=boundary_.end();it++)
	{
		inBoundary[trianglePosition[(*it)->getParent()->getId()]][(*it)->getId()] = true;
	}
	double minusonethird = -1/3.0;
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<3;j++)
		{
			if( !inBoundary[i][j] )
			{ 
				int pos = trianglePosition[triangles_[i]->getEdge(j)->getAdjacent()->getParent()->getId()];
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

void BoundaryMeasure::UpdateHistogram()
{
	int bSize = boundary_.size();
	int tSize = triangles_.size();

	std::vector<int> sample;
	int numSamples = std::min(2*bSize,tSize);
	triangulation_->RandomSample(0,tSize-1,numSamples,sample);
	for(int n=1;n <= max_delta_ && n <= bSize;n++)
	{
		for(int t=0;t<numSamples;t++)
		{
			double measure = 0.0;
			for(int i=0;i<n;i++)
			{
				measure += invLaplacian_(sample[t],i);
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
				measure += invLaplacian_(sample[t],(i+n)%bSize) - invLaplacian_(sample[t],i);
			}
		}
		measurements_[n-1] += numSamples;
	}
}

std::string BoundaryMeasure::OutputData() const 
{ 
	std::ostringstream stream;
	stream << std::fixed << "peelingmeasurement -> {minlogmeasure -> ";
	stream << minlogmeasure_ << ", maxlogmeasure -> " << maxlogmeasure_;
	stream << ",measurebins -> " << measure_bins_ << ", measurements ->";
	PrintToStream(stream, measurements_.begin(), measurements_.end());
	stream << ", histogram -> ";
	PrintToStream2D(stream, measure_histogram_.begin(), measure_histogram_.end());
	stream << "}";
	return stream.str();
}