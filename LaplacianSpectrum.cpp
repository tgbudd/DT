#include "LaplacianSpectrum.h"
#include "LinearAlgebra.h"
#include "Embedding.h"

LaplacianSpectrum::LaplacianSpectrum(const Triangulation * const triangulation) : triangulation_(triangulation),
	ev_bins_(1000), ev_max_(10.0), ev_min_(0.0), measurements_(0)
{
	eigenvalue_histogram_.resize(ev_bins_,0);
}


LaplacianSpectrum::~LaplacianSpectrum(void)
{
}

void LaplacianSpectrum::Measure()
{
	LaplacianMatrix laplacian(triangulation_);
	laplacian.ComputeEigenvalues();

	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		double eigenvalue = laplacian.GetEigenvalue(i);
		int val = (int)(ev_bins_ * (eigenvalue - ev_min_)/(ev_max_-ev_min_));
		if( val >=0 && val < ev_bins_ )
		{
			eigenvalue_histogram_[val]++;
		}
	}
	measurements_++;
}


std::string LaplacianSpectrum::OutputData() const
{
	std::ostringstream stream;
	stream << "laplacianspectrum -> {measurements -> " << measurements_ << ", eigenvaluehistogram -> ";
	PrintToStream(stream,eigenvalue_histogram_.begin(),eigenvalue_histogram_.end());
	stream << "}";
	return stream.str();
}