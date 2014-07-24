#include "BallSizeDistribution.h"

BallSizeDistribution::BallSizeDistribution(Embedding * const embedding, const Triangulation * const triangulation,
		double absmaxlogradius, int logradiusbins, int samples) :
	embedding_(embedding),
	triangulation_(triangulation),
	log_dist_(-absmaxlogradius,0.0,logradiusbins,false),
	samples_(samples)
{
}

BallSizeDistribution::~BallSizeDistribution()
{
}

void BallSizeDistribution::AddReferenceRadius(double epsilon)
{
	ref_epsilons_.push_back(epsilon);
}
void BallSizeDistribution::AddVolumeFraction(double delta)
{
	deltas_.push_back(delta);
}
void BallSizeDistribution::InitializeHistograms()
{
	std::sort(ref_epsilons_.begin(),ref_epsilons_.end());
	std::sort(deltas_.begin(),deltas_.end());
	ball_size_ = std::vector<std::vector<Histogram<int> > >(ref_epsilons_.size(),
		std::vector<Histogram<int> >(deltas_.size(), 
		Histogram<int>(0,log_dist_.GetBins()+1,log_dist_.GetBins()+1,false)));
	unrooted_ball_size_ = ball_size_;
}

void BallSizeDistribution::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		modulus_ = embedding_->CalculateModuli();

		points_.resize( triangulation_->NumberOfTriangles() );
		for(int i=0,endi=points_.size();i<endi;++i)
		{
			points_[i] = embedding_->GetCentroid(triangulation_->getTriangle(i));
		}
		for(int i=0;i<samples_;i++)
		{
			Vector2D center = points_[triangulation_->RandomInteger(0,points_.size()-1)];
			ComputeDistanceHistogram(points_,center,ball_size_);
		}
		for(int i=0;i<samples_;i++)
		{
			Vector2D center = {triangulation_->RandomReal(0.0,1.0),triangulation_->RandomReal(0.0,1.0)};
			ComputeDistanceHistogram(points_,center,unrooted_ball_size_);
		}
	}
}

void BallSizeDistribution::ComputeDistanceHistogram(const std::vector<Vector2D> & x, const Vector2D & c, std::vector<std::vector<Histogram<int> > > & histograms)
{
	log_dist_.Reset();
	for(int i=0,endi=x.size();i<endi;++i)
	{
		double distsq = std::max(1.0e-20,DistanceSquared(SubtractVectors2D(x[i],c)));
		double logdist = 0.5 * std::log(distsq);
		log_dist_.Insert(logdist);
	}

	for(int i=0,endi=ref_epsilons_.size();i<endi;++i)
	{
		double logepsilon = std::log(ref_epsilons_[i]);
		int referencevolume = log_dist_.Cumulative(logepsilon);

		std::vector<int> volumes;
		volumes.reserve(deltas_.size());
		BOOST_FOREACH(double delta,deltas_)
		{
			volumes.push_back( static_cast<int>(delta * referencevolume) );
		}
		std::vector<int> logeps_bin = log_dist_.CumulToBin(volumes);

		for(int j=0,endj=logeps_bin.size();j<endj;++j)
		{
			histograms[i][j].Insert(logeps_bin[j]);
		}
	}
}

double BallSizeDistribution::DistanceSquared(Vector2D x) const
{
	// this is the exact distance when x close to zero or when modulus_.first is small,
	// but otherwise may overestimate the distance for large x 
	x[0] = properfmod(x[0] + 0.5,1.0)-0.5;
	x[1] = properfmod(x[1] + 0.5,1.0)-0.5;
	return NormSquaredTransformedByModulus(x,modulus_);
}

std::string BallSizeDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << "ballsizedistribution -> {";
	stream << "samples -> " << samples_;
	stream << ", absmaxlogradius -> " << -log_dist_.GetMin();
	stream << ", logradiusbins -> " << log_dist_.GetBins();
	stream << ", refepsilons -> ";
	PrintToStream(stream,ref_epsilons_.begin(),ref_epsilons_.end());
	stream << ", deltas -> {";
	for(int i=0,endi=deltas_.size();i<endi;i++)
	{
		stream << (i>0?",":"") << "Exp[" << std::log(deltas_[i]) << "]";
	}
	stream << "}";
	stream << ", ballsize -> {";
	for(int i=0,endi=ref_epsilons_.size();i<endi;++i)
	{
		stream << (i>0?",":"") << "{";
		for(int j=0,endj=deltas_.size();j<endj;++j)
		{
			stream << (j>0?",":"");
			ball_size_[i][j].PrintTo(stream);
		}
		stream << "}";
	}
	stream << "}";
	stream << ", unrootedballsize -> {";
	for(int i=0,endi=ref_epsilons_.size();i<endi;++i)
	{
		stream << (i>0?",":"") << "{";
		for(int j=0,endj=deltas_.size();j<endj;++j)
		{
			stream << (j>0?",":"");
			unrooted_ball_size_[i][j].PrintTo(stream);
		}
		stream << "}";
	}
	stream << "}}";
	return stream.str();
}