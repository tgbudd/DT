#include <algorithm>
#include <queue>
#include <fstream>

#include "ConformalDistribution.h"

ConformalDistribution::ConformalDistribution(Embedding * const embedding, const Triangulation * const triangulation) :
	embedding_(embedding),
	measurements_(0),
	minlograd_(-40.0),
	maxlograd_(0.0),
	bins_(1600),
	triangulation_(triangulation),
	measure_euclidean_ball_size_(false),
	max_ball_triangles_(200),
	ball_minlograd_(-25.0),
	ball_maxlograd_(0.0),
	ball_bins_(125)
{
	histogram_.resize(bins_,0);
	for(int i=1,endi=max_ball_triangles_;i<=endi;i++)
	{
		record_ball_sizes_.push_back(i);
	}
}

ConformalDistribution::~ConformalDistribution()
{
}

void ConformalDistribution::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		MeasureRadiusDistribution();

		if( measure_euclidean_ball_size_ )
		{
			MeasureEuclideanBallSize();
		}

		if( save_embedding_ )
		{
			embedding_->SaveEmbedding(embedding_file_);
		}
	}
}

void ConformalDistribution::MeasureRadiusDistribution()
{
	std::vector<double> radius;
		
	if( embedding_->GetRadii(radius) )
	{
		for(int i=0,endi=radius.size();i<endi;i++)
		{
			if( radius[i] > 0.0 )
			{
				double logradius = std::log(radius[i]);
				if( logradius >= minlograd_ && logradius < maxlograd_ )
				{
					histogram_[static_cast<int>((logradius - minlograd_)/(maxlograd_-minlograd_)*bins_)]++;
				}
			}
		}
		measurements_++;
	}
}

void ConformalDistribution::setMeasureEuclideanBallSize(bool measure)
{
	measure_euclidean_ball_size_ = measure;
}

void ConformalDistribution::MeasureEuclideanBallSize()
{
	BOOST_ASSERT( triangulation_ != NULL );
	if( ball_histogram_.empty() )
	{
		ball_histogram_.resize(record_ball_sizes_.size(),std::vector<int>(ball_bins_,0));
		point_histogram_.resize(ball_bins_,std::vector<int>(record_ball_sizes_.size(),0));
	}
	triangle_visited_.ResetAndResize(triangulation_->NumberOfTriangles());
	vertex_visited_.ResetAndResize(triangulation_->NumberOfVertices());
	vertex_position_.resize(triangulation_->NumberOfVertices());
	modulus_ = embedding_->CalculateModuli();
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;i++)
	{
		MeasureEuclideanBallSize(triangulation_->getTriangle(i));
	}
}

void ConformalDistribution::setMaxBallTriangles(int n)
{
	max_ball_triangles_ = n;
	record_ball_sizes_.clear();
	for(int i=1,endi=max_ball_triangles_;i<=endi;i++)
	{
		record_ball_sizes_.push_back(i);
	}
}

void ConformalDistribution::setRecordBallSizesExponentially(int maxsize, double factor)
{
	record_ball_sizes_.clear();
	for(double size=1.1;size<maxsize;size *=factor)
	{
		if( record_ball_sizes_.empty() || static_cast<int>(size) > record_ball_sizes_.back() )
		{
			record_ball_sizes_.push_back(static_cast<int>(size));
		}
	}
	if( (maxsize - record_ball_sizes_.back()) > 0.3*factor*record_ball_sizes_.back() )
	{
		record_ball_sizes_.push_back( maxsize );
	}
	max_ball_triangles_ = record_ball_sizes_.back();
}

void ConformalDistribution::MeasureEuclideanBallSize(Triangle * triangle)
{
	std::priority_queue<double> smallestSqDistances;
	std::priority_queue<vertexNode> queue;
	vertex_visited_.Reset();
	triangle_visited_.Reset();
	for(int i=0;i<3;i++)
	{
		Edge * edge = triangle->getEdge(i);
		vertexNode v;
		v.vertex = edge->getOpposite();
		v.position = AddScaledVectors2D(1/3.0,embedding_->getForm(edge->getNext()),-1/3.0,embedding_->getForm(edge->getPrevious()));
		v.square_distance = NormSquaredTransformedByModulus(v.position,modulus_);
		queue.push(v);
		vertex_visited_.Set(v.vertex->getId());
	}
	triangle_visited_.Set(triangle->getId());

	while( !queue.empty() && (smallestSqDistances.size() < max_ball_triangles_ || smallestSqDistances.top() > queue.top().square_distance) )
	{
		vertexNode v = queue.top();
		queue.pop();

		Edge * edge = v.vertex->getParent()->getPrevious();
		do
		{
			Triangle * t = edge->getParent();
			if( !triangle_visited_.isSet(t->getId()) )
			{
				Vector2D y = AddVectors2D(v.position, AddScaledVectors2D( 1/3.0, embedding_->getForm(edge), -1/3.0, embedding_->getForm(edge->getPrevious())));
				double d = NormSquaredTransformedByModulus(y,modulus_);
				if( smallestSqDistances.size() < max_ball_triangles_ )
				{
					smallestSqDistances.push(d);
				} else if( smallestSqDistances.top() > d )
				{
					smallestSqDistances.pop();
					smallestSqDistances.push(d);
				}
				triangle_visited_.Set(t->getId());
			}

			vertexNode v2(v);
			v2.vertex = edge->getPrevious()->getOpposite();
			if( !vertex_visited_.isSet(v2.vertex->getId()) )
			{
				vertex_visited_.Set(v2.vertex->getId());
				v2.position = AddVectors2D(v2.position,embedding_->getForm(edge));
				v2.square_distance = NormSquaredTransformedByModulus(v2.position,modulus_);
				queue.push(v2);
			}

			edge = edge->getAdjacent()->getNext();
		} while( edge != v.vertex->getParent()->getPrevious() );
	}

	BOOST_ASSERT(smallestSqDistances.size() == max_ball_triangles_ || static_cast<int>(smallestSqDistances.size())+1 == triangulation_->NumberOfTriangles());

	std::vector<int> record(max_ball_triangles_,-1);
	for(int i=0,endi=record_ball_sizes_.size();i<endi;i++)
	{
		record[record_ball_sizes_[i]-1] = i;
	}

	int index = max_ball_triangles_;
	int currentbin = ball_bins_-1;
	double logepsilon = ball_maxlograd_;
	while( !smallestSqDistances.empty() )
	{
		index--;
		double logradius = 0.5*std::log(smallestSqDistances.top());
		smallestSqDistances.pop();
		if( record[index] >=0 && logradius >= ball_minlograd_ && logradius < ball_maxlograd_ )
		{
			ball_histogram_[record[index]][static_cast<int>((logradius - ball_minlograd_)/(ball_maxlograd_-ball_minlograd_)*ball_bins_)]++;
		}

		while( currentbin >=0 && logepsilon > logradius )
		{
			if( record[index] >= 0 )
			{
				point_histogram_[currentbin][record[index]]++;
			}
			currentbin--;
			logepsilon -= (ball_maxlograd_-ball_minlograd_)/ball_bins_;
		}
	}
}

void ConformalDistribution::setMethodString(const std::string & method)
{
	method_ = method;
}

std::string ConformalDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "conformaldistribution -> {method -> \"" << method_ << "\", minlograd -> " << minlograd_;
	stream << ", maxlograd -> " << maxlograd_;
	stream << ", bins -> " << bins_;
	stream << ", measurements -> " <<  measurements_;
	stream << ", histogram -> ";
	PrintToStream(stream,histogram_.begin(),histogram_.end());
	if( measure_euclidean_ball_size_ )
	{
		stream << ", ballminlograd -> " << ball_minlograd_;
		stream << ", ballmaxlograd -> " << ball_maxlograd_;
		stream << ", ballbins -> " << ball_bins_;
		stream << ", ballsizes -> ";
		PrintToStream(stream,record_ball_sizes_.begin(),record_ball_sizes_.end());
		stream << ", ballhistogram -> ";
		PrintToStream2D(stream,ball_histogram_.begin(),ball_histogram_.end());
		stream << ", pointhistogram -> ";
		PrintToStream2D(stream,point_histogram_.begin(),point_histogram_.end());
	}
	stream << "}";
	return stream.str();
}

void ConformalDistribution::setSaveEmbedding(bool save, std::string file)
{
	save_embedding_ = save;
	embedding_file_ = file;
}

