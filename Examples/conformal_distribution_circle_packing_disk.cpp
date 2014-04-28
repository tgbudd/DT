#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <queue>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "LinearAlgebra.h"
#include "DiskCirclePacking.h"

class DiskConformalDistribution : public Observable
{
public:
	DiskConformalDistribution(Triangulation * triangulation, CMinusTwoBuilder * builder);
	void Measure();
	void setRecordBallSizesExponentially(int maxsize, double factor);
	std::string OutputData() const;
private:
	void SaveDisk();
	const Edge * chooseEdge();
	void RecordDistances(int v);
	Triangulation * triangulation_;
	CMinusTwoBuilder * builder_;
	DiskCirclePacking diskcirclepacking_;


	std::list<const Edge*> boundary_;
	std::vector<std::pair<Vertex*,Vector2D> > coor_;
	std::vector<double> tmp_dist_;

	std::vector<int> record_;
	unsigned int max_ball_triangles_;
	std::vector<int> record_ball_sizes_;
	std::vector<std::vector<int> > ball_histogram_;
	double ball_minlograd_, ball_maxlograd_;
	int ball_bins_;
	std::vector<std::vector<int> > point_histogram_;
	int measurements_;
};

DiskConformalDistribution::DiskConformalDistribution(Triangulation * triangulation, CMinusTwoBuilder * builder) :
	triangulation_(triangulation),
	builder_(builder),
	diskcirclepacking_(triangulation),
	measurements_(0),
	max_ball_triangles_(40),
	ball_minlograd_(-20.0),
	ball_maxlograd_(8.0),
	ball_bins_(140)
{
	for(int i=1,endi=max_ball_triangles_;i<=endi;i++)
	{
		record_ball_sizes_.push_back(i);
	}
	record_.resize(max_ball_triangles_,-1);
	for(int i=0,endi=record_ball_sizes_.size();i<endi;i++)
	{
		record_[record_ball_sizes_[i]-1] = i;
	}
}

void DiskConformalDistribution::setRecordBallSizesExponentially(int maxsize, double factor)
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
	record_.resize(max_ball_triangles_,-1);
	for(int i=0,endi=record_ball_sizes_.size();i<endi;i++)
	{
		record_[record_ball_sizes_[i]-1] = i;
	}
}

void DiskConformalDistribution::Measure()
{
	if( ball_histogram_.empty() )
	{
		ball_histogram_.resize(record_ball_sizes_.size(),std::vector<int>(ball_bins_,0));
		point_histogram_.resize(ball_bins_,std::vector<int>(record_ball_sizes_.size(),0));
	}

	builder_->getDiskBoundary(boundary_);
	//const Edge * edge = chooseEdge();
	if( diskcirclepacking_.FindEmbedding(boundary_))
	{
		std::cout << "Embedding found.\n";
		coor_.clear();
		diskcirclepacking_.getHyperbolicCoordinates(coor_);
		
		tmp_dist_.resize(coor_.size());
		for(int i=0,endi=coor_.size();i<endi;i++)
		{
			if( Norm2D(coor_[i].second) < 0.999 )
			{
				RecordDistances(i);
				measurements_++;
			}
		}

	}
}

void DiskConformalDistribution::RecordDistances(int v)
{
	for(int i=0,endi=coor_.size();i<endi;i++)
	{
		tmp_dist_[i] = PoincareDistance(coor_[i].second,coor_[v].second);
	}
	std::sort(tmp_dist_.begin(),tmp_dist_.end());

	int currentbin = ball_bins_-1;
	double logepsilon = ball_maxlograd_;
	for(int i=max_ball_triangles_-1;i>=0;i--)
	{
		if( i >= static_cast<int>(tmp_dist_.size()) )
			continue;

		double logradius = std::log(tmp_dist_[i+1]);
		if( record_[i] >=0 && logradius >= ball_minlograd_ && logradius < ball_maxlograd_ )
		{
			ball_histogram_[record_[i]][static_cast<int>((logradius - ball_minlograd_)/(ball_maxlograd_-ball_minlograd_)*ball_bins_)]++;
		}

		while( currentbin >=0 && logepsilon > logradius )
		{
			if( record_[i] >= 0 )
			{
				point_histogram_[currentbin][record_[i]]++;
			}
			currentbin--;
			logepsilon -= (ball_maxlograd_-ball_minlograd_)/ball_bins_;
		}
	}
}

const Edge * DiskConformalDistribution::chooseEdge()
{
	std::vector<int> distance(triangulation_->NumberOfVertices(),-1);
	std::queue<Vertex *> q;
	for(std::list<const Edge *>::iterator it = boundary_.begin();it != boundary_.end();it++)
	{
		Vertex * v = (*it)->getNext()->getOpposite();
		q.push(v);
		distance[v->getId()] = 0;
	}
	int maxdist = 0;
	Vertex * furthestV;
	Vertex * centerVertex = boundary_.front()->getAdjacent()->getOpposite();
	while( !q.empty() )
	{
		Vertex * v = q.front();
		q.pop();

		Edge * edge = v->getParent()->getPrevious();
		do {
			Vertex * v2 = edge->getPrevious()->getOpposite();
			if( distance[v2->getId()] == -1 )
			{
				q.push(v2);
				distance[v2->getId()] = distance[v->getId()]+1;
				if( distance[v2->getId()] > maxdist && v2 != centerVertex )
				{
					furthestV = v2;
					maxdist = distance[v2->getId()];
				}
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != v->getParent()->getPrevious() );
	}

	return furthestV->getParent()->getPrevious();
}

/*void DiskConformalDistribution::SaveDisk()
{
	std::vector<int> vertexIds(triangulation_->NumberOfVertices(),-1);
	for(int i=0,endi=circles_.size();i<endi;i++)
	{
		vertexIds[circles_[i].first->getId()] = i;
	}
	std::vector<Triangle *> triangles = diskcirclepacking_.getDiskTriangles();
	
	std::ofstream file("disk.txt");
	file << std::fixed << std::setprecision(10) << "{triangles -> {{";
	for(int i=0,endi=triangles.size();i<endi;i++)
	{
		file << (i>0?"},{":"");
		for(int j=0;j<3;j++)
		{
			file << (j>0?",":"") << vertexIds[triangles[i]->getEdge(j)->getOpposite()->getId()];
		}
	}
	file << "}}, boundary -> {";
	for(std::list<const Edge *>::const_iterator it = boundary_.begin();it!= boundary_.end();it++)
	{
		file << (it == boundary_.begin()?"":",") << vertexIds[(*it)->getNext()->getOpposite()->getId()];
	}
	file << "}, circles -> {";
	for(int i=0,endi=circles_.size();i<endi;i++)
	{
		file << (i>0?",":"") << "{{" << circles_[i].second.first[0] << "," << circles_[i].second.first[1] << "}," << circles_[i].second.second << "}";
	}
	file << "}}\n";
}*/


std::string DiskConformalDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "diskconformaldistribution -> {";
	stream << " measurements -> " <<  measurements_;
	stream << ", ballminlograd -> " << ball_minlograd_;
	stream << ", ballmaxlograd -> " << ball_maxlograd_;
	stream << ", ballbins -> " << ball_bins_;
	stream << ", ballsizes -> ";
	PrintToStream(stream,record_ball_sizes_.begin(),record_ball_sizes_.end());
	stream << ", ballhistogram -> ";
	PrintToStream2D(stream,ball_histogram_.begin(),ball_histogram_.end());
	stream << ", pointhistogram -> ";
	PrintToStream2D(stream,point_histogram_.begin(),point_histogram_.end());
	stream << "}";
	return stream.str();
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	int b = param.Read<int>("boundary length");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n,b);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	DiskConformalDistribution conf(&triangulation,&builder);

	Simulation simulation( &triangulation, 0, secondsperoutput, output );
	simulation.AddObservable( &conf, 1 );

	simulation.Run();

	return 0;
}