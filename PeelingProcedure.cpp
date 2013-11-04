#include "boost/utility.hpp"
#include <set>

#include "PeelingProcedure.h"
#include "TriangulationProperties.h"
#include "BabyUniverseDetector.h"

PeelingProcedure::PeelingProcedure(const Triangulation * const triangulation)
	: triangulation_(triangulation), 
	  babyuniversedetector_(triangulation),
	  walk_time_(50000),
	  random_walk_measurements_(0),
	  random_walk_measurements_mother_(0),
	  max_trajectories_(1),
	  baby_walk_(true)
{
	in_mother_universe_.resize(triangulation_->NumberOfTriangles(),false);
	in_frontier_.resize(triangulation_->NumberOfVertices(),false);
	vertex_in_mother_universe_.resize(triangulation_->NumberOfVertices(),false);
	distance_square_.resize(walk_time_,0);
	distance_square_mother_.resize(walk_time_,0);
	distance_square_mother_effective_.resize(walk_time_,0);
	waiting_times_.resize(1000,0);
}


PeelingProcedure::~PeelingProcedure(void)
{
}


void PeelingProcedure::Measure(void)
{
	DoPeeling();
	//RandomWalk( false );
	//RandomWalk( true );
}


std::string PeelingProcedure::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "peeling -> {trajectories -> ";
	PrintToStream2D(stream,frontier_size_trajectory_.begin(),frontier_size_trajectory_.end());
	stream << ", volume -> ";
	PrintToStream2D(stream,volume_trajectory_.begin(),volume_trajectory_.end());
	stream << ", distance -> ";
	PrintToStream2D(stream,distance_trajectory_.begin(),distance_trajectory_.end());
	stream << ", waitingtimes -> ";
	PrintToStream(stream,waiting_times_.begin(),waiting_times_.end());
	stream << "/" << random_walk_measurements_ << ", distancesquare -> ";
	PrintToStream(stream,distance_square_.begin(),distance_square_.end());
	stream << "/" << random_walk_measurements_<< ", distancesquaremother -> ";
	PrintToStream(stream,distance_square_mother_.begin(),distance_square_mother_.end());
	stream << "/" << random_walk_measurements_mother_ << ", distancesquaremothereffective -> ";
	PrintToStream(stream,distance_square_mother_effective_.begin(),distance_square_mother_effective_.end());
	stream << "/" << random_walk_measurements_;
	stream << ", babyuniverses -> ";
	PrintToStream2D(stream,baby_universes_.begin(),baby_universes_.end());
	stream << ", babytime -> {";
	for(int i=0;i<static_cast<int>(baby_return_time_.size());i++)
	{
		if( i > 0 )
		{
			stream << ", ";
		}
		PrintToStream2D(stream,baby_return_time_[i].begin(),baby_return_time_[i].end());
	}
	stream << "}}";
	return stream.str();
}

void PeelingProcedure::DoPeeling()
{
	Triangle * startTriangle;
	do {
		 startTriangle = triangulation_->getRandomTriangle();
	} while( startTriangle->getEdge(0)->getOpposite() == startTriangle->getEdge(1)->getOpposite()
		|| startTriangle->getEdge(1)->getOpposite() == startTriangle->getEdge(2)->getOpposite()
		|| startTriangle->getEdge(0)->getOpposite() == startTriangle->getEdge(2)->getOpposite() );

	start_triangle_ = startTriangle;
	Edge * oppositeEdge = startTriangle->getEdge(triangulation_->RandomInteger(0,2));
	Vertex * startVertex = oppositeEdge->getOpposite();
	start_vertex_ = startVertex;
	frontier_.clear();
	frontier_.push_back(oppositeEdge);
	frontier_.push_back(oppositeEdge->getNext());
	frontier_.push_back(oppositeEdge->getPrevious());
	frontier_size_=3;
	volume_within_frontier_=1;

	std::fill(in_mother_universe_.begin(),in_mother_universe_.end(),false);
	in_mother_universe_[startTriangle->getId()] = true;
	std::fill(in_frontier_.begin(),in_frontier_.end(),false);
	in_frontier_[oppositeEdge->getOpposite()->getId()] = true;
	in_frontier_[oppositeEdge->getNext()->getOpposite()->getId()] = true;
	in_frontier_[oppositeEdge->getPrevious()->getOpposite()->getId()] = true;
	std::fill(vertex_in_mother_universe_.begin(),vertex_in_mother_universe_.end(),false);
	vertex_in_mother_universe_[oppositeEdge->getOpposite()->getId()] = true;
	vertex_in_mother_universe_[oppositeEdge->getNext()->getOpposite()->getId()] = true;
	vertex_in_mother_universe_[oppositeEdge->getPrevious()->getOpposite()->getId()] = true;

	properties::VertexDistanceList(triangulation_,startVertex,distance_);

	Vertex * finalVertex = ChooseFinalVertex();
	properties::VertexDistanceList(triangulation_,finalVertex,distance_to_final_);

	int steps = 0;
	bool recordTrajectory = false;
	if( static_cast<int>(frontier_size_trajectory_.size()) < max_trajectories_ )
	{
		frontier_size_trajectory_.push_back(std::vector<int>());
		volume_trajectory_.push_back(std::vector<int>());
		distance_trajectory_.push_back(std::vector<int>());
		recordTrajectory = true;
	}
	while( !in_frontier_[finalVertex->getId()] )
	{
		if( PeelingStep() )
		{
			// Fronteir has a figure-8 form and therefore a baby universe has been encountered.
			ProcessBabyUniverse();
		}
		volume_within_frontier_++;
		steps++;
		if( recordTrajectory )
		{
			frontier_size_trajectory_.back().push_back(frontier_size_);
			volume_trajectory_.back().push_back(volume_within_frontier_);
			distance_trajectory_.back().push_back(distance_[(*frontier_.begin())->getNext()->getOpposite()->getId()]);
		}
	}
}

bool PeelingProcedure::PeelingStep()
{
	Edge * adjEdge = frontier_.back()->getAdjacent();
	Triangle * newTriangle = adjEdge->getParent();

	BOOST_ASSERT( !in_mother_universe_[newTriangle->getId()] );
	in_mother_universe_[newTriangle->getId()] = true;
	
	if( frontier_.size() > 1 && adjEdge->getNext() == (*(boost::prior(boost::prior(frontier_.end()))))->getAdjacent() )
	{
		frontier_.pop_back();
		frontier_.pop_back();
		frontier_.push_back(adjEdge->getPrevious());
		in_frontier_[adjEdge->getPrevious()->getOpposite()->getId()] = false;
		frontier_size_--;
	} else if( frontier_.front()->getAdjacent() == adjEdge->getPrevious() )
	{
		frontier_.pop_back();
		frontier_.pop_front();
		frontier_.push_back(adjEdge->getNext());
		in_frontier_[adjEdge->getNext()->getOpposite()->getId()] = false;
		frontier_size_--;
	} else
	{
		frontier_.pop_back();
		frontier_.push_back(adjEdge->getNext());
		frontier_.push_front(adjEdge->getPrevious());
		Vertex * newVertex = adjEdge->getOpposite();
		frontier_size_++;
		if( !in_frontier_[newVertex->getId()] )
		{
			in_frontier_[newVertex->getId()] = true;
			vertex_in_mother_universe_[newVertex->getId()] = true;
		} else
		{
			return true;
		}
	}

	return false;
}

void PeelingProcedure::ProcessBabyUniverse()
{
	Vertex * saddleVertex = frontier_.front()->getNext()->getOpposite();
	std::list<const Edge*>::iterator loopBegin;
	int firstLoopSize=0;
	for(std::list<const Edge*>::iterator it=boost::next(frontier_.begin());it!=frontier_.end();it++)
	{
		firstLoopSize++;
		if( (*it)->getNext()->getOpposite() == saddleVertex )
		{
			loopBegin = it;
			break;
		}
	}

	int minDistanceToFinalVertex = 1000000;
	for(std::list<const Edge*>::iterator it=frontier_.begin();it!=loopBegin;it++)
	{
		minDistanceToFinalVertex = std::min(distance_to_final_[(*it)->getNext()->getOpposite()->getId()],minDistanceToFinalVertex);
	}
	bool firstLoopContainsFinalVertex = false; 
	for(std::list<const Edge*>::iterator it=frontier_.begin();it!=loopBegin;it++)
	{
		Vertex * vertex = (*it)->getNext()->getOpposite();
		if( distance_to_final_[vertex->getId()] == minDistanceToFinalVertex )
		{
			Edge * endEdge;
			if( it == frontier_.begin() )
			{
				endEdge = (*boost::prior(loopBegin))->getAdjacent();
			} else
			{
				endEdge = (*boost::prior(it))->getAdjacent();
			}
			Edge * edge = (*it)->getAdjacent()->getNext();
			while( edge != endEdge )
			{
				if( distance_to_final_[edge->getPrevious()->getOpposite()->getId()] == minDistanceToFinalVertex - 1 )
				{
					firstLoopContainsFinalVertex  = true;
					break;
				}
				edge = edge->getAdjacent()->getNext();
			}
			if( firstLoopContainsFinalVertex )
			{
				break;
			}
		}
	}
	if( firstLoopContainsFinalVertex )
	{
		int volume = babyuniversedetector_.VolumeEnclosed(loopBegin,frontier_.end(),false);
		DoMeasurementOnBabyUniverse(loopBegin,frontier_.end(),volume);
		for(std::list<const Edge*>::iterator it=boost::next(loopBegin);it!=frontier_.end();it++)
		{
			in_frontier_[(*it)->getNext()->getOpposite()->getId()] = false;
		}
		volume_within_frontier_ += volume;	
		frontier_.erase(loopBegin,frontier_.end());
		frontier_.push_back(frontier_.front());
		frontier_.pop_front();
		frontier_size_ = firstLoopSize;

	} else
	{
		int volume = babyuniversedetector_.VolumeEnclosed(frontier_.begin(),loopBegin,false);
		DoMeasurementOnBabyUniverse(frontier_.begin(),loopBegin,volume);
		for(std::list<const Edge*>::iterator it=boost::next(frontier_.begin());it!=loopBegin;it++)
		{
			in_frontier_[(*it)->getNext()->getOpposite()->getId()] = false;
		}
		volume_within_frontier_ += volume;	
		frontier_.erase(frontier_.begin(),loopBegin);
		frontier_size_ -= firstLoopSize;
	}
}

Vertex * PeelingProcedure::ChooseFinalVertex()
{
	// Choose a random vertex of maximal distance.
	int maxdist = 0;
	int VerticesAtMaxDist = 0;
	Vertex * vertex;
	for(int i=0,endi=triangulation_->NumberOfVertices();i<endi;i++)
	{
		if( distance_[i] > maxdist )
		{
			maxdist = distance_[i];
			VerticesAtMaxDist = 1;
			vertex = triangulation_->getVertex(i);
		} else if( distance_[i] == maxdist )
		{
			VerticesAtMaxDist++;
			if( 0 == triangulation_->RandomInteger(0,VerticesAtMaxDist-1) )
			{
				vertex = triangulation_->getVertex(i);
			}
		}
	}
	return vertex;
}

bool PeelingProcedure::FrontierIsSimpleClosedCurve()
{
	std::set<Vertex*> vertices;

	for(std::list<const Edge*>::iterator it=frontier_.begin();it!=frontier_.end();it++)
	{
		std::list<const Edge*>::iterator it2 = boost::next(it);
		if( it2==frontier_.end() )
		{
			it2 = frontier_.begin();
		}
		if( (*it)->getPrevious()->getOpposite() != (*it2)->getNext()->getOpposite() )
		{
			return false;
		}
		if( !vertices.insert((*it)->getPrevious()->getOpposite()).second )
		{
			return false;
		}
		if( !in_frontier_[(*it)->getPrevious()->getOpposite()->getId()] )
		{
			return false;
		}
		if( in_mother_universe_[(*it)->getAdjacent()->getParent()->getId()] )
		{
			return false;
		}
	}
	return true;
}

bool PeelingProcedure::FinalVertexBeyondFrontier()
{
	return FinalVertexInLoop(frontier_);
}

bool PeelingProcedure::FinalVertexInLoop(const std::list<const Edge *> & boundary)
{
	BabyUniverseDetector babyuniverse(triangulation_);
	std::vector<Triangle *> triangles;
	babyuniverse.EnclosedTriangles(boundary,triangles,false);
	int min=100000;
	Triangle * finalParent;
	for(int i=0,endi=distance_to_final_.size();i<endi;i++)
	{
		if( distance_to_final_[i] == 0 )
		{
			finalParent = triangulation_->getVertex(i)->getParent()->getParent();
		}
	}
	if( std::find(triangles.begin(),triangles.end(),finalParent) == triangles.end() )
	{
		return false;
	}
	return true;
}

void PeelingProcedure::DualRandomWalk(bool StayInMotherUniverse)
{
	Triangle * triangle = start_triangle_;

	int lastTimeInMotherUniverse = -1;
	int lastDistance = 0;
	for(int t=0;t<walk_time_;t++)
	{
		int minDistance = 1000000;
		for(int i=0;i<3;i++)
		{
			minDistance = std::min(minDistance,distance_[triangle->getEdge(i)->getOpposite()->getId()]);
		}
		if( StayInMotherUniverse )
		{
			distance_square_mother_[t] += minDistance * minDistance;
		} else
		{
			distance_square_mother_effective_[t] += lastDistance * lastDistance;
			distance_square_[t] += minDistance * minDistance;
		}

		Edge * edge = triangle->getEdge(triangulation_->RandomInteger(0,2));
		Triangle * nextTriangle = edge->getAdjacent()->getParent();

		if( !StayInMotherUniverse )
		{
			if( !in_mother_universe_[nextTriangle->getId()] && in_mother_universe_[triangle->getId()] )
			{
				lastTimeInMotherUniverse = t-1;
				lastDistance = minDistance;
			}
			if( in_mother_universe_[nextTriangle->getId()] && !in_mother_universe_[triangle->getId()] )
			{
				int waitingTime = t - lastTimeInMotherUniverse - 1;
				if( waitingTime <= static_cast<int>(waiting_times_.size()) )
				{
					waiting_times_[waitingTime-1]++;
				}
			}
		}

		if( !StayInMotherUniverse || in_mother_universe_[nextTriangle->getId()] )
		{
			triangle = nextTriangle;
		}

	}

	if( StayInMotherUniverse )
	{
		random_walk_measurements_mother_++;
	} else
	{
		random_walk_measurements_++;
	}
}

void PeelingProcedure::RandomWalk(bool StayInMotherUniverse)
{
	std::vector<int> degree;
	properties::DegreeList(triangulation_,degree);

	Vertex * vertex = start_vertex_;

	int lastTimeInMotherUniverse = -1;
	int lastDistance = 0;
	for(int t=0;t<walk_time_;t++)
	{
		int distance = distance_[vertex->getId()];

		if( StayInMotherUniverse )
		{
			distance_square_mother_[t] += distance * distance;
		} else
		{
			distance_square_mother_effective_[t] += lastDistance * lastDistance;
			distance_square_[t] += distance * distance;
		}

		Vertex * nextVertex = RandomNeighbour(vertex);

		if( !StayInMotherUniverse )
		{
			if( !vertex_in_mother_universe_[nextVertex->getId()] && vertex_in_mother_universe_[vertex->getId()] )
			{
				lastTimeInMotherUniverse = t-1;
				lastDistance = distance;
			}
			if( vertex_in_mother_universe_[nextVertex->getId()] && !vertex_in_mother_universe_[vertex->getId()] )
			{
				int waitingTime = t - lastTimeInMotherUniverse - 1;
				if( waitingTime <= static_cast<int>(waiting_times_.size()) )
				{
					waiting_times_[waitingTime-1]++;
				}
			}
		}

		if( !StayInMotherUniverse || vertex_in_mother_universe_[nextVertex->getId()] )
		{
			vertex = nextVertex;
		}

	}

	if( StayInMotherUniverse )
	{
		random_walk_measurements_mother_++;
	} else
	{
		random_walk_measurements_++;
	}
}

Vertex * PeelingProcedure::RandomNeighbour(Vertex * v)
{
	Edge * edge = v->getParent()->getPrevious();
	int degree = 0;
	do {
		degree++;
		edge = edge->getAdjacent()->getNext();
	} while( edge != v->getParent()->getPrevious() );
	int randomedge = triangulation_->RandomInteger(0,degree-1);
	degree = 0;
	do {
		if( degree == randomedge )
		{
			break;
		}
		degree++;
		edge = edge->getAdjacent()->getNext();
	} while( edge != v->getParent()->getPrevious() );
	return edge->getPrevious()->getOpposite();
}

void PeelingProcedure::DoMeasurementOnBabyUniverse(const std::list<const Edge*>::const_iterator & begin, const std::list<const Edge*>::const_iterator & end,int volume)
{
	if( baby_walk_ )
	{
		baby_walk_samples_ = 100;
		if( volume == -1 )
		{
			volume = babyuniversedetector_.VolumeEnclosed(begin,end,false);
		}
		int length = std::distance(begin,end);

		if( baby_return_time_.empty() )
		{
			baby_return_time_.resize(60,std::vector<std::vector<int> >(25,std::vector<int>(1000,0)));
			baby_universes_.resize(60,std::vector<int>(25,0));
		}
		baby_universes_[std::min(volume/10,59)][std::min(length/2,24)]++;
		for(int i=0;i<baby_walk_samples_;i++)
		{
			std::list<const Edge*>::const_iterator startEdge = std::next(begin,triangulation_->RandomInteger(0,length-1));
			const Edge* prevEdge = *std::prev( ( startEdge == begin ? end : startEdge ) );
			const Edge * edge = *startEdge;
			Vertex * startVertex = edge->getPrevious()->getOpposite();
			int index = 1;
			do {
				edge = edge->getAdjacent()->getNext();
				if( triangulation_->RandomInteger(0,index)==0 )
				{
					startVertex = edge->getPrevious()->getOpposite();
				}
				index++;
			} while( edge != prevEdge->getAdjacent() );

			int waitingtime = DoBabyWalk(startVertex,1000);
			if( waitingtime >= 0 )
			{
				baby_return_time_[std::min(volume/10,59)][std::min(length/2,24)][waitingtime]++;
			}
		}
	}
}

int PeelingProcedure::DoBabyWalk(Vertex * startVertex, int maxTime)
{
	for(int i=0;i<maxTime;i++)
	{
		if( in_frontier_[startVertex->getId()] )
		{
			return i;
		}
		startVertex = RandomNeighbour(startVertex);
	}
	return -1;
}