#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <map>

#include "ThetaModel.h"


ThetaModel::ThetaModel(Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis, int PiInUnits) : triangulation_(triangulation), dualcohomologybasis_(dualcohomologybasis), pi_in_units_(PiInUnits), use_cos_power_(false)
{
}


ThetaModel::~ThetaModel(void)
{
}

void ThetaModel::Initialize()
{
	// For now we have to assume that the lattice is regular since we have 
	// no algorithm to determine an admissable set of theta's for an arbitrary triangulation.
	for(int i=0;i<triangulation_->NumberOfVertices();i++)
	{
		BOOST_ASSERT(triangulation_->getVertex(i)->getDegree() == 6);
	}

	BOOST_ASSERT( pi_in_units_ % 3 == 0 );
	boost::array<int,3> equilateral = {pi_in_units_/3,pi_in_units_/3,pi_in_units_/3};
	theta_.assign(triangulation_->NumberOfTriangles(),equilateral);
	
}

void ThetaModel::TryThetaMove(int n)
{
	int SuccesfulMoves = 0;
	while( SuccesfulMoves < n)
	{
		if( TryThetaMove() )
			SuccesfulMoves++;
	}
}

void ThetaModel::DoSweep()
{
	int SuccesfulMoves = 0;
	int NumberOfEdges = 3*triangulation_->NumberOfTriangles()/2;
	while( SuccesfulMoves < NumberOfEdges)
	{
		if( TryThetaMove() )
			SuccesfulMoves++;
	}
}

bool ThetaModel::TryThetaMove()
{
	return TryThetaMove(triangulation_->getRandomEdge());
}

/// <image url="$(SolutionDir)\images\ThetaMove.png" scale="1.1" />

bool ThetaModel::TryThetaMove(Edge * moveEdge)
{
	// In order to deal with the move and its mirror in one go
	// we define the function pointers for getNext() and getPrevious()
	Edge * (ThetaModel::*next)(Edge *) = &ThetaModel::nextEdge;
	Edge * (ThetaModel::*previous)(Edge *) = &ThetaModel::previousEdge;
	bool mirror = (triangulation_->RandomInteger(0,1) == 0);
	if( mirror )
	{
		std::swap(next,previous);
	} 
	Edge * kiteEdge = (this->*next)(moveEdge)->getAdjacent();

	if( (this->*previous)(kiteEdge)->getAdjacent() == moveEdge || (this->*next)(kiteEdge)->getAdjacent() == (this->*previous)(moveEdge) )
	{
		// The kite is adjacent to itself, which is not allowed.
		return false;
	}
	
	int moveEdgeTheta = getTheta(moveEdge);
	
	int dThetaMin = 1 + std::max(
		std::max( -getTheta((this->*previous)(moveEdge)), -getTheta((this->*next)(kiteEdge)) ),
		std::max(  moveEdgeTheta - pi_in_units_,  getTheta((this->*previous)(kiteEdge))-pi_in_units_) );
	int dThetaMax = -1 + std::min(
		std::min( pi_in_units_ - getTheta((this->*previous)(moveEdge)), pi_in_units_ - getTheta((this->*next)(kiteEdge))),
		getTheta((this->*previous)(kiteEdge)));

	if( dThetaMax > moveEdgeTheta )
	{
		if( (this->*previous)(moveEdge->getAdjacent())->getAdjacent() == (this->*previous)(kiteEdge) )
		{
			dThetaMax = moveEdgeTheta - 1;
		} else
		{
			// Extend range to allow for a flip move on moveEdge
			dThetaMax = -1 + std::min( getTheta((this->*previous)(moveEdge->getAdjacent())) + moveEdgeTheta , 
				std::min( getTheta((this->*previous)(kiteEdge)) , pi_in_units_ - getTheta((this->*next)(kiteEdge)) ) );
		}
	}
	// Choose dTheta uniformly at random
	int dTheta = triangulation_->RandomInteger(dThetaMin, dThetaMax);

	if( dTheta == 0 )
	{
		// nothing will change 
		return true;
	}
	if( dTheta == moveEdgeTheta )
	{
		return false;
	}

	if( use_cos_power_ )
	{
		// Boltzmann factor is product over edges of |cos(theta)|^(cos_power_).
		double BoltzmannChange = 1.0;

		BoltzmannChange *= cosine(getTheta((this->*previous)(kiteEdge))-dTheta) / cosine(getTheta((this->*previous)(kiteEdge)));
		BoltzmannChange *= cosine(getTheta((this->*next)(kiteEdge))+dTheta) / cosine(getTheta((this->*next)(kiteEdge)));
		BoltzmannChange *= cosine(getTheta((this->*previous)(moveEdge))-dTheta) / cosine(getTheta((this->*previous)(moveEdge)));
		if( dTheta < moveEdgeTheta )
		{
			BoltzmannChange *= cosine(moveEdgeTheta-dTheta) /cosine(moveEdgeTheta);
		} else
		{
			BoltzmannChange *= cosine( getTheta((this->*previous)(moveEdge->getAdjacent())) - dTheta + moveEdgeTheta ) / cosine( getTheta((this->*previous)(moveEdge->getAdjacent())));
			BoltzmannChange *= cosine(dTheta-moveEdgeTheta)/cosine(moveEdgeTheta);
		}
		BoltzmannChange = std::pow(std::fabs(BoltzmannChange),cos_power_);

		if(	!triangulation_->SucceedWithProbability( BoltzmannChange ) )
		{
			return false;
		}
	}

	if( dTheta < 0 )
	{
		// Construct the integral of omega which is needed by TestCutCondition to distinguish contractible
		// paths from non-contractible ones.
		IntForm2D integral = {0,0};
		for(int i=0;i<2;i++)
		{
			integral[i] -= dualcohomologybasis_->getOmega((this->*next)(kiteEdge),i);
			integral[i] += dualcohomologybasis_->getOmega(kiteEdge,i);
			integral[i] += dualcohomologybasis_->getOmega((this->*previous)(moveEdge),i);
		}
		if( !TestCutCondition( (this->*previous)(moveEdge)->getAdjacent(),
							   (this->*next)(kiteEdge)->getAdjacent(), 
							   integral, 
							   2*pi_in_units_ - getTheta((this->*previous)(moveEdge)) - getTheta(kiteEdge) - getTheta((this->*next)(kiteEdge)) - 2*dTheta ) )
		{
			return false;
		}
	} else
	{
		IntForm2D integral = {0,0};
		for(int i=0;i<2;i++)
		{
			integral[i] -= dualcohomologybasis_->getOmega((this->*previous)(kiteEdge),i);
			integral[i] += dualcohomologybasis_->getOmega(kiteEdge,i);
			integral[i] += dualcohomologybasis_->getOmega(moveEdge,i);
		}

		if( !TestCutCondition( moveEdge->getAdjacent(),
							   (this->*previous)(kiteEdge)->getAdjacent(),
							   integral,
							   2*pi_in_units_ - moveEdgeTheta - getTheta(kiteEdge) - getTheta((this->*previous)(kiteEdge)) + 2*std::min(dTheta,moveEdgeTheta)) )
		{
			return false;
		}
		if( dTheta > moveEdgeTheta )
		{
			for(int i=0;i<2;i++)
			{
				integral[i] += dualcohomologybasis_->getOmega((this->*previous)(moveEdge->getAdjacent()),i);
			}
			if( !TestCutCondition( (this->*previous)(moveEdge->getAdjacent())->getAdjacent(),
								   (this->*previous)(kiteEdge)->getAdjacent(),
							       integral,
								   2*pi_in_units_ - moveEdgeTheta - getTheta(kiteEdge) - getTheta((this->*previous)(kiteEdge)) - getTheta((this->*previous)(moveEdge->getAdjacent())) + 2*dTheta) )
			{
			   return false;
			}
		}
	}

	// The kite move with change dTheta is possible as far as the thetas are concerned 

	Edge * otherEdge = moveEdge->getAdjacent();
	int otherEdgePreviousTheta = getTheta((this->*previous)(otherEdge));
	if( dTheta > moveEdgeTheta )
	{
		// perform a flip on moveEdge
		if( !triangulation_->TryFlipMove(moveEdge) )
		{
			return false;
		}
	}

	int moveEdgeNextTheta = getTheta((this->*next)(moveEdge));
	int moveEdgePrevTheta = getTheta((this->*previous)(moveEdge));
	setTheta((this->*previous)(kiteEdge),getTheta((this->*previous)(kiteEdge))-dTheta);
	setTheta((this->*next)(kiteEdge),getTheta((this->*next)(kiteEdge))+dTheta);

	if( dTheta < moveEdgeTheta )
	{
		setTheta((this->*previous)(moveEdge),moveEdgePrevTheta + dTheta);
		setTheta(moveEdge,moveEdgeTheta - dTheta);
	} else
	{
		if( !mirror )
		{
			setTheta(otherEdge,moveEdgeTheta+moveEdgePrevTheta);
			setTheta(moveEdge,otherEdgePreviousTheta-dTheta+moveEdgeTheta);
			setTheta(moveEdge->getPrevious(),dTheta - moveEdgeTheta);
		} else
		{
			setTheta(moveEdge->getNext(),moveEdgePrevTheta+moveEdgeTheta);
			setTheta(otherEdge->getNext(),otherEdgePreviousTheta-dTheta+moveEdgeTheta);
			setTheta(moveEdge->getPrevious(),dTheta - moveEdgeTheta);
			setTheta(otherEdge,getTheta(kiteEdge));
			setTheta(moveEdge,getTheta(moveEdge->getAdjacent()));
		}
	}

	/*// TESTS
	BOOST_ASSERT(getTheta(kiteEdge) == getTheta(kiteEdge->getAdjacent()));
	BOOST_ASSERT(getTheta(moveEdge) == getTheta(moveEdge->getAdjacent()));
	BOOST_ASSERT(getTheta(moveEdge->getPrevious()) == getTheta(moveEdge->getPrevious()->getAdjacent()));
	
	BOOST_ASSERT(TestVertexSum(kiteEdge->getNext()->getOpposite()));
	BOOST_ASSERT(TestVertexSum(kiteEdge->getPrevious()->getOpposite()));
	BOOST_ASSERT(TestVertexSum(moveEdge->getNext()->getOpposite()));
	BOOST_ASSERT( TestVertexSum(mirror ? moveEdge->getPrevious()->getOpposite() : otherEdge->getPrevious()->getOpposite()) );
	
	BOOST_ASSERT( TestCutCondition(kiteEdge) == TestCutCondition(kiteEdge->getAdjacent()) );

	BOOST_ASSERT( TestCutCondition(kiteEdge) );
	BOOST_ASSERT( TestCutCondition(moveEdge) );
	BOOST_ASSERT( TestCutCondition(moveEdge->getNext()) );
	BOOST_ASSERT( TestCutCondition(moveEdge->getPrevious()) );
	BOOST_ASSERT( TestCutCondition(otherEdge) );
	BOOST_ASSERT( TestCutCondition(otherEdge->getNext()) );
	BOOST_ASSERT( TestCutCondition(otherEdge->getPrevious()) );
	BOOST_ASSERT( TestCutCondition(kiteEdge->getNext()) );
	BOOST_ASSERT( TestCutCondition(kiteEdge->getPrevious()) );*/

	return true;
}

Edge * ThetaModel::nextEdge(Edge * edge)
{
	return edge->getNext();
}
Edge * ThetaModel::previousEdge(Edge * edge)
{
	return edge->getPrevious();
}

bool ThetaModel::TestVertexSum(Vertex * vertex)
{
	int totalTheta = 0;
	Edge * initialEdge = vertex->getParent()->getPrevious();
	Edge * edge = initialEdge;
	do
	{
		totalTheta += getTheta(edge);
		edge = edge->getPrevious()->getAdjacent();
	} while( edge != initialEdge );

	BOOST_ASSERT( totalTheta == 2*pi_in_units_ );

	return totalTheta == 2*pi_in_units_;
}


bool ThetaModel::TestCutCondition(Edge * fromEdge, Edge * toEdge, const IntForm2D & integral, int totalTheta) const
{
	// Return false if there exists a contractible path in the dual graph from
	// fromEdge->getParent() to toEdge->getParent() avoiding the dual edge fromEdge
	// and having total theta smaller than totalTheta.

	// Perform a Dijkstra algorithm on the (universal covering of) the dual graph	
	if( distance_.empty() )
	{
		distance_.resize(triangulation_->NumberOfTriangles());
	}

	std::priority_queue<triangleNode> q;
	triangleNode startNode;
	startNode.triangle = fromEdge->getParent();
	startNode.integral = integral;
	startNode.distance = 0;
	q.push(startNode);
	distance_[startNode.triangle->getId()].insert(std::pair<IntForm2D,int>(startNode.integral,startNode.distance));

	while( !q.empty() )
	{
		triangleNode node = q.top();
		q.pop();

		if( node.distance > distance_[node.triangle->getId()][node.integral] )
			continue;

		for(int i=0;i<3;i++)
		{
			Edge * edge = node.triangle->getEdge(i);
			if( edge == fromEdge )
				continue;

			triangleNode nextNode(node);			
			nextNode.triangle = edge->getAdjacent()->getParent();
			nextNode.distance += getTheta(edge);
			if( nextNode.distance <= totalTheta )	// we only need to search nodes that are closer than totalTheta from the start
			{
				nextNode.integral = AddForms( nextNode.integral, dualcohomologybasis_->getOmega(edge) );
				if( nextNode.triangle == toEdge->getParent() && FormIsZero(nextNode.integral) )
				{
					// found a path with total theta smaller than (or equal to) totalTheta
					CleanUpDistance(startNode.triangle);
					return false;
				}
				std::pair<std::map<IntForm2D,int>::iterator,bool> returnValue = distance_[nextNode.triangle->getId()].insert(std::pair<IntForm2D,int>(nextNode.integral,nextNode.distance));
				if( returnValue.second )
				{
					q.push(nextNode);
				}else
				{
					if( nextNode.distance < returnValue.first->second )
					{
						returnValue.first->second = nextNode.distance;
						q.push(nextNode);
					}
				}
			}
		}
	}
	CleanUpDistance(startNode.triangle);
	return true;
}

void ThetaModel::CleanUpDistance(Triangle * triangle) const
{
	// The non-empty maps in distance_ make up a connected subset of the triangulation, therefore to
	// clean up perform a breadth-first search.

	std::queue<Triangle *> q;
	q.push(triangle);

	while( !q.empty() )
	{
		Triangle * currentTriangle = q.front();
		q.pop();

		distance_[currentTriangle->getId()].clear();

		for(int i=0;i<3;i++)
		{
			Triangle * nbrTriangle = currentTriangle->getEdge(i)->getAdjacent()->getParent();
			if( !distance_[nbrTriangle->getId()].empty() )
			{
				q.push(nbrTriangle);
			}
		}
	}
}

bool ThetaModel::TestAllCutConditions() const
{
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);
			if( edge->getAdjacent()->getParent() < triangle )
				continue;
			if( !TestCutCondition(edge) )
				return false;
		}
	}
	return true;
}

bool ThetaModel::TestCutCondition(Edge * edge) const
{
	IntForm2D integral = dualcohomologybasis_->getOmega(edge);
	return TestCutCondition(edge->getAdjacent(),edge,integral,-1 + 2*pi_in_units_ - getTheta(edge));
}

std::string ThetaModel::ExportState() const
{
	std::ostringstream stream;
	stream << "piInUnits -> " << pi_in_units_ << ", theta -> ";
	PrintToStream2D(stream, theta_.begin(), theta_.end() );
	return stream.str();
}

