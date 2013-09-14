#include <algorithm>
#include <math.h>
#include <vector>
#include <queue>
#include <map>

#include "ThetaModel.h"


ThetaModel::ThetaModel(Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis, int PiInUnits) : triangulation_(triangulation), dualcohomologybasis_(dualcohomologybasis), pi_in_units_(PiInUnits)
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
			integral[i] += dualcohomologybasis_->getOmega((this->*next)(moveEdge),i);
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

	/* TESTS
	BOOST_ASSERT(getTheta(kiteEdge) == getTheta(kiteEdge->getAdjacent()));
	BOOST_ASSERT(getTheta(moveEdge) == getTheta(moveEdge->getAdjacent()));
	BOOST_ASSERT(getTheta(moveEdge->getPrevious()) == getTheta(moveEdge->getPrevious()->getAdjacent()));


	BOOST_ASSERT(TestVertexSum(kiteEdge->getNext()->getOpposite()));
	BOOST_ASSERT(TestVertexSum(kiteEdge->getPrevious()->getOpposite()));
	BOOST_ASSERT(TestVertexSum(moveEdge->getNext()->getOpposite()));
	BOOST_ASSERT( TestVertexSum(mirror ? moveEdge->getPrevious()->getOpposite() : otherEdge->getPrevious()->getOpposite()) );
	*/
	
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


struct triangleNode
{
	Triangle* triangle;
	IntForm2D integral;
};
bool operator<(const triangleNode &leftNode, const triangleNode &rightNode) {
	if (leftNode.triangle != rightNode.triangle) return leftNode.triangle < rightNode.triangle;
	if(leftNode.integral[0] != rightNode.integral[0] ) return leftNode.integral[0] < rightNode.integral[0];
	if(leftNode.integral[1] != rightNode.integral[1] ) return leftNode.integral[1] < rightNode.integral[1];
	return false;
}
bool operator<(const std::pair<triangleNode,int> &leftNode, const std::pair<triangleNode,int> &rightNode) {
	if (leftNode.second != rightNode.second) return leftNode.second > rightNode.second;
	if (leftNode.first.triangle != rightNode.first.triangle) return leftNode.first.triangle > rightNode.first.triangle;
	if(leftNode.first.integral[0] != rightNode.first.integral[0] ) return leftNode.first.integral[0] > rightNode.first.integral[0];
	if(leftNode.first.integral[1] != rightNode.first.integral[1] ) return leftNode.first.integral[1] > rightNode.first.integral[1];
	return false;
}

bool ThetaModel::TestCutCondition(Edge * fromEdge, Edge * toEdge, const IntForm2D & integral, int totalTheta) const
{
	// Return false if there exists a contractible path in the dual graph from
	// fromEdge->getParent() to toEdge->getParent() avoiding the dual edge fromEdge
	// and having total theta smaller than totalTheta.

	// Perform a Dijkstra algorithm on the (universal covering of) the dual graph	
	std::priority_queue<std::pair<triangleNode,int> > q;
	std::map<triangleNode,int> visited;
	triangleNode startNode;
	startNode.triangle = fromEdge->getParent();
	startNode.integral = integral;
	q.push(std::pair<triangleNode,int>(startNode,0));
	visited.insert(std::pair<triangleNode,int>(startNode,0));

	while( !q.empty() )
	{
		std::pair<triangleNode,int> node = q.top();
		q.pop();

		for(int i=0;i<3;i++)
		{
			Edge * edge = node.first.triangle->getEdge(i);
			if( edge == fromEdge )
				continue;

			std::pair<triangleNode,int> nextNode(node);			
			nextNode.first.triangle = edge->getAdjacent()->getParent();
			nextNode.second += getTheta(edge);
			if( nextNode.second < totalTheta )	// we only need to search nodes that are closer than totalTheta from the start
			{
				nextNode.first.integral[0] += dualcohomologybasis_->getOmega(edge,0);
				nextNode.first.integral[1] += dualcohomologybasis_->getOmega(edge,1);
				if( nextNode.first.triangle == toEdge->getParent() && nextNode.first.integral[0] == 0 && nextNode.first.integral[1] == 0 )
				{
					// found a path with total theta smaller than totalTheta
					return false;
				}
				std::pair<std::map<triangleNode,int>::iterator, bool> returnValue = visited.insert(nextNode);
				if( returnValue.second )
				{
					q.push(nextNode);
				}else
				{
					if( nextNode.second < returnValue.first->second )
					{
						returnValue.first->second = nextNode.second;
					}
				}
			}
		}
	}
	return true;
}