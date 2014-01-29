#include <queue>

#include "DiskCirclePacking.h"
#include "TriangulationProperties.h"

DiskCirclePacking::DiskCirclePacking(Triangulation * const triangulation)
	: triangulation_(triangulation),
	  babyuniversedetector_(triangulation)
{
	max_iterations_ = 5000;
	epsilon_ = 1.0e-5;
	delta_ = 0.04;
}


double DiskCirclePacking::Angle(Edge * edge)
{
	double radius = radius_[edge->getOpposite()->getId()];
	double nbrRadius1 = radius_[edge->getPrevious()->getOpposite()->getId()];
	double nbrRadius2 = radius_[edge->getNext()->getOpposite()->getId()];
	double x = radius * (1.0-nbrRadius1) / (1.0-nbrRadius1 * radius) * (1.0-nbrRadius2) / (1.0-nbrRadius2 * radius);
	BOOST_ASSERT( x >= 0.0 && x <= 1.0 );
	return 2.0 * std::asin(std::sqrt( x ) );
}


double DiskCirclePacking::AngleSum(int vertexId)
{
	Vertex * vertex = triangulation_->getVertex(vertexId);
	Edge * edge = vertex->getParent();
	double angle = 0.0;
	do {
		angle += Angle(edge);
		edge = edge->getNext()->getAdjacent()->getNext();
	} while( edge != vertex->getParent() );
	return angle;
}


bool DiskCirclePacking::FindDiskRadii(const std::list<const Edge*> & boundary)
{
	// Finds the radii of the disks in hyperbolic plane subject to the boundary conditions
	// that the vertices on the boundary have infinite (hyperbolic radius).
	// Notice that radius_[i] corresponds to exp(-2*r) where r is the hyperbolic radius.

	if( static_cast<int>(radius_.size()) != triangulation_->NumberOfVertices() )
	{
		radius_.resize(triangulation_->NumberOfVertices(),0.0);
	}
	double almostZero = 0.00001;

	disk_triangles_.clear();
	babyuniversedetector_.EnclosedTriangles(boundary,disk_triangles_,true);

	double startRadius = std::exp(-2.0 * std::sqrt(1.0/disk_triangles_.size()) );
	
	std::vector<bool> onBoundary(triangulation_->NumberOfVertices(),false);
	for(std::list<const Edge*>::const_iterator it = boundary.begin();it != boundary.end();it++)
	{
		onBoundary[(*it)->getNext()->getOpposite()->getId()] = true;
	}

	disk_vertices_.clear();
	std::vector<bool> indisk(triangulation_->NumberOfVertices(),false);
	for(int i=0,endi=disk_triangles_.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			Vertex * v = disk_triangles_[i]->getEdge(j)->getOpposite();
			if( onBoundary[v->getId()] )
			{
				radius_[v->getId()] = almostZero;
			}else if( !indisk[v->getId()] )
			{
				disk_vertices_.push_back(v);
				indisk[v->getId()] = true;
				if( radius_[v->getId()] <= 1.5 * almostZero || radius_[v->getId()] >= 1.0 )
				{
					radius_[v->getId()] = startRadius;
				}
			}
		}
	}

	std::vector<double> radius2(triangulation_->NumberOfVertices());

	properties::DegreeList(triangulation_,degree_);

	double c = epsilon_ + 1;
	double lambda = -1.0;
	bool flag = false;
	int steps = 0;
	while(c > epsilon_ )
	{
		double c0 = c;
		c = 0.0;
		double lambda0 = lambda;
		bool flag0 = flag;

		for(int i=0,endi=disk_vertices_.size();i<endi;i++)
		{
			// See Colling, Stephenson (2003): p. 244
			int vertexId = disk_vertices_[i]->getId();
			double theta = AngleSum(vertexId);
			double sintheta = std::sin(theta/(2*degree_[vertexId]));
			double sintarget = std::sin(PI/degree_[vertexId]);
			double vhat = (sintheta - std::sqrt(radius_[vertexId]))/(sintheta*radius_[vertexId] - std::sqrt(radius_[vertexId]));
			if( vhat < 0.0 )
			{
				vhat = 0.0;
			}
			double t = 2.0 * sintarget / (std::sqrt((1.0-vhat)*(1.0-vhat)+4.0*sintarget*sintarget*vhat)+1.0-vhat);
			radius2[vertexId] = t*t;
			c += (theta - 2.0 * PI)*(theta - 2.0 * PI);
		}
		c = std::sqrt(c);
		lambda = c / c0;
		flag = true;
		if( flag0 && lambda < 1.0 )
		{
			c *= lambda;
			if( std::abs( lambda - lambda0 ) < delta_ )
			{
				lambda = lambda / (1.0 - lambda);
			}
			double maxlambda = 1000.0;
			for(int i=0,endi=disk_vertices_.size();i<endi;i++)
			{
				int vertexId = disk_vertices_[i]->getId();
				if( radius2[vertexId] + maxlambda*(radius2[vertexId]-radius_[vertexId]) < 0.0 )
				{
					maxlambda = - radius2[vertexId] / (radius2[vertexId]-radius_[vertexId]);
				}
				if( radius2[vertexId] + maxlambda*(radius2[vertexId]-radius_[vertexId]) > 1.0 )
				{
					maxlambda = (1.0 - radius2[vertexId]) / (radius2[vertexId]-radius_[vertexId]);
				}
			}
			lambda = std::min(lambda, 0.5*maxlambda );
			for(int i=0,endi=disk_vertices_.size();i<endi;i++)
			{
				int vertexId = disk_vertices_[i]->getId();
				radius_[vertexId] = radius2[vertexId] + lambda * (radius2[vertexId]-radius_[vertexId]);
			}
			flag = false;
		} else
		{
			for(int i=0,endi=disk_vertices_.size();i<endi;i++)
			{
				int vertexId = disk_vertices_[i]->getId();
				radius_[vertexId] = radius2[vertexId];
			}
		}
		steps++;
		if( steps > max_iterations_ )
		{
			return false;
		}
	}

	return true;
}

bool DiskCirclePacking::FindEmbedding(const std::list<const Edge*> & boundary,  const Edge * centerEdge)
{
	if( !FindDiskRadii( boundary ) )
	{
		return false;
	}
	if( !DiskLayout(boundary,centerEdge) )
	{
		return false;
	}
	return true;
}

void DiskCirclePacking::getCircles(std::vector<std::pair<Vertex*,std::pair<Vector2D,double> > > & circles)
{
	circles.clear();
	for(int i=0,endi=circle_order_.size();i<endi;i++)
	{
		circles.push_back(std::pair<Vertex*,std::pair<Vector2D,double> >(triangulation_->getVertex(circle_order_[i]),std::pair<Vector2D,double>(coordinate_[circle_order_[i]],euclidean_radius_[circle_order_[i]])));
	}
}

bool DiskCirclePacking::DiskLayout(const std::list<const Edge*> & boundary, const Edge * centerEdge)
{
	Vertex * center = centerEdge->getNext()->getOpposite();
	
	circle_fixed_.resize(triangulation_->NumberOfVertices());
	std::fill(circle_fixed_.begin(),circle_fixed_.end(),false);
	euclidean_radius_.resize(triangulation_->NumberOfVertices());
	coordinate_.resize(triangulation_->NumberOfVertices());
	hyp_coordinate_.resize(triangulation_->NumberOfVertices());

	coordinate_[center->getId()] = MakeVector2D(0.0,0.0);
	hyp_coordinate_[center->getId()] = MakeVector2D(0.0,0.0);
	euclidean_radius_[center->getId()] = (1.0 - std::sqrt(radius_[center->getId()]))/(1.0 + std::sqrt(radius_[center->getId()]));
	circle_fixed_[center->getId()] = true;

	Vertex * endpoint = centerEdge->getPrevious()->getOpposite();
	double y0 = (1.0 - std::sqrt(radius_[center->getId()])*radius_[endpoint->getId()])/(1.0 + std::sqrt(radius_[center->getId()])*radius_[endpoint->getId()]);
	euclidean_radius_[endpoint->getId()] = 0.5*(y0 - euclidean_radius_[center->getId()]);
	coordinate_[endpoint->getId()] = MakeVector2D(0.0,euclidean_radius_[center->getId()]+euclidean_radius_[endpoint->getId()]);
	hyp_coordinate_[endpoint->getId()] = MakeVector2D(0.0,(1.0 - std::sqrt(radius_[center->getId()]*radius_[endpoint->getId()]))/(1.0 + std::sqrt(radius_[center->getId()]*radius_[endpoint->getId()])));
	circle_fixed_[endpoint->getId()] = true;
	
	std::vector<bool> visited(triangulation_->NumberOfTriangles(),true);
	for(int i=0,endi=disk_triangles_.size();i<endi;i++)
	{
		visited[disk_triangles_[i]->getId()] = false;
	}
	std::queue<const Edge* > q;
	q.push(centerEdge);
	visited[centerEdge->getParent()->getId()] = true;
	q.push(centerEdge->getAdjacent());
	visited[centerEdge->getAdjacent()->getParent()->getId()] = true;

	circle_order_.clear();
	circle_order_.push_back(center->getId());
	circle_order_.push_back(endpoint->getId());

	while( !q.empty() )
	{
		const Edge * edge = q.front();
		q.pop();

		int newvertId = edge->getOpposite()->getId();
		if( !circle_fixed_[newvertId] )
		{
			// assign radius and coordinate to vertex
			double theta = Angle(edge->getNext());
			int pivotId = edge->getNext()->getOpposite()->getId();
			Vector2D z = RotateVector2D(MakeVector2D((1.0-std::sqrt(radius_[pivotId]*radius_[newvertId]))/(1.0+std::sqrt(radius_[pivotId]*radius_[newvertId])),0.0),theta);
			hyp_coordinate_[newvertId] = InverseMobiusTransform(z,hyp_coordinate_[pivotId],hyp_coordinate_[edge->getPrevious()->getOpposite()->getId()]);
			
			std::pair<Vector2D,double> circle = EuclideanCircle(hyp_coordinate_[newvertId],radius_[newvertId]);
			coordinate_[newvertId] = circle.first;
			euclidean_radius_[newvertId] = circle.second;
			circle_fixed_[edge->getOpposite()->getId()] = true;
			circle_order_.push_back(edge->getOpposite()->getId());
		}

		const Edge * nextedge = edge;
		for(int i=0;i<2;i++)
		{
			nextedge = nextedge->getNext();
			if( !visited[nextedge->getAdjacent()->getParent()->getId()] )
			{
				visited[nextedge->getAdjacent()->getParent()->getId()] = true;
				q.push(nextedge->getAdjacent());
			}
		}
	}
	return true;
}