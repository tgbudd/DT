#include <queue>

#include "CirclePacking.h"
#include "TriangulationProperties.h"

CirclePacking::CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis)
	: triangulation_(triangulation), Embedding(triangulation,cohomologybasis), babyuniversedetector_(triangulation)
{
	max_iterations_ = 5000;
	epsilon_ = 1.0e-5;
	delta_ = 0.04;
}

double CirclePacking::Angle(Edge * edge)
{
	double radius = radius_[edge->getOpposite()->getId()];
	double nbrRadius1 = radius_[edge->getPrevious()->getOpposite()->getId()];
	double nbrRadius2 = radius_[edge->getNext()->getOpposite()->getId()];
	if( nbrRadius1 >= 0.0 || nbrRadius2 >= 0.0 )
	{
		double x = nbrRadius1 / (nbrRadius1 + radius) * nbrRadius2 / (nbrRadius2 + radius);
		BOOST_ASSERT( x >= 0.0 && x < 1.0 );
		return 2.0 * std::asin(std::sqrt( x ) );
	}
	else 
		return 0.0;
}

double CirclePacking::AngleSum(int vertexId)
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

bool CirclePacking::FindRadii()
{
	// given the theta, determine the radii of the circles in the circle pattern

	if( static_cast<int>(radius_.size()) != triangulation_->NumberOfVertices() )
	{
		radius_.resize(triangulation_->NumberOfVertices());
	}
	std::fill(radius_.begin(),radius_.end(),1.0);
	double startRadius = std::sqrt(1.0/triangulation_->NumberOfVertices());
	std::vector<double> radius2(triangulation_->NumberOfVertices(),startRadius);

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

		//double totradiussq = 0.0;
		for(int i=0,endi=radius_.size();i<endi;i++)
		{
			double theta = AngleSum(i);
			double sintheta = std::sin(theta/(2*degree_[i]));
			double sintarget = std::sin(PI/degree_[i]);
			radius2[i] = (1.0 - sintarget)/sintarget * sintheta / (1.0 - sintheta) * radius_[i];
			//totradiussq += radius2[i]*radius2[i];
			c += (theta - 2.0 * PI)*(theta - 2.0 * PI);
		}
		//double normalize = 1.0/std::sqrt(totradiussq);
		//for(int i=0,endi=radius_.size();i<endi;i++)
		//{
		//	radius2[i] *= normalize;
		//}
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
			for(int i=0,endi=radius_.size();i<endi;i++)
			{
				if( radius2[i] + maxlambda*(radius2[i]-radius_[i]) < 0.0 )
				{
					maxlambda = - radius2[i] / (radius2[i]-radius_[i]);
				}
			}
			lambda = std::min(lambda, 0.5*maxlambda );
			for(int i=0,endi=radius_.size();i<endi;i++)
			{
				radius_[i] = radius2[i] + lambda * (radius2[i]-radius_[i]);
			}
			flag = false;
		} else
		{
			std::copy(radius2.begin(),radius2.end(),radius_.begin());
		}
		steps++;
		if( steps > max_iterations_ )
		{
			return false;
		}
	}

	double area=0.0;
	for(int i=0,endi=radius_.size();i<endi;i++)
	{
		area += PI*radius_[i]*radius_[i];
	}
	double norm = 1.0/std::sqrt(area);
	for(int i=0,endi=radius_.size();i<endi;i++)
	{
		radius_[i] *= norm;
	}
	return true;
}

bool CirclePacking::FindEdgeMeasure()
{
	if( !FindRadii() )
	{
		return false;
	}

	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(j);
			Edge * adjEdge = edge->getAdjacent();
			double measure = 0.5 / std::tan( Angle(edge) ) + 0.5 / std::tan( Angle(adjEdge) );
			setEdgeMeasure(i,j,measure);
		}
	}

	return true;
}

void CirclePacking::CalculateEdgeMeasure(std::vector<boost::array<double,3> > & measure)
{
	if( static_cast<int>(measure.size()) !=  triangulation_->NumberOfTriangles() )
	{
		measure.resize(triangulation_->NumberOfTriangles());
	}
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(j);
			Edge * adjEdge = edge->getAdjacent();
			double angle1 = Angle(edge), angle2 = Angle(adjEdge);
			if( angle1 > 1.0e-7 && angle2 > 1.0e-7 )
			{
				measure[i][j] = 0.5 / std::tan( Angle(edge) ) + 0.5 / std::tan( Angle(adjEdge) );
			} else
			{
				measure[i][j] = 1.0;
			}
		}
	}
}

bool CirclePacking::GetRadii(std::vector<double> & radii)
{
	radii.clear();
	radii.reserve(triangulation_->NumberOfVertices());

	radii = radius_;

	return true;
}

bool CirclePacking::FindDiskRadii(const std::list<const Edge*> & boundary)
{
	if( static_cast<int>(radius_.size()) != triangulation_->NumberOfVertices() )
	{
		radius_.resize(triangulation_->NumberOfVertices(),-1.0);
	}
	std::fill(radius_.begin(),radius_.end(),-1.0);

	disk_triangles_.clear();
	babyuniversedetector_.EnclosedTriangles(boundary,disk_triangles_,true);
	double startRadius = std::sqrt(1.0/disk_triangles_.size());

	disk_vertices_.clear();
	for(int i=0,endi=disk_triangles_.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			Vertex * v = disk_triangles_[i]->getEdge(j)->getOpposite();
			if( radius_[v->getId()] < 0.0 )
			{
				disk_vertices_.push_back(v);
				radius_[v->getId()] = startRadius;
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
			int vertexId = disk_vertices_[i]->getId();
			double theta = AngleSum(vertexId);
			double sintheta = std::sin(theta/(2*degree_[vertexId]));
			double sintarget = std::sin(PI/degree_[vertexId]);
			radius2[vertexId] = (1.0 - sintarget)/sintarget * sintheta / (1.0 - sintheta) * radius_[vertexId];
			if( radius2[vertexId] > 0.5 )
			{
				int hy=0;
			}
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
				if( radius2[vertexId] + maxlambda*(radius2[vertexId]-radius_[vertexId]) >= 0.5 )
				{
					maxlambda = (0.5 - radius2[vertexId]) / (radius2[vertexId]-radius_[vertexId]);
				}
			}
			lambda = std::min(lambda, 0.2*maxlambda );
			for(int i=0,endi=disk_vertices_.size();i<endi;i++)
			{
				int vertexId = disk_vertices_[i]->getId();
				radius_[vertexId] = radius2[vertexId] + lambda * (radius2[vertexId]-radius_[vertexId]);
				if( radius_[vertexId] > 0.5 )
				{
					int hy=0;
				}
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

	/// TEST ////
	std::vector<double> anglesum;
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		anglesum.push_back(AngleSum(vertexId));
	}

	return true;
}

bool CirclePacking::DiskLayout(const std::list<const Edge*> & boundary)
{
	if( static_cast<int>(disk_coordinate_.size()) != triangulation_->NumberOfVertices() )
	{
		disk_coordinate_.resize(triangulation_->NumberOfVertices());
	}
	std::vector<bool> coordinate_fixed(triangulation_->NumberOfVertices(),false);

	double phi = 0.0;
	for(std::list<const Edge*>::const_iterator it = boundary.begin();it!=boundary.end();it++)
	{
		int start = (*it)->getNext()->getOpposite()->getId();
		int end = (*it)->getPrevious()->getOpposite()->getId();
		if( !coordinate_fixed[start] )
		{
			disk_coordinate_[end] = MakeVector2D(1.0 - radius_[end],0.0);
		} else
		{
			double a = 1.0-radius_[start], b = 1.0-radius_[end], c = radius_[start]+radius_[end];
			double theta = std::acos((a*a + b*b - c*c)/(2.0*a*b));
			phi += theta;
			disk_coordinate_[end][0] = a * std::cos(phi);
			disk_coordinate_[end][0] = a * std::sin(phi);
		}
		coordinate_fixed[end] = true;
	}

	std::vector<const Vertex*> free_vertices;
	std::vector<const Vertex*> fixed_vertices;
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		if( !coordinate_fixed[disk_vertices_[i]->getId()] )
		{
			free_vertices.push_back(disk_vertices_[i]);
		} else
		{
			fixed_vertices.push_back(disk_vertices_[i]);
		}
	}

	CalculateEdgeMeasure(disk_edge_measure_);
	DomainLaplacianMatrix lapl(triangulation_,disk_edge_measure_,free_vertices,fixed_vertices);

	std::vector<double> x(triangulation_->NumberOfVertices());
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			x[j] = disk_coordinate_[j][i];
		}
		if( !lapl.FindHarmonic(x,1e-6) )
		{
			return false;
		}
		for(int j=0;j<triangulation_->NumberOfVertices();j++)
		{
			disk_coordinate_[j][i] = x[j];
		}
	}
	return true;
}

bool CirclePacking::FindDiskEmbedding(const std::list<const Edge*> & boundary,  const Edge * centerEdge)
{
	if( !FindDiskRadii( boundary ) )
	{
		return false;
	}
	if( !DiskLayout( boundary ) )
	{
		return false;
	}
//	MobiusTransformation(centerEdge);

	return true;
}

void CirclePacking::MobiusTransformation(const Edge * centerEdge)
{
	Vector2D start = disk_coordinate_[centerEdge->getNext()->getOpposite()->getId()];
	Vector2D end = disk_coordinate_[centerEdge->getPrevious()->getOpposite()->getId()];
	double startradius = radius_[centerEdge->getNext()->getOpposite()->getId()];
	double endradius = radius_[centerEdge->getPrevious()->getOpposite()->getId()];
	double angle = 0.5*PI - VectorAngle( SubtractVectors2D( end, start ) );
	Vector2D center = AddScaledVectors2D( endradius/(startradius+endradius), start, startradius/(startradius+endradius), end);
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		std::pair<Vector2D,double> newcircle = CircleMobius(disk_coordinate_[i],radius_[i],center,angle);
		disk_coordinate_[i] = newcircle.first;
		radius_[i] = newcircle.second;
	}
}

std::pair<Vector2D,double> CirclePacking::CircleMobius(Vector2D c, double radius, Vector2D z0, double angle)
{
	double phi = triangulation_->RandomReal(0.0,2.0*PI);
	boost::array<Vector2D,3> v;
	for(int i=0;i<3;i++)
	{
		v[i] = c;
		v[i][0] += std::cos(phi) * radius;
		v[i][1] += std::sin(phi) * radius;
		v[i] = MobiusTransform(v[i],z0,angle);
	}
	Vector2D c2 = CenterOfCircle(v);
	return std::pair<Vector2D,double>(c2,Norm2D(SubtractVectors2D(c2,v[0])));
}

void CirclePacking::getCircles(std::vector<std::pair<Vector2D,double> > & circles)
{
	circles.clear();
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		circles.push_back(std::pair<Vector2D,double>(disk_coordinate_[disk_vertices_[i]->getId()],radius_[disk_vertices_[i]->getId()]));
	}
}