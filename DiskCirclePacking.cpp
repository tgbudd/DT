#include <queue>
#include <fstream>

//#define LOG_CIRCLE_PACKING

#include "DiskCirclePacking.h"
#include "TriangulationProperties.h"

DiskCirclePacking::DiskCirclePacking(Triangulation * const triangulation)
	: triangulation_(triangulation),
	  babyuniversedetector_(triangulation)
{
	use_one_minus_radius_ = false;
	max_iterations_ = 5000;
	epsilon_ = 1.0e-6;
	delta_ = 0.001;
	//delta_ = 0.001;
	reset_radius_ = true;
}


double DiskCirclePacking::Angle(Edge * edge)
{
	double radius = radius_[edge->getOpposite()->getId()];
	double nbrRadius1 = radius_[edge->getPrevious()->getOpposite()->getId()];
	double nbrRadius2 = radius_[edge->getNext()->getOpposite()->getId()];
	double x;
	if( use_one_minus_radius_ )
	{
		x = (1.0-radius) * nbrRadius1 / (nbrRadius1 * (1.0-radius) + radius) * nbrRadius2 / (nbrRadius2 * (1.0-radius) + radius);
	} else
	{
		x = radius * (1.0-nbrRadius1) / (1.0-nbrRadius1 * radius) * (1.0-nbrRadius2) / (1.0-nbrRadius2 * radius);
	}
	BOOST_ASSERT( x >= 0.0 && x <= 1.0 );
	return 2.0 * std::asin(std::sqrt( x ) );
}

double DiskCirclePacking::Angle(Edge * edge, const std::vector<int> & vertexPos, const std::vector<double> & radius)
{
	double radius0 = radius[vertexPos[edge->getOpposite()->getId()]];
	double radius1 = radius[vertexPos[edge->getNext()->getOpposite()->getId()]];
	double radius2 = radius[vertexPos[edge->getPrevious()->getOpposite()->getId()]];
	return std::acos(1.0 - 2.0 * radius1 * radius2 / ((radius0 + radius1)*(radius0+radius2)));
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


	double startRadius = std::exp(-2.0 * std::sqrt(1.0/disk_triangles_.size()) );
	
	std::vector<bool> onBoundary(triangulation_->NumberOfVertices(),false);
	for(std::list<const Edge*>::const_iterator it = boundary.begin();it != boundary.end();it++)
	{
		onBoundary[(*it)->getNext()->getOpposite()->getId()] = true;
	}

#ifdef LOG_CIRCLE_PACKING
	std::vector<double> radius3(radius_);

lab:
	std::cin >> delta_;

	std::ofstream file("log.txt");
	file << std::fixed << "{0.0,{";
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		file << (i>0?",":"") << radius_[i];
	}	
	file << "},{";
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		file << (i>0?",":"") << AngleSum(vertexId);
	}
	file << "}}\n";
#endif

	disk_vertices_.clear();
	std::vector<bool> indisk(triangulation_->NumberOfVertices(),false);
	for(int i=0,endi=disk_triangles_.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			Vertex * v = disk_triangles_[i]->getEdge(j)->getOpposite();
			if( onBoundary[v->getId()] )
			{
				radius_[v->getId()] = (use_one_minus_radius_? 1.0-almostZero : almostZero);
			}else if( !indisk[v->getId()] )
			{
				disk_vertices_.push_back(v);
				indisk[v->getId()] = true;
				if( reset_radius_ || radius_[v->getId()] <= 1.2*almostZero || radius_[v->getId()] >= 1.0 )
				{
					radius_[v->getId()] = (use_one_minus_radius_? 1.0-startRadius:startRadius);
				}
			}
		}
	}

	std::vector<double> radius2(triangulation_->NumberOfVertices());

	properties::DegreeList(triangulation_,degree_);

	//// TEST /////
	std::vector<Vertex *> order;
	std::queue<Vertex *> q;
	std::vector<bool> visited(triangulation_->NumberOfVertices(),false);
	for(std::list<const Edge*>::const_iterator it = boundary.begin();it != boundary.end();it++)
	{
		q.push((*it)->getPrevious()->getOpposite());
		visited[(*it)->getPrevious()->getOpposite()->getId()] = true;
	}
	while( !q.empty() )
	{
		Vertex * v = q.front();
		q.pop();
		
		Edge * edge = v->getParent()->getPrevious();
		do {
			Vertex * nbr = edge->getPrevious()->getOpposite();
			if( indisk[nbr->getId()] && !visited[nbr->getId()] )
			{
				order.push_back(nbr);
				visited[nbr->getId()] = true;
				q.push(nbr);
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != v->getParent()->getPrevious());
	}
	//////////////
	

	double c = epsilon_ + 1;
	double lambda = -1.0;
	bool flag = false;
	steps_ = 0;
	while(c > epsilon_ )
	{
		double c0 = c;
		c = 0.0;
		double lambda0 = lambda;
		bool flag0 = flag;
		
		/// TEST
		//std::copy(radius_.begin(),radius_.end(),radius2.begin());

		//for(int i=0,endi=disk_vertices_.size();i<endi;i++)
		for( std::vector<Vertex *>::iterator it = order.begin();it!=order.end();it++)
		{
			// See Colling, Stephenson (2003): p. 244
			//int vertexId = disk_vertices_[i]->getId();
			int vertexId = (*it)->getId();
			double theta = AngleSum(vertexId);
			double sintheta = std::sin(theta/(2*degree_[vertexId]));
			double sintarget = std::sin(PI/degree_[vertexId]);
			if( use_one_minus_radius_ )
			{
				double oneminvhat = sintheta * radius_[vertexId] / (std::sqrt(1.0-radius_[vertexId]) - sintheta * (1.0-radius_[vertexId]));
				if( oneminvhat > 1.0 )
				{
					oneminvhat = 1.0;
				}
				radius2[vertexId] = oneminvhat * (-oneminvhat + 2.0 * (oneminvhat-1.0) * sintarget * sintarget + std::sqrt(oneminvhat * oneminvhat + 4.0 * (1.0-oneminvhat) * sintarget * sintarget)) 
					/ ( 2.0 * (1.0-oneminvhat) * (1.0-oneminvhat) * sintarget * sintarget );
				//radius_[vertexId] = oneminvhat * (-oneminvhat + 2.0 * (oneminvhat-1.0) * sintarget * sintarget + std::sqrt(oneminvhat * oneminvhat + 4.0 * (1.0-oneminvhat) * sintarget * sintarget)) 
				//	/ ( 2.0 * (1.0-oneminvhat) * (1.0-oneminvhat) * sintarget * sintarget );			
			} else
			{
				double vhat = (sintheta - std::sqrt(radius_[vertexId]))/(sintheta*radius_[vertexId] - std::sqrt(radius_[vertexId]));
				if( vhat < 0.0 )
				{
					vhat = 0.0;
				}
				double t = 2.0 * sintarget / (std::sqrt((1.0-vhat)*(1.0-vhat)+4.0*sintarget*sintarget*vhat)+1.0-vhat);
				radius2[vertexId] = t*t;
				//radius_[vertexId] = t*t;
			}
			c += (theta - 2.0 * PI)*(theta - 2.0 * PI);
		}
		//std::swap(radius_,radius2);
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

#ifdef LOG_CIRCLE_PACKING
	file << "{0.0,{";
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		file << (i>0?",":"") << radius_[vertexId];
	}	
	file << "},{";
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		file << (i>0?",":"") << radius2[vertexId];
	}	
	file << "},{";
	for(int i=0,endi=disk_vertices_.size();i<endi;i++)
	{
		int vertexId = disk_vertices_[i]->getId();
		file << (i>0?",":"") << AngleSum(vertexId);
	}
	file << "}}\n";
#endif
		steps_++;
		if( steps_ > max_iterations_ )
		{
			return false;
		}
	}
#ifdef LOG_CIRCLE_PACKING
file.close();
radius_ = radius3;
goto lab;

#endif
	return true;
}

bool DiskCirclePacking::FindEmbedding(const std::list<const Edge*> & boundary,  const Edge * centerEdge)
{
	boundary_ = boundary;
	center_edge_ = centerEdge;

	disk_triangles_.clear();
	babyuniversedetector_.EnclosedTriangles(boundary,disk_triangles_,true);

/*	if( !FindDiskRadii( boundary ) )
	{
		return false;
	}
	if( !DiskLayout(boundary,centerEdge) )
	{
		return false;
	}
	*/
	FindFlatEmbedding();

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

void DiskCirclePacking::getHyperbolicCoordinates(std::vector<std::pair<Vertex*,Vector2D> > & coor) const
{
	coor.clear();
	for(int i=0,endi=circle_order_.size();i<endi;i++)
	{
		coor.push_back(std::pair<Vertex*,Vector2D >(triangulation_->getVertex(circle_order_[i]),hyp_coordinate_[circle_order_[i]]));
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

void DiskCirclePacking::getBoundaryPositions(std::vector<double> & angles)
{
	angles.clear();
	for(std::list<const Edge*>::const_iterator it = boundary_.begin(); it!=boundary_.end();it++)
	{
		angles.push_back( VectorAngle(hyp_coordinate_[(*it)->getNext()->getOpposite()->getId()]) );
	}
}

double DiskCirclePacking::getCenterRadius() const
{
	return euclidean_radius_[center_edge_->getNext()->getOpposite()->getId()];
}

bool DiskCirclePacking::FindFlatEmbedding()
{
	std::vector<double> radius;
	std::vector<Vector2D> coordinate;
	std::vector<Vertex*> diskVertices;
	std::vector<int> diskPosition(triangulation_->NumberOfTriangles(),-1);
	int bLength = boundary_.size();

	for(std::list<const Edge*>::iterator it=boundary_.begin();it!=boundary_.end();it++)
	{
		Vertex * v = (*it)->getNext()->getOpposite();
		diskPosition[v->getId()] = diskVertices.size();
		diskVertices.push_back(v);
	}
	for(int i=0,endi=disk_triangles_.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			Vertex * v = disk_triangles_[i]->getEdge(j)->getOpposite();
			if( diskPosition[v->getId()] == -1 )
			{
				diskPosition[v->getId()] = diskVertices.size();
				diskVertices.push_back(v);
			}
		}
	}

	radius.resize(diskVertices.size(),1.0);
	coordinate.resize(diskVertices.size(),MakeVector2D(0,0));
	int nIntVertices = static_cast<int>(diskVertices.size())-bLength;

	///
	/*std::ofstream file("log.txt");
	file << std::fixed << "{";
	for(int i=0,endi=diskVertices.size();i<endi;i++)
	{
		Vertex * v = diskVertices[i];
		Edge * edge = v->getParent()->getPrevious();
		file << (i>0?",":"") << "{";
		bool first=true;
		do
		{
			int nbrPos = diskPosition[edge->getPrevious()->getOpposite()->getId()];
			if( nbrPos >= 0 )
			{
				file << (first?"":",") << nbrPos;
				first=false;
			}
			edge = edge->getAdjacent()->getNext();
		} while( edge != v->getParent()->getPrevious() );
		file << "}";
	}
	file << "}\n";
	*////

	int steps = 0;
	double totalerror = 1000.0;
	while( totalerror > 0.0001 )
	{
		steps++;
		scaleRadii(radius,bLength);
		LayoutBoundary(radius,coordinate,bLength);

		lapl_rules_.clear();
		boundary_rules_.clear();
		for(int i=bLength,endi=diskVertices.size();i<endi;i++)
		{
			Vertex * v = diskVertices[i];
			double totalWeight = 0.0;

			Edge * edge = v->getParent()->getPrevious();
			do
			{
				double weight = 0.5/std::tan(Angle(edge,diskPosition,radius)) + 0.5/std::tan(Angle(edge->getAdjacent(),diskPosition,radius));
				totalWeight += weight;

				int nbrPos = diskPosition[edge->getPrevious()->getOpposite()->getId()];
				if( nbrPos < bLength )
				{
					boundary_rules_.push_back(Eigen::Triplet<double>(i-bLength,nbrPos,weight));
				}else
				{
					lapl_rules_.push_back(Eigen::Triplet<double>(i-bLength,nbrPos-bLength,weight));
				}
				edge = edge->getAdjacent()->getNext();
			} while( edge != v->getParent()->getPrevious() );
			lapl_rules_.push_back(Eigen::Triplet<double>(i-bLength,i-bLength,-totalWeight));
		}

		Eigen::SparseMatrix<double> laplMatrix(nIntVertices,nIntVertices);
		laplMatrix.setFromTriplets(lapl_rules_.begin(),lapl_rules_.end());

		Eigen::SparseMatrix<double> bMatrix(nIntVertices,bLength);
		bMatrix.setFromTriplets(boundary_rules_.begin(),boundary_rules_.end());

		Eigen::MatrixXd solveVec(bLength,2);
		for(int i=0;i<bLength;i++)
		{
			solveVec(i,0) = -coordinate[i][0];
			solveVec(i,1) = -coordinate[i][1];
		}

		Eigen::MatrixXd solveVec2 = bMatrix * solveVec;

		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > laplChol(laplMatrix);
		if( laplChol.info() != Eigen::Success )
		{
			return false;
		}

		Eigen::MatrixXd newCoor = laplChol.solve(solveVec2);
		if( laplChol.info() != Eigen::Success )
		{
			return false;
		}

		for(int i=0;i<nIntVertices;i++)
		{
			coordinate[i+bLength][0] = newCoor(i,0);
			coordinate[i+bLength][1] = newCoor(i,1);
		}

		//PrintToStream(file,radius.begin(),radius.end());
		//file << "\n";

		totalerror = 0.0;
		double maxradius = 0.0;
		for(int i=0,endi=diskVertices.size();i<endi;i++)
		{
			Vertex * v = diskVertices[i];

			double totalangle = 0.0;
			double totalarea = 0.0;
			Edge * edge = v->getParent();
			do
			{
				int v1Pos = diskPosition[edge->getNext()->getOpposite()->getId()];
				int v2Pos = diskPosition[edge->getPrevious()->getOpposite()->getId()];
				if( v1Pos >= 0 && v2Pos >= 0 )
				{
					double len0 = Norm2D(SubtractVectors2D(coordinate[v1Pos],coordinate[v2Pos]));
					double len1 = Norm2D(SubtractVectors2D(coordinate[v2Pos],coordinate[i]));
					double len2 = Norm2D(SubtractVectors2D(coordinate[v1Pos],coordinate[i]));
					double angle = std::acos((len1*len1+len2*len2-len0*len0)/(2.0*len1*len2));
					totalangle += angle;
					totalarea += 0.5 * angle * 0.25 * (len1 + len2 - len0) * (len1 + len2 - len0);
				}
				edge = edge->getNext()->getAdjacent()->getNext();
			} while( edge != v->getParent() );
			
			BOOST_ASSERT( i < bLength || std::fabs(totalangle - 2.0*PI) < 0.001 );

			double oldradius = radius[i];

			//radius[i] = std::sqrt(2.0*totalarea / totalangle );
			radius[i] = 0.2*radius[i] + 0.8*std::sqrt(2.0*totalarea / totalangle );
			BOOST_ASSERT( radius[i] > 0.0 );

			if( radius[i] > maxradius )
			{
				maxradius = radius[i];
			}

			totalerror += std::fabs(oldradius-radius[i]);
		}

		////
		/*PrintToStream(file,radius.begin(),radius.end());
		file << "\n";
		PrintToStream2D(file,coordinate.begin(),coordinate.end());
		file << "\n";
		*////

		//std::cout << totalerror << " " << maxradius << "\n";
	}
	//file.close();

	return true;
}

void DiskCirclePacking::LayoutBoundary(std::vector<double> & radius, std::vector<Vector2D> & coordinate, int bLength)
{
	double angle = 0.0;
	for(int i=0;i<bLength;i++)
	{
		coordinate[i][0] = (1.0 - radius[i]) * std::cos(angle);
		coordinate[i][1] = (1.0 - radius[i]) * std::sin(angle);

		angle += std::acos(1.0 - 2.0 * radius[i] * radius[(i+1)%bLength] / (1.0 - radius[i])/(1.0-radius[(i+1)%bLength]));
	}
}	

void DiskCirclePacking::scaleRadii(std::vector<double> & radius, int bLength)
{
	double r = 0.0;
	for(int i=0;i<bLength;i++)
	{
		if( radius[i] > r )
		{
			r = radius[i];
		}
	}
	r *= 3.0;

	double maxError = 1.0e-9;
	double angle;
	do
	{
		angle = boundaryAngle(radius,bLength,r);
		double derivative = boundaryAngleDerivative(radius,bLength,r);

		r = r / (1.0 + (angle - 2.0*PI) / (r * derivative));
	} while( std::abs(angle - 2.0*PI) > maxError );

	double scale = 1.0/r;

	BOOST_ASSERT( scale > 0.0 );

	for(int i=0,endi=radius.size();i<endi;i++)
	{
		radius[i] *= scale;
	}
}

double DiskCirclePacking::boundaryAngle(std::vector<double> & radius, int bLength, double r)
{
	double angle = 0.0;
	for(int i=0;i<bLength;i++)
	{
		int nexti = (i+1)%bLength;
		angle += std::acos(1.0 - 2.0 * radius[i] * radius[nexti] / (r - radius[i])/(r-radius[nexti]));
	}
	return angle;
}

double DiskCirclePacking::boundaryAngleDerivative(std::vector<double> & radius, int bLength, double r)
{
	double derivative = 0.0;
	for(int i=0;i<bLength;i++)
	{
		int nexti = (i+1)%bLength;
		derivative += radius[i] * radius[nexti] * (radius[i] + radius[nexti] - 2.0*r) / ( (r - radius[i])*(r-radius[nexti])*std::sqrt(r*radius[i]*radius[nexti]*(r-radius[i]-radius[nexti])));
	}
	return derivative;
}