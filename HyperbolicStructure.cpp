#include "HyperbolicStructure.h"


HyperbolicStructure::HyperbolicStructure(const Triangulation * const triangulation, ThetaModel * const thetamodel, CirclePattern * const circlepattern, DualCohomologyBasis * const dualcohomologybasis)
	: triangulation_(triangulation),
	  thetamodel_(thetamodel),
	  circlepattern_(circlepattern),
	  dualcohomologybasis_(dualcohomologybasis),
	  maxtheta_(2.05*PI),
	  sqlength_(0.0,25.0,250),
	  legs_(0,1),
	  angle_surplus_(0.0,0.05*PI,100),
	  leg_shear_(0.0,12.5,125),
	  log_sqlength_shear_ratio_(-2.0,10.0,120),
	  cusp_size_(-2.0,16.0,180),
	  num_log_(100),
	  use_log_(true),
	  log_fresh_(true),
	  maxhyplength_(0.5)
{
	turn_right_mat_(0,0) = 0;
	turn_right_mat_(0,1) = -1;
	turn_right_mat_(1,0) = 1;
	turn_right_mat_(1,1) = 1;
}

HyperbolicStructure::~HyperbolicStructure()
{
}

void HyperbolicStructure::setMaxTheta(double max)
{
	maxtheta_=max;
}

void HyperbolicStructure::setMaxLength(double length)
{
	maxhyplength_=length;
}

void HyperbolicStructure::Measure()
{
	if( circlepattern_->IsUpToDate() || circlepattern_->MakeUpToDate() )
	{
		/*int maxThetaInUnits = static_cast<int>(maxtheta_ * thetamodel_->getPiInUnits()/PI);
		std::list<std::pair<int,std::list<const Edge*> > > paths;
		thetamodel_->FindAllShortCurves(maxThetaInUnits,paths);

		for(std::list<std::pair<int,std::list<const Edge*> > >::iterator pathIt = paths.begin(); pathIt!= paths.end();pathIt++)
		{
			double length = HyperbolicLength(pathIt->second);
			sqlength_.Insert(length*length);

			int legs = NumberOfLegs(pathIt->second);
			legs_.Insert(legs);

			angle_surplus_.Insert( (static_cast<double>(pathIt->first)/thetamodel_->getPiInUnits() - 2.0)*PI );
			
			if( use_log_ )
			{
				log_ << (log_fresh_?"":",") << "{" << length << ", " << (static_cast<double>(pathIt->first)/thetamodel_->getPiInUnits() - 2.0)*PI << ",{";
				log_fresh_ = false;
			}
			MeasureLegShears(pathIt->second, length);


			if( use_log_ )
			{
				log_ << "}}";
				num_log_--;
				if( num_log_ == 0 )
				{
					use_log_ = false;
				}
			}
		}*/
		FindShortHyperbolicCurves();
		std::pair<double,double> mod = circlepattern_->CalculateModuli();
		Complex modulus(mod.first,mod.second);

		if( !paths_.empty() )
		{
			if( use_log_ )
			{
				log_ << (log_fresh_?"":",") << "{";
				log_fresh_ = false;
			}
			for(std::list<std::pair<double,std::list<const Edge*> > >::iterator pathIt = paths_.begin(); pathIt!= paths_.end();pathIt++)
			{
				double length = pathIt->first;
				std::cout << length << "\n";
				sqlength_.Insert(length*length);

				int legs = NumberOfLegs(pathIt->second);
				legs_.Insert(legs);

				MeasureLegShears(pathIt->second, length);

				if( use_log_ )
				{
					log_ << (pathIt == paths_.begin()?"":",") << "{" << std::fixed << length << "," << legs << "," ;
					for(std::list<const Edge*>::iterator it = pathIt->second.begin();it!=pathIt->second.end();it++)
					{
						log_ << (it == pathIt->second.begin()?"":",") << "{" << (*it)->getParent()->getId() << "," << (*it)->getId() << "}";
						log_ << ",{" << (*it)->getAdjacent()->getParent()->getId() << "," << (*it)->getAdjacent()->getId() << "}";
					}
					log_ << "}";
				}
			}
			if( use_log_ )
			{
				log_ << "}";
				num_log_--;
				if( num_log_ == 0 )
				{
					use_log_ = false;
				}
			}

		}
		MeasureCuspSize(modulus);
	}
}

void HyperbolicStructure::MeasureCuspSize(Complex modulus)
{
	for(int i=0,endi=triangulation_->NumberOfVertices();i<endi;i++)
	{
		double cusp = CuspSize(triangulation_->getVertex(i),modulus);
		cusp_size_.Insert(cusp);
	}
}

double HyperbolicStructure::CuspSize(const Vertex * vertex, Complex modulus) const
{
	const Edge * edge = vertex->getParent()->getNext();

	std::list<const Edge*> path;
	const Edge * startEdge = edge;
	do {
		path.push_back(edge);
		edge = edge->getAdjacent()->getPrevious();
	} while( edge != startEdge );

	SL2Mat mat = Holonomy(path,true);
	
	boost::array<double,3> loglengths;
	for(int i=0;i<3;i++)
	{
		loglengths[i] = 0.5 * std::log(NormSquaredTransformedByModulus(circlepattern_->getForm(edge), std::make_pair(modulus.real(),modulus.imag()) ));
		edge = edge->getNext();
	}

	double cusp = - loglengths[0] - loglengths[1] + loglengths[2] + std::log(-mat(1,0));

	return cusp;
}

int HyperbolicStructure::FindShortHyperbolicCurves()
{
	paths_.clear();
	visited_.ResetAndResize(triangulation_->NumberOfTriangles());
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;i++)
	{
		const Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			const Edge * edge = triangle->getEdge(j);
			std::list<const Edge *> path;
			path.push_back(edge);
			PathInfo info;
			info.holonomy = Eigen::Matrix2d::Identity();
			info.has_turned_left = false;
			info.has_turned_right = false;
			info.omega[0] = 0;
			info.omega[1] = 0;
			visited_.Reset();
			SearchPath(path,info,std::cosh(maxhyplength_));
		}
	}
	return paths_.size();
}

void HyperbolicStructure::SearchPath(std::list<const Edge *> & path, PathInfo & info, double maxcoshlength)
{
	if( path.back()->getAdjacent()->getParent()->getId() <= path.front()->getParent()->getId() &&
		( path.back()->getAdjacent()->getParent() != path.front()->getParent() || 
		  path.back()->getAdjacent()->getId() < path.front()->getId() ) )
	{
		// to make sure that each path is recorder only once, we require the starting edge to
		// have ``smallest'' ID.
		return;
	}
	double tmp = 2.0 * info.holonomy(0,0) * info.holonomy(1,1) - 1.0;
	if( 2.0 * info.holonomy(0,0) * info.holonomy(1,1) - 1.0 > maxcoshlength )
	{
		// the lower bound on the hyperbolic length already exceeds the maximum length
		return;
	}

	if( path.size() > 1 && path.back() == path.front() && FormIsZero(info.omega) )
	{
		if( info.has_turned_left && info.has_turned_right )
		{
			double coshlength = 0.5 * info.holonomy.trace() * info.holonomy.trace() - 1.0;
			if( coshlength < maxcoshlength )
			{
				double length = std::log( coshlength + std::sqrt(coshlength *coshlength-1.0));
				// have found a short path
				paths_.push_back( std::pair<double,std::list<const Edge *> >(length,path));
				paths_.back().second.pop_back();
			}
		} 
		return;
	}

	visited_.Set(path.back()->getAdjacent()->getParent()->getId());

	SL2Mat shearMat = ShearMatrix(path.back());

	PathInfo info2(info);	
	info2.holonomy = shearMat * info.holonomy;
	IntForm2D form = dualcohomologybasis_->getOmega(path.back());
	info2.omega = AddForms(info.omega,form);

	for(int i=0;i<2;i++)
	{
		const Edge * edge;
		info2.holonomy = turn_right_mat_ * info2.holonomy;
		if( i==0 )
		{
			info2.has_turned_right = true;
			edge = path.back()->getAdjacent()->getNext();
		} else
		{
			info2.has_turned_right = info.has_turned_right;
			info2.has_turned_left = true;
			edge = path.back()->getAdjacent()->getPrevious();
		}

		if( !visited_.isSet(edge->getAdjacent()->getParent()->getId()) || (edge == path.front() && FormIsZero(info2.omega)) )
		{
			path.push_back(edge);
			SearchPath(path,info2,maxcoshlength);
			path.pop_back();
		}
	}
	visited_.Set(path.back()->getAdjacent()->getParent()->getId(),false);
}

HyperbolicStructure::SL2Mat HyperbolicStructure::ShearMatrix(const Edge * edge) const
{
	double shear = circlepattern_->getShear(edge);
	double exphalfshear = std::exp(0.5*shear);
	SL2Mat shearMat = Eigen::Matrix2d::Zero();
	shearMat(0,1) = exphalfshear;
	shearMat(1,0) = -1.0/exphalfshear;
	return shearMat;
}

HyperbolicStructure::SL2Mat HyperbolicStructure::Holonomy(const std::list<const Edge*> & path, bool startWithTurnRight) const
{
	SL2Mat mat = Eigen::Matrix2d::Identity();
	bool lastleft = false;
	for(std::list<const Edge*>::const_iterator edgeIt = path.begin();edgeIt != path.end();edgeIt++)
	{
		SL2Mat shearMat = ShearMatrix(*edgeIt);
		mat = shearMat * mat;
		mat = turn_right_mat_ * mat;

		std::list<const Edge*>::const_iterator nextIt = boost::next(edgeIt);
		if( nextIt == path.end() )
		{
			nextIt = path.begin();
		}
		lastleft = false;
		if( (*nextIt) == (*edgeIt)->getAdjacent()->getPrevious() )
		{
			mat = turn_right_mat_ * mat;
			lastleft = true;
		}
	}
	if( startWithTurnRight )
	{
		if( lastleft )
		{
			mat = turn_right_mat_ * mat * turn_right_mat_ * turn_right_mat_;
		} else
		{
			mat = turn_right_mat_ * turn_right_mat_ * mat * turn_right_mat_;
		}
	}
	return mat;
}

double HyperbolicStructure::HyperbolicLength(const std::list<const Edge*> & path) const
{
	SL2Mat mat = Holonomy(path);
	double trace = mat.trace();
	return 2.0*std::log(trace/2.0 + std::sqrt(trace*trace/4.0-1.0));
}

void HyperbolicStructure::GeodesicCoordinates(const std::list<const Edge*> & path, std::list<boost::array<Vector2D,2> > & coor, Complex modulus) const
{
	std::list<const Edge*> path2 = path;
	coor.clear();
	for(int i=0,endi = path.size();i<endi;i++)
	{
		path2.push_back(path2.front());
		path2.pop_front();
		SL2Mat mat = Holonomy(path2,true);
		mat = turn_right_mat_ * mat * turn_right_mat_ * turn_right_mat_;
		// mat is now of the form mat = R.X(i).R^?.X(i-1)....X(0).R^?.X(n).R^?....X(i+1).R^?.R.R

		boost::array<double,2> x = FixedPoints(mat);
		if( x[1] > -1.0 && x[1] < 0.0 )
		{
			std::swap(x[0],x[1]);
		}
		BOOST_ASSERT( x[0] > -1.0 && x[0] < 0.0 );
		BOOST_ASSERT( x[1] < -1.0 || x[1] > 0.0 );

		// determine the intersections of the geodesic with the ideal triangle (-1,0,infinity)
		boost::array<Vector2D,2> intersections;
		intersections[0][0] = x[0] * x[1] / (1.0 + x[0] + x[1]);
		intersections[0][1] = std::sqrt( -x[0] * (1.0 + x[0]) * x[1] * (1.0 + x[1]) ) / std::abs( 1.0 + x[0] +x[1] );
		if( x[1] < -1.0 )
		{
			intersections[1][0] = -1.0;
			intersections[1][1] = std::sqrt(-(1.0 + x[0]) * (1.0 + x[1]) );
		} else
		{
			intersections[1][0] = 0.0;
			intersections[1][1] = std::sqrt( - x[0] * x[1] );
		}

		const Edge * base = path2.back()->getAdjacent();
		std::pair<double,double> moduluspair(modulus.real(),modulus.imag());
		boost::array<Complex,3> tricoor;
		tricoor[0] = ToComplex(TransformByModulus(circlepattern_->getCoordinate(base->getNext()->getOpposite()),moduluspair));
		tricoor[1] = tricoor[0] + ToComplex(TransformByModulus(circlepattern_->getForm(base),moduluspair));
		tricoor[2] = tricoor[1] + ToComplex(TransformByModulus(circlepattern_->getForm(base->getNext()),moduluspair));

		boost::array<Vector2D,2> endpoints;
		for(int j=0;j<2;j++)
		{
			Complex p = PoincareDiskToKleinDisk(MapToTriangle(ToComplex(intersections[j]),tricoor),tricoor);
			endpoints[j][0] = p.real() - modulus.real()*p.imag()/modulus.imag();
			endpoints[j][1] = p.imag()/modulus.imag();
		}
		coor.push_back(endpoints);
	}
}

void HyperbolicStructure::RetrieveCurves(std::list<Curve> & curves) const
{
	curves.clear();
	if( paths_.empty() )
	{
		return;
	}
	std::pair<double,double> moduluspair = circlepattern_->CalculateModuli();
	Complex modulus(moduluspair.first,moduluspair.second);
	for(std::list<std::pair<double,std::list<const Edge*> > >::const_iterator pathIt = paths_.begin(); pathIt!= paths_.end();pathIt++)
	{
		curves.push_back(Curve());
		curves.back().length = HyperbolicLength(pathIt->second);
		GeodesicCoordinates(pathIt->second,curves.back().segments,modulus);
	}
}

boost::array<double,2> HyperbolicStructure::FixedPoints(const SL2Mat & mat)
{
	double r = std::sqrt(mat.trace()*mat.trace()-4.0)/(2.0 * std::abs(mat(1,0)) );
	double c = (mat(0,0)-mat(1,1))/(2.0*mat(1,0));
	boost::array<double,2> x = {c-r,c+r};
	return x;
}

HyperbolicStructure::Complex HyperbolicStructure::ToComplex(Vector2D v)
{
	return Complex(v[0],v[1]);
}

Vector2D HyperbolicStructure::ToVector2D(Complex x)
{
	return MakeVector2D(x.real(),x.imag());
}

HyperbolicStructure::Complex HyperbolicStructure::MapToTriangle(Complex x, const boost::array<Complex,3> & v)
{
	// mobius transformation that maps (-1,0,infinity) -> (v[0],v[1],v[2])
	return ( x * v[2] * (v[0]-v[1]) - v[1] * (v[2] - v[0]) ) / (x * (v[0] - v[1]) - (v[2] - v[0]) );
}

HyperbolicStructure::Complex HyperbolicStructure::PoincareDiskToKleinDisk(Complex x, const boost::array<Complex,3> & v)
{
	// map the Poincare disk corresponding to the circumscribed circle of (v[0],v[1],v[2]) to the corresponding Klein model
	boost::array<Vector2D,3> vec = {ToVector2D(v[0]),ToVector2D(v[1]),ToVector2D(v[2])};
	Complex c = ToComplex(CenterOfCircle(vec));
	return 2.0 / (1.0 + std::norm(x-c)/std::norm(v[0]-c)) * (x-c) + c;
}

std::string HyperbolicStructure::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "hyperbolicstructure -> {maxtheta -> " << maxtheta_;
	stream << ", anglesurplus -> ";
	angle_surplus_.PrintTo(stream);
	stream << ", sqlength -> ";
	sqlength_.PrintTo(stream);
	stream << ", legshear -> ";
	leg_shear_.PrintTo(stream);
	stream << ", legs -> ";
	legs_.PrintTo(stream);
	stream << ", logsqlengthshearratio -> ";
	log_sqlength_shear_ratio_.PrintTo(stream);
	stream << ", cuspsize -> ";
	cusp_size_.PrintTo(stream);
	stream << ", log -> {" << log_.str() << "}";
	stream << "}";
	return stream.str();
}

int HyperbolicStructure::NumberOfLegs(const std::list<const Edge*> & path) const
{
	bool right = (path.back()->getAdjacent()->getNext() == path.front());
	int legs = 0;
	for( std::list<const Edge*>::const_iterator edgeIt = path.begin();edgeIt!=path.end();edgeIt++)
	{
		std::list<const Edge*>::const_iterator nextIt = boost::next(edgeIt);
		if( nextIt == path.end() )
		{
			nextIt = path.begin();
		}
		bool nextRight = ((*edgeIt)->getAdjacent()->getNext() == (*nextIt));
		if( right != nextRight )
		{
			legs++;
			right = nextRight;
		}
	}
	BOOST_ASSERT( legs%2 == 0 );
	return legs/2;
}

void HyperbolicStructure::MeasureLegShears(const std::list<const Edge*> & path, double length)
{
	bool right = (path.back()->getAdjacent()->getNext() == path.front());
	bool first = true;
	for( std::list<const Edge*>::const_iterator edgeIt = path.begin();edgeIt!=path.end();edgeIt++)
	{
		std::list<const Edge*>::const_iterator nextIt = boost::next(edgeIt);
		if( nextIt == path.end() )
		{
			nextIt = path.begin();
		}
		bool nextRight = ((*edgeIt)->getAdjacent()->getNext() == (*nextIt));
		if( right != nextRight )
		{
			double shear = std::fabs(circlepattern_->getShear(*edgeIt));
			leg_shear_.Insert(shear);
			log_sqlength_shear_ratio_.Insert( shear + 2.0*std::log(length/2.0) );
			/*if( use_log_ )
			{
				log_ << std::fixed << (first?"":",") << shear;
			}*/
			first = false;
			right = nextRight;
		}
	}
}