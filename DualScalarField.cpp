#include "DualScalarField.h"
#include "Triangle.h"

DualScalarField::DualScalarField(Triangulation * const triangulation) 
	: triangulation_(triangulation), massive_(false), mass_squared_(0.0), babyuniversedetector_(triangulation)
{
	field_.resize(triangulation->NumberOfTriangles(),0.0);	
}

DualScalarField::DualScalarField(Triangulation * const triangulation, double massSquared) 
	: triangulation_(triangulation), massive_(true), mass_squared_(massSquared), babyuniversedetector_(triangulation)
{
	field_.resize(triangulation->NumberOfTriangles(),0.0);
	if( std::fabs(mass_squared_) < 1.0e-12 )
	{
		massive_ = false;
		mass_squared_ = 0.0;
	}
}

DualScalarField::~DualScalarField(void)
{
}

void DualScalarField::Initialize()
{
	field_.resize(triangulation_->NumberOfTriangles());	
	std::fill(field_.begin(),field_.end(),0.0);
	DoSweep();
}

double DualScalarField::getField(Triangle * triangle) const
{
	return getField(triangle->getId());
}

double DualScalarField::getField(int triangle) const
{
	return field_[triangle];
}

void DualScalarField::setField(Triangle * triangle, const double & field)
{
	setField(triangle->getId(),field);
}

void DualScalarField::setField(int triangle, const double & field)
{
	field_[triangle] = field;
}

void DualScalarField::addToField(Triangle * triangle, const double & field)
{
	addToField(triangle->getId(),field);
}

void DualScalarField::addToField(int triangle, const double & field)
{
	field_[triangle] += field;
}

void DualScalarField::DoSweep()
{
	int sweepsize = ( UsesCustomSweepSize() ? CustomSweepSize() : triangulation_->NumberOfTriangles() );
	for(int i=0;i<sweepsize;i++)
	{
		if( triangulation_->CalculateGenus() !=0 || triangulation_->SucceedWithProbability(0.7) )
		{
			DoMove();
		} else
		{
			TryBabyUniverseMove();
		}
	}
}

void DualScalarField::DoMove()
{
	Triangle * triangle = triangulation_->getRandomTriangle();

	double mean = 0.0;
	double neighbours = 0.0;
	for(int i=0;i<3;i++)
	{
		Triangle * nbr = triangle->getEdge(i)->getAdjacent()->getParent();
		if( nbr != triangle )
		{
			mean += getField(nbr);
			neighbours += 1.0;
		}
	}
	mean = mean / ( neighbours + mass_squared_ );
	double sigma = 1.0 / std::sqrt(2.0 * (neighbours + mass_squared_ ));
	
	setField(triangle, triangulation_->RandomNormal(mean,sigma));
}

void DualScalarField::TryBabyUniverseMove()
{
	Edge * edge = triangulation_->getRandomEdge();
	if( edge->getNext()->getOpposite() == edge->getPrevious()->getOpposite() )
	{
		std::list<const Edge*> boundary( 1, edge );
		std::pair<int,bool> volume = babyuniversedetector_.VolumeEnclosed(boundary);
		std::vector<Triangle *> triangles;
		triangles.reserve( volume.first );
		babyuniversedetector_.EnclosedTriangles( boundary, triangles, volume.second );

		double sigma = 1.0 / std::sqrt( 2.0 * ( 1.0 + mass_squared_ * volume.first ) );
		double mean = ( volume.second ? 1.0 : -1.0 ) * (getField(edge->getAdjacent()->getParent()) - getField(edge->getParent()));
		if( massive_ )
		{
			for(std::vector<Triangle *>::iterator it=triangles.begin();it != triangles.end();++it)
			{
				mean -= mass_squared_ * getField(*it);
			}
			mean /= 1.0 + mass_squared_ * volume.first;
		}
		double fieldchange = triangulation_->RandomNormal( mean, sigma );
		BOOST_ASSERT( fieldchange < 50.0 );
		for(std::vector<Triangle *>::iterator it=triangles.begin();it != triangles.end();++it)
		{
			addToField(*it,fieldchange);
		}
	}
}

double DualScalarField::BoltzmannChangeUnderFlipMove(const Edge * const edge) const
{
	double daction = 2.0 * ( getField( edge->getParent()) - getField( edge->getAdjacent()->getParent() ) )
		* ( getField( edge->getPrevious()->getAdjacent()->getParent() ) - getField( edge->getAdjacent()->getPrevious()->getAdjacent()->getParent() ) );

	return std::exp(-daction);
}

double DualScalarField::BoltzmannChangeUnderGeneralMove(const std::vector<boost::array<Triangle *,2> > & toBeDeleted, const std::vector<boost::array<Triangle *,2> > & toBeAdded ) const
{
	double daction = 0.0;
	for(std::vector<boost::array<Triangle *,2> >::const_iterator it = toBeDeleted.begin(); it!=toBeDeleted.end();++it)
	{
		daction -= math::square( getField( (*it)[0] ) - getField( (*it)[1] ) );
	}
	for(std::vector<boost::array<Triangle *,2> >::const_iterator it = toBeAdded.begin(); it!=toBeAdded.end();++it)
	{
		daction += math::square( getField( (*it)[0] ) - getField( (*it)[1] ) );
	}
	return std::exp(-daction);
}

double DualScalarField::CentralCharge() const {
	return 1.0;
}

double DualScalarField::AverageField() const {
	double total=0.0;
	for(int i=0,end=field_.size();i<end;i++)
	{
		total += field_[i];
	}
	return total/field_.size();
}

std::string DualScalarField::ConfigurationData() const {
	std::ostringstream stream;
	stream << std::fixed << "{type -> \"dualscalarfield\", centralcharge -> 1, massive -> ";
	stream << (massive_? "true" : "false" ) << ", masssquared -> " << mass_squared_ << "}";
	return stream.str();
}


