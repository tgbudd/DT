#include "DualScalarField.h"
#include "Triangle.h"

DualScalarField::DualScalarField(Triangulation * const triangulation) 
	: triangulation_(triangulation), massive_(false), mass_squared_(0.0)
{
	field_.resize(triangulation->NumberOfTriangles(),0.0);	
}

DualScalarField::DualScalarField(Triangulation * const triangulation, double massSquared) 
	: triangulation_(triangulation), massive_(true), mass_squared_(massSquared)
{
	field_.resize(triangulation->NumberOfTriangles(),0.0);	
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

void DualScalarField::DoSweep()
{
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		DoMove();
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

std::string DualScalarField::ConfigurationData() const {
	std::ostringstream stream;
	stream << std::fixed << "{type -> \"dualscalarfield\", centralcharge -> 1, massive -> ";
	stream << (massive_? "true" : "false" ) << ", masssquared -> " << mass_squared_ << "}";
	return stream.str();
}