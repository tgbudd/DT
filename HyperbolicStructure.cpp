#include "HyperbolicStructure.h"
#include <Eigen/dense>

HyperbolicStructure::HyperbolicStructure(const Triangulation * const triangulation, ThetaModel * const thetamodel, CirclePattern * const circlepattern)
	: triangulation_(triangulation),
	  thetamodel_(thetamodel),
	  circlepattern_(circlepattern),
	  maxtheta_(2.05*PI),
	  sqlength_max_(25.0),
	  sqlength_bins_(250),
	  angle_surplus_bins_(100),
	  leg_shear_min_(0.0),
	  leg_shear_max_(12.5),
	  leg_shear_bins_(125),
	  measurements_(0),
	  num_log_(3000),
	  use_log_(true),
	  log_fresh_(true)
{
	sqlength_histogram_.resize(sqlength_bins_,0);
	angle_surplus_max_ = maxtheta_ - 2.0*PI;
	angle_surplus_histogram_.resize(angle_surplus_bins_,0);
	leg_shear_histogram_.resize(leg_shear_bins_,0);
}

HyperbolicStructure::~HyperbolicStructure()
{
}

void HyperbolicStructure::setMaxTheta(double max)
{
	maxtheta_=max;
	angle_surplus_max_ = maxtheta_ - 2.0*PI;
}

void HyperbolicStructure::Measure()
{
	if( circlepattern_->IsUpToDate() || circlepattern_->MakeUpToDate() )
	{
		int maxThetaInUnits = static_cast<int>(maxtheta_ * thetamodel_->getPiInUnits()/PI);
		std::list<std::pair<int,std::list<const Edge*> > > paths;
		thetamodel_->FindAllShortCurves(maxThetaInUnits,paths);

		for(std::list<std::pair<int,std::list<const Edge*> > >::iterator pathIt = paths.begin(); pathIt!= paths.end();pathIt++)
		{
			double length = HyperbolicLength(pathIt->second);
			int sqLengthBin = static_cast<int>(sqlength_bins_ *length*length/sqlength_max_);
			if( sqLengthBin >= 0 && sqLengthBin < sqlength_bins_ )
			{
				sqlength_histogram_[sqLengthBin]++;
			}

			int legs = NumberOfLegs(pathIt->second);
			if( legs >= static_cast<int>(legs_histogram_.size()) )
			{
				legs_histogram_.resize(legs+1);
			}
			legs_histogram_[legs]++;

			int anglesurplusBin = static_cast<int>(angle_surplus_bins_ * (static_cast<double>(pathIt->first)/thetamodel_->getPiInUnits() - 2.0)*PI/angle_surplus_max_);
			if( anglesurplusBin >= 0 && anglesurplusBin < angle_surplus_bins_ )
			{
				angle_surplus_histogram_[anglesurplusBin]++;
			}

			if( use_log_ )
			{
				log_ << (log_fresh_?"":",") << "{" << length << ", " << (static_cast<double>(pathIt->first)/thetamodel_->getPiInUnits() - 2.0)*PI << ",{";
				log_fresh_ = false;
			}
			MeasureLegShears(pathIt->second);
			if( use_log_ )
			{
				log_ << "}}";
				num_log_--;
				if( num_log_ == 0 )
				{
					use_log_ = false;
				}
			}
			/*std::cout << std::fixed << static_cast<double>(pathIt->first)/thetamodel_->getPiInUnits() << "*PI , " << length << ", (";
			for(std::list<const Edge*>::const_iterator it = pathIt->second.begin(); it != pathIt->second.end();it++)
			{
				std::cout << (it == pathIt->second.begin()?"":",") << circlepattern_->getShear(*it);
			}
			std::cout << ")\n";*/
		}
		measurements_++;
	}
}

double HyperbolicStructure::HyperbolicLength(const std::list<const Edge*> & path) const
{
	Eigen::Matrix2d mat = Eigen::Matrix2d::Identity();
	Eigen::Matrix2d turnRightMat;
	turnRightMat(0,0) = 0;
	turnRightMat(0,1) = -1;
	turnRightMat(1,0) = 1;
	turnRightMat(1,1) = 1;
	for(std::list<const Edge*>::const_iterator edgeIt = path.begin();edgeIt != path.end();edgeIt++)
	{
		double shear = circlepattern_->getShear(*edgeIt);
		double exphalfshear = std::exp(0.5*shear);
		Eigen::Matrix2d shearMat = Eigen::Matrix2d::Zero();
		shearMat(0,1) = exphalfshear;
		shearMat(1,0) = -1.0/exphalfshear;

		mat = shearMat * mat;
		mat = turnRightMat * mat;

		std::list<const Edge*>::const_iterator nextIt = boost::next(edgeIt);
		if( nextIt == path.end() )
		{
			nextIt = path.begin();
		}
		if( (*nextIt) == (*edgeIt)->getAdjacent()->getPrevious() )
		{
			mat = turnRightMat * mat;
		}
	}

	double trace = mat.trace();
	return 2.0*std::log(trace/2.0 + std::sqrt(trace*trace/4.0-1.0));
}

std::string HyperbolicStructure::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "hyperbolicstructure -> {maxtheta -> " << maxtheta_;
	stream << ", anglemax -> " << angle_surplus_max_;
	stream << ", anglebins -> " << angle_surplus_bins_;
	stream << ", anglehistogram -> ";
	PrintToStream(stream, angle_surplus_histogram_.begin(),angle_surplus_histogram_.end());
	stream << ", legshistogram -> ";
	PrintToStream(stream, legs_histogram_.begin(),legs_histogram_.end());
	stream << ", sqlengthmax -> " << sqlength_max_;
	stream << ", sqlengthbins -> " << sqlength_bins_;
	stream << ", sqlengthhistogram -> ";
	PrintToStream(stream,sqlength_histogram_.begin(),sqlength_histogram_.end());
	stream << ", legshearmin -> " << leg_shear_min_;
	stream << ", legshearmax -> " << leg_shear_max_;
	stream << ", legshearbins -> " << leg_shear_bins_;
	stream << ", legshearhistogram -> ";
	PrintToStream(stream,leg_shear_histogram_.begin(),leg_shear_histogram_.end());
	stream << ", measurements -> " << measurements_;
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

void HyperbolicStructure::MeasureLegShears(const std::list<const Edge*> & path)
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
			int shearBin = static_cast<int>(leg_shear_bins_*(shear - leg_shear_min_)/(leg_shear_max_-leg_shear_min_));
			if( shearBin >= 0 && shearBin < leg_shear_bins_ )
			{
				leg_shear_histogram_[shearBin]++;
			}
			if( use_log_ )
			{
				log_ << std::fixed << (first?"":",") << shear;
			}
			first = false;
			right = nextRight;
		}
	}
}