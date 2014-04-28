#ifndef THETA_MODEL_H
#define THETA_MODEL_H

#include <vector>
#include <map>

#include "Triangulation.h"
#include "DominantMatter.h"
#include "DualCohomologyBasis.h"

class ThetaModel :
	public DominantMatter
{
public:
	ThetaModel(Triangulation * const triangulation, const DualCohomologyBasis * const dualcohomologybasis, int PiInUnits = 1800);
	~ThetaModel(void);

	void Initialize();
	void DoSweep();

	void setTheta(Edge * const & edge, int theta) {
		BOOST_ASSERT( theta > 0 && theta < 2 * pi_in_units_ );
		theta_[edge->getParent()->getId()][edge->getId()] = theta;
		theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] = theta;
	}
	void addToTheta(Edge * const & edge, int theta) {
		BOOST_ASSERT( theta_[edge->getParent()->getId()][edge->getId()] == theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] );
		int firsttheta = (theta_[edge->getParent()->getId()][edge->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0 && firsttheta < 2 * pi_in_units_ );
		firsttheta = (theta_[edge->getAdjacent()->getParent()->getId()][edge->getAdjacent()->getId()] += theta);
		BOOST_ASSERT( firsttheta > 0 && firsttheta < 2 * pi_in_units_ );
	}
	int getTheta(const Edge * const & edge) const {
		return theta_[edge->getParent()->getId()][edge->getId()]; 
	}
	int getTheta(int triangleId, int edgeId) const {
		return theta_[triangleId][edgeId]; 
	}
	double getRealTheta(Edge * const & edge ) const {
		return PI * (double)(getTheta(edge))/pi_in_units_;
	}
	double getRealTheta(int triangleId, int edgeId) const {
		return PI * (double)(getTheta(triangleId,edgeId))/pi_in_units_;
	}
	int getSize() const {
		return static_cast<int>(theta_.size());
	}
	void setCosinePower(double CosinePower)
	{
		cos_power_ = CosinePower;
		if( cos_power_ < -1e-8 || cos_power_ > 1e-8 )
		{
			use_cos_power_ = true;
		}
	}

	bool TryThetaMove();
	void TryThetaMove(int n);
	bool TestAllCutConditions() const; // for debugging

	std::string ExportState() const;

	void FindAllShortCurves(int maxtheta, std::list<std::pair<int,std::list<const Edge*> > > & paths) const; 
	int getPiInUnits() const {
		return pi_in_units_;
	}
private:
	double cosine(int theta) {
		return std::cos( static_cast<double>(theta * PI/pi_in_units_) );
	}

	std::vector<boost::array<int,3> > theta_;
	const int pi_in_units_;	// use integer angles where pi_in_units corresponds to an angle of pi
	bool use_cos_power_;
	double cos_power_;

	Triangulation * const triangulation_;
	const DualCohomologyBasis * const dualcohomologybasis_;

	bool TryThetaMove(Edge * edge);

	void CleanUpDistance(Triangle * triangle) const;
	bool TestCutCondition(Edge *, Edge *, const IntForm2D &, int theta) const;
	bool TestCutConditionOld(Edge *, Edge *, const IntForm2D &, int theta) const;
	bool TestCutCondition(Edge *) const;

	void FindShortCurve(Edge *, int maxtheta, std::pair<int,std::list<const Edge*> > & path) const;


	Edge * previousEdge(Edge * edge);  // functions of which pointers are used in the thetamove
	Edge * nextEdge(Edge * edge);

	bool TestVertexSum(Vertex *);

	struct triangleNode
	{
		Triangle* triangle;
		IntForm2D integral;
		int distance;
		bool operator<(const triangleNode & node) const { return distance > node.distance; }
	};
	mutable std::vector<std::map<IntForm2D,int> > distance_; // used in the Dijkstra algorithm of TestCutCondition
	void RetrievePath(const Edge * lastEdge, const triangleNode & end, const triangleNode & begin, std::pair<int,std::list<const Edge*> > & path) const;

};

#endif
