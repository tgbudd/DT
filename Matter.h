#ifndef DT_MATTER_H
#define DT_MATTER_H

#include <vector>
#include <string>

#include "Triangulation.h"
#include "Decoration.h"
#include "Edge.h"

class Matter : public Decoration
{
public:
	Matter() : use_custom_sweep_size_(false) {}
	~Matter(void) {}
	virtual bool IsFlipMoveAllowed(const Edge * const) {
		return true;
	}
	virtual void Initialize() = 0;
	virtual double BoltzmannChangeUnderFlipMove(const Edge * const ) const { 
		return 1.0; 
	}
	virtual double BoltzmannChangeUnderGeneralMove(const std::vector<boost::array<Triangle *,2> > & toBeDeleted, const std::vector<boost::array<Triangle *,2> > & toBeAdded ) const { 
		return 1.0; 
	}
	virtual void UpdateAfterFlipMove(const Edge * const) {}
	virtual void DoSweep() = 0;

	virtual double CentralCharge() const {
		return 0.0;
	}
	virtual std::string ConfigurationData() const {
		return "{}";
	}
	void SetCustomSweepSize(int n) {
		use_custom_sweep_size_ = (n>0);
		custom_sweep_size_ = n;
	}
	bool UsesCustomSweepSize() const {
		return use_custom_sweep_size_;
	}
	int CustomSweepSize() const {
		return custom_sweep_size_;
	}
private:
	bool use_custom_sweep_size_;	// Only used in special circumstances where we want to do 
	int custom_sweep_size_;		
};

#endif
