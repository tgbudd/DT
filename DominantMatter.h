#ifndef DOMINANT_MATTER_H
#define DOMINANT_MATTER_H

#include <string>

class DominantMatter
{
public:
	DominantMatter() : use_custom_sweep_size_(false) {}
	~DominantMatter() {}
	virtual void DoSweep() = 0;
	virtual std::string ExportState() const { return ""; }

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
