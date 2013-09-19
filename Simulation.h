#ifndef DT_SIMULATION_H
#define DT_SIMULATION_H

#include <map>
#include <ctime>

#include "Triangulation.h"
#include "Observable.h"

class Simulation
{
public:
	Simulation(Triangulation * const & triangulation, int ThermalizationSweeps, int SweepsPerOutput );
	~Simulation(void);

	void AddObservable(Observable * observable, int SweepsPerMeasurement);
	void Run();
	void Output();

	void AddConfigurationInfo(std::string info);
private:
	void setStartTime();

	Triangulation * const triangulation_;
	int thermalization_sweeps_;
	int sweeps_per_output_;

	int sweeps_;
	time_t start_time_;

	std::string identifier_;
	std::string process_id_;
	std::string directory_;
	std::string prefix_;

	struct MeasurementInfo {
		int SweepsPerMeasurement;
	};

	std::map<Observable *,MeasurementInfo> observables_;
	std::vector<std::string> configuration_info_;
};

#endif

