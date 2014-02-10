#ifdef _WIN32
#include <Windows.h>
#define GET_PROCESS_ID() GetCurrentProcessId()
#else
#include <sys/types.h>
#include <unistd.h>
#define GET_PROCESS_ID() getpid()
#endif

#include <fstream>
#include <sstream>

#include "utilities.h"
#include "Simulation.h"


Simulation::Simulation(Triangulation * const & triangulation, int ThermalizationSweeps, int SecondsPerOutput, bool output)
	: triangulation_(triangulation), thermalization_sweeps_(ThermalizationSweeps), seconds_per_output_(SecondsPerOutput), output_(output)
{
	char s[] = "aaaaaa";
	for(int i=0;i<6;i++)
	{
		s[i] = (char)triangulation_->RandomInteger((int)'a',(int)'z');
	}
	identifier_ = std::string(s);
	directory_ = "";
	prefix_ = "data-";
	std::ostringstream stream;
	stream << GET_PROCESS_ID();
	process_id_ = stream.str();
}


Simulation::~Simulation(void)
{
}

std::string Simulation::GetIdentifier() const
{
	return identifier_;
}

void Simulation::AddObservable(Observable * observable, int SweepsPerMeasurement)
{
	MeasurementInfo info;
	info.SweepsPerMeasurement = SweepsPerMeasurement;
	observables_.insert(std::pair<Observable *,MeasurementInfo>(observable,info));
}

void Simulation::Run()
{
	setStartTime();

	if( output_ )
	{
		std::cout << "Start thermalization.\n";
	}

	triangulation_->DoSweep(thermalization_sweeps_);

	time_t last_output = time(NULL);

	if( output_ )
	{
		std::cout << "Finished thermalization.\n";
	}

	for(sweeps_=1;true;sweeps_++)
	{
		triangulation_->DoSweep();


		for(std::map<Observable*,MeasurementInfo>::iterator observable=observables_.begin();observable!=observables_.end();observable++)
		{
			if( sweeps_ % observable->second.SweepsPerMeasurement == 0 )
			{
				if( output_ )
				{
					std::cout << "Measurement after " << sweeps_ << " sweeps.\n";
				}				
				observable->first->Measure();
			}
		}
		if( true )
		{
			time_t current_time = time(NULL);
			if( static_cast<int>(difftime(current_time,last_output)) > seconds_per_output_ )
			{
				Output();
				last_output = current_time;
			}
		}
	}

}

void Simulation::setStartTime()
{
	start_time_ = time(NULL);
}

void Simulation::SetDirectory(const std::string & dir)
{
	directory_ = dir;
}

void Simulation::Output()
{
	time_t currentTime = time(NULL);

	std::ostringstream filename;
	filename << directory_ << prefix_ << process_id_ << "-" << identifier_ << ".txt";
	std::ofstream file(filename.str().c_str());

	file << std::fixed;
	
	file << "{ identifier -> \"" << identifier_ << "\"";
	file << ", processid -> \"" << process_id_ << "\"";
	file << ", starttime -> \"" << RemoveNewline(ctime( &start_time_)) << "\"";
	file << ", outputtime -> \"" << RemoveNewline(ctime( &currentTime)) << "\"";
	file << ", runtimeinhours -> " << difftime(currentTime, start_time_)/3600.0;
	file << ", sweeps -> " << sweeps_;
	file << ", thermalizationsweeps -> " << thermalization_sweeps_;
	file << ", secondsperoutput -> " << seconds_per_output_;

	file << ", configuration -> ";

	std::ostringstream config;
	PrintToStream(config,configuration_info_.begin(),configuration_info_.end());
	file << config.str();

	file << ", " << triangulation_->OutputData();
	for(std::map<Observable*,MeasurementInfo>::iterator observable=observables_.begin();observable!=observables_.end();observable++)
	{
		file << ", " << observable->first->OutputData();
	}
	file << "}\n";

	if( output_ )
	{
		std::cout << "Output written to " << filename.str() << ".\n";
	}				

}

void Simulation::AddConfigurationInfo(std::string info)
{
	configuration_info_.push_back(info);
}
