#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "DualCohomologyBasis.h"
#include "ThetaModel.h"
#include "CirclePattern.h"
#include "BitmapDrawer.h"
#include "HarmonicEmbedding.h"

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, ThetaModel * thetamodel );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void Save();
	void UpdateEdgeShade();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;
	ThetaModel * thetamodel_;

	std::string info_;
	std::string prefix_;
	int counter_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, ThetaModel * thetamodel ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  thetamodel_(thetamodel),
	  bitmap_(1280,720,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  counter_(0)
{
	textdrawer_.setColor(160,0,30);
}

void Snapshot::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		UpdateEdgeShade();
		
		bitmap_.Clear();

		bitmap_.SetPeriodicDomain(embedding_->CalculateModuli(),0.24);
		//tridrawer_.DrawShading(bitmap_);
		bitmap_.setPenWidth(4);
		bitmap_.setPenColor(0,0,0);
		tridrawer_.Draw(bitmap_);
		std::ostringstream text;
		text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
		text << info_;
		textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

		Save();
	}
}

void Snapshot::Save()
{
	std::ostringstream os;
	os << prefix_ << std::setw( 4 ) << std::setfill( '0' ) << counter_++ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
}

void Snapshot::UpdateEdgeShade()
{
	for(int i=0,end=triangulation_->NumberOfTriangles();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			tridrawer_.SetEdgeShade( i, j, 1.0 - thetamodel_->getRealTheta(i,j) / PI );
		}
	}
}

int main(int argc, char* argv[])
{
	bool output = true;
	int w,h;
	int thermalizationSweeps, MeasurementSweeps, OutputSweeps;
	double cosinepower;
	if( argc < 7 )
	{
		std::cout << "width = ";
		std::cin >> w;
		std::cout << "height = ";
		std::cin >> h;
		std::cout << "power = ";
		std::cin >> cosinepower;
		std::cout << "thermalization sweeps = ";
		std::cin >> thermalizationSweeps;
		std::cout << "measurement sweeps = ";
		std::cin >> MeasurementSweeps;
		std::cout << "output sweeps = ";
		std::cin >> OutputSweeps;
	}else
	{
		std::istringstream is(argv[1]);
		is >> w;
		std::istringstream is2(argv[2]);
		is2 >> h;
		std::istringstream is3(argv[3]);
		is3 >> cosinepower;
		std::istringstream is4(argv[4]);
		is3 >> thermalizationSweeps;
		std::istringstream is5(argv[5]);
		is4 >> MeasurementSweeps;
		std::istringstream is6(argv[6]);
		is5 >> OutputSweeps;
		output = false;
	}

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	
	DualCohomologyBasis dualcohomologybasis( &triangulation );
	dualcohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &dualcohomologybasis );

	ThetaModel thetamodel( &triangulation, &dualcohomologybasis, 18000 );
	triangulation.setDominantMatter( &thetamodel );
	thetamodel.Initialize();
	thetamodel.setCosinePower( cosinepower );

	Simulation simulation( &triangulation, thermalizationSweeps, OutputSweeps );
	std::ostringstream config;
	config << "cosinepower -> " << std::fixed << cosinepower;
	simulation.AddConfigurationInfo(config.str());

	CohomologyBasis cohom( &triangulation, dualcohomologybasis );
	cohom.SetMakeUpToDateVia( &dualcohomologybasis );

	CirclePattern circlepattern( &triangulation, &cohom, &thetamodel );
	HarmonicEmbedding harm(  &triangulation, &cohom );

	Snapshot snapshot( &triangulation, &circlepattern, &thetamodel );
	std::ostringstream prefix;
	prefix << "./output/snapshot-" << cosinepower << "-";
	snapshot.SetPrefix( prefix.str() );
	std::ostringstream info;
	info << "Edge weight: |cos(theta)|^" << cosinepower;
	snapshot.SetInfo( info.str() );

	Snapshot snapshot2( &triangulation, &harm, &thetamodel );
	std::ostringstream prefix2;
	prefix2 << "./output/snapshot-harm-" << cosinepower << "-";
	snapshot2.SetPrefix( prefix2.str() );
	std::ostringstream info2;
	info2 << "Edge weight: |cos(theta)|^" << cosinepower;
	snapshot2.SetInfo( info2.str() );

	simulation.AddObservable( &snapshot, MeasurementSweeps );
	simulation.AddObservable( &snapshot2, MeasurementSweeps );

	simulation.Run();

	return 0;
}