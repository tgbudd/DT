#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "DualCohomologyBasis.h"
#include "ThetaModel.h"
#include "CirclePattern.h"
#include "BitmapDrawer.h"

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
		//UpdateEdgeShade();
		
		bitmap_.Clear();

		bitmap_.SetPeriodicDomain(embedding_->CalculateModuli(),0.24);
		tridrawer_.DrawShading(bitmap_);
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
	ParameterStream param(argc,argv);

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	double maxtheta = param.Read<double>("maxtheta (x PI)");
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	
	DualCohomologyBasis dualcohomologybasis( &triangulation );
	dualcohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &dualcohomologybasis );

	ThetaModel thetamodel( &triangulation, &dualcohomologybasis, 18000 );
	triangulation.setDominantMatter( &thetamodel );
	thetamodel.Initialize();

	Simulation simulation( &triangulation, thermalizationSweeps, 1000000 );

	CohomologyBasis cohom( &triangulation, dualcohomologybasis );
	cohom.SetMakeUpToDateViaReinitialization(true);

	CirclePattern circlepattern( &triangulation, &cohom, &thetamodel );

	Snapshot snapshot( &triangulation, &circlepattern, &thetamodel );
	std::ostringstream prefix;
	prefix << "./output/snap-";
	snapshot.SetPrefix( prefix.str() );

	simulation.AddObservable( &snapshot, measurementSweeps );

	simulation.Run();

	return 0;
}