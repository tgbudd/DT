#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "DualCohomologyBasis.h"
#include "ThetaModel.h"
#include "CirclePattern.h"
#include "BitmapDrawer.h"
#include "utilities.h"
#include "HyperbolicStructure.h"

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, ThetaModel * thetamodel, HyperbolicStructure * hyperbolicstructure );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void Save();
	void DrawCurves();
	std::pair<Vector2D,double> DetermineRange();

	std::list<HyperbolicStructure::Curve> curves_;

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;
	ThetaModel * thetamodel_;
	HyperbolicStructure * hyperbolicstructure_;

	std::string info_;
	std::string prefix_;
	int counter_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, ThetaModel * thetamodel, HyperbolicStructure * hyperbolicstructure ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  thetamodel_(thetamodel),
	  hyperbolicstructure_(hyperbolicstructure),
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
		int num = hyperbolicstructure_->FindShortHyperbolicCurves();	
		hyperbolicstructure_->RetrieveCurves(curves_);

		if( num == 0 )
		{
			return;
		}

		bitmap_.Clear();

		std::pair<Vector2D,double> range = DetermineRange();
		std::pair<double,double> modulus = embedding_->CalculateModuli();
		bitmap_.SetPeriodicDomain(modulus, std::min(1/(2.5 * range.second * range.second),500.0), 0.5, 0.5, range.first[0], range.first[1] );

		tridrawer_.DrawShading(bitmap_);
		bitmap_.setPenWidth(4);
		bitmap_.setPenColor(0,0,0);
		tridrawer_.Draw(bitmap_);

		bitmap_.setPenWidth(6);
		DrawCurves();

		std::ostringstream text;
		text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
		text << info_;
		textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

		Save();
	}
}

std::pair<Vector2D,double> Snapshot::DetermineRange()
{
	HyperbolicStructure::Curve & curve = curves_.front();
	Vector2D v = curve.segments.front()[0];
	Vector2D mean = {0,0};
	for(std::list<boost::array<Vector2D,2> >::iterator segment = curve.segments.begin();segment != curve.segments.end();segment++)
	{
		v = AddVectors2D(v,SubtractVectors2D((*segment)[1],(*segment)[0]));
		mean = AddVectors2D(mean,v);
	}
	mean[0] /= curve.segments.size();
	mean[1] /= curve.segments.size();
	double maxdist = 0.0;
	for(std::list<boost::array<Vector2D,2> >::iterator segment = curve.segments.begin();segment != curve.segments.end();segment++)
	{
		v = AddVectors2D(v,SubtractVectors2D((*segment)[1],(*segment)[0]));
		double dist = Norm2D(SubtractVectors2D(mean,v));
		if( dist > maxdist )
		{
			maxdist = dist;
		}
	}

	return std::pair<Vector2D,double>(mean,maxdist);
}

void Snapshot::DrawCurves()
{
	ColorScheme scheme(ColorScheme::BLUE_GREEN_YELLOW);

	for(std::list<HyperbolicStructure::Curve>::iterator curve = curves_.begin();curve!=curves_.end();curve++)
	{
		boost::array<unsigned char,3> color = scheme.getColor(1.0-std::min(curve->length/1.5,1.0));
		bitmap_.setPenColor(color[0],color[1],color[2]);

		for(std::list<boost::array<Vector2D,2> >::iterator segment = curve->segments.begin();segment != curve->segments.end();segment++)
		{
			bitmap_.domainLineSegment((*segment)[0][0],(*segment)[0][1],(*segment)[1][0],(*segment)[1][1]);
		}
	}
}

void Snapshot::Save()
{
	std::ostringstream os;
	os << prefix_ << std::setw( 4 ) << std::setfill( '0' ) << counter_++ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
}


int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("snapshot sweeps");
	double maxlength = param.Read<double>("maxlength");
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

	HyperbolicStructure hyp( &triangulation, &thetamodel, &circlepattern, &dualcohomologybasis);
	hyp.setMaxLength(maxlength);

	Snapshot snapshot( &triangulation, &circlepattern, &thetamodel, &hyp );
	std::ostringstream prefix;
	prefix << "./output/snap-geo-";
	snapshot.SetPrefix( prefix.str() );

	simulation.AddObservable( &snapshot, measurementSweeps );

	simulation.Run();

	return 0;
}