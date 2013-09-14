#include <iostream>
#include <sstream>
#include <iomanip>

#include "triangulation.h"
#include "DualCohomologyBasis.h"
#include "CohomologyBasis.h"
#include "ThetaModel.h"
#include "BabyUniverseDistribution.h"

#include "HarmonicEmbedding.h"
#include "CirclePattern.h"
#include "BitmapDrawer.h"

int main(int argc, char* argv[])
{
	bool output = true;
	int w,h;
	int thermalizationSweeps, MeasurementSweeps;
	if( argc < 5 )
	{
		std::cout << "width = ";
		std::cin >> w;
		std::cout << "height = ";
		std::cin >> h;
		std::cout << "thermalization sweeps = ";
		std::cin >> thermalizationSweeps;
		std::cout << "measurement sweeps = ";
		std::cin >> MeasurementSweeps;
	}else
	{
		std::istringstream is(argv[1]);
		is >> w;
		std::istringstream is2(argv[2]);
		is2 >> h;
		std::istringstream is3(argv[3]);
		is3 >> thermalizationSweeps;
		std::istringstream is4(argv[4]);
		is4 >> MeasurementSweeps;
		output = false;
	}

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	
	DualCohomologyBasis dualcohomologybasis( &triangulation );
	dualcohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &dualcohomologybasis );

	ThetaModel thetamodel( &triangulation, &dualcohomologybasis );
	triangulation.setDominantMatter( &thetamodel );
	thetamodel.Initialize();

	CohomologyBasis cohomologybasis(dualcohomologybasis);
	BabyUniverseDistribution babyuniv( &triangulation, &cohomologybasis, 3);

	triangulation.DoSweep(thermalizationSweeps);

	int measurements=0;

	BitmapDrawer bitmap(600,400,2);
	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );
	CirclePattern circlepattern( &triangulation, &cohomologybasis, &thetamodel );
	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	TriangulationDrawer tridrawer2( &triangulation, &circlepattern );

	int i=0;
	
	while(true)
	{
		cohomologybasis.SetToDualOf(dualcohomologybasis);
		cohomologybasis.Simplify();
		dualcohomologybasis.SetToDualOf(cohomologybasis);
		babyuniv.Measure();

		std::list<std::list<const Edge *> > paths;
		babyuniv.FindMinbuNecksOfLength3(paths);
		harmonicembedding.FindEmbedding();
		bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.3);
		bitmap.Clear();
		bitmap.setPenColor(100,100,250);
		BOOST_FOREACH(const std::list<const Edge *> & path, paths)
		{
			std::vector<Vector2D> coor;
			Vector2D currentCoor = harmonicembedding.getCoordinate(path.front()->getNext()->getOpposite());
			BOOST_FOREACH(const Edge * edge, path )
			{
				coor.push_back( currentCoor );
				currentCoor = AddVectors2D( currentCoor, harmonicembedding.getForm(edge) );
			}
			bitmap.domainPolygon(coor);
		}

		bitmap.setPenWidth(4);
		bitmap.setPenColor(0,0,0);
		tridrawer.Draw(bitmap);
		circlepattern.FindEmbedding();
		bitmap.SetPeriodicDomain(circlepattern.CalculateModuli(),0.3);
		bitmap.setPenColor(200,50,50);
		tridrawer2.Draw(bitmap);

		std::ostringstream os;
		os << "D:\\temp\\dt\\output\\test-" << std::setw( 4 ) << std::setfill( '0' ) << i++ << ".bmp";
		bitmap.SaveImage(os.str());

	
		measurements++;
		if( measurements % 1 == 0 )
		{
			std::cout << babyuniv.OutputData() << "\n\n";
		}

		triangulation.DoSweep( MeasurementSweeps );
	}


	return 0;
}