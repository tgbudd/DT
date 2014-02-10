#include <iostream>
#include <sstream>
#include <iomanip>


#define LOG_CIRCLE_PACKING

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "DiskCirclePacking.h"
#include "ConnectivityRestrictor.h"
#include "PeelingProcedure.h"
#include "BitmapDrawer.h"
#include "CMinusTwoBuilder.h"

class PeelingMeasurement : public Observable {
public:
	PeelingMeasurement( Triangulation * triangulation );
	void Measure();
	std::string OutputData() { return ""; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }

private:
	void MeasureBoundary();
	void DrawBitmap();
	Triangle * SelectTriangle();
	void SaveBitmap();
	void SaveDisk();

	Triangulation * triangulation_;
	DiskCirclePacking diskcirclepacking_;
	PeelingProcedure peelingprocedure_;

	std::list<const Edge *> boundary_;
	std::vector<std::pair<Vertex*,std::pair<Vector2D,double> > > circles_;

	BitmapDrawer bitmap_;
	int fullwidth_, fullheight_;
	TextDrawer textdrawer_;
	DiscreteColorScheme scheme_;

	std::string prefix_;
	int counter_;
};

PeelingMeasurement::PeelingMeasurement(Triangulation * triangulation ) 
	: triangulation_(triangulation), 
	  diskcirclepacking_(triangulation),
	  peelingprocedure_(triangulation),
	  fullwidth_(6400),
	  fullheight_(6400),
	  bitmap_(1600,1600,4),
	  counter_(0),
  	  textdrawer_("lucida",2),
	  scheme_(DiscreteColorScheme::COLORS8)
{
	textdrawer_.setColor(160,0,30);
}

void PeelingMeasurement::Measure()
{
	Triangle * triangle = SelectTriangle();
		
	peelingprocedure_.PreparePeelingOnSphere( triangle );
	Edge * centerEdge = peelingprocedure_.final_vertex_->getParent()->getPrevious();
	do {
		if( peelingprocedure_.volume_within_frontier_ > triangulation_->NumberOfTriangles()/2)
		{
			boundary_ = peelingprocedure_.frontier_;
			boundary_.reverse();
			for(std::list<const Edge *>::iterator it = boundary_.begin(); it != boundary_.end(); it++)
			{
				(*it) = (*it)->getAdjacent();
			}
			if( diskcirclepacking_.FindEmbedding(boundary_,centerEdge) )
			{
				std::cout << diskcirclepacking_.getStepsUsed() << " steps\n";
				circles_.clear();
				diskcirclepacking_.getCircles(circles_);

				MeasureBoundary();
				//DrawBitmap();
				//SaveBitmap();
				//SaveDisk();
				//int hy=0;
				//std::cout << "hy: ";
				//std::cin >> hy;
			} else
			{
				std::cout << "Failed to find embedding.\n";
			}
			break;
		}
	} while( !peelingprocedure_.DoPeelingStepOnSphere() );
}

void PeelingMeasurement::MeasureBoundary()
{
	std::vector<double> angles;
	diskcirclepacking_.getBoundaryPositions(angles);
	double time = -std::log(diskcirclepacking_.getCenterRadius());

	std::ofstream file("measure.txt",std::ios::app);
	file << std::fixed << "{volume -> " << triangulation_->NumberOfTriangles() << ", circles -> " << circles_.size() << ", time -> " << time << ", distance -> " << peelingprocedure_.distance_[boundary_.front()->getNext()->getOpposite()->getId()] << ", angles -> ";
	std::ostringstream os;
	PrintToStream(os,angles.begin(),angles.end());
	file << os.str() << "}\n";
}

Triangle * PeelingMeasurement::SelectTriangle() 
{
	return triangulation_->getRandomTriangle();
}

void PeelingMeasurement::DrawBitmap()
{
	bitmap_.Clear();
	bitmap_.setPenColor(200,200,200);
	bitmap_.Disk(std::pair<int,int>(fullwidth_/2,fullheight_/2),static_cast<int>(0.48*fullheight_));
	bitmap_.setPenColor(20,20,200);
	double centerradius;
	int peelpos=-1,peelpos2=-1;
	for(int i=0,endi=circles_.size();i<endi;i++)
	{
		boost::array<unsigned char,3> c = scheme_.getColor(peelingprocedure_.distance_[circles_[i].first->getId()]);
		if( circles_[i].first == peelingprocedure_.final_vertex_ )
		{
			bitmap_.setPenColor(10,10,10);
			centerradius = circles_[i].second.second;
		} else
		{
			bitmap_.setPenColor(c[0],c[1],c[2]);				
		}
		std::pair<int,int> x;
		x.first = static_cast<int>(fullwidth_/2 + fullheight_ * 0.48 * circles_[i].second.first[0] );
		x.second = static_cast<int>(fullheight_/2 + fullheight_ * 0.48 * circles_[i].second.first[1] );
		bitmap_.Disk(x, static_cast<int>(fullheight_ * 0.48 * circles_[i].second.second) );
		if( peelingprocedure_.in_frontier_[circles_[i].first->getId()] )
		{
			std::pair<int,int> y1,y2;
			if( circles_[i].first == (boundary_.front()->getNext()->getOpposite()) )
			{
				peelpos = i;
			} else if( circles_[i].first == (boundary_.front()->getPrevious()->getOpposite()) )
			{
				peelpos2 = i;
			}
			bitmap_.setPenColor(80,80,80);
			bitmap_.setPenWidth(6);
			y1.first = static_cast<int>(fullwidth_/2 + fullheight_ * 0.48 * circles_[i].second.first[0]/Norm2D(circles_[i].second.first) );
			y1.second = static_cast<int>(fullheight_/2 + fullheight_ * 0.48 * circles_[i].second.first[1]/Norm2D(circles_[i].second.first) );
			y2.first = static_cast<int>(fullwidth_/2 + fullheight_ * 0.495 * circles_[i].second.first[0]/Norm2D(circles_[i].second.first) );
			y2.second = static_cast<int>(fullheight_/2 + fullheight_ * 0.495 * circles_[i].second.first[1]/Norm2D(circles_[i].second.first) );
			bitmap_.LineSegment(y1,y2);
		}
	}

	if( peelpos >= 0 )
	{
		bitmap_.setPenColor(180,0,70);
		//bitmap_.setPenWidth(8);
		std::pair<int,int> y1;
		Vector2D v = Normalize(AddVectors2D(Normalize(circles_[peelpos].second.first),Normalize(circles_[peelpos2].second.first)));
		y1.first = static_cast<int>(fullwidth_/2 + fullheight_ * 0.488 * v[0] );
		y1.second = static_cast<int>(fullheight_/2 + fullheight_ * 0.488 * v[1] );
		bitmap_.Disk(y1,static_cast<int>(fullheight_ * 0.007));
	}
	std::ostringstream text;
	text << "Circles: " << circles_.size() << "\n";
	text << std::fixed << "-log(r): " << -std::log(centerradius) << "\n";
	text << "steps: " << diskcirclepacking_.getStepsUsed() << "\n";
	text << "triangles: " << triangulation_->NumberOfTriangles() - peelingprocedure_.volume_within_frontier_ << "/" << triangulation_->NumberOfTriangles() << "\n";
	text << "boundary: " << boundary_.size() << "\n";
	text << "distance: " << peelingprocedure_.distance_[boundary_.front()->getNext()->getOpposite()->getId()] << "/" << peelingprocedure_.distance_[peelingprocedure_.final_vertex_->getId()] << "\n";
	text << "model: c = -2";
	textdrawer_.DrawText(text.str(),bitmap_,30,30,0.);
}

void PeelingMeasurement::SaveBitmap()
{
	std::ostringstream os;
	os << prefix_ << "-" << std::setw( 4 ) << std::setfill( '0' ) << counter_ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
	counter_++;
}

void PeelingMeasurement::SaveDisk()
{
	std::vector<int> vertexIds(triangulation_->NumberOfVertices(),-1);
	for(int i=0,endi=circles_.size();i<endi;i++)
	{
		vertexIds[circles_[i].first->getId()] = i;
	}
	std::vector<Triangle *> triangles = diskcirclepacking_.getDiskTriangles();
	
	std::ofstream file("disk.txt");
	file << std::fixed << "{triangles -> {{";
	for(int i=0,endi=triangles.size();i<endi;i++)
	{
		file << (i>0?"},{":"");
		for(int j=0;j<3;j++)
		{
			file << (j>0?",":"") << vertexIds[triangles[i]->getEdge(j)->getOpposite()->getId()];
		}
	}
	file << "}}, boundary -> {";
	for(std::list<const Edge *>::const_iterator it = boundary_.begin();it!= boundary_.end();it++)
	{
		file << (it == boundary_.begin()?"":",") << vertexIds[(*it)->getNext()->getOpposite()->getId()];
	}
	file << "}, circles -> {";
	for(int i=0,endi=circles_.size();i<endi;i++)
	{
		file << (i>0?",":"") << "{{" << circles_[i].second.first[0] << "," << circles_[i].second.first[1] << "}," << circles_[i].second.second << "}";
	}
	file << "}}\n";
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	//int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	//int	measurementSweeps = param.Read<int>("measurement sweeps");
	bool output = param.UserInput();

	Triangulation triangulation;
	//triangulation.LoadSphericalBySubdivision(n);
	CMinusTwoBuilder builder(&triangulation,0,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, 0, 10000, output );

	//ConnectivityRestrictor connect( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	//triangulation.AddMatter( &connect );

	PeelingMeasurement peelingmeasurement( &triangulation );

	std::ostringstream prefix;
	prefix << "D:/temp/output/circ-" << simulation.GetIdentifier() << "-";
	peelingmeasurement.SetPrefix( prefix.str() );
	
	simulation.AddObservable( &peelingmeasurement, 1 );

	simulation.Run();
	return 0;
}