#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "ConnectivityRestrictor.h"
#include "PeelingProcedure.h"
#include "CMinusTwoBuilder.h"
#include "LinearAlgebra.h"

class PeelingMeasurement : public Observable {
public:
	PeelingMeasurement( Triangulation * triangulation );
	void Measure();
	std::string OutputData();

private:
	void MeasureBoundary();
	Triangle * SelectTriangle();

	Triangulation * triangulation_;
	PeelingProcedure peelingprocedure_;
	BabyUniverseDetector babyuniversedetector_;

	std::list<const Edge *> boundary_;
	std::vector<double> measure_;
	std::vector<boost::array<int,3> > in_boundary_;
	std::vector<int> walk_end_;

	int random_walks_;
	void DoRandomWalk();
	int max_diffusion_steps_;
	boost::array<std::vector<double>,2> field_;
	void DoDiffusion();

	void InvertMatrix();
};

PeelingMeasurement::PeelingMeasurement(Triangulation * triangulation ) 
	: triangulation_(triangulation), 
	  peelingprocedure_(triangulation),
	  random_walks_(100),
	  max_diffusion_steps_(100000),
	  babyuniversedetector_(triangulation)
{
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
			MeasureBoundary();

			break;
		}
	} while( !peelingprocedure_.DoPeelingStepOnSphere() );
}

void PeelingMeasurement::MeasureBoundary()
{
	measure_.resize(boundary_.size());
	std::fill(measure_.begin(),measure_.end(),0.0);
	walk_end_.resize(boundary_.size());
	std::fill(walk_end_.begin(),walk_end_.end(),0);
	boost::array<int,3> allminusone = {-1,-1,-1};
	in_boundary_.resize(triangulation_->NumberOfTriangles());
	std::fill(in_boundary_.begin(),in_boundary_.end(),allminusone);
	int index=0;
	for(std::list<const Edge *>::const_iterator it=boundary_.begin();it!=boundary_.end();it++,index++)
	{
		in_boundary_[(*it)->getParent()->getId()][(*it)->getId()] = index;
	}

	InvertMatrix();


/*	for(int i=0;i<random_walks_;i++)
	{
		DoRandomWalk();
	}
	DoDiffusion();
	*/

	double tot=0.0;
	for(int i=0,endi=measure_.size();i<endi;i++)
	{
		tot+=measure_[i];
	}
	std::cout << "Total measure = " << std::fixed << tot << "\n";
}

void PeelingMeasurement::DoRandomWalk()
{
	Vertex * vertex = peelingprocedure_.final_vertex_;
	int degree = vertex->getDegree();
	int direction = triangulation_->RandomInteger(0,degree-1);
	Edge * edge = vertex->getParent();
	for(int i=0;i<direction;i++)
	{
		edge = edge->getNext()->getAdjacent()->getNext();
	}
	Triangle * triangle = edge->getParent();

	int inboundary = -1;
	while( inboundary == -1 )
	{
		int nbr = triangulation_->RandomInteger(0,2);
		inboundary = in_boundary_[triangle->getId()][nbr];
		triangle = triangle->getEdge(nbr)->getAdjacent()->getParent();
	}

	walk_end_[inboundary]++;
}

void PeelingMeasurement::DoDiffusion()
{
	field_[0].resize(triangulation_->NumberOfTriangles());
	field_[1].resize(triangulation_->NumberOfTriangles());
	std::fill(field_[0].begin(),field_[0].end(),0.0);
	std::fill(field_[1].begin(),field_[1].end(),0.0);

	Vertex * vertex = peelingprocedure_.final_vertex_;
	double initial = 1.0/vertex->getDegree();
	Edge * edge = vertex->getParent();
	do
	{
		field_[0][edge->getParent()->getId()] = initial;
		edge = edge->getNext()->getAdjacent()->getNext();
	} while( edge != vertex->getParent() );

	std::vector<Triangle *> triangles;
	babyuniversedetector_.EnclosedTriangles(boundary_,triangles,true);

	int currentField = 0, nextField = 1;
	double totalBoundaryMeasure = 0.0;
	int steps = 0;
	std::cout << "start - ";
	double prevsave=totalBoundaryMeasure;
	while( totalBoundaryMeasure < 0.999 && steps < max_diffusion_steps_ )
	{
		int index=0;
		for(std::list<const Edge *>::const_iterator it=boundary_.begin();it!=boundary_.end();it++,index++)
		{
			double toboundary = field_[currentField][(*it)->getParent()->getId()]/3.0;
			measure_[index] += toboundary;
			totalBoundaryMeasure += toboundary;
		}
		for(std::vector<Triangle *>::const_iterator it=triangles.begin();it!=triangles.end();it++)
		{
			int id = (*it)->getId();
			double newField=0.0;
			for(int j=0;j<3;j++)
			{
				newField += field_[currentField][(*it)->getEdge(j)->getAdjacent()->getParent()->getId()];
			}
			field_[nextField][id] = newField/3.0;
		}
		std::swap(currentField,nextField);
		steps++;

		//////
		/*if( totalBoundaryMeasure - prevsave > 0.1 || (totalBoundaryMeasure > 0.95 && totalBoundaryMeasure - prevsave > 0.001 ) )
		{
			std::ofstream file("measure.txt",std::ios::app);
			file << "{" << std::fixed << totalBoundaryMeasure << "," << steps << ",";
			PrintToStream(file,measure_.begin(),measure_.end());
			file << "}\n";
			prevsave = totalBoundaryMeasure;
		}*/
	}

}

void PeelingMeasurement::InvertMatrix()
{
	std::vector<Triangle *> triangles;
	babyuniversedetector_.EnclosedTriangles(boundary_,triangles,true);

	linearalgebra::PositiveDefiniteDenseMatrix matrix(triangles.size());

	std::vector<int> trianglePosition(triangulation_->NumberOfTriangles(),-1);
	for(int i=0,endi=triangles.size();i<endi;i++)
	{
		trianglePosition[triangles[i]->getId()] = i;
	}
	double minusonethird = -1/3.0;
	for(int i=0,endi=triangles.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			int pos = trianglePosition[triangles[i]->getEdge(j)->getAdjacent()->getParent()->getId()];
			if( pos >= 0 )
			{
				matrix.Set(i,pos,minusonethird);
			}
		}
		matrix.Set(i,i,1.0);
	}

	if( matrix.ComputeInverse() )
	{
		std::cout << "Inverse computed.\n";
		int index=0;
		for(std::list<const Edge *>::const_iterator it=boundary_.begin();it!=boundary_.end();it++,index++)
		{
			measure_[index] = 0.0;
				
			Vertex * vertex = peelingprocedure_.final_vertex_;
			double factor = 1.0/vertex->getDegree()/3.0;
			Edge * edge = vertex->getParent();
			do
			{
				measure_[index] += factor * matrix.Get(trianglePosition[(*it)->getParent()->getId()],trianglePosition[edge->getParent()->getId()]);
				edge = edge->getNext()->getAdjacent()->getNext();
			} while( edge != vertex->getParent() );
		}
	} else
	{
		std::cout << "Inverse failed.\n";
	}
}

Triangle * PeelingMeasurement::SelectTriangle() 
{
	return triangulation_->getRandomTriangle();
}

std::string PeelingMeasurement::OutputData() 
{ 
	return ""; 
}


int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, 0, 10000, output );

	PeelingMeasurement peelingmeasurement( &triangulation );

	simulation.AddObservable( &peelingmeasurement, 1 );

	simulation.Run();
	return 0;
}