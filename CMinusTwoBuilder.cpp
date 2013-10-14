#include <stack>

#include "CMinusTwoBuilder.h"
#include "Triangle.h"
#include "Edge.h"

CMinusTwoBuilder::CMinusTwoBuilder(Triangulation * const triangulation, int Genus, int NumberOfTriangles)
	: triangulation_(triangulation), genus_(Genus), n_triangles_(NumberOfTriangles)
{
	BOOST_ASSERT( genus_ == 0 || genus_ == 1 );
	BOOST_ASSERT( n_triangles_ >= 2 && n_triangles_ % 2 == 0 );
	tree_list_.resize(2*n_triangles_ + 1,0);
}


CMinusTwoBuilder::~CMinusTwoBuilder(void)
{
}

void CMinusTwoBuilder::DoSweep()
{
	triangulation_->Clear();
	RandomDiskTriangulation();
	RandomBoundaryMatching();
	BOOST_ASSERT( genus_ == 0 );
	if( genus_ == 0 )
	{
		ApplyBoundaryMatching(matching_);
	}

	triangulation_->DetermineVertices();
}

void CMinusTwoBuilder::RandomDiskTriangulation()
{
	// use algorithm R in section 7.2.1.6 from Knuth TAOCP to generate a random binary tree in list form
	tree_list_[0] = 0;
	for(int n=0;n<n_triangles_;++n)
	{
		int x = triangulation_->RandomInteger( 0, 4 * n + 1 );
		tree_list_[2*n + 2 - (x%2)] = 2*(n+1);
		tree_list_[2*n + 1 + (x%2)] = tree_list_[x/2];
		tree_list_[x/2] = 2*(n+1)-1;
	}
	for(int i=0,end=tree_list_.size();i<end;i++)
	{
		if( tree_list_[i]%2 == 0 )
			tree_list_[i] = n_triangles_ + 1 + tree_list_[i]/2;
		else
			tree_list_[i] = (tree_list_[i]-1)/2;
	}
	// Now tree_list_[0] is the root simplex-id and is connected to boundary link n_triangles_.
	// The children of simplex i are tree_list_[2*i+1] and tree_list_[2*i+2].
	// Values >= n_triangles_ correspond to boundary links.

	for(int i=0;i<n_triangles_;++i)
	{
		triangulation_->NewTriangle();
	}
	for(int i=0;i<n_triangles_;++i)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<2;j++)
		{
			if( tree_list_[2*i+j+1] < n_triangles_ )
			{
				triangle->getEdge(j)->bindAdjacent(triangulation_->getTriangle(tree_list_[2*i+j+1])->getEdge(2));
			}
		}
	}
	
	// Determine the boundary of the triangulation
	std::stack<Edge *> frontier;

	Triangle * rootTriangle = triangulation_->getTriangle(tree_list_[0]);
	frontier.push( rootTriangle->getEdge(1) );
	frontier.push( rootTriangle->getEdge(0) );
	frontier.push( rootTriangle->getEdge(2) );

	boundary_.clear();
	while( !frontier.empty() )
	{
		Edge * edge = frontier.top();
		frontier.pop();

		if( edge->getAdjacent() == NULL )
		{
			boundary_.push_back(edge);
		} else
		{
			frontier.push( edge->getAdjacent()->getPrevious() );
			frontier.push( edge->getAdjacent()->getNext() );
		}
	}

	BOOST_ASSERT( boundary_.size() % 2 == 0 );

}

void CMinusTwoBuilder::MatchingToGenusOneMatching()
{

}

void CMinusTwoBuilder::RandomBoundaryMatching()
{
	// Use algorithm W from section 7.2.1.6 of Knuth TAOCP to generate
	// a random planar matching of the boundary edges of a polygon.

	int n = (n_triangles_+2)/2;
	int p=n,q=n,m=1;

	matching_.clear();
	std::stack<int> levels;

	while( q>0 )
	{
		int x = triangulation_->RandomInteger(0,(q+p)*(q-p+1)-1);
		if( x < (q+1)*(q-p) )
		{
			q--;
			matching_.push_back(std::pair<int,int>(levels.top()-1,m-1));
			levels.pop();
			m++;
		}else
		{
			p--;
			levels.push(m);
			m++;
		}
	}

}

void CMinusTwoBuilder::ApplyBoundaryMatching(const std::vector<std::pair<int,int> > & matching)
{
	for(std::vector<std::pair<int,int> >::const_iterator it = matching.begin(); it != matching.end(); ++it)
	{
		boundary_[it->first]->bindAdjacent(boundary_[it->second]);
	}
}

double CMinusTwoBuilder::CentralCharge() const
{
	return -2.0;
}

std::string CMinusTwoBuilder::ConfigurationData() const
{
	std::ostringstream stream;
	stream << std::fixed << "{type -> \"cminustwobuilder\", centralcharge -> -2 }";
	return stream.str();
}

void CMinusTwoBuilder::getSpanningTree(std::vector<boost::array<bool,3> > & intree) const
{
	if( intree.size() != n_triangles_ )
	{
		intree.resize( n_triangles_ );
	}
	boost::array<bool,3> allfalse = {false,false,false};
	std::fill(intree.begin(),intree.end(),allfalse);
	BOOST_FOREACH(Edge * edge, boundary_)
	{
		intree[edge->getParent()->getId()][edge->getId()] = true;
	}
}

