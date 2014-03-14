#include <stack>
#include <queue>

#include "boost/assert.hpp"

#include "CMinusTwoBuilder.h"
#include "CohomologyBasis.h"
#include "BabyUniverseRemover.h"
#include "Triangle.h"
#include "Edge.h"

LogFactorialTable::LogFactorialTable(int max) : max_(max), veryNegative_(-1.0e100)
{
	logfactorial_.reserve(max+1);
	logfactorial_.push_back(0.0);
	logfactorial_.push_back(0.0);
	for(int i=2;i<=max;i++)
	{
		logfactorial_.push_back(logfactorial_.back() + std::log(static_cast<double>(i)));
	}
}

double LogFactorialTable::LogFactorial(int n) const
{
	BOOST_ASSERT( n>=0 && n<=max_ );
	return logfactorial_[n];
}

double LogFactorialTable::LogBinomial(int n, int k) const
{
	BOOST_ASSERT( n >= k && k >= 0 && n <= max_ );
	return LogFactorial(n) - LogFactorial(k) - LogFactorial(n-k);
}

double LogFactorialTable::LogBinomialDifference(int n, int k1, int k2) const
{
	// gives log( binomial(n,k1)-binomial(n,k2) )
	double logb1 = LogBinomial(n,k1);
	if( k2 < 0 )
	{
		return logb1;
	}
	double logb2 = LogBinomial(n,k2);
	return logb1 + std::log(1.0-std::exp(logb2-logb1));
}

double LogFactorialTable::LogNumberWalksMark(int u, int d, int q) const
{
	if( u <= 0 || d <= q )
	{
		return veryNegative_;
	} 
	return LogBinomialDifference(d+u, (d-u>=q ? u-1 : d-q-1), u-q-2);
}

double LogFactorialTable::LogNumberWalks(int u, int d) const
{
	return std::log(static_cast<double>(d-u+1)/static_cast<double>(d+1)) + LogBinomial(d+u,d);
}

CMinusTwoBuilder::CMinusTwoBuilder(Triangulation * const triangulation, int Genus, int NumberOfTriangles, int boundaryLength)
	: triangulation_(triangulation), 
	genus_(Genus), 
	n_triangles_(NumberOfTriangles),
	remove_baby_universes_(false),
	table_( (boundaryLength == 0 ? 1 : NumberOfTriangles+2 ) ),
	boundary_length_(boundaryLength)
{
	BOOST_ASSERT( genus_ == 0 || genus_ == 1 );
	BOOST_ASSERT( n_triangles_ >= 2 && n_triangles_ % 2 == 0 );
	tree_list_.resize(2*n_triangles_ + 1,0);
	make_disk_ = (boundaryLength > 0);
}


CMinusTwoBuilder::~CMinusTwoBuilder(void)
{
}

void CMinusTwoBuilder::DoSweep()
{
	triangulation_->Clear();
	RandomDiskTriangulation();
	if( make_disk_ )
	{
		RandomBoundaryMatchingDisk();
	} else
	{
		RandomBoundaryMatching();
	}
	ApplyBoundaryMatching(matching_);

	if( remove_baby_universes_ && genus_ == 0 && make_disk_ )
	{
		InsertDiskAtBoundary();

	}

	triangulation_->DetermineVertices();

	BOOST_ASSERT( triangulation_->NumberOfTriangles() == 2*triangulation_->NumberOfVertices() - 4 );

	if( genus_ == 1 )
	{
		GenusZeroToGenusOneMatching();
		ApplyBoundaryMatching( matching_ );
		triangulation_->DetermineVertices();
		BOOST_ASSERT( triangulation_->NumberOfTriangles() == 2*triangulation_->NumberOfVertices() );
	}


	if( remove_baby_universes_ )
	{
		if( genus_ == 0 )
		{
			BabyUniverseRemover remover(triangulation_);
			if( make_disk_ )
			{
				// make sure to keep the inserted disk (to which the newest triangle belongs)
				remover.RemoveBabyUniverses(triangulation_->getTriangle(triangulation_->NumberOfTriangles()-1));
			} else
			{
				remover.RemoveBabyUniverses();
			}
		} else
		{
			CohomologyBasis cohom(triangulation_);
			cohom.SetMakeUpToDateViaReinitialization(true);
			cohom.MakeUpToDate();
			BabyUniverseRemover remover(triangulation_,&cohom);
			remover.RemoveBabyUniverses();
		}
	}
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

void CMinusTwoBuilder::GenusZeroToGenusOneMatching()
{
	// select three distinct vertices at random
	boost::array<Vertex*,3> vertices;
	for(int i=0;i<3;i++)
	{
		int v = triangulation_->RandomInteger(0,triangulation_->NumberOfVertices()-1-i);
		if( i == 1 && v >= vertices[0]->getId() )
		{
			v++;
		} else if( i == 2 && v >= std::min(vertices[0]->getId(),vertices[1]->getId()) )
		{
			v++;
			if( v >= std::max(vertices[0]->getId(),vertices[1]->getId()) )
			{
				v++;
			}
		}
		vertices[i] = triangulation_->getVertex(v);
	}
	BOOST_ASSERT( vertices[0] != vertices[1] && vertices[1] != vertices[2] && vertices[2] != vertices[0] );


	// For each edge we need to know what its position in the boundary is (if applicable).
	/*boost::array<int,3> allMinusOne = {-1,-1,-1};
	std::vector<boost::array<int,3> > edgeToBoundary(n_triangles_,allMinusOne);
	for(int i=0,endi=boundary_.size();i<endi;i++)
	{
		edgeToBoundary[boundary_[i]->getParent()->getId()][boundary_[i]->getId()] = i;
	}*/

	// Select the edge entering vertex v[i] which appears first in boundary_.
	boost::array<int,3> minimalEdgeIndices = {-1,-1,-1};
	for(int i=0,endi=boundary_.size();i<endi;i++)
	{
		for(int j=0;j<3;j++)
		{
			if( minimalEdgeIndices[j] == -1 && boundary_[i]->getPrevious()->getOpposite() == vertices[j] )
			{
				minimalEdgeIndices[j] = i;
			}
		}
	}

	// Order according to the position in boundary
	std::sort(minimalEdgeIndices.begin(),minimalEdgeIndices.end());

	// Define a relabeling of the boundary (see Chapuy, "A new combinatorial identity...", lemma 1)
	std::vector<int> label(boundary_.size());
	for(int i=0,endi=boundary_.size();i<endi;i++)
	{
		if(i <= minimalEdgeIndices[0] || i > minimalEdgeIndices[2]) 
			label[i]=i;
		else if(i <= minimalEdgeIndices[1])
			label[i] = i+(minimalEdgeIndices[2]-minimalEdgeIndices[1]);
		else 
			label[i] = i-(minimalEdgeIndices[1]-minimalEdgeIndices[0]);
	}

	for(std::vector<std::pair<int,int> >::iterator it = matching_.begin();it!=matching_.end();it++)
	{
		it->first = label[it->first];
		it->second = label[it->second];
	}
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

int CMinusTwoBuilder::ChooseRandomWithLogWeights(const std::vector<double> & w) const
{
	double max=-1e-10;
	for(int i=0,endi=w.size();i<endi;i++)
	{
		if( w[i] > max )
		{
			max = w[i];
		}
	}
	std::vector<double> expw;
	expw.reserve(w.size());
	double total = 0.0;
	for(int i=0,endi=w.size();i<endi;i++)
	{
		expw.push_back(std::exp(w[i]-max));
		total += expw[i];
	}
	double x = triangulation_->RandomReal(0.0,total);
	total = 0.0;
	for(int i=0,endi=w.size()-1;i<endi;i++)
	{
		total += expw[i];
		if( total > x )
		{
			return i;
		}
	}
	return w.size()-1;
}

void CMinusTwoBuilder::RandomBoundaryMatchingDisk()
{
	int TotalSteps = n_triangles_+2;
	int upStepsLeft = TotalSteps/2;
	int downStepsLeft = TotalSteps/2;
	int markLevel = boundary_length_/2-1;
	bool marked = false;
	boundary_positions_.resize(boundary_length_/2);

	BOOST_ASSERT( TotalSteps > 2*markLevel );

	matching_.clear();
	std::stack<int> levels;

	for( int currentStep = 0; currentStep < TotalSteps; currentStep++)
	{
		bool goDown = false;
		bool markHere = false;
		if( marked )
		{
			int x = triangulation_->RandomInteger(0,(downStepsLeft+upStepsLeft)*(downStepsLeft-upStepsLeft+1)-1);
			goDown = x < (downStepsLeft+1)*(downStepsLeft-upStepsLeft);
		} else if( downStepsLeft != upStepsLeft )
		{
			int m = downStepsLeft+upStepsLeft-1;
			std::vector<double> logWeights;
			logWeights.push_back( table_.LogNumberWalksMark(upStepsLeft-1,downStepsLeft,markLevel) );
			logWeights.push_back( table_.LogNumberWalksMark(upStepsLeft,downStepsLeft-1,markLevel) );
			if( markLevel == downStepsLeft-upStepsLeft )
			{
				logWeights.push_back( table_.LogNumberWalks(upStepsLeft-1,downStepsLeft) );
			}

			int p = ChooseRandomWithLogWeights(logWeights);
			goDown = (p==1);
			markHere = (p==2);
		}
		if( goDown )
		{
			matching_.push_back(std::pair<int,int>(levels.top(),currentStep));
			levels.pop();
			downStepsLeft--;
		}else
		{
			if( !marked && downStepsLeft-upStepsLeft < boundary_length_/2)
			{
				boundary_positions_[downStepsLeft-upStepsLeft] = currentStep;
			}
			if( markHere )
			{
				marked = true;
			}
			levels.push(currentStep);
			upStepsLeft--;
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

void CMinusTwoBuilder::getDiskBoundary(std::list<const Edge *> & boundary) const
{
	boundary.clear();
	if( remove_baby_universes_ )
	{
		for(std::list<const Edge *>::const_iterator it = reverse_boundary_.begin();it!=reverse_boundary_.end();it++)
		{
			boundary.push_front((*it)->getAdjacent());
		}
	} else
	{
		for(int i=0,endi=boundary_positions_.size();i<endi;i++)
		{
			boundary.push_back(boundary_[boundary_positions_[i]]);
		}
		for(int i=boundary_positions_.size()-1;i>=0;i--)
		{
			boundary.push_back(boundary_[boundary_positions_[i]]->getAdjacent());
		}
	}
}
void CMinusTwoBuilder::getDiskBoundary(std::list<Edge *> & boundary) const
{
	boundary.clear();
	for(int i=0,endi=boundary_positions_.size();i<endi;i++)
	{
		boundary.push_back(boundary_[boundary_positions_[i]]);
	}
	for(int i=boundary_positions_.size()-1;i>=0;i--)
	{
		boundary.push_back(boundary_[boundary_positions_[i]]->getAdjacent());
	}
}

double CMinusTwoBuilder::CentralCharge() const
{
	return -2.0;
}

std::string CMinusTwoBuilder::ConfigurationData() const
{
	std::ostringstream stream;
	stream << std::fixed << "{type -> \"cminustwobuilder\", centralcharge -> -2, numberoftriangles -> " << n_triangles_;
	if( boundary_length_ > 0 )
	{
		stream << ", boundarylength -> " << boundary_length_;
	}
	stream << "}";
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

void CMinusTwoBuilder::setRemoveBabyUniverses(bool remove)
{
	remove_baby_universes_ = remove;
}

void CMinusTwoBuilder::getPath(const Vertex * v0, const Vertex * v1, std::list<const Edge *> & path) const
{
	std::vector<std::list<const Edge *> > vertToEdge(triangulation_->NumberOfVertices());
	BOOST_FOREACH(const Edge * edge, boundary_)
	{
		vertToEdge[edge->getNext()->getOpposite()->getId()].push_back(edge);
	}
	std::queue<const Vertex *> q;
	q.push(v0);
	std::vector<const Edge *> pathEdge(triangulation_->NumberOfVertices(),NULL);
	bool found = false;
	while(!found )
	{
		const Vertex * v = q.front();
		q.pop();
		BOOST_FOREACH(const Edge * edge, vertToEdge[v->getId()])
		{
			Vertex * nextv = edge->getPrevious()->getOpposite();
			if( nextv != v0 && pathEdge[nextv->getId()] == NULL )
			{
				pathEdge[nextv->getId()] = edge;
				q.push(nextv);
				if( nextv == v1 )
				{
					found = true;
					break;
				}
			}
		}
	}

	path.clear();
	const Vertex * v = v1;
	while( v != v0 )
	{
		path.push_front(pathEdge[v->getId()]);
		v = pathEdge[v->getId()]->getNext()->getOpposite();
	}
}

void CMinusTwoBuilder::InsertDiskAtBoundary()
{
	std::list<Edge *> boundary;
	getDiskBoundary(boundary);
	std::vector<Triangle *> triangles;
	reverse_boundary_.clear();
	for(int i=0,endi=boundary.size();i<endi;i++)
	{
		triangles.push_back(triangulation_->NewTriangle());
	}
	int i=0;
	for(std::list<Edge *>::iterator it=boundary.begin();it!=boundary.end();it++,i++)
	{
		triangles[i]->getEdge(0)->bindAdjacent(*it);
		triangles[i]->getEdge(2)->bindAdjacent(triangles[(i+1)%triangles.size()]->getEdge(1));
		reverse_boundary_.push_front(triangles[i]->getEdge(0));
	}
}