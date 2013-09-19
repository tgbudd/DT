#include "BitmapDrawer.h"

ColorScheme::ColorScheme( Scheme scheme ) : scheme_(scheme)
{
	if( scheme_ == TEMPERATURE_MAP )
	{
		double data[13][4] = {{0., 0.178927, 0.305394, 0.933501}, {0.0833333, 0.308746, 0.441842, 0.940894}, {0.166667, 0.453318, 0.567063, 0.950106}, {0.25, 0.642359, 0.720535, 0.964988}, {0.333333, 0.819984, 0.859297, 0.982692}, {0.416667, 0.935699, 0.951565, 0.993729}, {0.5, 0.984192, 0.987731, 0.911643}, {0.583333, 0.995282, 0.992317, 0.727853}, {0.666667, 0.992503, 0.986373, 0.425376}, {0.75, 0.955963, 0.863115, 0.283425}, {0.833333, 0.904227, 0.657999, 0.241797}, {0.916667, 0.858405, 0.449932, 0.203562}, {1., 0.817319, 0.134127, 0.164218}};
		schemedata_.resize(13);
		for(int i=0;i<13;i++)
		{
			schemedata_[i].first = data[i][0];
			schemedata_[i].second[0] = data[i][1];
			schemedata_[i].second[1] = data[i][2];
			schemedata_[i].second[2] = data[i][3];
		}
	} else if( scheme_ == BLUE_GREEN_YELLOW )
	{
		double data[7][4] = {{0., 0.122103, 0.00901808, 0.39826}, {0.166667, 0.0839935, 0.279645, 0.510102}, {0.333333, 0.097699, 0.498132, 0.548165}, {0.5, 0.175507, 0.652273, 0.528496}, {0.666667, 0.329526, 0.762208, 0.474596}, {0.833333, 0.571909, 0.839991, 0.408102}, {1., 0.914809, 0.897673, 0.350652}};
		schemedata_.resize(7);
		for(int i=0;i<7;i++)
		{
			schemedata_[i].first = data[i][0];
			schemedata_[i].second[0] = data[i][1];
			schemedata_[i].second[1] = data[i][2];
			schemedata_[i].second[2] = data[i][3];
		}
	}
}

boost::array<unsigned char,3> ColorScheme::getColor( double x ) const
{
	boost::array<unsigned char,3> color;
	for(int i=1;i<static_cast<int>(schemedata_.size());i++)
	{
		if( schemedata_[i].first > x )
		{
			for(int j=0;j<3;j++)
			{
				color[j] = (unsigned char)( 256.0 * ( (x-schemedata_[i-1].first) * schemedata_[i].second[j] + (schemedata_[i].first-x) * schemedata_[i-1].second[j] ) / (schemedata_[i].first - schemedata_[i-1].first) );
			}
			break;
		}
	}
	return color;
}

void TriangulationDrawer::Draw(BitmapDrawer & drawer)
{
	ColorScheme scheme(ColorScheme::BLUE_GREEN_YELLOW);

	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangle->getEdge(j);

			if( edge->getAdjacent()->getParent()->getId() < edge->getParent()->getId() )
				continue;

			Vector2D start = embedding_->getCoordinate(edge->getNext()->getOpposite());
			Vector2D form = embedding_->getForm(edge);

			if( !edge_shade_.empty() )
			{
				boost::array<unsigned char,3> color = scheme.getColor( edge_shade_[i][j] );
				drawer.setPenColor(color[0],color[1],color[2]);
			}
			drawer.domainLineSegment(start[0],start[1],start[0] + form[0],start[1] + form[1]);
		}
	}
}

void TriangulationDrawer::DrawShading(BitmapDrawer & drawer)
{
	ColorScheme scheme(ColorScheme::TEMPERATURE_MAP);

	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		std::vector<Vector2D> polygon;
		polygon.push_back( embedding_->getCoordinate(triangle->getEdge(0)->getOpposite()) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,2)) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,0)) );

		double area = SignedPolygonArea( polygon );
		BOOST_ASSERT( area > -1.0e-8 );
		double minlogarea = std::max(0.0,std::min(8.0,-log(area)));
		boost::array<unsigned char,3> color = scheme.getColor( (minlogarea+1.0) / 10.0 );

		drawer.setPenColor(color[0],color[1],color[2]);
		drawer.domainPolygon( polygon );
	}
}

void ShortestLoopDrawer::Draw(BitmapDrawer & drawer)
{
	DrawPath(drawer,shortestloop_->getShortestLoop());
}

void ShortestLoopDrawer::DrawGenerators(BitmapDrawer & drawer)
{
	std::vector<std::list<Edge*> > generators = shortestloop_->getGenerators();
	for( std::vector<std::list<Edge*> >::iterator pathIt = generators.begin(); pathIt != generators.end(); pathIt++)
	{
		DrawPath(drawer,*pathIt);
	}
}	

void ShortestLoopDrawer::DrawPath(BitmapDrawer & drawer, const std::list<Edge*> & path )
{
	for( std::list<Edge*>::const_iterator edgeIt = path.begin(); edgeIt != path.end(); edgeIt++)
	{
		Vector2D start = embedding_->getCoordinate((*edgeIt)->getNext()->getOpposite());
		Vector2D form = embedding_->getForm(*edgeIt);
		drawer.domainLineSegment(start[0],start[1],start[0] + form[0],start[1] + form[1]);
	}
}