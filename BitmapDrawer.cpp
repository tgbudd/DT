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

DiscreteColorScheme::DiscreteColorScheme( DiscreteColorScheme::Scheme scheme ) : scheme_(scheme)
{
	if( scheme_ == COLORS8 )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.611765, 0.298039, 0.294118}, {1., 0.854902, 0.827451}, {0.858824, 0.560784, 0.458824}, {0.654902, 0.298039, 0.0823529}, {0.870588, 0.92549, 0.215686}, {0.501961, 0.686275, 0.0352941}, {0.215686, 0.364706, 0.164706}, {0.243137, 0.317647, 0.309804}, {0.686275, 0.192157, 0.247059}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	} else if( scheme_ == COLORS2 )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.901176, 0.30549, 0.30549}, {1., 0.486667, 0.3}, {1., 0.618431, 0.478431}, {0.637647, 0.431765, 0.324706}, {1., 0.914902, 0.654118}, {0.791373, 0.821569, 0.876471}, {0.44549, 0.467451, 0.525098}, {0.360392, 0.398824, 0.500392}, {0.536078, 0.582745, 0.736471}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	} else if( scheme_ == COLORS3 )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.3, 0.3, 0.3}, {0.997255, 0.552549, 0.319216}, {0.997255, 0.991765, 0.324706}, {0.678824, 0.799608, 0.319216}, {0.401569, 0.604706, 0.56902}, {0.30549, 0.656863, 0.950588}, {0.407059, 0.379608, 0.643137}, {0.629412, 0.483922, 0.70902}, {0.923137, 0.308235, 0.643137}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	} else if( scheme_ == COLORS6 )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.736471, 0.736471, 0.736471}, {0.857255, 0.642353, 0.617255}, {0.676863, 0.642353, 0.625098}, {0.916863, 0.893333, 0.844706}, {0.924706, 0.838431, 0.61098}, {0.642353, 0.701961, 0.664314}, {0.650196, 0.723922, 0.77098}, {0.631373, 0.658039, 0.727059}, {0.863529, 0.841569, 0.880784}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	} else if( scheme_ == COLORS10 )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.818824, 0.409412, 0.4}, {0.952941, 0.696471, 0.658824}, {0.962353, 0.776471, 0.501176}, {0.995294, 0.889412, 0.694118}, {0.835294, 0.88, 0.442353}, {0.590588, 0.694118, 0.447059}, {0.503529, 0.616471, 0.442353}, {0.616471, 0.644706, 0.72}, {0.536471, 0.543529, 0.670588}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	} else if( scheme_ == CMYK )
	{
		schemedata_.resize(9);
		double data[9][3] = {{0.510508, 0.776344, 0.931191}, {0.608986, 0.786259, 0.929747}, {0.723023, 0.662397, 0.846162}, {0.817232, 0.597492, 0.77304}, {0.890185, 0.630436, 0.712476}, {0.937868, 0.74166, 0.663667}, {0.961741, 0.852738, 0.637882}, {0.940713, 0.894623, 0.658068}, {0.864761, 0.849713, 0.696532}};
		for(int i=0;i<9;i++)
		{
			for(int j=0;j<3;j++)
			{
				schemedata_[i][j] = data[i][j];
			}
		}
	}
}

boost::array<unsigned char,3> DiscreteColorScheme::getColor( int i ) const
{
	boost::array<unsigned char,3> color;
	int index = ( i%schemedata_.size() + schemedata_.size())%schemedata_.size();
	for(int j=0;j<3;j++)
	{
		color[j] = static_cast<unsigned char>( 255.9 * schemedata_[index][j] );
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

void TriangulationDrawer::DrawShading(BitmapDrawer & drawer, ColorScheme::Scheme colorscheme)
{
	ColorScheme scheme(colorscheme);

	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		std::vector<Vector2D> polygon;
		polygon.push_back( embedding_->getCoordinate(triangle->getEdge(0)->getOpposite()) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,2)) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,0)) );

		boost::array<unsigned char,3> color;
		if( triangle_shade_.empty() )
		{
			double area = SignedPolygonArea( polygon );
			BOOST_ASSERT( area > -1.0e-8 );
			double minlogarea = std::max(0.0,std::min(8.0,-log(area)));
			boost::array<unsigned char,3> color = scheme.getColor( (minlogarea+1.0) / 10.0 );
		} else
		{
			color = scheme.getColor( triangle_shade_[i] );
		}
		drawer.setPenColor(color[0],color[1],color[2]);
		drawer.domainPolygon( polygon );
	}
}

void TriangulationDrawer::DrawShading(BitmapDrawer & drawer, DiscreteColorScheme::Scheme colorscheme)
{
	DiscreteColorScheme scheme(colorscheme);

	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Triangle * triangle = triangulation_->getTriangle(i);
		std::vector<Vector2D> polygon;
		polygon.push_back( embedding_->getCoordinate(triangle->getEdge(0)->getOpposite()) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,2)) );
		polygon.push_back( AddVectors2D(polygon.back(),embedding_->getForm(i,0)) );

		boost::array<unsigned char,3> color;
		if( !triangle_color_index_.empty() )
		{
			color = scheme.getColor( triangle_color_index_[i] );
		}
		drawer.setPenColor(color[0],color[1],color[2]);
		drawer.domainPolygon( polygon );
	}
}

void TriangulationDrawer::SetEdgeShade(int triangle, int edge, double shade)
{
	if( edge_shade_.empty() )
	{
		boost::array<double,3> zero = {0.0,0.0,0.0};
		edge_shade_.resize(triangulation_->NumberOfTriangles(),zero);
	}
	edge_shade_[triangle][edge] = shade;
}
void TriangulationDrawer::SetTriangleShade(int triangle, double shade)
{
	if( triangle_shade_.empty() )
	{
		triangle_shade_.resize(triangulation_->NumberOfTriangles(),0.0);
	}
	triangle_shade_[triangle] = shade;
}
void TriangulationDrawer::SetTriangleColorIndex(int triangle, int index)
{
	if( triangle_shade_.empty() )
	{
		triangle_color_index_.resize(triangulation_->NumberOfTriangles(),0);
	}
	triangle_color_index_[triangle] = index;
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

void SpanningTreeDrawer::Draw(BitmapDrawer & drawer)
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		for(int j=0;j<3;j++)
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(j);
			if( spanningtree_->InSpanningTree(edge) )
			{
				Vector2D start = embedding_->getCoordinate(edge->getNext()->getOpposite());
				Vector2D form = embedding_->getForm(edge);
				drawer.domainLineSegment(start[0],start[1],start[0] + form[0],start[1] + form[1]);
			}
		}
	}
}

