#pragma once

#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <map>

#include "triangulation.h"
#include "Decoration.h"
#include "Embedding.h"
#include "ShortestLoop.h"
#include "bitmap_image.h"

class BitmapDrawer
{
public:
	BitmapDrawer(int width, int height, int antialiasing = 1) : width_(width), height_(height), 
		fullwidth_(width * antialiasing), fullheight_(height * antialiasing),
		antialiasing_(antialiasing),
		image_(antialiasing * width, antialiasing * height), draw_(image_) 
	{
		BOOST_ASSERT( antialiasing == 1 || antialiasing == 2 || antialiasing == 4 || antialiasing == 8 );
		// start with a white background
		image_.clear(255);	
	}

	void Clear()
	{
		image_.clear(255);
	}
	void SetPeriodicDomain(const std::pair<double,double> & modulus, double scale )
	{
		SetPeriodicDomain(modulus,scale,fullwidth_/2,fullheight_/2);
	}
	void SetPeriodicDomain(const std::pair<double,double> & modulus, double scale, int offsetX, int offsetY )
	{
		modulus_ = modulus;
		scale_ = scale;
		offsetX_ = offsetX;
		offsetY_ = offsetY;
		absolutescale_ = (int)(scale_ * fullwidth_);
		draw_.periodicity = true;
		draw_.periodicboundary = false;
		draw_.period[0][0] = absolutescale_;
		draw_.period[0][1] = 0;
		draw_.period[1][0] = (int)(modulus_.first * scale_ * fullwidth_);
		draw_.period[1][1] = (int)(modulus_.second * scale_ * fullwidth_);
	}
	~BitmapDrawer() {}

	double getScale() const {
		return scale_;
	}
	void setScale(const double & scale) {
		scale_ = scale;
	}
	int getWidth() const {
		return width_;
	}
	void setWidth(const int & width) {
		width_ = width;
	}
	int getHeight() const {
		return height_;
	}
	void setHeight(const int & height) {
		height_ = height;
	}
	void setPixel(int x, int y, int r, int g, int b)
	{
		image_.set_pixel(x,y,r,g,b);
	}
	void setPenWidth(int w)
	{
		draw_.pen_width(w);
	}
	void setPenColor(int r, int g, int b)
	{
		draw_.pen_color(r,g,b);
	}
	void domainLineSegment(double x1, double y1, double x2, double y2)
	{
		LineSegment(Transform(x1,y1),Transform(x2,y2));
	}
	void domainPoint(double x, double y)
	{
		Point(Transform(x,y));
	}
	void domainPolygon(const std::vector<Vector2D> & polygon)
	{
		std::vector<std::pair<int,int> > transPolygon;
		for(int i=0;i<polygon.size();i++)
		{
			transPolygon.push_back(Transform(polygon[i][0],polygon[i][1]));
		}
		Polygon(transPolygon);
	}
	void LineSegment(std::pair<int,int> p1, std::pair<int,int> p2)
	{
		draw_.line_segment(p1.first,p1.second,p2.first,p2.second);
	}
	void Point(std::pair<int,int> p)
	{
		draw_.plot_pen_pixel(p.first,p.second);
	}
	void Polygon(const std::vector<std::pair<int,int> > & polygon)
	{
		std::vector<std::pair<double,double> > dPolygon;
		for(int i=0;i<polygon.size();i++)
		{
			dPolygon.push_back(std::pair<double,double>(polygon[i].first,polygon[i].second));
		}
		draw_.polygon(dPolygon);
	}
	void TransparentBox(int x1, int y1, int x2, int y2, double opacity)
	{
		for(int x=x1;x<=x2;x++)
		{
			for(int y=y1;y<=y2;y++)
			{
				unsigned char r,g,b;
				image_.get_pixel(x,y,r,g,b);
				image_.set_pixel(x,y,r+ (255-r)*opacity,g+(255-g)*opacity,b+(255-b)*opacity);
			}
		}
	}
	void SaveImage(std::string filename)
	{
		if( antialiasing_ > 1 )
		{
			int subpixels = 2;
			bitmap_image image1;
			image_.subsample(image1);
			while( subpixels < antialiasing_ )
			{
				bitmap_image image2;
				image1.subsample(image2);
				image1 = image2;
				subpixels *= 2;
			}
			image1.save_image(filename.c_str());
		} else
		{
			image_.save_image(filename.c_str());
		}
	}

private:

	std::pair<int,int> Transform(double x, double y) const
	{
		return std::pair<int,int>((int)(scale_ * fullwidth_ * ((x-0.5) + modulus_.first * (y-0.5))) + offsetX_,(int)(scale_ * fullwidth_ * modulus_.second * (y-0.5)) + offsetY_);
	}

	

	bitmap_image image_;
	image_drawer draw_;
	std::pair<double,double> modulus_;
	double scale_;		// the width of the domain compared to the width of the bitmap
	int absolutescale_;
	int offsetX_, offsetY_;
	int width_, height_;
	int fullwidth_, fullheight_;
	int antialiasing_;


};

class ComponentDrawer
{
public:
	ComponentDrawer() {}
	~ComponentDrawer() {}

	virtual void Draw(BitmapDrawer & drawer) = 0;
};

class TriangulationDrawer : public ComponentDrawer
{
public:
	TriangulationDrawer(Triangulation * const triangulation, const Embedding * const embedding) 
		: triangulation_(triangulation), embedding_(embedding) {}

	void Draw(BitmapDrawer & drawer);
private:
	Triangulation * const triangulation_;
	const Embedding * const embedding_;
};


class ShortestLoopDrawer : public ComponentDrawer
{
public:
	ShortestLoopDrawer(const ShortestLoop * const shortestloop, const Embedding * const embedding) 
		: shortestloop_(shortestloop), embedding_(embedding) {}

	void Draw(BitmapDrawer & drawer);
	void DrawGenerators(BitmapDrawer & drawer);
	void DrawPath(BitmapDrawer & drawer, const std::list<Edge*> & path );
private:
	const ShortestLoop * const shortestloop_;
	const Embedding * const embedding_;
};

class FundamentalDomainDrawer : public ComponentDrawer
{
public:
	FundamentalDomainDrawer() {}
	void Draw(BitmapDrawer & drawer)
	{
		drawer.domainLineSegment(0.0,0.0,1.0,0.0);
		drawer.domainLineSegment(0.0,0.0,0.0,1.0);
	}
};

class HotEdgesDrawer : public ComponentDrawer, public Decoration
{
public:
	HotEdgesDrawer(const Embedding * const embedding) : embedding_(embedding) {}

	void Initialize()
	{
		hotedges_.clear();
	}
	void UpdateAfterFlipMove(const Edge * const edge)
	{
		hotedges_.push_front(edge->getPrevious());
		if( hotedges_.size() > numberofedges_ )
		{
			hotedges_.pop_back();
		}
	}
	void setNumberOfEdges(int n)
	{
		numberofedges_ = n;
	}
	void setColor(double r, double g, double b)
	{
		r_ = r;
		g_ = g;
		b_ = b;
	}
	void Draw(BitmapDrawer & drawer)
	{
		int index = hotedges_.size()-1;
		for(std::list<const Edge *>::reverse_iterator edgeIt=hotedges_.rbegin();edgeIt!=hotedges_.rend();edgeIt++)
		{
			drawer.setPenColor( (int)(256.0*(r_+(index/((double)numberofedges_)*(1-r_)))),
								(int)(256.0*(g_+(index/((double)numberofedges_)*(1-g_)))),
								(int)(256.0*(b_+(index/((double)numberofedges_)*(1-b_)))) );

			std::vector<Vector2D> polygon;
			Vector2D startPoint = embedding_->getCoordinate((*edgeIt)->getNext()->getOpposite());
			polygon.push_back(startPoint);
			polygon.push_back(AddVectors2D(startPoint,embedding_->getForm((*edgeIt)->getAdjacent()->getNext())));
			polygon.push_back(AddVectors2D(startPoint,embedding_->getForm((*edgeIt))));
			polygon.push_back(SubtractVectors2D(startPoint,embedding_->getForm((*edgeIt)->getPrevious())));
			drawer.domainPolygon(polygon);
			index--;	
		}
	}
private:
	std::list<const Edge *> hotedges_;
	int numberofedges_;
	double r_, g_, b_;
	const Embedding * const embedding_;
};

class TextDrawer {
public:
	TextDrawer(std::string fontfile)
	{
		LoadFont(fontfile);
	}
	TextDrawer() {}
	~TextDrawer() {}
	
	void LoadFont(std::string fontfile)
	{
		std::ostringstream osbitmap, osdata;
		osbitmap << fontfile << ".bmp";
		osdata << fontfile << ".dat";
		font_bitmap_ = bitmap_image(osbitmap.str());
		std::ifstream datafile(osdata.str());
		while( !datafile.eof() )
		{
			int code;
			Character character;
			datafile >> code >> character.x >> character.y >> character.width >> character.height;
			characters_.insert(std::pair<char,Character>((char)code,character));
		}
	}
	void DrawText(const std::string & text, BitmapDrawer & drawer, int startX, int startY, double backgroundOpacity)
	{
		if( backgroundOpacity > 0.0001 )
		{
			std::pair<int,int> size = CalculateSize(text);
			drawer.TransparentBox(startX - 12, startY - 12, startX + size.first + 12, startY + size.second + 12, backgroundOpacity );
		}

		int x = startX;
		int y = startY;
		for(std::string::const_iterator charIt = text.begin(); charIt != text.end(); charIt ++)
		{
			if( *charIt == '\n' )
			{
				y += characters_['A'].height;
				x = startX;
			} else
			{
				std::map<char,Character>::iterator mapIt = characters_.find(*charIt);
				if( mapIt != characters_.end() )
				{
					DrawCharacter(mapIt->second,x,y,drawer);
					x += mapIt->second.width;
				}
			}
		}
	}
	void setColor(unsigned char r, unsigned char g, unsigned char b)
	{
		red_=r;
		green_=g;
		blue_=b;
	}
private:
	struct Character {
		int x, y, width, height;
	};
	void DrawCharacter(const Character & character, int x, int y, BitmapDrawer & drawer)
	{
		for(int y0=0;y0<character.height;y0++)
		{
			for(int x0=0;x0<character.width;x0++)
			{
				unsigned char red, green, blue;
				font_bitmap_.get_pixel(x0+character.x,y0+character.y,red,green,blue);
				if( red != 255 || green != 255 || blue != 255 )
				{
					drawer.setPixel(x0+x,y0+y,red_,green_,blue_);
				}
			}
		}
	}
	std::pair<int,int> CalculateSize(const std::string & text)
	{
		int x = 0;
		int y = 0;
		int width = 0;
		for(std::string::const_iterator charIt = text.begin(); charIt != text.end(); charIt ++)
		{
			if( *charIt == '\n' )
			{
				y += characters_['A'].height;
				x = 0;
			} else
			{
				std::map<char,Character>::iterator mapIt = characters_.find(*charIt);
				if( mapIt != characters_.end() )
				{
					x += mapIt->second.width;
					if( x > width )
					{
						width = x;
					}
				}
			}
		}
		y += characters_['A'].height;
		return std::pair<int,int>(width,y);
	}
	unsigned char red_, green_, blue_;
	std::map<char,Character> characters_;
	bitmap_image font_bitmap_;
};