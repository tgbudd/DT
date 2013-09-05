#pragma once

#include <string>

#include "triangulation.h"
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
	void LineSegment(std::pair<int,int> p1, std::pair<int,int> p2)
	{
		draw_.line_segment(p1.first,p1.second,p2.first,p2.second);
	}
	void Point(std::pair<int,int> p)
	{
		draw_.plot_pen_pixel(p.first,p.second);
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
