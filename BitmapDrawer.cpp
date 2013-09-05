#include "BitmapDrawer.h"

void TriangulationDrawer::Draw(BitmapDrawer & drawer)
{
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

			drawer.domainLineSegment(start[0],start[1],start[0] - form[0],start[1] - form[1]);
		}
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
		drawer.domainLineSegment(start[0],start[1],start[0] - form[0],start[1] - form[1]);
	}
}