#ifndef CONFORMAL_DISTRIBUTION_H
#define CONFORMAL_DISTRIBUTION_H

#include "Observable.h"
#include "Embedding.h"
#include "utilities.h"

class ConformalDistribution :
	public Observable
{
public:
	ConformalDistribution(Embedding * const embedding, const Triangulation * const triangulation=NULL);
	~ConformalDistribution(void);

	void Measure();
	std::string OutputData() const;
	void setMethodString(const std::string & method);
	void setSaveEmbedding(bool save, std::string file);
	void setMeasureEuclideanBallSize(bool measure);
	void setRecordBallSizesExponentially(int maxsize, double factor);
	void setMaxBallTriangles(int n);
private:
	void MeasureRadiusDistribution();
	void MeasureEuclideanBallSize();
	void MeasureEuclideanBallSize(Triangle * triangle);
	int measurements_;
	Embedding * const embedding_;
	const Triangulation * const triangulation_;

	struct vertexNode
	{
		Vertex * vertex;
		double square_distance;
		Vector2D position;
		bool operator<(const vertexNode & node) const { return square_distance > node.square_distance; }
	};

	std::pair<double,double> modulus_;
	ReusableFlag triangle_visited_, vertex_visited_;
	std::vector<Vector2D> vertex_position_;
	std::vector<double> vertex_distance_;
	bool measure_euclidean_ball_size_;
	unsigned int max_ball_triangles_;
	std::vector<int> record_ball_sizes_;
	std::vector<std::vector<int> > ball_histogram_;
	double ball_minlograd_, ball_maxlograd_;
	int ball_bins_;
	std::vector<std::vector<int> > point_histogram_;



	bool save_embedding_;
	std::string embedding_file_;
	void SaveEmbedding();

	std::vector<int> histogram_;
	double minlograd_, maxlograd_;
	int bins_;
	std::string method_;
};

#endif