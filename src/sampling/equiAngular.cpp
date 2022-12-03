#include "../sampling.h"

// equi-angular distance sampling strategy from Kulla et al. (Importance sampling of area lights in participating media, 2011)
void sampleEquiAngular(SamplingData samplingData, double u, Vec& sample, double& pdf) {
	double theta_o = atan(-samplingData.dist_lineOrigin_closestPoint / samplingData.dist_interestPoint_closestPoint);
	double theta_e = atan((samplingData.distance_line - samplingData.dist_lineOrigin_closestPoint) / samplingData.dist_interestPoint_closestPoint);
	double dist_sample_closestPoint = samplingData.dist_interestPoint_closestPoint * tan((1 - u) * theta_o + u * theta_e);
	pdf = samplingData.dist_interestPoint_closestPoint / ((theta_e - theta_o) * (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	sample = samplingData.closestPoint + samplingData.direction_line * dist_sample_closestPoint;
}