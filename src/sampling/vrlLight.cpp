#include "../sampling.h"

// sampling strategy from Virtual Ray Lights
void sampleVRLLight(SamplingData samplingData, double u, Vec& sample, double& pdf) {
	if (samplingData.theta_directions == 0.0) {
		pdf = 1.0;
		sample = samplingData.closestPoint;
	}
	double theta_o = asinh(-samplingData.dist_lineOrigin_closestPoint / samplingData.dist_both_lines * sin(samplingData.theta_directions));
	double theta_e = asinh((samplingData.distance_line - samplingData.dist_lineOrigin_closestPoint) / samplingData.dist_both_lines * sin(samplingData.theta_directions));
	double dist_sample_closestPoint = samplingData.dist_both_lines * sinh((1 - u) * theta_o + u * theta_e) / sin(samplingData.theta_directions);
	pdf = sin(samplingData.theta_directions) / ((theta_e - theta_o) * sqrt(samplingData.dist_both_lines * samplingData.dist_both_lines + dist_sample_closestPoint * dist_sample_closestPoint * sin(samplingData.theta_directions) * sin(samplingData.theta_directions)));
	sample = samplingData.closestPoint + samplingData.direction_line * dist_sample_closestPoint;
}