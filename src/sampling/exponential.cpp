#include "../sampling.h"

void sampleExponential(SamplingData samplingData, double u, Vec& sample, double& pdf) {
	// remap u to account for finite max distance
	double minU = exp(-(samplingData.medium->sigma_t.x) * samplingData.distance_line);
	double a = u * (1.0 - minU) + minU;

	// sample with pdf proportional to exp(-sig*d)
	double x = -log(a) / (samplingData.medium->sigma_t.x);
	pdf = (samplingData.medium->sigma_t.x) * a / (1.0 - minU);

	sample = samplingData.lineStart + samplingData.direction_line * x;
}