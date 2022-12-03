#include "../sampling.h"

void sampleUniform(SamplingData samplingData, double u, Vec& sample, double& pdf) {
	pdf = 1.0 / samplingData.distance_line;
	sample = samplingData.lineStart + samplingData.direction_line * samplingData.distance_line * u;
}