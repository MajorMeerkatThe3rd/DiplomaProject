#include "../sampling.h"
#include "../helperFunctions.h"

const int M = 10;
struct PiecewiseSegments {
	double k[M - 1], d[M - 1];
	double cdf[M - 1];
};

double computeSegment(double theta_s, double theta_e, double fs_s, double fs_e, double& k, double& d, double& cdf) {
	k = (fs_e - fs_s) / (theta_e - theta_s);
	d = fs_e - k * theta_e;
	double integral_end = (k * theta_e * theta_e) / 2.0 + d * theta_e;
	double integral_start = (k * theta_s * theta_s) / 2.0 + d * theta_s;
	cdf = integral_end - integral_start;
	return cdf;
}

// advance sampling proportional to the product of the phase functions and the inverse squared distance
void sampleVRLAnisotropic(SamplingData samplingData, double u, Vec& sample, double& pdf) {
	double theta[M], fs[M];
	PiecewiseSegments segments;

	Line cameraRay = Line(samplingData.lineEnd, samplingData.direction_line * -1, samplingData.lineStart, 0, samplingData.distance_line);

	// construct all variables to find the peak
	Vec a = (cameraRay.origin - samplingData.interestPoint).norm();
	Vec b = (cameraRay.end - samplingData.interestPoint).norm();
	Vec c = (a % b).norm();
	// the projection of the light ray on the plane that is spanned by a and b
	Vec e = ((c % samplingData.direction_interestPoint) % c).norm();

	samplingData.dist_lineOrigin_closestPoint = cameraRay.distance - samplingData.dist_lineOrigin_closestPoint;

	// evaluate the angle to the start and the product of the phase functions with this angle
	double theta_1 = atan(-samplingData.dist_lineOrigin_closestPoint / samplingData.dist_interestPoint_closestPoint);
	theta[0] = theta_1;
	fs[0] = pfHenyeyGreenstein(cameraRay.direction, a, samplingData.medium->g) * pfHenyeyGreenstein(samplingData.direction_interestPoint, a, samplingData.medium->g);

	// evaluate the angle to the end and the product of the phase functions with this angle
	double theta_M = atan((cameraRay.distance - samplingData.dist_lineOrigin_closestPoint) / samplingData.dist_interestPoint_closestPoint);
	theta[M - 1] = theta_M;
	fs[M - 1] = pfHenyeyGreenstein(cameraRay.direction, b, samplingData.medium->g) * pfHenyeyGreenstein(samplingData.direction_interestPoint, b, samplingData.medium->g);

	double theta_peak;
	bool inRange = false;

	// check if the peak lies on the arc
	if ((a % e).dot((a % b)) >= 0 && (b % e).dot((b % a)) >= 0) {
		// calculate the angle to the peak
		theta_peak = arccos(a.dot(e)) + theta_1;
		inRange = true;
	}
	else {
		// if the peak is not on the arc, search for the negative peak
		e = e * -1;
		if ((a % e).dot((a % b)) >= 0 && (b % e).dot((b % a)) >= 0) {
			// calculate the angle to the peak
			theta_peak = arccos(a.dot(e)) + theta_1;
			inRange = true;
		}
	}
	int j_peak;
	double cdf = 0.0;

	// case separation for if the peak lies within the arc
	if (inRange) {
		// find the index of the peak
		// this formula form the paper has been adjusted to fit the current indexing
		j_peak = floor((theta_peak - theta_1) / (theta_M - theta_1) * (M - 1) + .5);
		j_peak = j_peak < 2 ? 2 : j_peak;

		// calculate the first sub interval
		for (int j = 2; j <= j_peak; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_peak - theta_1) / 2.0 * (1 - cos((PI * (j - 1)) / (j_peak - 1)));
			theta[j - 1] = theta_j + theta_1;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, a, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(cameraRay.direction, dir, samplingData.medium->g) * pfHenyeyGreenstein(samplingData.direction_interestPoint, dir, samplingData.medium->g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}

		// calculate the second sub interval
		for (int j = j_peak + 1; j <= M - 1; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_M - theta_peak) / 2.0 * (1 - cos((PI * (j - j_peak)) / (M - j_peak)));
			theta[j - 1] = theta_j + theta_peak;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, e, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(cameraRay.direction, dir, samplingData.medium->g) * pfHenyeyGreenstein(samplingData.direction_interestPoint, dir, samplingData.medium->g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}
		// add the last function from the last sampled angle to the end
		cdf += computeSegment(theta[M - 2], theta[M - 1], fs[M - 2], fs[M - 1], segments.k[M - 2], segments.d[M - 2], segments.cdf[M - 2]);

	}
	else {
		// the peak was not in the interval, therefore distribute the angle in just one interval
		for (int j = 2; j <= M - 1; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_M - theta_1) / 2.0 * (1 - cos((PI * (j - 1)) / (M - 1)));
			theta[j - 1] = theta_j + theta_1;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, a, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(cameraRay.direction, dir, samplingData.medium->g) * pfHenyeyGreenstein(samplingData.direction_interestPoint, dir, samplingData.medium->g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}
		// add the last function from the last sampled angle to the end
		cdf += computeSegment(theta[M - 2], theta[M - 1], fs[M - 2], fs[M - 1], segments.k[M - 2], segments.d[M - 2], segments.cdf[M - 2]);
	}

	// adjist the xi by the cdf of the piecewise function to be able to find the right function for the evalutaion
	double xi_adj = u * cdf;
	double theta_sample = -1.0;
	double cdf_comb = 0.0;

	// loop through all piecewise segments
	for (int i = 0; i < M - 1; i++) {
		// see if the sampled value is in the interval of the cdf for this current segment
		if (cdf_comb < xi_adj && xi_adj <= cdf_comb + segments.cdf[i]) {
			// for a segment that is constant (d = 0) the formula can be simplified
			// this is always the case in isotropic media, but can also occur in anisotropic media
			if (-eps < segments.k[i] && segments.k[i] < eps) {
				// get the sampled angle with the inverse cdf
				theta_sample = (theta[i] * segments.d[i] + xi_adj - cdf_comb) / segments.d[i];
				// calcualte the cdf
				pdf = segments.d[i] / cdf;
			}
			else {
				// get the sampled angle with the inverse cdf
				double sqrt_term = segments.k[i] * (theta[i] * theta[i] * segments.k[i] + 2.0 * xi_adj - 2.0 * cdf_comb) + 2.0 * theta[i] * segments.d[i] * segments.k[i] + segments.d[i] * segments.d[i];
				theta_sample = (sqrt(sqrt_term) - segments.d[i]) / segments.k[i];
				// calcualte the cdf
				pdf = (segments.k[i] * theta_sample + segments.d[i]) / cdf;
			}
			break;
		}
		cdf_comb += segments.cdf[i];
	}

	// calculate the actual sample
	double dist_sample_closestPoint = samplingData.dist_interestPoint_closestPoint * tan(theta_sample);
	sample = samplingData.closestPoint + cameraRay.direction * dist_sample_closestPoint;
	//sample = theta_sample;

	// in a case where the pdf is exactly 1, the pdf is equivalent to the equi-angular case
	if (1.0 - eps < pdf && pdf < 1.0 + eps) {
		pdf = samplingData.dist_interestPoint_closestPoint / ((theta_M - theta_1) * (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	}
	else {
		// the pdf is multiplied by the already adjusted part for the inverse squared distance
		pdf = pdf * (samplingData.dist_interestPoint_closestPoint / (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	}
}