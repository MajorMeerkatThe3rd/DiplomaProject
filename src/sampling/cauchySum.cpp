#include "../sampling.h"
#include "../helperFunctions.h"

namespace cauchySum {
	double s_camera = 0.0, m_camera = 0.0;

	// pdf of Cauchy
	double cauchyPDF(double x, double s, double m, double multiplier = 1) {
		return (1.0 / (PI * s * (1.0 + ((x - m) / s) * ((x - m) / s)))) / multiplier;
	}

	// cdf of Cauchy
	double cauchyCDF(double x, double s, double m, double multiplier = 1, double start = 0) {
		return (0.5 + (1.0 / PI) * (atan((x - m) / s))) / multiplier - start;
	}

	// inverse cdf of product of Cauchy distributions
	double cauchyCDFMixtureInv(double u, double s1, double s2, double m1, double m2) {
		double v_norm = tan((2 * u - 1) * PI);
		double p_norm = (1 / v_norm) * (-m1 * v_norm - m2 * v_norm + s1 + s2);
		double q_norm = (1 / v_norm) * (m1 * m2 * v_norm - s1 * s2 * v_norm - s2 * m1 - s1 * m2);
		return u < 0.5 ? -(p_norm / 2) - sqrt((p_norm / 2) * (p_norm / 2) - q_norm) : u > 0.5 ? -(p_norm / 2) + sqrt((p_norm / 2) * (p_norm / 2) - q_norm) : (m1 * s2 + m2 * s1) / (s1 + s2);
	}

	double computeSWithMax(double max) {
		return 1 / (PI * max);
	}
}

// as the medium does not change in out test cases, we can make these calculations beforehand
void initCauchySum(Scene scene) {
	double cameraMax = scene.medium->g > 0 ? pfHenyeyGreenstein(1, scene.medium->g) : pfHenyeyGreenstein(-1, scene.medium->g);
	cauchySum::s_camera = cauchySum::computeSWithMax(cameraMax);
	cauchySum::m_camera = scene.medium->g > 0 ? PI / 2 : -PI / 2;  // the two rays are always perpendicular, dot is always 0, acos always PI/2
}

void sampleCauchySum(SamplingData samplingData, double random, Vec& sample, double& pdf) {
	Line cameraRay = Line(samplingData.lineEnd, samplingData.direction_line * -1, samplingData.lineStart, 0, samplingData.distance_line);
	samplingData.dist_lineOrigin_closestPoint = cameraRay.distance - samplingData.dist_lineOrigin_closestPoint;

	double theta_o = atan(-samplingData.dist_lineOrigin_closestPoint / samplingData.dist_interestPoint_closestPoint);
	double theta_e = atan((cameraRay.distance - samplingData.dist_lineOrigin_closestPoint) / samplingData.dist_interestPoint_closestPoint);

	// construct all variables to find the peak
	Vec a = (cameraRay.origin - samplingData.interestPoint).norm();
	Vec c = (a % samplingData.direction_interestPoint_closestPoint).norm();
	// the projection of the light ray on the plane that is spanned by a and direction_interestPoint_closestPoint
	Vec e = ((c % samplingData.direction_interestPoint) % c).norm();

	// create the scale for both, the light and the camera functions
	double maxDotLight = samplingData.direction_interestPoint.dot(e);
	double lightMax = samplingData.medium->g > 0 ? pfHenyeyGreenstein(maxDotLight, samplingData.medium->g) : pfHenyeyGreenstein(-maxDotLight, samplingData.medium->g);
	double s_light = cauchySum::computeSWithMax(lightMax);
	double s_camera = cauchySum::s_camera;

	// create the location for both, the light and the camera functions
	double lightDirectionDot = e.dot(samplingData.direction_interestPoint_closestPoint);
	double m_light = samplingData.medium->g > 0 ? -arccos(lightDirectionDot) : PI - arccos(lightDirectionDot);
	m_light = e.dot(cameraRay.direction) > 0.0 ? -m_light : m_light;
	double m_camera = cauchySum::m_camera;

	// calculate the CDFs to adjust the random value
	double cdf_theta_o = 0.5 * cauchySum::cauchyCDF(theta_o, s_light, m_light) + 0.5 * cauchySum::cauchyCDF(theta_o, s_camera, m_camera);
	double cdf_theta_e = 0.5 * cauchySum::cauchyCDF(theta_e, s_light, m_light) + 0.5 * cauchySum::cauchyCDF(theta_e, s_camera, m_camera);
	double u = random * (cdf_theta_e - cdf_theta_o) + cdf_theta_o;

	// sample via the inverse CDF
	double theta_sample = cauchySum::cauchyCDFMixtureInv(u, s_camera, s_light, m_camera, m_light);
	double dist_sample_closestPoint = samplingData.dist_interestPoint_closestPoint * tan(theta_sample);

	sample = samplingData.closestPoint + cameraRay.direction * dist_sample_closestPoint;

	pdf = 0.5 * cauchySum::cauchyPDF(theta_sample, s_camera, m_camera, (cdf_theta_e - cdf_theta_o)) + 0.5 * cauchySum::cauchyPDF(theta_sample, s_light, m_light, (cdf_theta_e - cdf_theta_o));
	pdf = pdf * (samplingData.dist_interestPoint_closestPoint / (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
}