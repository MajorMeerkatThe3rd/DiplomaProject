#include "../../sampling.h"
#include "../../helperFunctions.h"

namespace hyperbolicSecantMM {
	double s_camera = 0.0, m_camera = 0.0;
	// pdf of hyperbolic secant
	double hyperbolicSecantPDF(double x, double s, double m, double multiplier = 1) {
		return (1.0 / s)* (1.0 / cosh((x-m)/s)) / multiplier;
	}

	// cdf of hyperbolic secant
	double hyperbolicSecantCDF(double x, double s, double m, double multiplier = 1, double start = 0) {
		return atan(sinh((x-m)/s)) / multiplier - start;
	}

	// inverse cdf of hyperbolic secant
	double hyperbolicSecantCDFinv(double u, double s, double m) {
		return asinh(tan(u))*s+m;
	}

	double computeSWithMax(double max) {
		return 1 / max;
	}
}

void initHyperbolicSecantMM(Scene scene) {
	double cameraMax = scene.medium->g > 0 ? pfHenyeyGreenstein(1, scene.medium->g) : pfHenyeyGreenstein(-1, scene.medium->g);
	hyperbolicSecantMM::s_camera = hyperbolicSecantMM::computeSWithMax(cameraMax);
	hyperbolicSecantMM::m_camera = scene.medium->g > 0 ? PI / 2 : -PI / 2;  // the two rays are always perpendicular, dot is always 0, acos always PI/2
}

void sampleHyperbolicSecantMM(SamplingData samplingData, double random, Vec& sample, double& pdf) {
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
	double s_light = hyperbolicSecantMM::computeSWithMax(lightMax);
	double s_camera = hyperbolicSecantMM::s_camera;

	// create the location for both, the light and the camera functions
	double lightDirectionDot = e.dot(samplingData.direction_interestPoint_closestPoint);
	double m_light = samplingData.medium->g > 0 ? -arccos(lightDirectionDot) : PI - arccos(lightDirectionDot);
	m_light = e.dot(cameraRay.direction) > 0.0 ? -m_light : m_light;
	double m_camera = hyperbolicSecantMM::m_camera;

	double cdf_light_theta_o = hyperbolicSecantMM::hyperbolicSecantCDF(theta_o, s_light, m_light);
	double cdf_light_theta_e = hyperbolicSecantMM::hyperbolicSecantCDF(theta_e, s_light, m_light);
	double cdf_camera_theta_o = hyperbolicSecantMM::hyperbolicSecantCDF(theta_o, s_camera, m_camera);
	double cdf_camera_theta_e = hyperbolicSecantMM::hyperbolicSecantCDF(theta_e, s_camera, m_camera);

	double theta_sample = 0.0;
	if (random < 0.5) {
		double u = (random / 0.5) * (cdf_light_theta_e - cdf_light_theta_o) + cdf_light_theta_o;
		theta_sample = hyperbolicSecantMM::hyperbolicSecantCDFinv(u, s_light, m_light);
	}
	else {
		double u = ((random - 0.5) / 0.5) * (cdf_camera_theta_e - cdf_camera_theta_o) + cdf_camera_theta_o;
		theta_sample = hyperbolicSecantMM::hyperbolicSecantCDFinv(u, s_camera, m_camera);
	}

	double dist_sample_closestPoint = samplingData.dist_interestPoint_closestPoint * tan(theta_sample);
	sample = samplingData.closestPoint + cameraRay.direction * dist_sample_closestPoint;

	pdf = 0.5 * hyperbolicSecantMM::hyperbolicSecantPDF(theta_sample, s_camera, m_camera, (cdf_camera_theta_e - cdf_camera_theta_o)) + 0.5 * hyperbolicSecantMM::hyperbolicSecantPDF(theta_sample, s_light, m_light, (cdf_light_theta_e - cdf_light_theta_o));
	pdf = pdf * (samplingData.dist_interestPoint_closestPoint / (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
}