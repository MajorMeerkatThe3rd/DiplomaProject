#include "../../sampling.h"
#include "../../helperFunctions.h"

namespace gaussianMM {
	double s_camera = 0.0, m_camera = 0.0;

	double erfinv(long double x) {

		if (x < -1 || x > 1) {
			return std::numeric_limits<long double>::quiet_NaN();
		}
		else if (x == 1.0) {
			return std::numeric_limits<long double>::infinity();
		}
		else if (x == -1.0) {
			return -std::numeric_limits<long double>::infinity();
		}

		const long double LN2 = 6.931471805599453094172321214581e-1L;

		const long double A0 = 1.1975323115670912564578e0L;
		const long double A1 = 4.7072688112383978012285e1L;
		const long double A2 = 6.9706266534389598238465e2L;
		const long double A3 = 4.8548868893843886794648e3L;
		const long double A4 = 1.6235862515167575384252e4L;
		const long double A5 = 2.3782041382114385731252e4L;
		const long double A6 = 1.1819493347062294404278e4L;
		const long double A7 = 8.8709406962545514830200e2L;

		const long double B0 = 1.0000000000000000000e0L;
		const long double B1 = 4.2313330701600911252e1L;
		const long double B2 = 6.8718700749205790830e2L;
		const long double B3 = 5.3941960214247511077e3L;
		const long double B4 = 2.1213794301586595867e4L;
		const long double B5 = 3.9307895800092710610e4L;
		const long double B6 = 2.8729085735721942674e4L;
		const long double B7 = 5.2264952788528545610e3L;

		const long double C0 = 1.42343711074968357734e0L;
		const long double C1 = 4.63033784615654529590e0L;
		const long double C2 = 5.76949722146069140550e0L;
		const long double C3 = 3.64784832476320460504e0L;
		const long double C4 = 1.27045825245236838258e0L;
		const long double C5 = 2.41780725177450611770e-1L;
		const long double C6 = 2.27238449892691845833e-2L;
		const long double C7 = 7.74545014278341407640e-4L;

		const long double D0 = 1.4142135623730950488016887e0L;
		const long double D1 = 2.9036514445419946173133295e0L;
		const long double D2 = 2.3707661626024532365971225e0L;
		const long double D3 = 9.7547832001787427186894837e-1L;
		const long double D4 = 2.0945065210512749128288442e-1L;
		const long double D5 = 2.1494160384252876777097297e-2L;
		const long double D6 = 7.7441459065157709165577218e-4L;
		const long double D7 = 1.4859850019840355905497876e-9L;

		const long double E0 = 6.65790464350110377720e0L;
		const long double E1 = 5.46378491116411436990e0L;
		const long double E2 = 1.78482653991729133580e0L;
		const long double E3 = 2.96560571828504891230e-1L;
		const long double E4 = 2.65321895265761230930e-2L;
		const long double E5 = 1.24266094738807843860e-3L;
		const long double E6 = 2.71155556874348757815e-5L;
		const long double E7 = 2.01033439929228813265e-7L;

		const long double F0 = 1.414213562373095048801689e0L;
		const long double F1 = 8.482908416595164588112026e-1L;
		const long double F2 = 1.936480946950659106176712e-1L;
		const long double F3 = 2.103693768272068968719679e-2L;
		const long double F4 = 1.112800997078859844711555e-3L;
		const long double F5 = 2.611088405080593625138020e-5L;
		const long double F6 = 2.010321207683943062279931e-7L;
		const long double F7 = 2.891024605872965461538222e-15L;

		long double abs_x = abs(x);

		if (abs_x <= 0.85L) {
			long double r = 0.180625L - 0.25L * x * x;
			long double num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0);
			long double den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + B0);
			return x * num / den;
		}

		long double r = sqrt(LN2 - log(1.0L - abs_x));

		long double num, den;
		if (r <= 5.0L) {
			r = r - 1.6L;
			num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0);
			den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + D0);
		}
		else {
			r = r - 5.0L;
			num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0);
			den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + F0);
		}

		if (x < 0L) {
			return -num / den;
		}
		else {
			return num / den;
		}
	}

	// pdf of Gaussian
	double gaussianPDF(double x, double s, double m, double multiplier = 1) {
		return multiplier == 0.0 ? 0.0 : (1.0 / (s * SQRT_2PI)) * exp(-0.5 * ((x - m) / s) * ((x - m) / s)) / multiplier;
	}

	// cdf of Gaussian
	double gaussianCDF(double x, double s, double m, double multiplier = 1, double start = 0) {
		return 0.5 * (1.0 + erf((x - m) / (s * SQRT_2))) / multiplier - start;
	}

	// inverse cdf of Gaussian
	double gaussianCDFinv(double u, double s, double m) {
		return s * SQRT_2* erfinv(2 * u - 1) + m;
	}

	double computeSWithMax(double max) {
		return 1 / (SQRT_2PI * max);
	}
}

void initGaussianMM(Scene scene) {
	double cameraMax = scene.medium->g > 0 ? pfHenyeyGreenstein(1, scene.medium->g) : pfHenyeyGreenstein(-1, scene.medium->g);
	gaussianMM::s_camera = gaussianMM::computeSWithMax(cameraMax);
	gaussianMM::m_camera = scene.medium->g > 0 ? PI / 2 : -PI / 2;  // the two rays are always perpendicular, dot is always 0, acos always PI/2
}

void sampleGaussianMM(SamplingData samplingData, double random, Vec& sample, double& pdf) {
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
	double s_light = gaussianMM::computeSWithMax(lightMax);
	double s_camera = gaussianMM::s_camera;

	// create the location for both, the light and the camera functions
	double lightDirectionDot = e.dot(samplingData.direction_interestPoint_closestPoint);
	double m_light = samplingData.medium->g > 0 ? -arccos(lightDirectionDot) : PI - arccos(lightDirectionDot);
	m_light = e.dot(cameraRay.direction) > 0.0 ? -m_light : m_light;
	double m_camera = gaussianMM::m_camera;

	double cdf_light_theta_o = gaussianMM::gaussianCDF(theta_o, s_light, m_light);
	double cdf_light_theta_e = gaussianMM::gaussianCDF(theta_e, s_light, m_light);
	double cdf_camera_theta_o = gaussianMM::gaussianCDF(theta_o, s_camera, m_camera);
	double cdf_camera_theta_e = gaussianMM::gaussianCDF(theta_e, s_camera, m_camera);

	double theta_sample = 0.0, u = 0.0;
	if (random < 0.5) {
		u = (random / 0.5) * (cdf_light_theta_e - cdf_light_theta_o) + cdf_light_theta_o;
		theta_sample = u > 1.0 - eps ? theta_e : u < eps ? theta_o : gaussianMM::gaussianCDFinv(u, s_light, m_light);
	}
	else {
		u = ((random - 0.5) / 0.5) * (cdf_camera_theta_e - cdf_camera_theta_o) + cdf_camera_theta_o;
		theta_sample = u > 1.0-eps ? theta_e : u < eps ? theta_o : gaussianMM::gaussianCDFinv(u, s_camera, m_camera);
	}

	double dist_sample_closestPoint = samplingData.dist_interestPoint_closestPoint * tan(theta_sample);
	sample = samplingData.closestPoint + cameraRay.direction * dist_sample_closestPoint;

	pdf = 0.5 * gaussianMM::gaussianPDF(theta_sample, s_camera, m_camera, (cdf_camera_theta_e - cdf_camera_theta_o)) + 0.5 * gaussianMM::gaussianPDF(theta_sample, s_light, m_light, (cdf_light_theta_e - cdf_light_theta_o));
	pdf = pdf * (samplingData.dist_interestPoint_closestPoint / (samplingData.dist_interestPoint_closestPoint * samplingData.dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
}