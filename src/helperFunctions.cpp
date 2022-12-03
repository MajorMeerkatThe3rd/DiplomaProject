#include "helperFunctions.h"

// given v1, set v2 and v3 so they form an orthonormal system
// (we assume v1 is already normalized)
void ons(const Vec& v1, Vec& v2, Vec& v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		double invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * invLen, 0.0f, v1.x * invLen);
	}
	else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		double invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0f, v1.z * invLen, -v1.y * invLen);
	}
	v3 = v1 % v2;
}


// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y, double fov) {
	double w = width;
	double h = height;
	double fovx = fov;
	double fovy = (h / w) * fovx;
	return Vec(((2 * x - w) / w) * tan(fovx),
		-((2 * y - h) / h) * tan(fovy),
		-1.0);
}

double arccos(double x) {
	return acos(x < -1.0 ? -1.0 : x > 1.0 ? 1.0 : x);
}

// rotate a Vector v around another Vector n
Vec rotateVec(Vec n, Vec v, double theta) {
	return v * cos(theta) + n % v * sin(theta);
}

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec sampleHemisphere(double u1, double u2) {
	double r = sqrt(1.0 - u1 * u1);
	double phi = 2 * PI * u2;
	return Vec(cos(phi) * r, sin(phi) * r, u1);
}

Vec sampleSphere(double u1, double u2) {
	const double cos_theta = u1 * 2.0 - 1.0; // remap to -1 .. 1
	const double r = sqrt(1.0 - cos_theta * cos_theta);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi) * r, sin(phi) * r, cos_theta);
}

Vec transmittance(Vec sigma_t, double dist) {
	Vec exp_term = sigma_t * dist * -1;
	return Vec(exp(exp_term.x), exp(exp_term.y), exp(exp_term.z));
}

Vec sampleHenyeyGreenstein(Vec incomingDir, double g, double u1, double u2) {
	double cos_theta;
	if (g == 0) {
		cos_theta = u1 * 2.0 - 1.0;
	}
	else {
		double g_spu = g * g;
		double squ_term = (1.0 - g_spu) / (1.0 - g + 2 * g * u1);
		cos_theta = (1.0 + g_spu - squ_term * squ_term) / (2.0 * g);
	}
	double sin_theta = sqrt(max(0.0, 1.0 - cos_theta * cos_theta));
	double phi = 2 * PI * u2;
	Vec v1, v2;
	ons(incomingDir, v1, v2);
	return v1 * sin_theta * cos(phi) + v2 * sin_theta * sin(phi) + incomingDir * cos_theta;
}

double pfHenyeyGreenstein(Vec v1, Vec v2, double g) {
	if (g == 0.0)
		return INV_4PI;
	double cos_theta = v1.dot(v2);
	double g_squ = g * g;
	double temp = 1.0 + g_squ - 2.0 * g * cos_theta;
	return INV_4PI * ((1.0 - g_squ) / (temp * sqrt(temp)));
}

double pfHenyeyGreenstein(double dot, double g) {
	if (g == 0.0)
		return INV_4PI;
	double g_squ = g * g;
	double temp = 1.0 + g_squ - 2.0 * g * dot;
	return INV_4PI * ((1.0 - g_squ) / (temp * sqrt(temp)));
}

void clear(Scene& scene, Vec** pix) {
	scene.clear();
	for (int i = 0; i < width; i++) {
		delete[] pix[i];
	}
	delete[] pix;
}

int color(double clr) {
	return max(min((int)(max(min(clr, 255.0), 0.0) * 255.0), 255), 0);
}

void saveImage(Vec** pix, std::string name, int spp) {
	FILE* f = fopen(std::string(name + ".ppm").c_str(), "w");
	fprintf(f, "P3\n%d %d\n%d\n ", width, height, 255);
	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {
			fprintf(f, "%d %d %d ", color(pix[col][row].x / (double)spp), color(pix[col][row].y / (double)spp), color(pix[col][row].z / (double)spp));
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
