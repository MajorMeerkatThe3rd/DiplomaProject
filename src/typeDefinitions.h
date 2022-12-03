#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>

#define PI 3.1415926536
#define INV_4PI 1.0/(4.0*PI)
#define SQRT_2 sqrt(2.0)
#define SQRT_2PI sqrt(2.0*PI)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define DIFFUSE 1
#define SPECULAR 2
#define REFRACTIVE 3

#define LIGHTPASS 0
#define RENDERPASS 1

#define LIGHT_BIDIR 0
#define LIGHT_DIR 1

#define SAMPLE_UNI 0
#define SAMPLE_EXP 1
#define SAMPLE_ORJ 2
#define SAMPLE_SIJ 3
#define SAMPLE_ADJ 4
#define SAMPLE_PRO_GAU 5
#define SAMPLE_SUM_CAU 6
#define SAMPLE_MM_CAU 7
#define SAMPLE_MM_GAU 8
#define SAMPLE_MM_HYS 9
#define SAMPLE_MM_LOG 10

extern int width, height;
const double inf = 1e9;
const double eps = 1e-6;
using namespace std;

struct RenderSettings {
	int scene = -1;
	int spp = 100, vrlps = 1;
	double timeComparison = INT_MAX;
	double sigma_a = -1, sigma_s = -1, g = -1, intensity = -1;
	int bounces = -1;
	int sampling = -1;
	int lightType = -1;

	double fov = PI / 4;

	bool surfaceRadiance = true;
	bool mediumRadiance = true;
	bool lightInMedium = true;
	bool cameraInMedium = true;
};

struct Vec {
	double x, y, z;
	Vec(double x0, double y0, double z0) { x = x0; y = y0; z = z0; }
	Vec(double x0) { x = x0; y = x0; z = x0; }
	Vec() { x = 0; y = 0; z = 0; }
	Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	Vec operator*(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
	Vec operator/(const Vec& b) const { return Vec(x / b.x, y / b.y, z / b.z); }
	Vec& norm() { return *this = *this * (1 / std::sqrt(x * x + y * y + z * z)); }
	double length() { return std::sqrt(x * x + y * y + z * z); }
	double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(const Vec& b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
	double dist(const Vec& b) const { return std::sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y) + (z - b.z) * (z - b.z)); }
	double average() const { return (x + y + z) / 3.0; }
};

// Rays have origin and direction.
// The direction vector should always be normalized.
struct Ray {
	Vec o, d;
	Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

// struct to store segments, photon power at the origin of the ray and the distances
struct Line {
public:
	Vec origin, direction, end, photonPower;
	double distance;

	Line(Vec o0, Vec d0, Vec e0, Vec power0, double dist0) :
		origin(o0), direction(d0), end(e0), photonPower(power0), distance(dist0) {}
	Line() {}
	void changeEnd(Vec e0, double dist0) {
		end = e0; distance = dist0;
	}
	Ray getRay() { return Ray(origin, direction); }
	void setPhotonPower(Vec p) { photonPower = p; }
};

// Objects have color, emission, type (diffuse, specular, refractive)
// All object should be intersectable and should be able to compute their surface normals.
class Obj {
public:
	Vec cl;
	double emission;
	int type;
	void setMat(Vec cl_ = 0, double emission_ = 0, int type_ = 0) { cl = cl_; emission = emission_; type = type_; }
	virtual double intersect(const Ray&) const = 0;
	virtual Vec normal(const Vec&) const = 0;
};

class Plane : public Obj {
public:
	Vec n;
	double d;
	Plane(double d_ = 0, Vec n_ = 0) {
		d = d_;
		n = n_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if (d0 != 0) {
			double t = -1 * (((n.dot(ray.o)) + d) / d0);
			return (t > eps) ? t : 0;
		}
		else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }
};

class RestrictedPlane : public Obj {
public:
	Vec n, min, max;
	double d;
	RestrictedPlane(double d_ = 0, Vec n_ = 0, Vec min_ = 0, Vec max_ = 0) {
		d = d_;
		n = n_;
		min = min_;
		max = max_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if (d0 != 0) {
			double t = -1 * (((n.dot(ray.o)) + d) / d0);

			Vec p;
			if (t > eps) {
				p = ray.o + ray.d * t;
			}
			else {
				return 0;
			}

			if (n.z != 0) {
				if (min.x < p.x && p.x < max.x &&
					min.y < p.y && p.y < max.y) {
					return t;
				}
				else {
					return 0;
				}
			}
			else if (n.y != 0) {
				if (min.x < p.x && p.x < max.x &&
					min.z < p.z && p.z < max.z) {
					return t;
				}
				else {
					return 0;
				}
			}
			else if (n.x != 0) {
				if (min.y < p.y && p.y < max.y &&
					min.z < p.z && p.z < max.z) {
					return t;
				}
				else {
					return 0;
				}
			}
			else {
				return 0;
			}

			return 0;
		}
		else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }

	Vec getRandomPoint(double u1, double u2) {
		if (n.z != 0) {
			double distanceX = max.x - min.x;
			double distanceY = max.y - min.y;
			return Vec(min.x + distanceX * u1, min.y + distanceY * u2, -n.z * d);
		}
		else if (n.y != 0) {
			double distanceX = max.x - min.x;
			double distanceZ = max.z - min.z;
			return Vec(min.x + distanceX * u1, -n.y * d, min.z + distanceZ * u2);
		}
		else if (n.x != 0) {
			double distanceY = max.y - min.y;
			double distanceZ = max.z - min.z;
			return Vec(-n.x * d, min.y + distanceY * u1, min.z + distanceZ * u2);
		}
		return Vec();
	}
};

class AABox : public Obj {
public:
	Vec min, max;
	AABox(Vec min_ = 0, Vec max_ = 0) {
		min = min_;
		max = max_;
	}

	//intersection routine from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
	double intersect(const Ray& ray) const {
		double tmin = (min.x - ray.o.x) / ray.d.x;
		double tmax = (max.x - ray.o.x) / ray.d.x;

		if (tmin > tmax) std::swap(tmin, tmax);

		double tymin = (min.y - ray.o.y) / ray.d.y;
		double tymax = (max.y - ray.o.y) / ray.d.y;

		if (tymin > tymax) std::swap(tymin, tymax);

		if ((tmin > tymax) || (tymin > tmax))
			return 0;

		if (tymin > tmin)
			tmin = tymin;

		if (tymax < tmax)
			tmax = tymax;

		double tzmin = (min.z - ray.o.z) / ray.d.z;
		double tzmax = (max.z - ray.o.z) / ray.d.z;

		if (tzmin > tzmax) std::swap(tzmin, tzmax);

		if ((tmin > tzmax) || (tzmin > tmax))
			return 0;

		if (tzmin > tmin)
			tmin = tzmin;

		if (tzmax < tmax)
			tmax = tzmax;

		if (tmin > tmax)
			std::swap(tmin, tmax);

		return tmin > eps ? tmin : tmax > eps ? tmax : 0;
	}

	Vec normal(const Vec& p0) const {
		if (p0.x >= min.x - eps && p0.x <= min.x + eps) return Vec(-1, 0, 0);
		if (p0.x >= max.x - eps && p0.x <= max.x + eps) return Vec(1, 0, 0);
		if (p0.y >= min.y - eps && p0.y <= min.y + eps) return Vec(0, -1, 0);
		if (p0.y >= max.y - eps && p0.y <= max.y + eps) return Vec(0, 1, 0);
		if (p0.z >= min.z - eps && p0.z <= min.z + eps) return Vec(0, 0, -1);
		if (p0.z >= max.z - eps && p0.z <= max.z + eps) return Vec(0, 0, 1);
	}
};

class Sphere : public Obj {
public:
	Vec c;
	double r;

	Sphere(double r_ = 0, Vec c_ = 0) { c = c_; r = r_; }
	double intersect(const Ray& ray) const {
		double b = ((ray.o - c) * 2).dot(ray.d);
		double c_ = (ray.o - c).dot((ray.o - c)) - (r * r);
		double disc = b * b - 4 * c_;
		if (disc < 0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b + disc;
		double sol2 = -b - disc;
		return (sol2 > eps) ? sol2 / 2 : ((sol1 > eps) ? sol1 / 2 : 0);
	}

	Vec normal(const Vec& p0) const {
		return (p0 - c).norm();
	}
};

class Intersection {
public:
	Intersection() { t = inf; object = nullptr; }
	Intersection(double t_, Obj* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Obj* object;
};

class Medium {
public:
	double g, k;
	Vec sigma_a, sigma_s, sigma_t;
	Obj* obj;

	Medium(Vec s_a, Vec s_s, double g0, Obj* o0)
		: sigma_a(s_a), sigma_s(s_s), g(g0), obj(o0) {
		sigma_t = sigma_a + sigma_s;
		k = 1.55 * g - 0.55 * g * g * g;
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		if (!obj) return closestIntersection;
		double t = obj->intersect(ray);
		if (t > eps&& t < closestIntersection.t) {
			closestIntersection.t = t;
			closestIntersection.object = obj;
		}
		return closestIntersection;
	}
};

class Pointlight {
public:
	Vec pos, intensity;
	Pointlight(Vec pos0, Vec intensity0)
		: pos(pos0), intensity(intensity0) {}
};

//multiple restricted plane objects, all have to have the same surface area 
class PlaneAreaLights {
	vector<RestrictedPlane*> planes;
public:
	void add(RestrictedPlane* plane) {
		planes.push_back(plane);
	}

	Vec getRandomPoint(double u1, double u2, double u3) {
		int size = planes.size();
		double invSize = 1.0 / (double)size;
		double prevSection = 0.0;

		for (int i = 0; i < size; i++) {
			double nextSection = (i + 1) * invSize;
			//if in it
			if (prevSection <= u1 && u1 <= nextSection) {
				return planes.at(i)->getRandomPoint(u2, u3);
			}
			prevSection = nextSection;
		}
		return Vec(1);
	}

	Vec getDirection() {
		return Vec(-.1, -1, .075).norm();
	}
};

struct SamplingData {
	SamplingData();
	SamplingData(Line l, Medium* m) {
		lineStart = l.origin;
		lineEnd = l.end;
		direction_line = l.direction;
		distance_line = l.distance;
		medium = m;
	}

	Vec lineStart, lineEnd, direction_line;
	Vec interestPoint, direction_interestPoint;
	Vec closestPoint, direction_interestPoint_closestPoint;
	double distance_line, dist_lineOrigin_closestPoint, dist_interestPoint_closestPoint, dist_both_lines, theta_directions;
	Medium* medium;
};

class Scene {
	vector<Obj*> objects;

public:
	Medium* medium;
	Pointlight* pointlight;
	PlaneAreaLights* areaLights;
	std::vector<Vec>* lightDirections;

	double refr_index;

	void add(Obj* object) {
		objects.push_back(object);
	}

	void add(Pointlight* p) {
		pointlight = p;
	}

	void add(Medium* m) {
		medium = m;
	}

	void add(PlaneAreaLights* a) {
		areaLights = a;
	}

	void add(std::vector<Vec>* d) {
		lightDirections = d;
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		// intersect all objects, one after the other
		for (auto iter = objects.begin(); iter != objects.end(); ++iter) {
			double t = (*iter)->intersect(ray);
			if (t > eps&& t < closestIntersection.t) {
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}
		return closestIntersection;
	}

	void clear() {
		for (int i = 0; i < objects.size(); i++) {
			delete objects.at(i);
		}
	}
};