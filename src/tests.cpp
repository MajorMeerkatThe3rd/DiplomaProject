#pragma once

#include "tests.h"

// adjusted version of: https://people.sc.fsu.edu/~jburkardt/cpp_src/sphere_fibonacci_grid/sphere_fibonacci_grid.html
// to generate equally spaces samples on a sphere
void fibonacciSamples(int amount, std::vector<Vec>& directions) {
	int i, j;
	double cphi, i_r8, ng_r8, r8_phi, sphi, theta;
	r8_phi = (1.0 + sqrt(5.0)) / 2.0;
	ng_r8 = (double)(amount);

	for (j = 0; j < amount; j++) {
		i_r8 = (double)(-amount + 1 + 2 * j);
		theta = 2.0 * PI * i_r8 / r8_phi;
		sphi = i_r8 / ng_r8;
		cphi = sqrt((ng_r8 + i_r8) * (ng_r8 - i_r8)) / ng_r8;
		double x = cphi * sin(theta);
		double y = cphi * cos(theta);
		double z = sphi;
		directions.push_back(Vec(x, y, z));
	}
}

void loadScene(Scene& scene, RenderSettings& s) {
	auto add = [&scene](Obj* s, Vec cl, double emission, int type) {
		s->setMat(cl, emission, type);
		scene.add(s);
	};

	scene.refr_index = 1.5;

	bool loadCornellBox = false;
	s.g = s.g == -1 ? 0.0 : s.g;
	s.bounces = s.bounces == -1 ? 1 : s.bounces;
	s.sampling = s.sampling == -1 ? SAMPLE_ADJ : s.sampling;

	if (s.scene == 0 || s.scene == -1) {
		/** vrl default scene settings */
		s.scene = 0;
		s.sigma_a = s.sigma_a == -1 ? .1 : s.sigma_a;
		s.sigma_s = s.sigma_s == -1 ? .01 : s.sigma_s;
		s.lightType = s.lightType == -1 ? LIGHT_DIR : s.lightType;
		s.intensity = s.intensity == -1 ? 50000 : s.intensity;

		Pointlight* pointlight = new Pointlight(Vec(0, 1.9, -3), Vec(s.intensity));
		scene.add(pointlight);
		Medium* medium = new Medium(s.sigma_a, s.sigma_s, s.g, NULL);
		scene.add(medium);
		// Position, normal, color, emission, type for planes*
		add(new Plane(2.55, Vec(1, 0, 0)), Vec(1, .5, 0), 0, 1); // Left plane
		add(new Plane(2.55, Vec(-1, 0, 0)), Vec(0, 1, 1), 0, 1); // Right plane
		add(new Plane(2.5, Vec(0, 1, 0)), Vec(.5, .5, .5), 0, 1); // Bottom plane
		add(new Plane(2.5, Vec(0, -1, 0)), Vec(.5, .5, .5), 0, 1); // Ceiling plane
		add(new Plane(4.5, Vec(0, 0, 1)), Vec(.5, .5, .5), 0, 1); // Back plane
		add(new Plane(0.5, Vec(0, 0, -1)), Vec(.5, .5, .5), 0, 1); // Front plane
		add(new AABox(Vec(0.25, -2.5, -3.75), Vec(1.75, -1, -2.75)), Vec(.5, .5, .5), 0, 1); // Front box
		add(new AABox(Vec(-1.75, -2.5, -4.5), Vec(-0.25, .5, -3.75)), Vec(.5, .5, .5), 0, 1); // Back box

	}
	else if (s.scene == 1) {
		/** infinite medium point light settings - bidirectional */
		s.sigma_a = s.sigma_a == -1 ? .1 : s.sigma_a;
		s.sigma_s = s.sigma_s == -1 ? .25 : s.sigma_s;
		s.lightType = s.lightType == -1 ? LIGHT_BIDIR : s.lightType;
		s.intensity = s.intensity == -1 ? 250 : s.intensity;

		Pointlight* pointlight = new Pointlight(Vec(0, 0, -2), Vec(s.intensity));
		scene.add(pointlight);
		Medium* medium = new Medium(s.sigma_a, s.sigma_s, s.g, NULL);
		scene.add(medium);
	}
	else if (s.scene == 2) {
		/** infinite medium point light settings - unidirectional */
		s.sigma_a = s.sigma_a == -1 ? .1 : s.sigma_a;
		s.sigma_s = s.sigma_s == -1 ? .25 : s.sigma_s;
		s.intensity = s.intensity == -1 ? 100 : s.intensity;

		Pointlight* pointlight = new Pointlight(Vec(0, 0, -2), Vec(s.intensity));
		scene.add(pointlight);
		Medium* medium = new Medium(s.sigma_a, s.sigma_s, s.g, NULL);
		scene.add(medium);
	}
	else if (s.scene == 3 || s.scene == 4) {
		/** sphere corner settings */
		s.sigma_a = s.sigma_a == -1 ? .65 : s.sigma_a;
		s.sigma_s = s.sigma_s == -1 ? .5 : s.sigma_s;
		s.intensity = s.intensity == -1 ? s.scene == 3 ? 250000 : 3500000 : s.intensity;

		add(new Plane(5.2, Vec(0.973562, -0.226953, -0.025893)), Vec(1, 129.0 / 255.0, 0), 0, 1); // Left plane
		add(new Plane(0.3, Vec(-0.973562, 0.226953, 0.025893)), Vec(0, 242.0 / 255.0, 1), 0, 1); // Right plane
		add(new Plane(0.4, Vec(-0.105381, 0.974352, 0.198829)), Vec(.5, .5, .5), 0, 1); // Bottom plane
		add(new Plane(5.1, Vec(0.105381, -0.974352, -0.198829)), Vec(.5, .5, .5), 0, 1); // Ceiling plane
		add(new Plane(1.5, Vec(-0.122788, 0.007596, 0.992404)), Vec(.5, .5, .5), 0, 1); // Back plane
		add(new Plane(4.0, Vec(0.122788, -0.007596, -0.992404)), Vec(.5, .5, .5), 0, 1); // Front plane

		add(new Sphere(0.25, Vec(0, 0, -.8)), Vec(1, 1, 1), 0, 1); // Left sphere
		Pointlight* pointlight = new Pointlight(Vec(-2.45, 2.35, 1.25), Vec(s.intensity));
		scene.add(pointlight);
		Medium* medium = new Medium(s.sigma_a, s.sigma_s, s.g, NULL);
		scene.add(medium);
		s.fov = PI / 6;
		if (s.scene == 3) {
			s.surfaceRadiance = false;
			s.mediumRadiance = true;
		}
		else {
			s.surfaceRadiance = true;
			s.mediumRadiance = false;
		}
	}
	else if (s.scene == 5) {
		/** small medium cube settings */
		s.sigma_a = s.sigma_a == -1 ? .5 : s.sigma_a;
		s.sigma_s = s.sigma_s == -1 ? .1 : s.sigma_s;
		s.intensity = s.intensity == -1 ? 500000 : s.intensity;

		Pointlight* pointlight = new Pointlight(Vec(0, 1.9, -3), Vec(s.intensity));
		scene.add(pointlight);
		double size = .75;
		Vec mid = Vec(0, -1, -3);
		Vec min = Vec(mid.x - size, mid.y - size, mid.z - size);
		Vec max = Vec(mid.x + size, mid.y + size, mid.z + size);
		AABox* mediumBox = new AABox(min, max);
		Medium* medium = new Medium(s.sigma_a, s.sigma_s, s.g, mediumBox);
		scene.add(medium);
		add(new Sphere(1.05, Vec(-0.75, -1.45, -4.4)), Vec(.5, .5, .5), 0, 2); // Middle sphere
		add(new Sphere(0.5, Vec(2.0, -2.05, -3.7)), Vec(.5, .5, .5), 0, 3); // Right sphere
		add(new Sphere(0.6, Vec(-1.75, -1.95, -3.1)), Vec(0, 0, 1), 0, 1); // Left sphere
		loadCornellBox = true;
		s.lightInMedium = false;
		s.cameraInMedium = false;
	}

	if (loadCornellBox) {
		// Position, normal, color, emission, type for planes*
		add(new Plane(2.75, Vec(1, 0, 0)), Vec(1, 0, 0), 0, 1); // Left plane
		add(new Plane(2.75, Vec(-1, 0, 0)), Vec(0, 1, 0), 0, 1); // Right plane
		add(new Plane(2.5, Vec(0, 1, 0)), Vec(.5, .5, .5), 0, 1); // Bottom plane
		add(new Plane(3.0, Vec(0, -1, 0)), Vec(.5, .5, .5), 0, 1); // Ceiling plane
		add(new Plane(5.5, Vec(0, 0, 1)), Vec(.5, .5, .5), 0, 1); // Back plane
		add(new Plane(0.5, Vec(0, 0, -1)), Vec(.5, .5, .5), 0, 1); // Front plane
	}

	// create an area light with four planes if needed
	if (s.lightType == LIGHT_DIR) {
		PlaneAreaLights *areaLights = new PlaneAreaLights();
		areaLights->add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.65, 0, -3.5), Vec(-.495, 0, -3.18)));
		areaLights->add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.425, 0, -3.5), Vec(-.27, 0, -3.18)));
		areaLights->add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.2, 0, -3.5), Vec(-.045, 0, -3.18)));
		areaLights->add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(.025, 0, -3.5), Vec(.18, 0, -3.18)));
		scene.add(areaLights);
	}

	// generate equally spaced direction
	std::vector<Vec>* directions = new std::vector<Vec>();
	fibonacciSamples(s.vrlps, *directions);
	scene.add(directions);
}
