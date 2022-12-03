// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
// virtual ray lights implementation by michael oppitz - e1227129@student.tuwien.ac.at
//
// g++ vrl.cpp -O3 -std=gnu++0x -fopenmp -static-libgcc -static-libstdc++
//
// This is an implementation of the "Virtual Ray Lights for Rendering Scenes with Participating Media"
// paper by Nov�k et al.
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <random>
#include <ctime>

#include "typeDefinitions.h"
#include "helperFunctions.h"
#include "sampling.h"
#include "tests.h"

using namespace std;


struct VecInt {
	int x, y, z;
	VecInt(int x0, int y0, int z0) { x = x0; y = y0; z = z0; }
	VecInt(int x0) { x = x0; y = x0; z = x0; }
	VecInt() { x = 0; y = 0; z = 0; }
};

// Helpers for random number generation
std::mt19937_64 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

int width, height;

void getInitialLightSample(SamplingData samplingData, int method, double random, Vec& sample, double& pdf) {
	if (method == SAMPLE_UNI) {
		sampleUniform(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_EXP) {
		sampleExponential(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_ORJ) {
		sampleUniform(samplingData, random, sample, pdf);
	}
	else {
		sampleVRLLight(samplingData, random, sample, pdf);
	}
}

void getSample(SamplingData samplingData, int method, double random, Vec& sample, double& pdf) {
	if (method == SAMPLE_UNI) {
		sampleUniform(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_EXP) {
		sampleExponential(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_ORJ) {
		sampleEquiAngular(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_SIJ) {
		sampleEquiAngular(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_ADJ) {
		sampleVRLAnisotropic(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_PRO_GAU) {
		sampleGaussianProduct(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_SUM_CAU) {
		sampleCauchySum(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_MM_CAU) {
		sampleCauchyMM(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_MM_GAU) {
		sampleGaussianMM(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_MM_HYS) {
		sampleHyperbolicSecantMM(samplingData, random, sample, pdf);
	}
	else if (method == SAMPLE_MM_LOG) {
		sampleLogisticMM(samplingData, random, sample, pdf);
	}
}

// calclualte the radiance from the medium
Vec mediumRadiance(const Scene& scene, std::vector<Line> lightRays, Line cameraRay, RenderSettings settings) {
	Vec contribution;
	SamplingData samplingData = SamplingData(cameraRay, scene.medium);

	for (auto iter = lightRays.begin(); iter != lightRays.end(); ++iter) {
		Line lightRay = *iter;
		SamplingData samplingDataVPL = SamplingData(lightRay, scene.medium);
		samplingData.direction_interestPoint = lightRay.direction;

		double pdf_light, pdf_camera;
		Vec lightSample, cameraSample;

		// find the points with the shortest distance between the two ray segments
		shortestPointsOnLineSegments(lightRay, cameraRay, samplingDataVPL, samplingData);
		samplingDataVPL.theta_directions = arccos((lightRay.direction).dot(cameraRay.direction * -1));

		getInitialLightSample(samplingDataVPL, settings.sampling, RND2, lightSample, pdf_light);
		samplingData.interestPoint = lightSample;

		// find the closest point on the camera ray to the light ray sample
		closestPointOnLine(samplingData);
		samplingData.direction_interestPoint_closestPoint = samplingData.closestPoint - samplingData.interestPoint;
		samplingData.dist_interestPoint_closestPoint = samplingData.direction_interestPoint_closestPoint.length();
		samplingData.direction_interestPoint_closestPoint = samplingData.direction_interestPoint_closestPoint / samplingData.dist_interestPoint_closestPoint;

		// get the sample for the camera ray
		getSample(samplingData, settings.sampling, RND2, cameraSample, pdf_camera);

		// calculate distance and direction between samples
		Vec dir_samples = (cameraSample - lightSample).norm();
		double dist_samples = cameraSample.dist(lightSample);

		// calculating visibility
		Intersection intersection = scene.intersect(Ray(lightSample, dir_samples));
		if (intersection) {
			if (intersection.t < dist_samples - eps) {
				continue;
			}
		}
		// phase functions with the angle as input (Henyey-Greenstein phase function)
		double fs_u = pfHenyeyGreenstein(cameraRay.direction * -1, dir_samples, scene.medium->g);
		double fs_v = pfHenyeyGreenstein(lightRay.direction, dir_samples, scene.medium->g);

		double dist_lightOrigin_lightSample = (lightRay.origin - lightSample).length();
		double dist_cameraOrigin_cameraSample = (cameraRay.origin - cameraSample).length();
		// transmittance functions
		Vec trans = transmittance(scene.medium->sigma_t, dist_cameraOrigin_cameraSample + dist_lightOrigin_lightSample + dist_samples);

		// inverse squared distance term
		double geom = 1.0 / (dist_samples * dist_samples);

		// putting the integrand of the formula together
		Vec integrand = lightRay.photonPower * scene.medium->sigma_s * scene.medium->sigma_s * trans * fs_u * fs_v * geom;
		contribution = (pdf_light * pdf_camera) == 0.0 ? 0.0 : contribution + (integrand / (pdf_light * pdf_camera)) / settings.vrlps;
	}

	return contribution;
}

// calculate the surface radiance
Vec surfaceRadiance(const Scene& scene, std::vector<Line> lightRays, Line cameraRay, Vec surfaceColor, Vec surfaceNormal, RenderSettings settings) {
	Vec contribution;
	for (auto iter = lightRays.begin(); iter != lightRays.end(); ++iter) {
		Line lightRay = *iter;
		SamplingData samplingData = SamplingData(lightRay, scene.medium);
		samplingData.interestPoint = cameraRay.end;
		samplingData.direction_interestPoint = cameraRay.direction;

		double pdf_light;
		Vec lightSample;

		// find the closest point on the light ray to the camera ray.e
		closestPointOnLine(samplingData);

		samplingData.direction_interestPoint_closestPoint = samplingData.closestPoint - cameraRay.end;
		samplingData.dist_interestPoint_closestPoint = samplingData.direction_interestPoint_closestPoint.length();
		samplingData.direction_interestPoint_closestPoint = samplingData.direction_interestPoint_closestPoint / samplingData.dist_interestPoint_closestPoint;

		getSample(samplingData, settings.sampling, RND2, lightSample, pdf_light);

		// calculate distance and direction between samples
		Vec dir_sample = (lightSample - cameraRay.end).norm();
		double dist_sample = lightSample.dist(cameraRay.end);

		// calculating visibility
		Intersection intersection = scene.intersect(Ray(cameraRay.end, dir_sample));
		if (intersection) {
			if (intersection.t < dist_sample - eps) {
				continue;
			}
		}

		// phase functions with the angle as input (Henyey-Greenstein phase function)
		double fs_v = pfHenyeyGreenstein(lightRay.direction * -1, dir_sample, scene.medium->g);

		// cosine-weight the BRDF
		Vec cosWeightedSurfaceColor = surfaceColor * abs(surfaceNormal.dot(dir_sample * -1));
		double dist_lightOrigin_lightSample = (lightRay.origin - lightSample).length();
		// transmittance functions
		Vec trans = transmittance(scene.medium->sigma_t, dist_sample + dist_lightOrigin_lightSample);

		// inverse squared distance term
		double geom = 1.0 / (dist_sample * dist_sample);

		// putting the integrand of the formula together
		Vec integrand = lightRay.photonPower * cosWeightedSurfaceColor * scene.medium->sigma_s * fs_v * trans * geom;

		contribution = contribution + (integrand / pdf_light) / settings.vrlps;
	}
	return contribution;
}

void trace(Line line, const Scene& scene, std::vector<Line>& lightRays, bool insideMedium, int depth, Vec& clr, RenderSettings settings, int pass) {

	double freeFlightDistance = 0.0;
	if (pass == LIGHTPASS) {
		if (depth >= settings.bounces)
			return;

		double freeFlightDistance = -log(RND2) / scene.medium->sigma_t.average();

		if (!insideMedium)
			freeFlightDistance = inf;

	}

	double rrFactor = 1.0;
	if (pass == RENDERPASS) {
		// Russian roulette: starting at depth 3, each recursive step will stop with a probability of 0.1
		if (depth >= 3) {
			const double rrStopProbability = 0.1;
			if (RND2 <= rrStopProbability) {
				return;
			}
			rrFactor = 1.0 / (1.0 - rrStopProbability);
		}
	}

	Intersection intersection = scene.intersect(line.getRay());
	Intersection intersectionMedium = scene.medium->intersect(line.getRay());

	bool mediumTransition = false;
	bool scatteringEvent = false;

	if (freeFlightDistance < intersection.t && freeFlightDistance < intersectionMedium.t)
		scatteringEvent = true;

	// a medium transition can only occur at an intersection with the medium where there is no intersection before this event
	if (intersectionMedium) {
		if (!intersection) {
			mediumTransition = true;
		}
		else if (intersectionMedium.t < intersection.t) {
			mediumTransition = true;
		}
	}

	Line original = line;

	// no intersection with geometry, no medium transition and outside the medium
	if (!intersection && !mediumTransition && !insideMedium) return;

	// no intersection with geometry, no medium transition and inside the medium
	if (!intersection && !mediumTransition && insideMedium) {
		line.changeEnd(line.origin + line.direction * inf, inf);
		if (pass == LIGHTPASS) {
			lightRays.push_back(line);
		}
		if (pass == RENDERPASS) {
			if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, line, settings) * rrFactor;
			return;
		}
	}

	// intersection, but no medium transition
	if (intersection && !mediumTransition) {
		// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
		Vec hp = line.origin + line.direction * intersection.t;
		Vec N = intersection.object->normal(hp);
		Vec surfaceClr = intersection.object->cl;
		int type = intersection.object->type;


		if (pass == LIGHTPASS) {
			line.changeEnd(hp, intersection.t);
			if (insideMedium) lightRays.push_back(line);
			if (depth + 1 >= settings.bounces) return;
		}

		if (pass == RENDERPASS) {
			// calculate VRL contribution
			if (insideMedium) {
				line.changeEnd(hp, intersection.t);
				if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, line, settings) * rrFactor;
			}
			else {
				line.changeEnd(hp, 0);
			}
		}

		// Diffuse BRDF
		if (type == DIFFUSE) {
			if (pass == LIGHTPASS) {
				// create new ray direction
				Vec rotX, rotY;
				ons(N, rotX, rotY);
				Vec sampledDir = sampleHemisphere(RND2, RND2);
				Vec rotatedDir;
				rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir);
				rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
				rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
				Vec direction = rotatedDir;	// already normalized

				// calculate the diffuse reflection probability
				double diff_reflect_prob = max(surfaceClr.x, max(surfaceClr.y, surfaceClr.z));
				// calculate the new photon photonPower
				Vec photonPower;
				photonPower.x = line.photonPower.x * (surfaceClr.x / diff_reflect_prob);
				photonPower.y = line.photonPower.y * (surfaceClr.y / diff_reflect_prob);
				photonPower.z = line.photonPower.z * (surfaceClr.z / diff_reflect_prob);
				line.setPhotonPower(photonPower);

				// multiply by transmittance
				if (insideMedium) {
					line.setPhotonPower(line.photonPower * transmittance(scene.medium->sigma_t, line.distance));
				}
				// trace new light ray

				line = Line(line.end, direction, Vec(), line.photonPower, 0);
				trace(line, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				// calculate surface radiance and multiply by transmittance and russion roulette factor
				if (settings.surfaceRadiance) clr = clr + surfaceRadiance(scene, lightRays, line, surfaceClr / 255.0, N, settings) * rrFactor;
			}
		}

		// Specular BRDF - this is a singularity in the rendering equation that follows
		// delta distribution, therefore we handle this case explicitly - one incoming
		// direction -> one outgoing direction, that is, the perfect reflection direction.
		if (type == SPECULAR) {
			Vec direction = (line.direction - N * (line.direction.dot(N) * 2)).norm();

			if (pass == LIGHTPASS) {
				// multiply by transmittance
				if (insideMedium) {
					line.setPhotonPower(line.photonPower * transmittance(scene.medium->sigma_t, line.distance));
				}
				// create new light ray
				line = Line(line.end, direction, Vec(), line.photonPower, 0);
				trace(line, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				line = Line(line.end, direction, Vec(), Vec(), 0);
				Vec temp;
				trace(line, scene, lightRays, insideMedium, depth + 1, temp, settings, RENDERPASS);
				clr = clr + temp;
			}
		}


		// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
		// to compute the outgoing reflection and refraction directions and probability weights.
		if (type == REFRACTIVE) {
			Vec direction;
			double n = scene.refr_index;
			double R0 = (1.0 - n) / (1.0 + n);
			R0 = R0 * R0;
			if (N.dot(line.direction) > 0) { // we're inside the medium
				N = N * -1;
				n = 1 / n;
			}
			n = 1 / n;
			double cost1 = (N.dot(line.direction)) * -1; // cosine of theta_1
			double cost2 = 1.0 - n * n * (1.0 - cost1 * cost1); // cosine of theta_2
			double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
			if (cost2 > 0 && RND2 > Rprob) { // refraction direction
				direction = ((line.direction * n) + (N * (n * cost1 - sqrt(cost2)))).norm();
			}
			else { // reflection direction
				direction = (line.direction + N * (cost1 * 2)).norm();
			}

			if (pass == LIGHTPASS) {
				// multiply by transmittance
				if (insideMedium) {
					line.setPhotonPower(line.photonPower * transmittance(scene.medium->sigma_t, line.distance));
				}
				// create new light ray
				line = Line(line.end, direction, Vec(), line.photonPower, 0);
				trace(line, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				line = Line(line.end, direction, Vec(), Vec(), 0);
				Vec temp;
				trace(line, scene, lightRays, insideMedium, depth + 1, temp, settings, RENDERPASS);
				clr = clr + temp;
			}
		}
	}
	else if (mediumTransition) {
		Vec mediumTransitionPoint = line.origin + line.direction * intersectionMedium.t;


		if (pass == LIGHTPASS) {
			line.changeEnd(mediumTransitionPoint, intersectionMedium.t);

			if (insideMedium) {
				lightRays.push_back(line);
				line.setPhotonPower(line.photonPower * transmittance(scene.medium->sigma_t, line.distance));
			}

			line = Line(line.end, line.direction, Vec(), line.photonPower, 0);
			trace(line, scene, lightRays, !insideMedium, depth, clr, settings, LIGHTPASS);
		}

		if (pass == RENDERPASS) {
			if (insideMedium) {
				line.changeEnd(mediumTransitionPoint, intersectionMedium.t);
				if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, line, settings) * rrFactor;
			}
			else {
				line.changeEnd(mediumTransitionPoint, 0);
			}
			line = Line(line.end, line.direction, Vec(), Vec(), 0);
			Vec temp;
			trace(line, scene, lightRays, !insideMedium, depth, temp, settings, RENDERPASS);
			clr = clr + temp;
		}
	}

	if (pass == LIGHTPASS && scatteringEvent) {
		if (depth >= settings.bounces - 1)
			return;

		line = original;

		Vec scatteringPoint = line.origin + line.direction * freeFlightDistance;
		line.changeEnd(scatteringPoint, freeFlightDistance);

		double scattering_prob = (scene.medium->sigma_s / scene.medium->sigma_t).average();
		if (RND2 > scattering_prob) return;

		line.setPhotonPower(line.photonPower * transmittance(scene.medium->sigma_t, line.distance));

		// create new ray direction
		Vec direction = sampleHenyeyGreenstein(line.direction, scene.medium->g, RND2, RND2);

		// create new light ray
		line = Line(line.end, direction, Vec(), line.photonPower, 0);
		trace(line, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
	}
}

double mse(VecInt** gtPix, Vec** evalPix, int spp) {
	double sum = 0.0;
	for (int col = 0; col < width; col++) {
		for (int row = 0; row < height; row++) {
			double resX = ((double)(gtPix[col][row].x - color((std::isinf(evalPix[col][row].x) || std::isnan(evalPix[col][row].x)) ? 0.0 : evalPix[col][row].x / (double)spp)));
			double resY = ((double)(gtPix[col][row].y - color((std::isinf(evalPix[col][row].y) || std::isnan(evalPix[col][row].y)) ? 0.0 : evalPix[col][row].y / (double)spp)));
			double resZ = ((double)(gtPix[col][row].z - color((std::isinf(evalPix[col][row].z) || std::isnan(evalPix[col][row].z)) ? 0.0 : evalPix[col][row].z / (double)spp)));

			sum +=
				(resX == 0.0 ? 0.0 : pow(resX, 2.0) +
					resY == 0.0 ? 0.0 : pow(resY, 2.0) +
					resZ == 0.0 ? 0.0 : pow(resZ, 2.0));
		}
	}
	double mseR = sum / (width * height * 3.0);
	return mseR;
}

double psnr(double mseResult) {
	return 20 * log10(255) - 10 * log10(mseResult);
}


vector<int> extractIntegerWords(string str)
{
	stringstream ss;
	ss << str;

	string temp;
	vector<int> allFound;
	int found;

	while (!ss.eof()) {
		ss >> temp;
		if (stringstream(temp) >> found)
			allFound.push_back(found);
		temp = "";
	}
	return allFound;
}

VecInt** readFile(string path) {

	ifstream input;
	input.open(path);

	if (!input) { // file couldn't be opened
		cerr << "Error: file could not be opened" << endl;
	}
	string line;
	int counter = 0;
	string found;
	VecInt** pix = new VecInt * [0];

	for (string line; getline(input, line); )
	{
		if (counter == 0 && line != "P3") {
			// first line has to be "P3"
			break;
		}

		if (counter == 1) {
			// second line is the size
			vector<int> dimensions = extractIntegerWords(line);
			if ((width != 0 && width != dimensions.at(0)) || (height != 0 && height != dimensions.at(1)))
				exit;

			width = dimensions.at(0);
			height = dimensions.at(1);

			pix = new VecInt * [width];
			for (int i = 0; i < width; i++)
				pix[i] = new VecInt[height];
		}

		if (counter == 2 && line != "255") {
			// third line is the max value, has to be 255 for us
			break;
		}

		if (counter > 2) {
			vector<int> imageLine = extractIntegerWords(line);
			int row = counter - 3;
			for (int i = 0; i < imageLine.size(); i += 3) {
				pix[i == 0 ? 0 : i / 3][row] = VecInt(imageLine[i], imageLine[i + 1], imageLine[i + 2]);
			}
		}

		counter++;
	}

	return pix;
}

string leadingZeros(string text, int places = 2) {
	while (text.length() < places)
		text = "0" + text;
	return text;
}

string msToHHMMSS(double time) {
	int hr = time / 3600000;
	time = time - 3600000 * hr;
	int min = time / 60000;
	time = time - 60000 * min;
	int sec = time / 1000;
	int ms = time - 1000 * sec;
	return leadingZeros(to_string(hr)) + ":" + leadingZeros(to_string(min)) + ":" + leadingZeros(to_string(sec)) + ":" + leadingZeros(to_string(ms), 3);
}

double render(RenderSettings& settings) {
	Scene scene;
	loadScene(scene, settings);

	if (settings.sampling == SAMPLE_PRO_GAU) {
		initGaussianProduct(scene);
	}
	else if (settings.sampling == SAMPLE_SUM_CAU) {
		initCauchySum(scene);
	}
	else if (settings.sampling == SAMPLE_MM_CAU) {
		initCauchyMM(scene);
	}
	else if (settings.sampling == SAMPLE_MM_GAU) {
		initGaussianMM(scene);
	}
	else if (settings.sampling == SAMPLE_MM_HYS) {
		initHyperbolicSecantMM(scene);
	}
	else if (settings.sampling == SAMPLE_MM_LOG) {
		initLogisticMM(scene);
	}

	width = 512;
	height = 512;

	Vec** pix = new Vec * [width];
	for (int i = 0; i < width; i++)
		pix[i] = new Vec[height];

	clock_t start = clock();

	int vrls = 0;
	for (int i = 0; i < settings.spp; i++) {
		std::vector<Line> lightRays;
		// for the number of vrl per sample, cast a light ray in and follow it through the scene
		for (int j = 0; j < settings.vrlps; j++) {
			Vec direction;
			if (settings.lightType == LIGHT_DIR) {
				scene.pointlight->pos = scene.areaLights->getRandomPoint(RND2, RND2, RND2);
				direction = scene.areaLights->getDirection();
			}
			else if (settings.lightType == LIGHT_BIDIR) {
				direction = scene.lightDirections->at(j);
			}
			else {
				direction = sampleSphere(RND2, RND2).norm();
			}

			Line line = Line(scene.pointlight->pos, direction, Vec(), scene.pointlight->intensity, 0);
			Vec temp;
			trace(line, scene, lightRays, settings.lightInMedium, 0, temp, settings, LIGHTPASS);
		}
			
		double currentTime = (double)(clock() - start);
		double loopTimeAverage = (currentTime / (double)(i + 1));
		double fullTime = loopTimeAverage * (double)settings.spp;

		fprintf(stdout, string("\r\r\r\rRendering: %f\tCurrent Time: " + msToHHMMSS(currentTime) + (settings.timeComparison != INT_MAX ? "" : "\tFull Time : " + msToHHMMSS(fullTime))).c_str(), i / (double)settings.spp * 100);
		vrls = vrls + lightRays.size();

#pragma omp parallel for schedule(dynamic)
		for (int col = 0; col < width; col++) {
			for (int row = 0; row < height; row++) {
				Vec color;
				Ray ray;
				ray.o = (Vec(0, 0, 0)); // rays start out from here
				Vec cam = camcr(col, row, settings.fov); // construct image plane coordinates
				cam.x = cam.x + RND / 700; // anti-aliasing for free
				cam.y = cam.y + RND / 700;
				ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
				Line line = Line(ray.o, ray.d, Vec(), Vec(), 0);

				trace(line, scene, lightRays, settings.cameraInMedium, 0, color, settings, RENDERPASS);
				pix[col][row] = pix[col][row] + color;
			}
		}

		double time = (double)(clock() - start) / CLOCKS_PER_SEC;

		if (time > settings.timeComparison) {
			settings.spp = i + 1;
			settings.timeComparison = time;
			break;
		}
	}

	clock_t end = clock();
	double time = (double)(end - start);
	fprintf(stdout, "\nTime: %f", time);

	string str_g = to_string(abs(settings.g)).substr(0, 4);
	str_g.erase(remove(str_g.begin(), str_g.end(), '.'), str_g.end());
	std::string name = "vrl_vrlps_" + to_string(settings.vrlps) + "_spp_" + to_string(settings.spp) + "_b_" + to_string(settings.bounces) + (settings.g > 0 ? "_g_+" : "_g_-") + str_g + "_t_" + to_string(settings.sampling);
	if (settings.timeComparison != INT_MAX) name = "vrl_timeComp_vrlps_" + to_string(settings.vrlps) + "_b_" + to_string(settings.bounces) + (settings.g > 0 ? "_g_+" : "_g_-") + str_g + "_t_" + to_string(settings.sampling);
	saveImage(pix, name, settings.spp);

	clear(scene, pix);
	return time;
}


void saveLog(std::string name, std::string log) {
	FILE* f = fopen(std::string(name + ".txt").c_str(), "w");
	fprintf(f, "%s", log.c_str());
	fclose(f);
}


int main(int argc, char* argv[]) {
	RenderSettings settings = RenderSettings();
	string arg1, arg2, help;

	for (int i = 1; i < argc; i += 2) {
		arg1 = argv[i];
		if (i + 1 >= argc && arg1 != "-h") break;

		if (arg1 != "-h")
			arg2 = argv[i + 1];

		if (arg1 == "-h") {
			help = "Usage: " + string(argv[0]) +
				" [ -s <scene_id> | -spp <samples_per_pixel> | -vrlps <vrl_per_sample> | -is <importance_sampling_strategy> | -b <bounces> | -sa <sigma_a> | -ss <sigma_s> | -g <henyey_greenstein> | -i <light_intensity> | -mr <medium_radiance> | -sr <surface_radiance> ]\n\n";
			help += "\t-s\tselects a scene (default 0):\n\t\t\t0\tcornell box with directional lights\n\t\t\t1\tbidirectional VRLs\n\t\t\t2\tunidirectional VRLs\n\t\t\t3\tmedium radiance test\n\t\t\t4\tsurface radiance test\n\t\t\t5\tmedia cube\n";
			help += "\t-spp\tsamples per pixel\n";
			help += "\t-vrlps\tVRLs per sample\n";
			help += "\t-is\tselects the importance sampling techniques for the light / camera rays (default 4):\n\t\t\t0\tuniform / uniform\n\t\t\t1\texponential / exponential\n\t\t\t2\tuniform / equi-angular\n\t\t\t3\tsimple joint distribution\n\t\t\t4\tadvanced joint distribution\n\t\t\t5\tproduct of Gaussian distributions\n\t\t\t6\tproduct of Cauchy distributions\n\t\t\t7\tMixture Model with Cauchy distributions\n\t\t\t8\tMixture Model with Gaussian distributions\n\t\t\t9\tMixture Model with Hyperbolic Secant distributions\n\t\t\t10\tMixture Model with Logistic distributions\n";
			help += "\t-b\tnumber of allowed light bounces\n";
			help += "\t-sa\tsigma_a\n";
			help += "\t-ss\tsigma_s\n";
			help += "\t-g\thenyey-greestein\n";
			help += "\t-i\tthe intensity of the light\n";
			help += "\t-mr\toption for medium radiance (default true): true / false\n";
			help += "\t-sr\toption for surface radiance (default true): true / false\n";
			fprintf(stdout, help.c_str());
			return 0;
		}
		else if (arg1 == "-s") {
			settings.scene = max(0, min(5, stoi(arg2)));
		}
		else if (arg1 == "-spp") {
			settings.spp = max(0, min(INT_MAX, stoi(arg2)));
		}
		else if (arg1 == "-vrlps") {
			settings.vrlps = max(0, min(INT_MAX, stoi(arg2)));
		}
		else if (arg1 == "-time") {
			settings.spp = max(0, min(INT_MAX, stoi(arg2)));
		}
		else if (arg1 == "-is") {
			settings.sampling = max(0, min(10, stoi(arg2)));
		}
		else if (arg1 == "-b") {
			settings.bounces = max(0, min(INT_MAX, stoi(arg2)));
		}
		else if (arg1 == "-sa") {
			settings.sigma_a = max(0.0, min(100.0, stod(arg2)));
		}
		else if (arg1 == "-ss") {
			settings.sigma_s = max(0.0, min(100.0, stod(arg2)));
		}
		else if (arg1 == "-g") {
			settings.g = max(-1.0, min(1.0, stod(arg2)));
		}
		else if (arg1 == "-i") {
			settings.intensity = max(0.0, min(inf, stod(arg2)));
		}
		else if (arg1 == "-mr") {
			settings.mediumRadiance = arg2 == "t";
		}
		else if (arg1 == "-sr") {
			settings.surfaceRadiance = arg2 == "t";
		}
	}

	if (settings.timeComparison != INT_MAX)
		settings.spp = INT_MAX;

	render(settings);
	getchar();
}