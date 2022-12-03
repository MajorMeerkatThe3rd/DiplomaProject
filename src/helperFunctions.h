#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include "typeDefinitions.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <filesystem>

// find the shortest point on a line to a point p
void closestPointOnLine(SamplingData& samplingData);

// shortest Distance between two line segments adapted from
// implementation from user Nick on Matlab (http://de.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments)
// adapts the algorithm found on dan sunday's website (http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment)
void shortestPointsOnLineSegments(Line line1, Line line2, SamplingData& samplingDataVPL, SamplingData& samplingData);

// given v1, set v2 and v3 so they form an orthonormal system
// (we assume v1 is already normalized)
void ons(const Vec& v1, Vec& v2, Vec& v3);

// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y, double fov);

double arccos(double x);

// rotate a Vector v around another Vector n
Vec rotateVec(Vec n, Vec v, double theta);

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec sampleHemisphere(double u1, double u2);

Vec sampleSphere(double u1, double u2);

Vec transmittance(Vec sigma_t, double dist);

Vec sampleHenyeyGreenstein(Vec incomingDir, double g, double u1, double u2);

double pfHenyeyGreenstein(Vec v1, Vec v2, double g);

double pfHenyeyGreenstein(double dot, double g);

void clear(Scene& scene, Vec** pix);

int color(double clr);

void saveImage(Vec** pix, std::string name, int spp);