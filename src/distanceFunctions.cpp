#include "helperFunctions.h"

// find the shortest point on a line to a point p
void closestPointOnLine(SamplingData& samplingData) {
	Vec dir_point_origin = samplingData.interestPoint - samplingData.lineStart;
	Vec dir_line = samplingData.direction_line * samplingData.distance_line;

	double length = samplingData.distance_line;
	samplingData.dist_lineOrigin_closestPoint = dir_point_origin.dot(dir_line) / (length * length);

	samplingData.closestPoint = samplingData.lineStart + dir_line * samplingData.dist_lineOrigin_closestPoint;
	samplingData.dist_lineOrigin_closestPoint = length * samplingData.dist_lineOrigin_closestPoint;
}

// shortest Distance between two line segments adapted from
// implementation from user Nick on Matlab (http://de.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments)
// adapts the algorithm found on dan sunday's website (http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment)
void shortestPointsOnLineSegments(Line line1, Line line2, SamplingData& samplingDataVPL, SamplingData& samplingData) {
	Vec dir_line1 = line1.direction * line1.distance;
	Vec dir_line2 = line2.direction * line2.distance;
	Vec dir_origin2_origin1 = line1.origin - line2.origin;

	double dot_dir_l1, dot_dir_l1_l2, dot_dir_l2, dot_dir_l1_o2o1, dot_dir_l2_o2o1, denominator, l1_denominator, l2_denominator, l1_numerator, l2_numerator, l1_fraction, l2_fraction;
	dot_dir_l1 = dir_line1.dot(dir_line1);
	dot_dir_l1_l2 = dir_line1.dot(dir_line2);
	dot_dir_l2 = dir_line2.dot(dir_line2);
	dot_dir_l1_o2o1 = dir_line1.dot(dir_origin2_origin1);
	dot_dir_l2_o2o1 = dir_line2.dot(dir_origin2_origin1);
	denominator = dot_dir_l1 * dot_dir_l2 - dot_dir_l1_l2 * dot_dir_l1_l2;

	if (denominator < eps) {
		l1_numerator = 0.0;
		l1_denominator = 1.0;
		l2_numerator = dot_dir_l2_o2o1;
		l2_denominator = dot_dir_l2;
	}
	else {
		l1_denominator = denominator;
		l2_denominator = denominator;
		l1_numerator = (dot_dir_l1_l2 * dot_dir_l2_o2o1 - dot_dir_l2 * dot_dir_l1_o2o1);
		l2_numerator = (dot_dir_l1 * dot_dir_l2_o2o1 - dot_dir_l1_l2 * dot_dir_l1_o2o1);
		if (l1_numerator < 0.0) {
			l1_numerator = 0.0;
			l2_numerator = dot_dir_l2_o2o1;
			l2_denominator = dot_dir_l2;
		}
		else if (l1_numerator > l1_denominator) {
			l1_numerator = l1_denominator;
			l2_numerator = dot_dir_l2_o2o1 + dot_dir_l1_l2;
			l2_denominator = dot_dir_l2;
		}
	}

	if (l2_numerator < 0.0) {
		l2_numerator = 0.0;
		if (-dot_dir_l1_o2o1 < 0.0) {
			l1_numerator = 0.0;
		}
		else if (-dot_dir_l1_o2o1 > dot_dir_l1) {
			l1_numerator = l1_denominator;
		}
		else {
			l1_numerator = -dot_dir_l1_o2o1;
			l1_denominator = dot_dir_l1;
		}
	}
	else if (l2_numerator > l2_denominator) {
		l2_numerator = l2_denominator;
		if ((-dot_dir_l1_o2o1 + dot_dir_l1_l2) < 0.0) {
			l1_numerator = 0;
		}
		else if ((-dot_dir_l1_o2o1 + dot_dir_l1_l2) > dot_dir_l1) {
			l1_numerator = l1_denominator;
		}
		else {
			l1_numerator = (-dot_dir_l1_o2o1 + dot_dir_l1_l2);
			l1_denominator = dot_dir_l1;
		}
	}

	if (abs(l1_numerator) < eps) {
		l1_fraction = 0.0;
	}
	else {
		l1_fraction = l1_numerator / l1_denominator;
	}
	if (abs(l2_numerator) < eps) {
		l2_fraction = 0.0;
	}
	else {
		l2_fraction = l2_numerator / l2_denominator;
	}

	Vec lineClosest = dir_origin2_origin1 + (dir_line1 * l1_fraction) - (dir_line2 * l2_fraction);

	samplingDataVPL.closestPoint = line1.origin + dir_line1 * l1_fraction;
	samplingDataVPL.dist_lineOrigin_closestPoint = l1_fraction * line1.distance;

	samplingData.closestPoint = line2.origin + dir_line2 * l2_fraction;
	samplingData.dist_lineOrigin_closestPoint = l2_fraction * line2.distance;

	samplingDataVPL.dist_both_lines = lineClosest.length();
}