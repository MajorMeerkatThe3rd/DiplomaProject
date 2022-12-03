#pragma once
#include "typeDefinitions.h"

// adjusted version of: https://people.sc.fsu.edu/~jburkardt/cpp_src/sphere_fibonacci_grid/sphere_fibonacci_grid.html
// to generate equally spaces samples on a sphere
void fibonacciSamples(int amount, std::vector<Vec>& directions);

void loadScene(Scene& scene, RenderSettings& s);