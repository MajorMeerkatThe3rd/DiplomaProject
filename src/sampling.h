#include "typeDefinitions.h"

void sampleUniform(SamplingData samplingData, double u, Vec& sample, double& pdf);
void sampleExponential(SamplingData samplingData, double u, Vec& sample, double& pdf);
void sampleEquiAngular(SamplingData samplingData, double u, Vec& sample, double& pdf);
void sampleVRLAnisotropic(SamplingData samplingData, double u, Vec& sample, double& pdf);
void sampleVRLLight(SamplingData samplingData, double u, Vec& sample, double& pdf);
void initGaussianProduct(Scene scene);
void sampleGaussianProduct(SamplingData samplingData, double random, Vec& sample, double& pdf);
void initCauchySum(Scene scene);
void sampleCauchySum(SamplingData samplingData, double random, Vec& sample, double& pdf);
void initCauchyMM(Scene scene);
void sampleCauchyMM(SamplingData samplingData, double random, Vec& sample, double& pdf);
void initGaussianMM(Scene scene);
void sampleGaussianMM(SamplingData samplingData, double random, Vec& sample, double& pdf);
void initHyperbolicSecantMM(Scene scene);
void sampleHyperbolicSecantMM(SamplingData samplingData, double random, Vec& sample, double& pdf);
void initLogisticMM(Scene scene);
void sampleLogisticMM(SamplingData samplingData, double random, Vec& sample, double& pdf);