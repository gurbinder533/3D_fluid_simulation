#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;

    timeStep = 0.001;
    NewtonMaxIters = 20;
    NewtonTolerance = 1e-8;

    clickMode = CM_ADDDENSITY;
    densityRadius = 2;
    densityMagnitude = 200;
    velocityRadius = 3;
    velocityMagnitude = 1000000;
    diffusionConstant = 0.03;
    viscosityFluid = 0.2;

    buoyancy = 4.0;

    source1 = true;
    source2 = false;
}
