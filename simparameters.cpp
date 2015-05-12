#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;

    timeStep = 0.001;

    densityMagnitude = 200;
    velocityMagnitude = 20;
    diffusionConstant = 0.02;
    viscosityFluid = 0.2;

    buoyancy = 4.0;

    source1 = true;
    source2 = false;
    source3 = false;
}
