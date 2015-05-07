#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = true;

    timeStep = 0.001;
    NewtonMaxIters = 20;
    NewtonTolerance = 1e-8;

    clickMode = CM_ADDDENSITY;
    densityRadius = 2;
    densityMagnitude = 1000;
    velocityRadius = 3;
    velocityMagnitude = 1000000;
    diffusionConstant = 0.03;
    viscosityFluid = 0.2;

    connector = CT_SPRING;

    springStiffness = 100;
    maxSpringStrain = 0.2;
    dampingStiffness = 1.0;

    particleMass = 1.0;
    particleFixed = false;
    maxSpringDist = 0.25;

    rodDensity = 2.0;
    rodStretchStiffness = 100.0;
    rodBendingStiffness = 0.05;
    rodSegments = 5;

    ropeDensity = 2.0;
    ropeBend = 0.01;
    ropeSegments = 5;
}
