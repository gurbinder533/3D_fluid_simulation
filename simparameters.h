#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters();

    enum ClickMode {CM_ADDVELOCITY, CM_ADDDENSITY};
    enum ConnectorType{ CT_SPRING, CT_RIGID_ROD, CT_FLEXIBLE_ROD, CT_ROPE};

    bool simRunning;

    double timeStep;
    double NewtonTolerance;
    int NewtonMaxIters;

    ClickMode clickMode;
    int densityRadius;
    double densityMagnitude;
    int velocityRadius;
    double velocityMagnitude;
    double diffusionConstant;
    double viscosityFluid;

    float buoyancy;

    bool source1;
    bool source2;
    bool source3;
    bool source4;

    bool velsource1;
    bool velsource2;
    bool velsource3;


};

#endif // SIMPARAMETERS_H
