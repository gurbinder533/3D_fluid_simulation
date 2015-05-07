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


    ConnectorType connector;

    double springStiffness;
    double maxSpringStrain;
    double dampingStiffness;

    double particleMass;
    bool particleFixed;
    double maxSpringDist;

    double rodDensity;
    double rodStretchStiffness;
    double rodBendingStiffness;
    int rodSegments;

    double ropeDensity;
    double ropeBend;
    double ropeSegments;
};

#endif // SIMPARAMETERS_H
