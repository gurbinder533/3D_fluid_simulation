#ifndef FLUID_H
#define FLUID_H

#include <Eigen/Core>
#include <vector>
#include <string>


class Fluid
{
public:
    Fluid();

    Eigen::MatrixXd fluidDensity, fluidDensityOld;
    Eigen::MatrixXd vy, vyOld;
    Eigen::MatrixXd vx, vxOld;

    int n;
    double sizeOfVoxel;
    bool debug;

    void zeroEverything();
    double getTotalDensity();
    void render();

};

#endif // FLUID_H
