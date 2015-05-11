#ifndef FLUID_H
#define FLUID_H

#include <Eigen/Core>
#include <vector>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//#define COFF(x,y,z) (((z)<<12) + ((y)<<6)+x)
#define COFF(x,y,z) (((z)<<10) + ((y)<<5)+x)
//#define COFF(x,y,z) (((z)*16*16) + ((y)*16)+x)

class Fluid
{
public:
    Fluid();

    Eigen::VectorXf fluidDensity3d, fluidDensity3dOld;
    Eigen::VectorXf vy3d, vy3dOld;
    Eigen::VectorXf vx3d, vx3dOld;
    Eigen::VectorXf vz3d, vz3dOld;



    int n;
    int size;
    double sizeOfVoxel;
    bool debug;
    int type;

    void zeroEverything();
    double getTotalDensity();
    void render();

};

#endif // FLUID_H
