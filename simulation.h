#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "simparameters.h"
#include <QGLWidget>
#include "fluid.h"


class SimParameters;

struct Plane
{
    Eigen::Vector3d pos;
    Eigen::Vector3d normal;
};

class Simulation
{
public:
    Simulation(const SimParameters &params);
    ~Simulation();

    void takeSimulationStep();

    // FLUID STUFF
///////////////////////////////////////////////
   void render();
   void fluidSimulationStep();
   void advection(int boundary, Eigen::VectorXf &d, Eigen::VectorXf &dOld, Eigen::VectorXf &uOld, Eigen::VectorXf &vOld, Eigen::VectorXf &zOld);
   void diffuse(int boundary, Eigen::VectorXf &d, Eigen::VectorXf &dOld, double diffFactor);
   void linearSolver(int b, Eigen::VectorXf &x, Eigen::VectorXf &xOld, double a, double c);
   void project(Eigen::VectorXf &x, Eigen::VectorXf &y,Eigen::VectorXf &z, Eigen::VectorXf &xOld, Eigen::VectorXf &yOld, Eigen::VectorXf &zOld);
   void swap(Eigen::VectorXf &left, Eigen::VectorXf &right);
   void addSource(Eigen::VectorXf &d, Eigen::VectorXf& dOld);
   void addVelocity(int sourceNo);
   void addDensity(int sourceNo);
   void setBoundry(int b, Eigen::VectorXf& m);
   void addBuoyancy();
   void smokeEffect();
///////////////////////////////////////////////


    void renderPlanes(bool transparent);
    void clearScene();

private:

    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;
    GLuint wallTex_;

    ////FLUID STUFF
    Fluid *fluid_;
    ////////////


};

#endif // SIMULATION_H
