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

typedef Eigen::Triplet<double> Tr;

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
    void initializeGL();

    // FLUID STUFF
///////////////////////////////////////////////
   void render();
   void fluidSimulationStep();
   void advection(int boundary, Eigen::MatrixXd &d, Eigen::MatrixXd &dOld, Eigen::MatrixXd &uOld, Eigen::MatrixXd &vOld);
   void diffuse(int boundary, Eigen::MatrixXd &d, Eigen::MatrixXd &dOld, double diffFactor);
   void linearSolver(int b, Eigen::MatrixXd &x, Eigen::MatrixXd &xOld, double a, double c);
   void project(Eigen::MatrixXd &x, Eigen::MatrixXd &y, Eigen::MatrixXd &xOld, Eigen::MatrixXd &yOld);
   void swap(Eigen::MatrixXd &left, Eigen::MatrixXd &right);
   void addSource(Eigen::MatrixXd &d, Eigen::MatrixXd& dOld);
   void addVelocity(double x, double y, double velX, double velY);
   void addDensity(double x, double y);
   void setBoundry(int b, Eigen::MatrixXd& m);
///////////////////////////////////////////////


    void renderPlanes(bool transparent);
    void clearScene();

private:
    void loadFloorTexture();
    void loadWallTexture();

    void renderPlane(const Plane &p, bool isFloor);


    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;
    GLuint wallTex_;

    ////FLUID STUFF
    Fluid *fluid_;
    void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);
    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H, Eigen::VectorXd &v);
    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
    Eigen::SparseMatrix<double> computeGradGTranspose(const Eigen::VectorXd &q);

    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void computeStepProject(Eigen::VectorXd &q, Eigen::VectorXd &oldq, Eigen::VectorXd &v);
    void computeLagrangeMultipliers(const Eigen::VectorXd &q, const Eigen::VectorXd &F, Eigen::VectorXd &v);


    ////////////


    std::vector<Plane> planes_;
};

#endif // SIMULATION_H
