#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include <QPainter>
#include "vectormath.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include "mesh.h"
#include <iostream>
const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0), wallTex_(0)
{
     fluid_ = new Fluid();
}

Simulation::~Simulation()
{
    clearScene();
}


void Simulation::renderPlanes(bool transparent)
{
    renderLock_.lock();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glDisable(GL_BLEND);
    glPushMatrix();

     glLineWidth(3.0);
        glBegin(GL_LINES);

        glColor3f(0.4, 0.4, 0.4);

        glVertex3f(-1, -1, -1);
        glVertex3f(-1, 1, -1);

        glVertex3f(-1, -1, -1);
        glVertex3f(1, -1, -1);

        glVertex3f(1, -1, -1);
        glVertex3f(1, 1, -1);

        glVertex3f(-1, 1, -1);
        glVertex3f(1, 1, -1);

        glVertex3f(-1, -1, 1);
        glVertex3f(-1, 1, 1);

        glVertex3f(-1, -1, 1);
        glVertex3f(1, -1, 1);

        glVertex3f(1, -1, 1);
        glVertex3f(1, 1, 1);

        glVertex3f(-1, 1, 1);
        glVertex3f(1, 1, 1);

        glVertex3f(-1, -1, -1);
        glVertex3f(-1, -1, 1);

        glVertex3f(-1, 1, -1);
        glVertex3f(-1, 1, 1);

        glVertex3f(1, -1, -1);
        glVertex3f(1, -1, 1);

        glVertex3f(1, 1, -1);
        glVertex3f(1, 1, 1);

        glEnd();

    glPopMatrix();
    glEnable(GL_BLEND);

    renderLock_.unlock();
}



void Simulation::render()
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;
    double baselinewidth = 0.5;
    int numcirclewedges = 20;

    renderLock_.lock();
    {
        fluid_->render();
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    fluidSimulationStep();
    time_ += params_.timeStep;
}

void Simulation::fluidSimulationStep()
{
    // Velocity Code
    addSource(fluid_->vx3d, fluid_->vx3dOld);
    addSource(fluid_->vy3d, fluid_->vy3dOld);
    addSource(fluid_->vz3d, fluid_->vz3dOld);
    addBuoyancy();
    smokeEffect();




/*******************ADVECTION*********************************************/
    this->swap(fluid_->vx3d, fluid_->vx3dOld);
    this->swap(fluid_->vy3d, fluid_->vy3dOld);
    this->swap(fluid_->vy3d, fluid_->vz3dOld);

    //TODO: we need to CHANGE THIS
    advection(1, fluid_->vx3d, fluid_->vx3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);
    advection(2, fluid_->vy3d, fluid_->vy3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);
    advection(3, fluid_->vz3d, fluid_->vz3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);

    project(fluid_->vx3d, fluid_->vy3d, fluid_->vz3d, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);

/*******************ADVECTION*********************************************/

/************************DIFFUSE******************************************/
    this->swap(fluid_->vx3d, fluid_->vx3dOld);
    this->swap(fluid_->vy3d, fluid_->vy3dOld);
    this->swap(fluid_->vz3d, fluid_->vz3dOld);

    diffuse(1, fluid_->vx3d, fluid_->vx3dOld, params_.viscosityFluid);
    diffuse(2, fluid_->vy3d, fluid_->vy3dOld, params_.viscosityFluid);
    diffuse(3, fluid_->vz3d, fluid_->vz3dOld, params_.viscosityFluid);

    project(fluid_->vx3d, fluid_->vy3d, fluid_->vz3d, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);
/************************DIFFUSE****************************************/


    //Density Code
    addSource(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);


    this->swap(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);
    diffuse(0, fluid_->fluidDensity3d, fluid_->fluidDensity3dOld, params_.diffusionConstant);


    this->swap(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);
    advection(0, fluid_->fluidDensity3d, fluid_->fluidDensity3dOld, fluid_->vx3d, fluid_->vy3d, fluid_->vz3d);




    fluid_->fluidDensity3dOld.setZero();
    fluid_->vx3dOld.setZero();
    fluid_->vy3dOld.setZero();
    fluid_->vz3dOld.setZero();


    cout << "fluid density: " << fluid_->getTotalDensity() << endl;
}

void Simulation::addSource(Eigen::VectorXf &d, Eigen::VectorXf &dOld)
{

    d += params_.timeStep * dOld;

}

void Simulation::diffuse(int boundary, VectorXf &d, VectorXf &dOld, double diffFactor)
{
    double a = params_.timeStep * diffFactor * fluid_->n * fluid_->n * fluid_->n;

    for(int loop = 0; loop < 20; loop++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                for(int k = 1; k <= fluid_->n; k++)
                {

                    d[COFF(i,j,k)] = (dOld[COFF(i,j,k)] + a * (d[COFF(i-1,j,k)] + d[COFF(i+1,j,k)] +
                                         d[COFF(i,j-1,k)] + d[COFF(i,j+1,k)] + d[COFF(i,j,k-1)] + d[COFF(i,j,k+1)]))/(1+6*a);

                }
            }
        }
        setBoundry(boundary, d);
    }
}

void Simulation::project(Eigen::VectorXf &x, Eigen::VectorXf &y, Eigen::VectorXf &z, Eigen::VectorXf &xOld, Eigen::VectorXf &yOld, Eigen::VectorXf &zOld)
{
    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            for(int k = 1; k <= fluid_->n; k++)
            {

                yOld[COFF(i,j,k)] = (x[COFF(i+1,j,k)] - x[COFF(i-1,j,k)]
                                    + y[COFF(i,j+1,k)] - y[COFF(i,j-1,k)]
                                    + z[COFF(i,j,k+1)] - z[COFF(i,j,k-1)]) * (-1)/(3*fluid_->n);


            }

        }
    }

    xOld.setZero();

    setBoundry(0,yOld);
    setBoundry(0,xOld);

    for(int loop = 0; loop < 20; loop++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                for(int k = 1; k <= fluid_->n; k++)
                {
                    xOld[COFF(i,j,k)] =  (xOld[COFF(i-1,j,k)] + xOld[COFF(i+1,j,k)]
                                         + xOld[COFF(i,j-1,k)] + xOld[COFF(i,j+1,k)]
                                         + xOld[COFF(i,j,k-1)] + xOld[COFF(i,j,k+1)]
                                         + yOld[COFF(i,j,k)]) / 6 ;
                }
            }
        }
        setBoundry(0, xOld);
    }


    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            for(int k = 1; k <= fluid_->n; k++)
            {
                x[COFF(i,j,k)] -= (1/3) * fluid_->n * (xOld[COFF(i+1,j,k)] - xOld[COFF(i-1,j,k)]);
                y[COFF(i,j,k)] -= (1/3) * fluid_->n * (yOld[COFF(i,j+1,k)] - zOld[COFF(i,j-1,k)]);
                x[COFF(i,j,k)] -= (1/3) * fluid_->n * (zOld[COFF(i,j,k+1)] - zOld[COFF(i,j,k-1)]);
            }
        }
    }
    setBoundry(1,x);
    setBoundry(2,y);
}

void Simulation::setBoundry(int b, VectorXf& x)
{


           int i, j;
           int N = fluid_->n;
           for (i=1; i<= fluid_->n; i++)
           {
                   for (j=1; j<= fluid_->n; j++) {
                           x[COFF(0,i,j)]    = (b==1) ? -x[COFF(1,i,j)] : x[COFF(1,i,j)];
                           x[COFF(N+1,i,j)]  = (b==1) ? -x[COFF(N,i,j)] : x[COFF(N,i,j)];
                           x[COFF(i,0,j)]    = (b==2) ? -x[COFF(i,1,j)] : x[COFF(i,1,j)];
                           x[COFF(i,N+1,j)]  = (b==2) ? -x[COFF(i,N,j)] : x[COFF(i,N,j)];
                           x[COFF(i,j,0)]    = (b==3) ? -x[COFF(i,j,1)] : x[COFF(i,j,1)];
                           x[COFF(i,j,N+1)]  = (b==3) ? -x[COFF(i,j,N)] : x[COFF(i,j,N)];
                   }
           }

           x[COFF(0,0,0)]       = (x[COFF(1,0,0)]    +x[COFF(0,1,0)]    +x[COFF(0,0,1)])    /3;
           x[COFF(0,N+1,0)]     = (x[COFF(1,N+1,0)]  +x[COFF(0,N,0)]    +x[COFF(0,N+1,1)])  /3;
           x[COFF(N+1,0,0)]     = (x[COFF(N,0,0)]    +x[COFF(N+1,1,0)]  +x[COFF(N+1,0,1)])  /3;
           x[COFF(N+1,N+1,0)]   = (x[COFF(N,N+1,0)]  +x[COFF(N+1,N,0)]  +x[COFF(N+1,N+1,1)])/3;
           x[COFF(0,0,N+1)]     = (x[COFF(1,0,N+1)]  +x[COFF(0,1,N+1)]  +x[COFF(0,0,N)])    /3;
           x[COFF(0,N+1,N+1)]   = (x[COFF(1,N+1,N+1)]+x[COFF(0,N,N+1)]  +x[COFF(0,N+1,N)])  /3;
           x[COFF(N+1,0,N+1)]   = (x[COFF(N,0,N+1)]  +x[COFF(N+1,1,N+1)]+x[COFF(N+1,0,N)])  /3;
           x[COFF(N+1,N+1,N+1)] = (x[COFF(N,N+1,N+1)]+x[COFF(N+1,N,N+1)]+x[COFF(N+1,N+1,N)])/3;
}

void Simulation::swap(VectorXf &left, VectorXf &right)
{
    VectorXf tmp(left);
    left = right;
    right = tmp;
}

void Simulation::linearSolver(int b, Eigen::VectorXf &x, Eigen::VectorXf &xOld, double a, double c)
{
    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                x.coeffRef(i,j) = (a * ( x.coeff(i-1, j) + x.coeff(i+1, j)
                                         + x.coeff(i,j-1) + x.coeff(i,j+1))
                                         + xOld.coeff(i,j)) / c ;
            }
        }
        setBoundry(b, x);
    }
}

void Simulation::advection(int boundry, VectorXf &d, VectorXf &dOld, VectorXf &xOld, VectorXf &yOld, VectorXf &zOld)
{
    double dt0 = params_.timeStep * fluid_->n;
    int n = fluid_->n;
    float v0, v1;

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            for(int k = 1; k <= n; k++)
            {
                double x = i - dt0 * xOld[COFF(i,j,k)];
                double y = j - dt0 * yOld[COFF(i,j,k)];
                double z = k - dt0 * zOld[COFF(i,j,k)];

                if(x > n + 0.5f)
                {
                    x = n + 0.5f;
                }
                if(x < 0.5f)
                {
                    x = 0.5f;
                }

                int i0 = (int) x;
                int i1 = i0 + 1;

                if(y > n + 0.5f)
                {
                    y = n + 0.5f;
                }
                if(y < 0.5f)
                {
                    y = 0.5f;
                }

                int j0 = (int) y;
                int j1 = j0 +1;



                if(z > n + 0.5f)
                {
                    z = n + 0.5f;
                }
                if(z < 0.5f)
                {
                    z = 0.5f;
                }

                int k0 = (int) z;
                int k1 = k0 + 1;


                double sx1 = x - i0;
                double sx0 = 1 - sx1;
                double sy1 = y - j0;
                double sy0 = 1 - sy1;
                double sz1 = z - k0;
                double sz0 = 1 - sz1;

                v0 = sx0 * (sy0 * dOld[COFF(i0, j0, k0)] + sy1 * dOld[COFF(i0,j1,k0)]) + sx1*(sy0*dOld[COFF(i1,j0,k0)] + sy1*dOld[COFF(i1,j1,k0)]);
                v1 = sx0 * (sy0 * dOld[COFF(i0, j0, k1)] + sy1 * dOld[COFF(i0,j1,k1)]) + sx1*(sy0*dOld[COFF(i1,j0,k1)] + sy1*dOld[COFF(i1,j1,k1)]);

                d[COFF(i,j,k)] = sz0*v0 + sz1*v1;
           }
        }
    }
    setBoundry(boundry, d);
}

void Simulation::addDensity(int sourceNo)
{
    if (sourceNo == 1)
    {
        fluid_->type = 0;
        int i = 1;
        int j = 1;
        int k = 1;
        cout << "HERE1 : " << COFF(i,j,k) << endl;
        fluid_->fluidDensity3dOld[COFF(i,j,k)] += params_.densityMagnitude;

    }
    if (sourceNo == 2)
    {
        fluid_->type = 1;
        int i = fluid_->n-1;
        int j = fluid_->n-1;
        int k = fluid_->n-1;
        cout << "HERE2 : " << COFF(i,j,k) << endl;
        fluid_->fluidDensity3dOld[COFF(i,j,k)] += params_.densityMagnitude;

    }
    if (sourceNo == 3)
    {
        fluid_->type = 1;
        int i = (fluid_->n-1)/2;
        int j = (fluid_->n-1)/2;
        int k = (fluid_->n-1)/2;
        cout << "HERE2 : " << COFF(i,j,k) << endl;
        fluid_->fluidDensity3dOld[COFF(i,j,k)] += params_.densityMagnitude;

    }
}

void Simulation::addVelocity(int sourceNo, double velX, double velY, double velZ)
{
    if(sourceNo == 1)
    {
        int i = (fluid_->n -1)/2;
        int j = (fluid_->n - 1)/2;
        int k = (fluid_->n -1)/2;

        int n = fluid_->n;
        //std::cout << "prev: velx " << fluid_->vx3dOld[COFF(i,j,k)] << std::endl;
        for(i = 0; i <= n; i++)
        {
           for(j = 0; j <= n; j++)
            {
            for(k = 0; k <= n; ++k)
            {
             fluid_->vy3dOld[COFF(i,j,k)] += params_.velocityMagnitude;
            }
            }
        }
        //std::cout << ": velz " << fluid_->vz3dOld[COFF(i,j,k)] << " params: " << params_.velocityMagnitude << std::endl;

    }


}



void Simulation::addBuoyancy()
{
    fluid_->vz3d += -fluid_->fluidDensity3d*params_.buoyancy*params_.timeStep;
}

void Simulation::smokeEffect()
{
    int index;
    float x,y,z;
    Eigen::VectorXf cx,cy,cz,cd;
    int size = (fluid_->n + 2)*(fluid_->n + 2)*(fluid_->n + 2);
    cx.resize(size);
    cy.resize(size);
    cz.resize(size);
    cd.resize(size);


    int n = fluid_->n;
    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            for(int k = 1; k <= n; k++)
            {
                index = COFF(i,j,k);
                x = cx[index] = (fluid_->vz3d[COFF(i,j+1,k)] - fluid_->vz3d[COFF(i,j-1,k)] )*0.5f -
                                 (fluid_->vy3d[COFF(i,j,k+1)] - fluid_->vy3d[COFF(i,j,k-1)] )*0.5f;

                y = cy[index] = (fluid_->vx3d[COFF(i,j,k+1)] - fluid_->vx3d[COFF(i,j,k-1)] )*0.5f -
                                 (fluid_->vz3d[COFF(i+1,j,k)] - fluid_->vz3d[COFF(i-1,j,k)] )*0.5f;

                z = cz[index] = (fluid_->vy3d[COFF(i+1,j,k)] - fluid_->vy3d[COFF(i-1,j,k)] )*0.5f -
                                 (fluid_->vx3d[COFF(i,j+1,k)] - fluid_->vx3d[COFF(i,j-1,k)] )*0.5f;

                cd[index] = sqrt(x*x + y*y + z*z);

            }

        }
    }

    float Nx, Ny, Nz, len;
    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            for(int k = 1; k <= n; k++)
            {
                index = COFF(i,j,k);
                 Nx = cd[COFF(i+1,j,k)] - cd[COFF(i-1,j,k)]*0.5f;
                 Nx = cd[COFF(i,j+1,k)] - cd[COFF(i,j-1,k)]*0.5f;
                 Nx = cd[COFF(i,j,k+1)] - cd[COFF(i,j,k-1)]*0.5f;

                 len = 1.0f/(sqrt(Nx*Nx + Ny*Ny + Nz*Nz) + 0.0000000001f);
                 Nx *= len;
                 Ny *= len;
                 Nz *= len;

                 fluid_->vx3d[index] += (Ny*cz[index] - Nz*cy[index]) * params_.timeStep;
                 fluid_->vy3d[index] += (Nz*cx[index] - Nx*cz[index]) * params_.timeStep;
                 fluid_->vz3d[index] += (Nx*cy[index] - Ny*cx[index]) * params_.timeStep;


            }
        }
    }
}


void Simulation::clearScene()
{
    renderLock_.lock();
    {
        fluid_->zeroEverything();
    }
    renderLock_.unlock();
}
