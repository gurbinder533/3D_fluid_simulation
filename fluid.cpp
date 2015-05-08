#include "fluid.h"
#include <QGLWidget>
#include <iostream>

using namespace Eigen;
using namespace std;

Fluid::Fluid()
{

    this->n = 100;
    this->fluidDensity.resize(this->n+1, this->n+1);
    this->fluidDensityOld.resize(this->n+1, this->n+1);
    this->vx.resize(this->n+1, this->n+1);
    this->vxOld.resize(this->n+1, this->n+1);
    this->vy.resize(this->n+1, this->n+1);
    this->vyOld.resize(this->n+1, this->n+1);

    fluidDensity.setZero();
    fluidDensityOld.setZero();
    vx.setZero();
    vy.setZero();
    vxOld.setZero();
    vyOld.setZero();

    this->sizeOfVoxel = 2.0/(this->n+1);
    this->debug = true;
}

void Fluid::render()
{

    std::cout << "FLUID RENDER\n";
    double xCell = -1;
    double yCell = 1;
    for(int i =0; i <= n; i++)
    {
        xCell = i * this->sizeOfVoxel - 1 ;
        for(int j = 0; j <= n; j++)
        {
            yCell = 1 - j * this->sizeOfVoxel;
            double dens = fluidDensity.coeff(i,j);
            //glColor3f(1-dens,1-dens,1);
            glColor3f(0,0,1);
            glBegin(GL_QUADS);
            {
                glVertex3f(xCell,yCell, 6.0);
                glVertex3f(xCell + this->sizeOfVoxel, yCell, 0.0);
                glVertex3f(xCell + this->sizeOfVoxel, yCell - this->sizeOfVoxel, 0.0);
                glVertex3f(xCell, yCell - this->sizeOfVoxel, 0.0);
            }
            glEnd();
        }
    }
}

double Fluid::getTotalDensity()
{
    double densTotal = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            densTotal += fluidDensity.coeff(i,j);
        }
    }
    return densTotal;
}

void Fluid::zeroEverything()
{
    fluidDensity.setZero();
    fluidDensityOld.setZero();
    vx.setZero();
    vy.setZero();
    vxOld.setZero();
    vyOld.setZero();
}
