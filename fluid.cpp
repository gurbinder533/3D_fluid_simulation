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

    double xCell = -1;
    double yCell = 1;
    for(int i =0; i <= n; i++)
    {
        xCell = i * this->sizeOfVoxel - 1 ;
        for(int j = 0; j <= n; j++)
        {
            yCell = 1 - j * this->sizeOfVoxel;
            double dens = fluidDensity.coeff(i,j);
            glColor3f(1-dens,1-dens,1);
            glBegin(GL_QUADS);
            {
                glVertex2f(xCell,yCell);
                glVertex2f(xCell + this->sizeOfVoxel, yCell);
                glVertex2f(xCell + this->sizeOfVoxel, yCell - this->sizeOfVoxel);
                glVertex2f(xCell, yCell - this->sizeOfVoxel);
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
