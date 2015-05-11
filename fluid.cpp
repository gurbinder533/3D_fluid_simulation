#include "fluid.h"
#include <QGLWidget>
#include <iostream>

using namespace Eigen;
using namespace std;

Fluid::Fluid()
{

    this->type = 1; // blue
    this->n = 30;
    this->size = (n+2)*(n+2)*(n+2);
    this->fluidDensity3d.resize(size);



    this->fluidDensity3dOld.resize(size);

    this->vx3d.resize(size);
    this->vx3dOld.resize(size);

    this->vy3d.resize(size);
    this->vy3dOld.resize(size);

    this->vz3d.resize(size);
    this->vz3dOld.resize(size);


    fluidDensity3d.setZero();
    fluidDensity3dOld.setZero();
    vx3d.setZero();
    vy3d.setZero();
    vz3d.setZero();
    vx3dOld.setZero();
    vy3dOld.setZero();
    vz3dOld.setZero();


    this->sizeOfVoxel = 2.0/(this->n+1);
    this->debug = true;

//     for(int i = 0; i < n; ++i)
//     {
//         for(int j = 0; j < n; ++j)
//         {
//             for(int k = 0; k < n; ++k)
//             {
//                //std::cout << COFF(i,j,0) << "\n";
//                this->fluidDensity3d[COFF(i,j,k)] = 255.0;
//             }
//         }

//     }

}

void Fluid::render()
{

    //std::cout << "FLUID RENDER\n";
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    double xCell = -1;
    double yCell = 1;
    double zCell = 1;
//    cout << "Render" << endl;
    for(int i =0; i <= n; i++)
    {
        xCell = i * this->sizeOfVoxel - 1 ;
        for(int j = 0; j <= n; j++)
        {
            yCell = 1 - j * this->sizeOfVoxel;
            for(int k = 0; k <= n; k++)
            {
                zCell = 1 - k * this->sizeOfVoxel;
                float dens = fluidDensity3d[COFF(i, j, k)];
//                cout << dens << endl;



                if (dens > 0)
                {
                    //glColor4f(255.0-dens, 255.0-dens, 255.0, 0.2);
                    if(type == 0) // blue
                    {
                        glColor4f(1.0-dens, 1.0-dens, 255.0, dens*20*n);
//                      glColor4f(0, 0, 255.0, 1);
                    }
                    else if(type == 1)
                    {
                       glColor4f(255.0, 1.0-dens, 1.0 - dens, dens*20*n);
                    }
                }
                else
                {
                    glColor4f(255.0, 255.0, 255.0, 0.0);
                }

                //glColor3f(1,0,0);
                //std::cout << " xcell " << xCell << " , ycell " << yCell << " , zCell " << zCell << "\n";
                //std::cout << " Voxel size : " << this->sizeOfVoxel << "\n";

                glBegin(GL_QUADS);
                {

                    glVertex3f(xCell,yCell, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell - this->sizeOfVoxel, zCell);
                    glVertex3f(xCell, yCell - this->sizeOfVoxel, zCell);

                    glVertex3f(xCell,yCell, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell, yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);

                    glVertex3f(xCell,yCell, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell , zCell - this->sizeOfVoxel);
                    glVertex3f(xCell, yCell, zCell - this->sizeOfVoxel);

                    glVertex3f(xCell,yCell - this->sizeOfVoxel, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell- this->sizeOfVoxel, zCell);
                    glVertex3f(xCell + this->sizeOfVoxel, yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell, yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);

                    glVertex3f(xCell,yCell , zCell);
                    glVertex3f(xCell , yCell- this->sizeOfVoxel, zCell);
                    glVertex3f(xCell , yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell, yCell, zCell - this->sizeOfVoxel);


                    glVertex3f(xCell +this->sizeOfVoxel,yCell , zCell);
                    glVertex3f(xCell  +this->sizeOfVoxel , yCell- this->sizeOfVoxel, zCell);
                    glVertex3f(xCell +this->sizeOfVoxel , yCell - this->sizeOfVoxel, zCell - this->sizeOfVoxel);
                    glVertex3f(xCell +this->sizeOfVoxel, yCell, zCell - this->sizeOfVoxel);
                }
                glEnd();
            }
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
            for(int k = 0; k < n; k++)
            {
               densTotal += fluidDensity3d[COFF(i,j,k)];
            }
        }
    }
    return densTotal;
}

void Fluid::zeroEverything()
{

/*
    std::fill(fluidDensity3d.begin(), fluidDensity3d.end(), 0);
    std::fill(fluidDensity3d.begin(), fluidDensity3d.end(), 0);

    std::fill(vx3d.begin(), vx3d.end(), 0);
    std::fill(vx3dOld.begin(), vx3dOld.end(), 0);

    std::fill(vy3d.begin(), vy3d.end(), 0);
    std::fill(vy3dOld.begin(), vy3dOld.end(), 0);

    std::fill(vz3d.begin(), vz3d.end(), 0);
    std::fill(vz3dOld.begin(), vz3dOld.end(), 0);
*/

    fluidDensity3d.setZero();
    fluidDensity3dOld.setZero();
    vx3d.setZero();
    vy3d.setZero();
    vz3d.setZero();
    vx3dOld.setZero();
    vy3dOld.setZero();
    vz3dOld.setZero();

}
