#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include <QPainter>
#include "SOIL.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include "mesh.h"
#include "quadprog/eiquadprog.hpp"
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

void Simulation::initializeGL()
{
   // std::cout << "intializeGL\n";
   // loadFloorTexture();
   // loadWallTexture();
}


void Simulation::loadFloorTexture()
{
    floorTex_ = SOIL_load_OGL_texture("resources/grid.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(floorTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}

void Simulation::loadWallTexture()
{
    wallTex_ = SOIL_load_OGL_texture("resources/wall.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(wallTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}


void Simulation::renderPlanes(bool transparent)
{
    renderLock_.lock();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glDisable(GL_BLEND);
    glPushMatrix();
    // construct the cub

     glLineWidth(3.0);
        glBegin(GL_LINES);

        glColor3f(1.0, 0.0, 0.0);

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


    /*
    glEnable(GL_CULL_FACE);
    if(transparent)
    {
        glCullFace(GL_FRONT);
        glColor4f(1.0, 1.0, 1.0, 0.5);
    }
    else
    {
        glCullFace(GL_BACK);
        glColor4f(1.0, 1.0, 1.0, 1.0);
    }

    for(vector<Plane>::iterator it = planes_.begin(); it != planes_.end(); ++it)
        renderPlane(*it, it == planes_.begin());

    glDisable(GL_CULL_FACE);
*/
    renderLock_.unlock();
}

void Simulation::renderPlane(const Plane &p, bool isFloor)
{
    /*
    if(isFloor && floorTex_)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else if(!isFloor && wallTex_)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else
        glColor3f(0.5, 0.5, 0.5);

    double texsize = 5.0;
    double gridsize = 1000.0;

    double texmax = gridsize/texsize;

    Vector3d tangent1 = VectorMath::perpToAxis(p.normal);
    Vector3d tangent2 = tangent1.cross(p.normal);

    Vector3d corner;
    */


    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, 6.0f);
    // construct the cube
        glBegin(GL_QUADS);

        glColor4f (  0.8,  0.8, 0.8, 0.2);
        glVertex3f(  0.5, -0.5, 0.5 );
        glVertex3f(  0.5,  0.5, 0.5 );
        glVertex3f( -0.5,  0.5, 0.5 );
        glVertex3f( -0.5, -0.5, 0.5 );

        glColor4f(  1.0,  0.0,  1.0, 0.2);
        glVertex3f( 0.5, -0.5, -0.5 );
        glVertex3f( 0.5,  0.5, -0.5 );
        glVertex3f( 0.5,  0.5,  0.5 );
        glVertex3f( 0.5, -0.5,  0.5 );

        glColor4f(   0.0,  1.0,  0.0, 0.2);
        glVertex3f( -0.5, -0.5,  0.5 );
        glVertex3f( -0.5,  0.5,  0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glColor4f(   0.0,  0.0,  1.0 , 0.2);
        glVertex3f(  0.5,  0.5,  0.5 );
        glVertex3f(  0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5,  0.5 );

        glColor4f(   1.0,  0.0,  0.0, 0.2);
        glVertex3f(  0.5, -0.5, -0.5 );
        glVertex3f(  0.5, -0.5,  0.5 );
        glVertex3f( -0.5, -0.5,  0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glColor4f(   1.0,  1.0, 0.0, 0.2);
        glVertex3f(  0.5, -0.5, -0.5 );
        glVertex3f(  0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glEnd();

    glPopMatrix();

  /*
   glBegin(GL_QUADS);
    {


        glTexCoord2f(texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);
    }
    glDisable(GL_TEXTURE_2D);
    glEnd();
    */
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
    addSource(fluid_->vx3d, fluid_->vx3dOld);
    addSource(fluid_->vy3d, fluid_->vy3dOld);
    addSource(fluid_->vz3d, fluid_->vz3dOld);
    // Velocity Code

    this->swap(fluid_->vx3d, fluid_->vx3dOld);
    this->swap(fluid_->vy3d, fluid_->vy3dOld);
    this->swap(fluid_->vy3d, fluid_->vz3dOld);

//    cout<<"Velocity X : \n"<<fluid_->vx<<endl;
//    cout<<"Velocity Y : \n"<<fluid_->vy<<endl;

    //TODO: we need to CHANGE THIS
    advection(1, fluid_->vx3d, fluid_->vx3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);
    advection(2, fluid_->vy3d, fluid_->vy3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);
    advection(2, fluid_->vz3d, fluid_->vz3dOld, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);

    project(fluid_->vx3d, fluid_->vy3d, fluid_->vz3d, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);

    this->swap(fluid_->vx3d, fluid_->vx3dOld);
    this->swap(fluid_->vy3d, fluid_->vy3dOld);
    this->swap(fluid_->vz3d, fluid_->vz3dOld);

    diffuse(1, fluid_->vx3d, fluid_->vx3dOld, params_.viscosityFluid);
    diffuse(2, fluid_->vy3d, fluid_->vy3dOld, params_.viscosityFluid);
    diffuse(2, fluid_->vz3d, fluid_->vz3dOld, params_.viscosityFluid);

    project(fluid_->vx3d, fluid_->vy3d, fluid_->vz3d, fluid_->vx3dOld, fluid_->vy3dOld, fluid_->vz3dOld);

    //Density Code
    addSource(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);
    this->swap(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);


    diffuse(0, fluid_->fluidDensity3d, fluid_->fluidDensity3dOld, params_.diffusionConstant);
    this->swap(fluid_->fluidDensity3d, fluid_->fluidDensity3dOld);

    advection(0, fluid_->fluidDensity3d, fluid_->fluidDensity3dOld, fluid_->vx3d, fluid_->vy3d, fluid_->vz3d);

 /*
    std::fill(fluid_->fluidDensity3dOld.begin(), fluid_->fluidDensity3dOld.end(), 0);
    std::fill(fluid_->vx3dOld.begin(), fluid_->vx3dOld.end(), 0);
    std::fill(fluid_->vy3dOld.begin(), fluid_->vy3dOld.end(), 0);
    std::fill(fluid_->vz3dOld.begin(), fluid_->vz3dOld.end(), 0);
 */



    fluid_->fluidDensity3dOld.setZero();
    fluid_->vx3dOld.setZero();
    fluid_->vy3dOld.setZero();
    fluid_->vz3dOld.setZero();


    //cout << "fluid density: " << fluid_->getTotalDensity() << endl;
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
                    /*
                    d.coeffRef(i,j) = (dOld.coeff(i,j) + a * (d.coeff(i-1,j) + d.coeff(i+1,j) +
                                                          d.coeff(i,j-1) + d.coeff(i,j+1)))/(1+4*a);
                     */
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

                /*
                yOld.coeffRef(i,j) = (x.coeff(i+1, j) - x.coeff(i-1, j)
                                    + y.coeff(i,j+1) - y.coeff(i,j-1)) * -0.5f / fluid_->n;

               */


            }

        }
    }

    xOld.setZero();
   // zOld.setZero();

    setBoundry(0,yOld);
    setBoundry(0,xOld);

    //setBoundry(0,zOld);


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

 /*
    for(int i = 1; i <= fluid_->n; i++)
    {
        m.coeffRef(0,i) =               (b==1) ? -m.coeff(1,i) : m.coeff(1,i);
        m.coeffRef(fluid_->n, i) =    (b==1) ? -m.coeff(fluid_->n-1,i) : m.coeff(fluid_->n-1,i);
        m.coeffRef(i,0) =               (b==2) ? -m.coeff(i,1) : m.coeff(i,1);
        m.coeffRef(i, fluid_->n) =    (b==2) ? -m.coeff(i, fluid_->n-1) : m.coeff(i,fluid_->n-1);
    }
    m.coeffRef(0,0) =                   0.5f *(m.coeff(1,0) + m.coeff(0,1));
    m.coeffRef(0,fluid_->n) =           0.5f *(m.coeff(1,fluid_->n) + m.coeff(0, fluid_->n-1));
    m.coeffRef(fluid_->n,0) =           0.5f *(m.coeff(fluid_->n-1,0) + m.coeff(fluid_->n,1));
    m.coeffRef(fluid_->n,fluid_->n) =   0.5f *(m.coeff(fluid_->n-1,fluid_->n) + m.coeff(fluid_->n,fluid_->n-1));
*/
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
                /*d.coeffRef(i,j) = s0 * (t0 * dOld.coeff(i0, j0) + t1 * dOld.coeff(i0,j1))
                        + s1 * (t0 * dOld.coeff(i1,j0) + t1 * dOld.coeff(i1, j1));
                        */
           }
        }
    }
    setBoundry(boundry, d);
}

/*
void Simulation::addDensity(double x, double y)
{
    int i = floor((x+1)/fluid_->sizeOfVoxel);
    int j = -1* floor((y-1)/fluid_->sizeOfVoxel);

    if(i >= 0 && i <fluid_->n && j >= 0 && j < fluid_->n)
    {
        fluid_->fluidDensityOld.coeffRef(i,j) += params_.densityMagnitude;
    }
    for(int r = 1; r <= params_.densityRadius; r++)
    {
        if(i-r >= 0 && i+r <fluid_->n && j-r >= 0 && j+r < fluid_->n)
        {
            fluid_->fluidDensityOld.coeffRef(i-r,j-r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i-r,j) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i-r,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j+r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i+r,j-r) += params_.densityMagnitude;
            fluid_->fluidDensityOld.coeffRef(i,j-r) += params_.densityMagnitude;
        }
    }
}

void Simulation::addVelocity(double x, double y, double velX, double velY)
{
    int i = floor((x+1)/fluid_->sizeOfVoxel);
    int j = -1* floor((y-1)/fluid_->sizeOfVoxel);
    if(i >= 0 && i <fluid_->n && j >= 0 && j < fluid_->n)
    {
        fluid_->vxOld.coeffRef(i,j) += velX;
        fluid_->vyOld.coeffRef(i,j) += velY;
    }
    for(int r = 1; r <= params_.velocityRadius; r++)
    {
        if(i-r >= 0 && i+r <fluid_->n && j-r >= 0 && j+r < fluid_->n)
        {
//            cout<<"Here : "<<velX<<" : "<<velY<<endl;
            fluid_->vxOld.coeffRef(i-r,j-r) += velX;
            fluid_->vxOld.coeffRef(i-r,j) += velX;
            fluid_->vxOld.coeffRef(i-r,j+r) += velX;
            fluid_->vxOld.coeffRef(i,j+r) += velX;
            fluid_->vxOld.coeffRef(i+r,j+r) += velX;
            fluid_->vxOld.coeffRef(i+r,j) += velX;
            fluid_->vxOld.coeffRef(i+r,j-r) += velX;
            fluid_->vxOld.coeffRef(i,j-r) += velX;

            fluid_->vyOld.coeffRef(i,j) += velY;
            fluid_->vyOld.coeffRef(i-r,j-r) += velY;
            fluid_->vyOld.coeffRef(i-r,j) += velY;
            fluid_->vyOld.coeffRef(i-r,j+r) += velY;
            fluid_->vyOld.coeffRef(i,j+r) += velY;
            fluid_->vyOld.coeffRef(i+r,j+r) += velY;
            fluid_->vyOld.coeffRef(i+r,j) += velY;
            fluid_->vyOld.coeffRef(i+r,j-r) += velY;
            fluid_->vyOld.coeffRef(i,j-r) += velY;
        }
    }
}
*/
void Simulation::clearScene()
{
    renderLock_.lock();
    {

        Plane groundPlane;
        groundPlane.pos << 0,0,0;
        groundPlane.normal << 0,0,1;
        planes_.push_back(groundPlane);


        fluid_->zeroEverything();
    }
    renderLock_.unlock();
}
