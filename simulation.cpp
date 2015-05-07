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
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, 6.0f);
    // construct the cube
        glBegin(GL_QUADS);

        glColor4f (  0.8,  0.8, 0.8, 0.1);
        glVertex3f(  0.5, -0.5, 0.5 );
        glVertex3f(  0.5,  0.5, 0.5 );
        glVertex3f( -0.5,  0.5, 0.5 );
        glVertex3f( -0.5, -0.5, 0.5 );

        glColor4f(  1.0,  0.0,  1.0, 0.1);
        glVertex3f( 0.5, -0.5, -0.5 );
        glVertex3f( 0.5,  0.5, -0.5 );
        glVertex3f( 0.5,  0.5,  0.5 );
        glVertex3f( 0.5, -0.5,  0.5 );

        glColor4f(   0.0,  1.0,  0.0, 0.1);
        glVertex3f( -0.5, -0.5,  0.5 );
        glVertex3f( -0.5,  0.5,  0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glColor4f(   0.0,  0.0,  1.0 , 0.1);
        glVertex3f(  0.5,  0.5,  0.5 );
        glVertex3f(  0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5,  0.5 );

        glColor4f(   1.0,  0.0,  0.0, 0.1);
        glVertex3f(  0.5, -0.5, -0.5 );
        glVertex3f(  0.5, -0.5,  0.5 );
        glVertex3f( -0.5, -0.5,  0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glColor4f(   1.0,  1.0, 0.0, 0.1);
        glVertex3f(  0.5, -0.5, -0.5 );
        glVertex3f(  0.5,  0.5, -0.5 );
        glVertex3f( -0.5,  0.5, -0.5 );
        glVertex3f( -0.5, -0.5, -0.5 );

        glEnd();

    glPopMatrix();



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
    addSource(fluid_->vx, fluid_->vxOld);
    addSource(fluid_->vy, fluid_->vyOld);
    // Velocity Code

    this->swap(fluid_->vx, fluid_->vxOld);
    this->swap(fluid_->vy, fluid_->vyOld);
//    cout<<"Velocity X : \n"<<fluid_->vx<<endl;
//    cout<<"Velocity Y : \n"<<fluid_->vy<<endl;

    advection(1, fluid_->vx, fluid_->vxOld, fluid_->vxOld, fluid_->vyOld );
    advection(2, fluid_->vy, fluid_->vyOld, fluid_->vxOld, fluid_->vyOld );

    project(fluid_->vx, fluid_->vy, fluid_->vxOld, fluid_->vyOld);

    this->swap(fluid_->vx, fluid_->vxOld);
    this->swap(fluid_->vy, fluid_->vyOld);

    diffuse(1, fluid_->vx, fluid_->vxOld, params_.viscosityFluid);
    diffuse(2, fluid_->vy, fluid_->vyOld, params_.viscosityFluid);

    project(fluid_->vx, fluid_->vy, fluid_->vxOld, fluid_->vyOld);

    //Density Code
    addSource(fluid_->fluidDensity, fluid_->fluidDensityOld);
    this->swap(fluid_->fluidDensity, fluid_->fluidDensityOld);


    diffuse(0, fluid_->fluidDensity, fluid_->fluidDensityOld, params_.diffusionConstant);
    this->swap(fluid_->fluidDensity, fluid_->fluidDensityOld);

    advection(0, fluid_->fluidDensity, fluid_->fluidDensityOld, fluid_->vx, fluid_->vy);
    fluid_->fluidDensityOld.setZero();
    fluid_->vxOld.setZero();
    fluid_->vyOld.setZero();

    cout << "fluid density: " << fluid_->getTotalDensity() << endl;
}

void Simulation::addSource(MatrixXd &d, MatrixXd &dOld)
{

    d += params_.timeStep * dOld;

}

void Simulation::diffuse(int boundary, MatrixXd &d, MatrixXd &dOld, double diffFactor)
{
    double a = params_.timeStep * diffFactor * fluid_->n * fluid_->n;

    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                d.coeffRef(i,j) = (dOld.coeff(i,j) + a * (d.coeff(i-1,j) + d.coeff(i+1,j) +
                                                          d.coeff(i,j-1) + d.coeff(i,j+1)))/(1+4*a);
            }
        }
        setBoundry(boundary, d);
    }
}

void Simulation::project(Eigen::MatrixXd &x, Eigen::MatrixXd &y, Eigen::MatrixXd &xOld, Eigen::MatrixXd &yOld)
{
    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            yOld.coeffRef(i,j) = (x.coeff(i+1, j) - x.coeff(i-1, j)
                                    + y.coeff(i,j+1) - y.coeff(i,j-1)) * -0.5f / fluid_->n;

        }
    }

    xOld.setZero();

    setBoundry(0,yOld);
    setBoundry(0,xOld);

    for(int k = 0; k < 20; k++)
    {
        for(int i  =1; i <= fluid_->n; i++)
        {
            for(int j = 1; j <= fluid_->n; j++)
            {
                xOld.coeffRef(i,j) =  ( xOld.coeff(i-1, j) + xOld.coeff(i+1, j)
                                         + xOld.coeff(i,j-1) + xOld.coeff(i,j+1)
                                         + yOld.coeff(i,j)) / 4 ;
            }
        }
        setBoundry(0, xOld);
    }

    for(int i = 1; i <= fluid_->n; i++)
    {
        for(int j = 1; j <= fluid_->n; j++)
        {
            x.coeffRef(i,j) -= 0.5f * fluid_->n * (xOld.coeff(i+1, j) - xOld.coeff(i-1, j));
            y.coeffRef(i,j) -= 0.5f * fluid_->n * (xOld.coeff(i, j+1) - xOld.coeff(i, j-1));
        }
    }
    setBoundry(1,x);
    setBoundry(2,y);
}

void Simulation::setBoundry(int b, MatrixXd& m)
{
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

}

void Simulation::swap(MatrixXd &left, MatrixXd &right)
{
    MatrixXd tmp(left);
    left = right;
    right = tmp;
}

void Simulation::linearSolver(int b, Eigen::MatrixXd &x, Eigen::MatrixXd &xOld, double a, double c)
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

void Simulation::advection(int boundry, MatrixXd &d, MatrixXd &dOld, MatrixXd &xOld, MatrixXd &yOld)
{
    double dt0 = params_.timeStep * fluid_->n;
    int n = fluid_->n;

    for(int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            double x = i - dt0 * xOld.coeff(i,j);
            double y = j - dt0 * yOld.coeff(i,j);

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

            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;

            d.coeffRef(i,j) = s0 * (t0 * dOld.coeff(i0, j0) + t1 * dOld.coeff(i0,j1))
                    + s1 * (t0 * dOld.coeff(i1,j0) + t1 * dOld.coeff(i1, j1));
        }
    }
    setBoundry(boundry, d);
}

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
