#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include "mesh.h"
#include "quadprog/eiquadprog.hpp"
#include <iostream>
const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0), wallTex_(0)
{

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


void Simulation::clearScene()
{
    renderLock_.lock();
    {

        Plane groundPlane;
        groundPlane.pos << 0,0,0;
        groundPlane.normal << 0,0,1;
        planes_.push_back(groundPlane);


    }
    renderLock_.unlock();
}
