#include "controller.h"
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <Eigen/Core>

using namespace Eigen;

Controller::Controller(int fps) : QThread(), mw_(NULL), fps_(fps)
{
}

Controller::~Controller()
{
    delete sim_;
}

void Controller::initialize(MainWindow *mw)
{
    mw_ = mw;
    sim_ = new Simulation(params_);
    dragOn = false;
    dragXOld = 0.0;
    dragYOld = 0.0;
}


void Controller::run()
{
    reset();
    connect(&simtimer_, SIGNAL(timeout()), this, SLOT(simTick()));
    simtimer_.start(1000/fps_);
    exec();
}

void Controller::reset()
{
    params_ = SimParameters();
    QMetaObject::invokeMethod(mw_, "setUIFromParameters", Q_ARG(SimParameters, params_));
    clearScene();
}

void Controller::clearScene()
{
    sim_->clearScene();
}

void Controller::updateParameters(SimParameters params)
{
    params_ = params;
}

void Controller::renderPlanes(bool transparent)
{
    sim_->renderPlanes(transparent);
}


void Controller::renderFluid()
{
    sim_->render();
}


void Controller::mouseClicked(double x, double y, double z, double dx, double dy, double dz)
{
    Vector3d pos(x,y,z);
    Vector3d dir(dx,dy,dz);
    std::cout << " clicked \n";
    //sim_->addRigidBody(pos, dir);
    //sim_->addFluid; // to be added.
}

void Controller::simTick()
{
    if(params_.simRunning)
       sim_->takeSimulationStep();
}

void Controller::mouseDrag(double x, double y)
{
    if (!this->dragOn)
    {
        this->dragXOld = x;
        this->dragYOld = y;
        return;
    }
    double velToAdd = params_.velocityMagnitude;
    double velX,velY;

    velX = (x - this->dragXOld) * velToAdd;
    velY = -(y - this->dragYOld) * velToAdd;

    this->dragXOld = x;
    this->dragYOld = y;

    if(params_.simRunning)
    {
        switch(params_.clickMode)
        {

            case SimParameters::CM_ADDVELOCITY:
                break;
            case SimParameters::CM_ADDDENSITY:
                break;
        }
    }
}

void Controller::resetDrag()
{
    this->dragXOld = 0;
    this->dragYOld = 0;
    this->dragOn = false;
//    cout<<"Mouse Released :"<<this->dragOn<<endl;
}

void Controller::leftMouseClicked(double x, double y)
{
    if (!this->dragOn)
    {
        this->dragOn = true;
    }
//    cout<<"Left Mouse Clicked :"<<this->dragOn<<endl;
}


void Controller::keyToAddVel(int i, double velX, double velY, double velZ)
{
    //std::cout << "V pressed" << std::endl;
    if(params_.velsource1)
        sim_->addVelocity(1);
    if(params_.velsource2)
        sim_->addVelocity(2);
    if(params_.velsource3)
        sim_->addVelocity(3);
}

void Controller::keyToAddFluid(int i)
{
    //std::cout << " Q pressed  " << i << " \n";
    //std::cout << std::endl;
    if(params_.source2)
        sim_->addDensity(1);
    if(params_.source1)
        sim_->addDensity(2);
    if(params_.source3)
        sim_->addDensity(3);
    if(params_.source4)
        sim_->addDensity(4);

}
