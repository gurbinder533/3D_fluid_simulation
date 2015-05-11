#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "simparameters.h"
#include <QThread>
#include <QTimer>
#include <iostream>

class MainWindow;
class Simulation;

class Controller : public QThread
{
    Q_OBJECT

public:
    Controller(int fps);
    virtual ~Controller();
    void initialize(MainWindow *mw);
    void initializeGL();
    void renderPlanes(bool transparent);
    void renderObjects();

public slots:
    void reset();
    void clearScene();
    void updateParameters(SimParameters params);
    void mouseClicked(double x, double y, double z, double dx, double dy, double dz);
    void leftMouseClicked(double x, double y);
    void keyToAddFluid(int);
    void keyToAddVel(int, double, double, double);
    void resetDrag();
    void mouseDrag(double x, double y);
    void simTick();

    void renderFluid();
protected:
    virtual void run();

private:
    MainWindow *mw_;
    Simulation *sim_;
    SimParameters params_;

    bool dragOn;
    double dragXOld;
    double dragYOld;
    int fps_;
    QTimer simtimer_;
};

#endif // CONTROLLER_H
