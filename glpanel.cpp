#include "glpanel.h"
#include "controller.h"
#include <QMouseEvent>

using namespace Eigen;

GLPanel::GLPanel(QWidget *parent) :
    QGLWidget(parent), c_(), translateDir_(0), rotator_(c_)
{
    cont_ = NULL;
    c_.setDefault3D();
    lightPos_.setZero();
    lightPos_[2] = 1.0;
}

void GLPanel::setController(Controller *cont)
{
    cont_ = cont;
}

void GLPanel::initializeGL()
{
    assert(cont_);
    glShadeModel(GL_SMOOTH);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClearDepth(1.0);
    glDisable(GL_DEPTH_TEST);
    glCullFace(GL_BACK);
    glCullFace(GL_FRONT);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

void GLPanel::resizeGL(int w, int h)
{
    c_.setPerpective(60.0, 1.0);
    c_.setViewport(w, h);
}

void GLPanel::paintGL()
{
    assert(cont_);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3f (0.0, 0.0, 0.0);

    c_.applyViewport();
    c_.applyProjection();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    c_.applyLookAt();

    GLfloat lightColor0[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat lightPosition[4] = { lightPos_[0], lightPos_[1], lightPos_[2], lightPos_[3] };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    glDisable(GL_LIGHTING);
    cont_->renderPlanes(false);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);


    glPushMatrix();
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(-1.0, -1.0);
    glDisable(GL_LIGHT0);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPopMatrix();

    cont_->renderFluid();
}

void GLPanel::scaleMousePos(int x, int y, double &scaledx, double &scaledy) const
{
    int w, h;
    c_.getViewport(w,h);
    scaledx = 2.0 * x / double(w-1) - 1.0;
    scaledy = 2.0 * (h - y -1) / double(h-1) - 1.0;
}

GLPanel::MouseAction GLPanel::deduceAction(QMouseEvent *event)
{
    if(event->buttons() & Qt::LeftButton)
        return MA_LAUNCH;
    else if(event->buttons() & Qt::RightButton)
        return MA_ROTATE;
    return MA_NONE;
}

void GLPanel::mousePressEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);

    MouseAction ma = deduceAction(event);
    switch(ma)
    {
    case MA_ROTATE:
    {
        rotator_.startRotation(pos);
        break;
    }
    default:
        break;
    }
}

void GLPanel::mouseMoveEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);
    rotator_.updateRotation(pos);
}

void GLPanel::mouseReleaseEvent(QMouseEvent *event)
{
    MouseAction ma = deduceAction(event);
    if(ma != MA_ROTATE)
    {
        int x = event->pos().x();
        int y = event->pos().y();
        Vector2d pos;
        scaleMousePos(x,y,pos[0],pos[1]);
        rotator_.stopRotation();
    }
}

void GLPanel::keyPressEvent(QKeyEvent *ke)
{
    switch(ke->key())
    {
    case 'W':
        translateDir_ |= TD_FWD;
        break;
    case 'A':
        translateDir_ |= TD_LEFT;
        break;
    case 'S':
        translateDir_ |= TD_BACK;
        break;
    case 'D':
        translateDir_ |= TD_RIGHT;
        break;
    case 'Q':
        addFluid_   |= WALL1;
        QMetaObject::invokeMethod(cont_, "keyToAddFluid", Q_ARG(int, 1));
        break;

    case 'V':
        addFluid_   |= WALL2;
        QMetaObject::invokeMethod(cont_, "keyToAddVel", Q_ARG(int, 1), Q_ARG(double, 10000), Q_ARG(double, 200), Q_ARG(double, 200));
        break;
    }
}

void GLPanel::tick()
{

    Vector3d right, up, center;
    c_.getSpanningSet(right, up, center);
    center[2] = 0;
    center.normalize();
    right[2] = 0;
    right.normalize();
    int fwdamt = 0;
    if(translateDir_ & TD_FWD)
        fwdamt++;
    if(translateDir_ & TD_BACK)
        fwdamt--;
    Vector3d dir = fwdamt*center;
    int rghtamt = 0;
    if(translateDir_ & TD_LEFT)
        rghtamt--;
    if(translateDir_ & TD_RIGHT)
        rghtamt++;
    dir += rghtamt*right;
    double dnorm = dir.norm();
    if(dnorm != 0)
        dir /= dnorm;
    c_.translateCenter(dir);
    c_.translateEye(dir);
}

void GLPanel::keyReleaseEvent(QKeyEvent *ke)
{
    switch(ke->key())
    {
    case 'W':
        translateDir_ &= ~TD_FWD;
        break;
    case 'A':
        translateDir_ &= ~TD_LEFT;
        break;
    case 'S':
        translateDir_ &= ~TD_BACK;
        break;
    case 'D':
        translateDir_ &= ~TD_RIGHT;
        break;

    case 'Q':
        addFluid_   &= ~WALL1;
        break;

    case 'V':
        addVelocity_   &= ~WALL2;
        break;
    }
}
