#-------------------------------------------------
#
# Project created by QtCreator 2014-09-03T16:42:53
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = fluids3d
TEMPLATE = app

INCLUDEPATH += $$_PRO_FILE_PWD_/eigen_3.2.2/
LIBS += -lGLU

SOURCES += main.cpp\
        mainwindow.cpp \
    glpanel.cpp \
    controller.cpp \
    simulation.cpp \
    simparameters.cpp \
    vectormath.cpp \
    camera.cpp \
    rotator.cpp \
    fluid.cpp

HEADERS  += mainwindow.h \
    glpanel.h \
    controller.h \
    simulation.h \
    simparameters.h \
    vectormath.h \
    camera.h \
    rotator.h \
    fluid.h

FORMS    += mainwindow.ui

libSOIL.commands = make
