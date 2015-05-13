/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include "glpanel.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionExit;
    QAction *actionReset;
    QAction *actionReset_Everything;
    QWidget *centralWidget;
    GLPanel *GLWidget;
    QFrame *parameterFrame;
    QGroupBox *SimulationOptions;
    QPushButton *Simulate;
    QLabel *timeStep;
    QLineEdit *timeStepEdit;
    QPushButton *clearScene;
    QLabel *label;
    QLineEdit *diffusionK;
    QLineEdit *viscosityK;
    QLabel *label_2;
    QLabel *De;
    QLabel *label_4;
    QLineEdit *densityMag;
    QLineEdit *velocityMag;
    QGroupBox *SimulationVar;
    QCheckBox *source1;
    QCheckBox *source2;
    QCheckBox *source3;
    QCheckBox *source4;
    QCheckBox *velSrc1;
    QCheckBox *velSrc2;
    QCheckBox *velSrc3;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1200, 800);
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionReset = new QAction(MainWindow);
        actionReset->setObjectName(QString::fromUtf8("actionReset"));
        actionReset_Everything = new QAction(MainWindow);
        actionReset_Everything->setObjectName(QString::fromUtf8("actionReset_Everything"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        GLWidget = new GLPanel(centralWidget);
        GLWidget->setObjectName(QString::fromUtf8("GLWidget"));
        GLWidget->setGeometry(QRect(10, 0, 731, 731));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(GLWidget->sizePolicy().hasHeightForWidth());
        GLWidget->setSizePolicy(sizePolicy);
        GLWidget->setFocusPolicy(Qt::StrongFocus);
        parameterFrame = new QFrame(centralWidget);
        parameterFrame->setObjectName(QString::fromUtf8("parameterFrame"));
        parameterFrame->setGeometry(QRect(749, -1, 441, 731));
        parameterFrame->setFrameShape(QFrame::StyledPanel);
        parameterFrame->setFrameShadow(QFrame::Raised);
        SimulationOptions = new QGroupBox(parameterFrame);
        SimulationOptions->setObjectName(QString::fromUtf8("SimulationOptions"));
        SimulationOptions->setGeometry(QRect(0, 10, 441, 301));
        Simulate = new QPushButton(SimulationOptions);
        Simulate->setObjectName(QString::fromUtf8("Simulate"));
        Simulate->setGeometry(QRect(10, 40, 181, 31));
        timeStep = new QLabel(SimulationOptions);
        timeStep->setObjectName(QString::fromUtf8("timeStep"));
        timeStep->setGeometry(QRect(20, 100, 66, 21));
        timeStepEdit = new QLineEdit(SimulationOptions);
        timeStepEdit->setObjectName(QString::fromUtf8("timeStepEdit"));
        timeStepEdit->setGeometry(QRect(99, 96, 113, 31));
        clearScene = new QPushButton(SimulationOptions);
        clearScene->setObjectName(QString::fromUtf8("clearScene"));
        clearScene->setGeometry(QRect(310, 40, 94, 29));
        label = new QLabel(SimulationOptions);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 150, 101, 21));
        diffusionK = new QLineEdit(SimulationOptions);
        diffusionK->setObjectName(QString::fromUtf8("diffusionK"));
        diffusionK->setGeometry(QRect(110, 143, 113, 31));
        viscosityK = new QLineEdit(SimulationOptions);
        viscosityK->setObjectName(QString::fromUtf8("viscosityK"));
        viscosityK->setGeometry(QRect(110, 187, 113, 31));
        label_2 = new QLabel(SimulationOptions);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(13, 190, 91, 21));
        De = new QLabel(SimulationOptions);
        De->setObjectName(QString::fromUtf8("De"));
        De->setGeometry(QRect(20, 230, 131, 21));
        label_4 = new QLabel(SimulationOptions);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(20, 265, 131, 21));
        densityMag = new QLineEdit(SimulationOptions);
        densityMag->setObjectName(QString::fromUtf8("densityMag"));
        densityMag->setGeometry(QRect(152, 226, 113, 31));
        velocityMag = new QLineEdit(SimulationOptions);
        velocityMag->setObjectName(QString::fromUtf8("velocityMag"));
        velocityMag->setGeometry(QRect(152, 260, 113, 31));
        SimulationVar = new QGroupBox(parameterFrame);
        SimulationVar->setObjectName(QString::fromUtf8("SimulationVar"));
        SimulationVar->setGeometry(QRect(10, 360, 421, 361));
        SimulationVar->setCheckable(false);
        source1 = new QCheckBox(SimulationVar);
        source1->setObjectName(QString::fromUtf8("source1"));
        source1->setGeometry(QRect(0, 50, 95, 26));
        source1->setAcceptDrops(false);
        source1->setChecked(true);
        source2 = new QCheckBox(SimulationVar);
        source2->setObjectName(QString::fromUtf8("source2"));
        source2->setGeometry(QRect(99, 50, 95, 26));
        source3 = new QCheckBox(SimulationVar);
        source3->setObjectName(QString::fromUtf8("source3"));
        source3->setGeometry(QRect(200, 50, 95, 26));
        source4 = new QCheckBox(SimulationVar);
        source4->setObjectName(QString::fromUtf8("source4"));
        source4->setGeometry(QRect(305, 52, 95, 26));
        velSrc1 = new QCheckBox(SimulationVar);
        velSrc1->setObjectName(QString::fromUtf8("velSrc1"));
        velSrc1->setGeometry(QRect(1, 114, 111, 26));
        velSrc1->setChecked(false);
        velSrc2 = new QCheckBox(SimulationVar);
        velSrc2->setObjectName(QString::fromUtf8("velSrc2"));
        velSrc2->setGeometry(QRect(121, 116, 111, 26));
        velSrc3 = new QCheckBox(SimulationVar);
        velSrc3->setObjectName(QString::fromUtf8("velSrc3"));
        velSrc3->setGeometry(QRect(257, 117, 111, 26));
        MainWindow->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Fluid3d", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionReset->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
        actionReset_Everything->setText(QApplication::translate("MainWindow", "Reset Everything", 0, QApplication::UnicodeUTF8));
        SimulationOptions->setTitle(QApplication::translate("MainWindow", "Simulation Options", 0, QApplication::UnicodeUTF8));
        Simulate->setText(QApplication::translate("MainWindow", "Start Simulation", 0, QApplication::UnicodeUTF8));
        timeStep->setText(QApplication::translate("MainWindow", "Time Step", 0, QApplication::UnicodeUTF8));
        timeStepEdit->setText(QApplication::translate("MainWindow", "0.001", 0, QApplication::UnicodeUTF8));
        clearScene->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Diffusion Coff", 0, QApplication::UnicodeUTF8));
        diffusionK->setText(QApplication::translate("MainWindow", "0.02", 0, QApplication::UnicodeUTF8));
        viscosityK->setText(QApplication::translate("MainWindow", "0.2", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Viscosity Coff", 0, QApplication::UnicodeUTF8));
        De->setText(QApplication::translate("MainWindow", "Density Magnitude", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Velocity Magnitude", 0, QApplication::UnicodeUTF8));
        densityMag->setText(QApplication::translate("MainWindow", "200", 0, QApplication::UnicodeUTF8));
        velocityMag->setText(QApplication::translate("MainWindow", "20", 0, QApplication::UnicodeUTF8));
        SimulationVar->setTitle(QApplication::translate("MainWindow", "Simulation Variables", 0, QApplication::UnicodeUTF8));
        source1->setText(QApplication::translate("MainWindow", "Source1", 0, QApplication::UnicodeUTF8));
        source2->setText(QApplication::translate("MainWindow", "Source2", 0, QApplication::UnicodeUTF8));
        source3->setText(QApplication::translate("MainWindow", "Source3", 0, QApplication::UnicodeUTF8));
        source4->setText(QApplication::translate("MainWindow", "Source4", 0, QApplication::UnicodeUTF8));
        velSrc1->setText(QApplication::translate("MainWindow", "VelocitySrc1", 0, QApplication::UnicodeUTF8));
        velSrc2->setText(QApplication::translate("MainWindow", "VelocitySrc2", 0, QApplication::UnicodeUTF8));
        velSrc3->setText(QApplication::translate("MainWindow", "VelocitySrc3", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
