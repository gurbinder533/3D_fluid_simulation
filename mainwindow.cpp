#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simparameters.h"
#include "controller.h"

MainWindow::MainWindow(Controller &cont, int fps, QWidget *parent) :
    QMainWindow(parent),
    cont_(cont),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->GLWidget->setController(&cont);
    simRunning_ = false;
    connect(&renderTimer_, SIGNAL(timeout()), this, SLOT(updateGL()));
    renderTimer_.start(1000/fps);
}

MainWindow::~MainWindow()
{
    delete ui;
}
/*
void MainWindow::on_actionExit_triggered()
{
    close();
}
*/
void MainWindow::setParametersFromUI()
{

    SimParameters params;

    params.simRunning = simRunning_;
    params.timeStep = ui->timeStepEdit->text().toDouble();

    params.source1 = ui->source1->isChecked();
    params.source2 = ui->source2->isChecked();
    params.source3 = ui->source3->isChecked();
    params.source4 = ui->source4->isChecked();

    params.velsource1 = ui->velSrc1->isChecked();
    params.velsource2 = ui->velSrc2->isChecked();
    params.velsource3 = ui->velSrc3->isChecked();

    params.diffusionConstant = ui->diffusionK->text().toDouble();
    params.viscosityFluid = ui->viscosityK->text().toDouble();
    params.densityMagnitude = ui->densityMag->text().toDouble();
    params.velocityMagnitude = ui->velocityMag->text().toDouble();

    setUIFromParameters(params);
    QMetaObject::invokeMethod(&cont_, "updateParameters", Q_ARG(SimParameters, params));
}

void MainWindow::setUIFromParameters(const SimParameters &params)
{

    if(params.simRunning)
    {
        ui->Simulate->setText(QString("Pause Simulation"));
        simRunning_ = true;
    }
    else
    {
        ui->Simulate->setText(QString("Start Simulation"));
        simRunning_ = false;
    }

    ui->timeStepEdit->setText(QString::number(params.timeStep));
}

void MainWindow::updateGL()
{

    {
        ui->GLWidget->tick();
    }
    ui->GLWidget->update();
}


void MainWindow::on_Simulate_clicked()
{
    simRunning_ = !simRunning_;
    setParametersFromUI();
}

void MainWindow::on_timeStepEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_clearScene_clicked()
{
    simRunning_ = false;
    setParametersFromUI();
    QMetaObject::invokeMethod(&cont_, "clearScene");
}

void MainWindow::on_source1_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_source2_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_source3_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_diffusionK_editingFinished()
{
   setParametersFromUI();
}


void MainWindow::on_viscosityK_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_densityMag_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_source4_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_velSrc1_clicked()
{
     setParametersFromUI();
}

void MainWindow::on_velSrc2_clicked()
{
     setParametersFromUI();
}

void MainWindow::on_velSrc3_clicked()
{
     setParametersFromUI();
}
