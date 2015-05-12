#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

class Controller;
struct SimParameters;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(Controller &cont, int fps, QWidget *parent = 0);
    ~MainWindow();   

public slots:
    void setUIFromParameters(const SimParameters &params);

private slots:
    void updateGL();

    void on_Simulate_clicked();

    void on_timeStepEdit_editingFinished();

    void on_clearScene_clicked();

    void on_source1_clicked();

    void on_source2_clicked();

    void on_source3_clicked();

    void on_diffusionK_editingFinished();

    void on_viscosityK_editingFinished();

    void on_densityMag_editingFinished();

private:
    Controller &cont_;
    Ui::MainWindow *ui;
    bool simRunning_;
    QTimer renderTimer_;

    void setParametersFromUI();
};

#endif // MAINWINDOW_H
