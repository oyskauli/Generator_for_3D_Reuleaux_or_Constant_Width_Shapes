#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "constantwidthgen.h"
#include <3rdparty/eigen-eigen-323c052e1731/Eigen/Core>
#include <QScrollBar>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    SetProgressTask(0);
    SetProgressSubTask(0);
    SetPage(0);
}

void MainWindow::SetDebugText(QString string)
{
    ui->debugText->setText(string);
}

void MainWindow::AddDebugText(QString string)
{
    this->SetDebugText(this->ui->debugText->toPlainText() + string);
    ui->debugText->verticalScrollBar()->setValue(ui->debugText->verticalScrollBar()->maximum());
}

void MainWindow::SetPage(int p)
{
    ui->stackedWidget->setCurrentIndex(p);
    QCoreApplication::processEvents();
}

void MainWindow::SetProgressTask(int p)
{
    ui->progress_Task->setValue(p);
    QCoreApplication::processEvents();
}

void MainWindow::SetProgressSubTask(int p)
{
    ui->progress_subTask->setValue(p);
    QCoreApplication::processEvents();
}



MainWindow::~MainWindow()
{
    delete ui;

}



void MainWindow::on_genShape_clicked()
{
    SetPage(1);

    QString toPrint;
    int npoints = 4000;
    int nmodels = 200;
    Eigen::Vector3d shift;
    shift << 0.00, 0.00, 0.008;
    ConstantWidthGen gen = ConstantWidthGen();
    gen.GenerateShapePoints(npoints, .0, .0, .0, this);
    gen.points_main[0] << 0, 0, 1;
    gen.points_oposite[0] << 0, 0, -1;
    //gen.import_points(npoints, .3, .3, .3, "C:/Users/OysteinAdm/Documents/untitled2.obj", this);
    for(int i = 0; i < nmodels; i++){
        auto outp = gen.lim_max_diam(this, 1e-4);

        //save to File
        auto filename = ui->outPath->text() + QString::number(i);
        gen.savePointsAsModel(outp, filename);
        gen.points_main[0] += shift;
        gen.points_oposite[0] += shift;
    }

}
