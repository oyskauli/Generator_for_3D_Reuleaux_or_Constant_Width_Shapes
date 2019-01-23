#ifndef CONSTANTWIDTHGEN_H
#define CONSTANTWIDTHGEN_H
#include <QVector>
#include <QVector3D>
#include "mainwindow.h"
#include <3rdparty/eigen-eigen-323c052e1731/Eigen/Core>

class ConstantWidthGen
{
public:
    ConstantWidthGen();
    //~ConstantWidthGen();
    void GenerateShapePoints(int npoints, double xvar, double yvar, double zvar, MainWindow *w, int start = 0);

    void import_points(int npoints, double xvar, double yvar, double zvar, QString filename, MainWindow *w);

    QVector<Eigen::Vector3d> keep_max_diam(MainWindow *w);

    QVector<Eigen::Vector3d> lim_max_diam(MainWindow *w, double tol = 1e-8);
    Eigen::Vector3d get_pair_shift(QVector<int> const &points, Eigen::Vector3d const a, Eigen::Vector3d const b, int i, int MAXIT);

    void savePointsAsModel(QVector<Eigen::Vector3d> &points, QString filename);
    Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, 1> points_main;
    Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, 1> points_oposite;
    Eigen::ArrayXd pitch;
    Eigen::ArrayXd ang;
    int npoints;
    double radius;
private:
   int imported_points = 0;
   bool isValidDiam(QVector<int> &points, Eigen::Vector3d a, Eigen::Vector3d b, double tol = 1e-8);
   void setElementsWSimilarAng(QVector<int> &satisfied,QVector<int> &sub_group, int c_point, double max_angle);
   QVector<Eigen::Vector3d> makeVectorOfPoints(QVector<int> &indexes);

   double totalAng(Eigen::Vector3d a, Eigen::Vector3d b);


};


#endif // CONSTANTWIDTHGEN_H
