#include "constantWidthGen.h"
#include <QtMath>
#include <QVector>
#include <QVector3D>
#include <QString>
#include <QFile>
#include <QTextStream>
#include <3rdparty/eigen-eigen-323c052e1731/Eigen/Core>
#include "random"


ConstantWidthGen::ConstantWidthGen()
{

}

/*ConstantWidthGen::~ConstantWidthGen()
{
    int rows = sizeof(points) / sizeof(points[0]);

    for (int i = 0; i < rows; i++) {
        delete[] points[i];
    }

    delete[] points;
}*/

void ConstantWidthGen::savePointsAsModel(QVector<Eigen::Vector3d> &points, QString filename){
    int nleft = points.length();
    QFile file(filename);
    file.remove();
    if(file.open(QIODevice::ReadWrite)){
        QTextStream stream(&file);
        for (int i = 0; i < nleft; i++) {
            stream << "v ";
            for (int j = 0; j < 3; j++) {
                stream << QString::number(points[i][j]) << " ";
            }
            stream << "\n";
        }
    }
}

void ConstantWidthGen::GenerateShapePoints(int npoints, double xvar, double yvar, double zvar, MainWindow *w, int start)
{
    //ui update helpers
    int ITERBEFOREUPDATE = 10;
    int itercount = 0;
    w->SetProgressSubTask(100);//no subtask here
    w->SetDebugText("Generating random point cloud...");

    double radius = 1;
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<double> ang_dist(0, 1);
    std::normal_distribution<double> pos_dist(0, 1);
    double r_ptch;
    double r_ang;

    auto points_main = Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, 1>(npoints);
    auto points_oposite = Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, 1>(npoints);
    auto pitch = Eigen::ArrayXd(npoints);
    auto ang = Eigen::ArrayXd(npoints);

    //transform to get approx. even dist. of points
    auto f = [](double v){
        return (asin(v))/M_PI;
    };

    //QVector<QVector3D> points;
    for (int i = start; i < npoints; i++)
    {
        //angle from 0 to 1
        r_ptch = f(ang_dist(e2));

        //angle from 0 to 2
        r_ang = ang_dist(e2)*2;
        //this only needs to cover upper half since the opposite point is also always added

        pitch[i] = r_ptch;
        ang[i] = r_ang;

        //convert to radians
        r_ptch = r_ptch*M_PI;
        r_ang = r_ang*M_PI;

        double tmp = cos(r_ptch);
        points_main[i] << sin(r_ptch), sin(r_ang)*tmp, cos(r_ang)*tmp;
        points_main[i] *= radius;
        points_oposite[i] << -1*points_main[i];

        Eigen::Vector3d loc(pos_dist(e2)*xvar, pos_dist(e2)*yvar, pos_dist(e2)*zvar);

        points_main[i] += loc;
        points_oposite[i] += loc;

        if(itercount < ITERBEFOREUPDATE){
            itercount++;
        }else{
            w->SetProgressTask((i*100)/npoints + 1);
        }

    }

    this->points_main = points_main;
    this->points_oposite = points_oposite;
    this->pitch = pitch;
    this->ang = ang;
    this->npoints = npoints;
    this->radius = radius;
}

void ConstantWidthGen::import_points(int npoints, double xvar, double yvar, double zvar,QString filename, MainWindow *w){
    QVector<Eigen::Vector3d> imported;
    Eigen::Vector3d tmp;
    QFile inputFile(filename);
    if (inputFile.open(QIODevice::ReadOnly))
    {
       QTextStream in(&inputFile);
       while (!in.atEnd())
       {
          QString line = in.readLine();
          if(line[0] == "v"){
              auto list = line.split(" ", QString::SkipEmptyParts);
              if(list.length() == 4){
                  tmp << list[1].toDouble(), list[2].toDouble(), list[3].toDouble();
                  imported.append(tmp);
              }else{
                  w->AddDebugText("Expexted 3 numbers got: " + line);
              }
          }
       }
       inputFile.close();

       this->imported_points = imported.length();
       int start = this->imported_points/2;
       this->GenerateShapePoints(npoints + start, xvar, yvar, zvar, w, start);

       for(int i = 0; i < start; i++){
           this->points_main[i] = imported[i*2];
           this->points_oposite[i] = imported[i*2 + 1];
       }

    }
}

QVector<Eigen::Vector3d> ConstantWidthGen::keep_max_diam(MainWindow *w){
    /*
    Generates a set of points valid under max diam. criteria
    If a pair arent by default valid they will be removed
    */
    //ui update helpers
    int ITERBEFOREUPDATE = 10;
    int itercount = 0;
    w->SetProgressSubTask(100);//no subtask here
    w->SetDebugText(QString("Finding points that already satify max diam criteria, ") + QString::number(std::pow(this->radius*2, 2)) + QString(", \n"));

    auto satisfied = QVector<int>(this->imported_points);
    if(this->imported_points == 0){
        satisfied.append(0);
    }else{
        for(int i = 0; i < this->imported_points; i++){
            satisfied[i] = i;
        }
    }

    for (int i = 1; i < this->npoints; i++) {
        /*
        max_len1 = 0;
        for(int k = 0; k < satisfied.length(); k++){
            //max len from i pair to k pair
            len1 = std::max((satisfied[k] - this->points_main[i]).squaredNorm(),
                            (satisfied[k] - this->points_oposite[i]).squaredNorm());

            max_len1 = std::max(max_len1, len1);
        }

        if(std::pow(this->radius*2, 2) > max_len1){
            satisfied.append(points_main[i]);
            satisfied.append(points_oposite[i]);
            w->AddDebugText("Found a pair!, ");
            w->AddDebugText(QString::number(max_len1));
            w->AddDebugText(", \n");
        }
        */

        if(isValidDiam(satisfied, points_main[i], points_oposite[i])){
            satisfied.append(i);
            w->AddDebugText("Found a pair!, \n");
        }

        if(itercount < ITERBEFOREUPDATE){
            itercount++;
        }else{
            w->SetProgressTask((i*100)/npoints + 1);
        }
    }
    return makeVectorOfPoints(satisfied);
}

QVector<Eigen::Vector3d> ConstantWidthGen::lim_max_diam(MainWindow *w, double tol){
/*
Generates a set of points valid under max diam. criteria
If a pair arent by default valid they will be shifted up to MAXITER times, if still not valid it will not be added to set
*/
    //ui update helpers
    int ITERBEFOREUPDATE = 10;
    int itercount = 0;

    w->SetDebugText("Converging points to set...\n");

    int MAXITER = 2000;
    int ITERPREVELOCITY = 100;//the first shifts tend to be large so velocity method can be divergent
    int MINSUBGROUP = 40;//minimum size of iteration poits for it to be used;

    //COULD ME MUCH LESS WITH SMART SELECTION OF POINTS, SELECTING ANGLES "AROUND" THE CURRENT ANGLE

    auto iteration_points = QVector<int>(4);

    auto satisfied = QVector<int>();//index of points that satisfy
    //auto debug_points = QVector<Eigen::Vector3d>();
    satisfied.append(0);

    bool valid;
    bool completed_iter;

    Eigen::Vector3d a;
    Eigen::Vector3d b;
    Eigen::Vector3d shift;
    Eigen::Vector3d velocity;

    double damping = 0.02;

    //loop over all points
    for (int i = 1; i < this->npoints; i++) {
        a = points_main[i];
        b = points_oposite[i];
        completed_iter = false;
/*
        //first run iteration on subset
        if(i > 100){
            //see if this pair is valid by default
            setElementsWSimilarAng(satisfied, iteration_points, i, M_PI/16);
            valid = isValidDiam(iteration_points, a, b, tol);

            if(!valid){
            //    debug_points.append(a);
            //    debug_points.append(b);
                //Shift without velocity, less likely to diverge
                for(int k = 0; k < ITERPREVELOCITY; k++){
                    shift = get_pair_shift(iteration_points, a, b, k, MAXITER);
                    a = a + shift;
                    b = b + shift;


                    if(shift.norm() > this->radius*5){
                        w->SetDebugText("worriyngly large shift (31)...\n");
                        break;
                    }
                    if(isValidDiam(iteration_points, a, b, tol)){
                    //    debug_points.clear();
                        break;
                    }
                    //debug_points.append(a);
                    //debug_points.append(b);
                }
                //last part is slow, simulate velocity to speed it up
                if(!completed_iter){
                    velocity << 0, 0, 0;
                    for(int k = ITERPREVELOCITY; k < MAXITER; k++){
                        shift = get_pair_shift(iteration_points, a, b, k, MAXITER);
                        velocity = (velocity + shift)*(1-damping);
                        a = a + velocity;
                        b = b + velocity;

                        if(i > 50){
                            setElementsWSimilarAng(satisfied, iteration_points, i);
                            auto pointsToPrint = makeVectorOfPoints(iteration_points);
                            Eigen::Vector3d t;
                            t << 2, 0, 0;
                            pointsToPrint.append(a + t);
                            pointsToPrint.append(b + t);
                            savePointsAsModel(pointsToPrint, "C:/Users/OysteinAdm/Desktop/debug.obj");
                            return makeVectorOfPoints(satisfied);
                        }
                        if(shift.norm() > this->radius*5){
                            w->SetDebugText(("worriyngly large shift (1)...\n"));
                            break;
                        }
                        if(isValidDiam(iteration_points, a, b, tol)){
                            completed_iter = true;
                            if(!isValidDiam(satisfied, a, b, tol)){
                                w->SetDebugText("NOT VALID AFTER SUB ITER\n(Exited and saved debug files)\n");
                                auto tr = makeVectorOfPoints(iteration_points);
                                Eigen::Vector3d m;
                                m << 2, 0, 0;
                                tr.append(a+m);
                                tr.append(b+m);
                                savePointsAsModel(tr, "C:/Users/OysteinAdm/Desktop/iterpoints.obj");
                                return makeVectorOfPoints(satisfied);

                            }else{
                                w->SetDebugText("VALID AFTER SUB ITER!!!\n");
                            }
                            //debug_points.clear();
                            break;
                        }
                        //debug_points.append(a);
                        //debug_points.append(b);
                    }
                }
            }
        }//----END SUBSET ITERATION---
*/
        //see if this pair is valid by default or after sub set iteration
        valid = isValidDiam(satisfied, a, b, tol);

        if(!valid){
        //    debug_points.append(a);
        //    debug_points.append(b);
            //Shift without velocity, less likely to diverge
            for(int k = 0; k < ITERPREVELOCITY; k++){
                shift = get_pair_shift(satisfied, a, b, k, MAXITER);
                a = a + shift;
                b = b + shift;


                if(shift.norm() > this->radius*5){
                    w->SetDebugText("worriyngly large shift (31)...\n");
                    break;
                }
                if(isValidDiam(satisfied, a, b, tol)){
                    completed_iter = true;
                //    debug_points.clear();
                    break;
                }
                //debug_points.append(a);
                //debug_points.append(b);
            }
            //last part is slow, simulate velocity to speed it up
            if(!completed_iter){
                velocity << 0, 0, 0;
                for(int k = ITERPREVELOCITY; k < MAXITER; k++){
                    shift = get_pair_shift(satisfied, a, b, k, MAXITER);
                    velocity = (velocity + shift)*(1-damping);
                    a = a + velocity;
                    b = b + velocity;

                    /*if(i > 50){
                        setElementsWSimilarAng(satisfied, iteration_points, i);
                        auto pointsToPrint = makeVectorOfPoints(iteration_points);
                        Eigen::Vector3d t;
                        t << 2, 0, 0;
                        pointsToPrint.append(a + t);
                        pointsToPrint.append(b + t);
                        savePointsAsModel(pointsToPrint, "C:/Users/OysteinAdm/Desktop/debug.obj");
                        return makeVectorOfPoints(satisfied);
                    }*/
                    if(shift.norm() > this->radius*5){
                        w->SetDebugText(("worriyngly large shift (1)...\n"));
                        break;
                    }
                    if(isValidDiam(satisfied, a, b, tol)){
                        completed_iter = true;
                        //debug_points.clear();
                        break;
                    }
                    //debug_points.append(a);
                    //debug_points.append(b);
                }
            }
        }
        if(valid || completed_iter){
            satisfied.append(i);
            this->points_main[i] = a;
            this->points_oposite[i] = b;

            if(completed_iter){
                w->AddDebugText("Converged\n");
            }
            //w->AddDebugText("Found a point,\n");
        }else{
            w->AddDebugText("Did not converge,\n");
        }

        if(itercount < ITERBEFOREUPDATE){
            itercount++;
        }else{
            w->SetProgressTask((i*100)/npoints + 1);
            w->SetDebugText("");
        }
    }
    //savePointsAsModel(debug_points, "C:/Users/OysteinAdm/Desktop/debug.obj");
    return makeVectorOfPoints(satisfied);
}

Eigen::Vector3d ConstantWidthGen::get_pair_shift(QVector<int> const &points, Eigen::Vector3d const a, Eigen::Vector3d const b, int i, int MAXIT){
/*
Determines how points a, b should be shifted to be valid in points
note, |a-b| assumed to be radius*2
i - iteration index
MAXIT - max iteration
prevPoints - previous positions of a
*/
    int npoints = points.length();

    //start at zero shift
    Eigen::Vector3d shift_total;
    Eigen::Vector3d shift_point;

    //average shift
    double len = 0;
    shift_total << 0, 0, 0;
    int c = 0;
    for (int i = 0; i < npoints; i++){
        len = (this->points_main[points[i]] - a).norm();
        if((len - this->radius*2) > 0){
            shift_total += (this->points_main[points[i]] - a).normalized()*(len-2);
            c++;
        }
        len = (this->points_main[points[i]] - b).norm();
        if((len - this->radius*2) > 0){
            shift_total += (this->points_main[points[i]] - b).normalized()*(len-2);
            c++;
        }


        len = (this->points_oposite[points[i]] - a).norm();
        if((len - this->radius*2) > 0){
            shift_total += (this->points_oposite[points[i]]- a).normalized()*(len-2);
            c++;
        }
        len = (this->points_oposite[points[i]] - b).norm();
        if((len - this->radius*2) > 0){
            shift_total += (this->points_oposite[points[i]] - b).normalized()*(len-2);
            c++;
        }

    }
    shift_total = shift_total/c;

/*
    //max Shift

    double len = 0;
    double tmpa = 1;
    double tmpb = 1;

    for (int i = 0; i < npoints; i++){
        tmpa = (points[i] - a).norm() - this->radius*2;
        tmpb = (points[i] - b).norm() - this->radius*2;
        if(tmpa > len){
            shift_total = (points[i] - a).normalized()*tmpa*(1+tmpa);
            len = tmpa;
        }
        if(tmpb > len){
            shift_total = (points[i] - b).normalized()*tmpb*(1+tmpb);
            len = tmpb;
        }

    }
*/
    return shift_total;
}

bool ConstantWidthGen::isValidDiam(QVector<int> &points, Eigen::Vector3d a, Eigen::Vector3d b, double tol){
/*
Determines if points a, b valid in points under the max diam. criteria
note, |a-b| assumed to be radius*2
*/
    double len;

    for(int k = 0; k < points.length(); k++){
        //max len from i pair to k pair
        len = std::max((this->points_main[points[k]] - a).squaredNorm(),
                        (this->points_main[points[k]] - b).squaredNorm());

        if((len - std::pow(this->radius*2, 2)) > tol) return false;

        len = std::max((this->points_oposite[points[k]] - a).squaredNorm(),
                        (this->points_oposite[points[k]] - b).squaredNorm());

        if((len - std::pow(this->radius*2, 2)) > tol) return false;
    }
    return true;
}

void ConstantWidthGen::setElementsWSimilarAng(QVector<int> &satisfied, QVector<int> &sub_group, int c_point, double max_angle){
/*
    double min_val =  cos(max_angle);
    Eigen::Vector3d c_vec = (this->points_main[c_point] - this->points_oposite[c_point]).normalized();
    Eigen::Vector3d vec;
    double val;
    for(int i : satisfied){
        vec = (this->points_main[i] - this->points_oposite[i]).normalized();
        val = c_vec.dot(vec);
        if(val > min_val){
            sub_group.append(i);
        }
    }
*/
    //Find projection matrix with normal equal to the vector between the pair c_point
    Eigen::Vector3d plane_norm = (this->points_main[c_point] - this->points_oposite[c_point]).normalized();

    Eigen::Matrix3d identity;
    identity.setIdentity();
    Eigen::Matrix3d project_plane = identity - plane_norm * plane_norm.transpose();

}

QVector<Eigen::Vector3d> ConstantWidthGen::makeVectorOfPoints(QVector<int> &indexes){
    auto out = QVector<Eigen::Vector3d>(indexes.length()*2);
    for(int i = 0; i < indexes.length(); i++){
        out[i*2] = this->points_main[indexes[i]];
        out[i*2 + 1] = this->points_oposite[indexes[i]];
    }
    return out;
}

double ConstantWidthGen::totalAng(Eigen::Vector3d a, Eigen::Vector3d b){
    return acos(a.normalized().dot(b.normalized()));
}
