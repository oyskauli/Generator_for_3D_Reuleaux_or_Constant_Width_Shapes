#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void SetDebugText(QString string);
    void AddDebugText(QString string);
    void SetPage(int p);
    void SetProgressTask(int p);
    void SetProgressSubTask(int p);

private slots:

    void on_genShape_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
