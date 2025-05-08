#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "DirichletSolverModel.hpp"
#include "DirichletWidget.hpp"
#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void setupTestProblem();
    void setupMainProblem();

    QTabWidget *m_tabWidget;
    DirichletSolverModel *m_testModel;
    DirichletSolverModel *m_mainModel;
    DirichletWidget *m_testWidget;
    DirichletWidget *m_mainWidget;
};
#endif // MAINWINDOW_H
