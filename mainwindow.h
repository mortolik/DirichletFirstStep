#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "Dirichle3StepModel.hpp"
#include "Dirichlet2step.hpp"
#include "DirichletSolverModel.hpp"
#include "DirichletSolverModel2.hpp"
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
    void setupFirstTestProblem();
    void setupFirstMainProblem();

    void setupSecond1TestProblem();
    void setupSecond1MainProblem();

    void setupSecond2TestProblem();
    void setupSecond2MainProblem();

    void setupThirdTestProblem();
    void setupThirdMainProblem();

    QTabWidget *m_tabWidget;
    QTabWidget *m_stepOneTabs;
    QTabWidget *m_stepTwo1Tabs;
    QTabWidget *m_stepTwo2Tabs;
    QTabWidget *m_stepThreeTabs;

    DirichletSolverModel *m_testModel;
    DirichletSolverModel *m_mainModel;

    DirichletSolverModel2 *m_test21Model;
    DirichletSolverModel2 *m_main21Model;

    Dirichlet2step *m_test22Model;
    Dirichlet2step *m_main22Model;

    Dirichle3StepModel *m_test3Model;
    Dirichle3StepModel *m_main3Model;

    DirichletWidget *m_testWidget;
    DirichletWidget *m_mainWidget;
    QDockWidget* m_tableDock {nullptr};

    void updateDockWidget(int index);
};
#endif // MAINWINDOW_H
