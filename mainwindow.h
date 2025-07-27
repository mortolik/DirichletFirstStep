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
    void setupFirstTestProblem();
    void setupFirstMainProblem();

    QTabWidget *m_tabWidget;
    QTabWidget *m_stepOneTabs;

    DirichletSolverModel *m_testModel;
    DirichletSolverModel *m_mainModel;

    DirichletWidget *m_testWidget;
    DirichletWidget *m_mainWidget;
    QDockWidget* m_tableDock {nullptr};

    void updateDockWidget();
};
#endif // MAINWINDOW_H
