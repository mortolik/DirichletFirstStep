#include "mainwindow.h"
#include "DirichletDisplayWidget.hpp"
#include "DirichletSolverModel.hpp"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    m_tabWidget(new QTabWidget(this)),
    m_testModel(new DirichletSolverModel(this)),
    m_mainModel(new DirichletSolverModel(this))
{
    setCentralWidget(m_tabWidget);
    resize(800, 600);
    setWindowTitle("Лабораторная работа №1 - Задача Дирихле (вариант 4)");

    setupTestProblem();
    setupMainProblem();
}

void MainWindow::setupTestProblem()
{
    int n = 50;
    int m = 50;
    double omega = 1.8;
    double eps = 1e-6;
    int maxIter = 10000;

    m_testModel->setup(n, m, omega, eps, maxIter);
    m_testModel->solveTestProblem();

    auto *display = new DirichletDisplayWidget(m_testModel, true, this);
    m_tabWidget->addTab(display, "Тестовая задача");
}

void MainWindow::setupMainProblem()
{
    int n = 50;
    int m = 50;
    double omega = 1.8;
    double eps = 1e-6;
    int maxIter = 10000;

    m_mainModel->setup(n, m, omega, eps, maxIter);
    m_mainModel->solveMainProblem();

    auto *display = new DirichletDisplayWidget(m_mainModel, false, this);
    m_tabWidget->addTab(display, "Основная задача");
}


MainWindow::~MainWindow()
{
}
