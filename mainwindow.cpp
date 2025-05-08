#include "mainwindow.h"
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

    // Создаём внутренний QTabWidget
    QTabWidget *testTabs = new QTabWidget;

    DirichletWidget *exactVsNumerical = new DirichletWidget(m_testModel, true, this);
    testTabs->addTab(exactVsNumerical, "u*, v");

    auto *errorWidget = new DirichletWidget(new QVector<QVector<double>>(m_testModel->error()), 1.0, 2.0, 2.0, 3.0, this);
    testTabs->addTab(errorWidget, "Ошибка |u* - v|");

    m_tabWidget->addTab(testTabs, "Тестовая задача");
}


void MainWindow::setupMainProblem()
{
    // Настройка основной задачи (вариант 4)
    int n = 50;
    int m = 50;
    double omega = 1.8;
    double eps = 1e-6;
    int maxIter = 10000;

    m_mainModel->setup(n, m, omega, eps, maxIter);
    m_mainModel->solveMainProblem();

    m_mainWidget = new DirichletWidget(m_mainModel, false, this);
    m_tabWidget->addTab(m_mainWidget, "Основная задача");
}

MainWindow::~MainWindow()
{
}
