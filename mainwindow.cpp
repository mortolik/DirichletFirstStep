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
    int n = 100;
    int m = 100;
    double omega = 1.8;
    double eps = 1e-6;
    int maxIter = 10000;

    m_testModel->setup(n, m, omega, eps, maxIter);
    m_testModel->solveTestProblem();

    // Создаём внутренний QTabWidget
    QTabWidget *testTabs = new QTabWidget;

    QString rep = m_testModel->reportString(true);
    auto *vTable = new DirichletDisplayWidget(m_testModel, true, this);
    testTabs->addTab(vTable, "v(x,y) + таблица");

    // Uncomment when need error widget
    // auto *errorWidget = new DirichletWidget(new QVector<QVector<double>>(m_testModel->error()), 1.0, 2.0, 2.0, 3.0, this);
    // testTabs->addTab(errorWidget, "Ошибка |u* - v|");

    m_tabWidget->addTab(testTabs, "Тестовая задача");
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

    // Вкладка 1: численное решение
    DirichletWidget *mainSol = new DirichletWidget(m_mainModel, false, this);

    // Вкладка 2: сравнение с сеткой 2n × 2m
    double eps2;
    QVector<QVector<double>> diff = m_mainModel->compareWithFinerGrid(n * 2, m * 2, eps2);

    QString accReport = m_mainModel->reportString(false, eps2);
    auto *accuracyWidget = new  DirichletWidget(m_mainModel, false, this);


    QTabWidget *mainTabs = new QTabWidget;
    mainTabs->addTab(mainSol, "Численное решение");
    mainTabs->addTab(accuracyWidget, QString("Разность с сеткой 2n×2m (ε₂ ≈ %1)").arg(eps2, 0, 'e', 2));

    m_tabWidget->addTab(mainTabs, "Основная задача");
}


MainWindow::~MainWindow()
{
}
