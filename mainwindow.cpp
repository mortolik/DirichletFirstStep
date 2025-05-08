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
    // Настройка тестовой задачи (вариант 4)
    int n = 50; // количество разбиений по x
    int m = 50; // количество разбиений по y
    double omega = 1.8; // параметр релаксации
    double eps = 1e-6; // точность
    int maxIter = 10000; // максимальное число итераций

    m_testModel->setup(n, m, omega, eps, maxIter);
    m_testModel->solveTestProblem();

    m_testWidget = new DirichletWidget(m_testModel);
    m_tabWidget->addTab(m_testWidget, "Тестовая задача");
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
    // Здесь нужно добавить решение основной задачи

    m_mainWidget = new DirichletWidget(m_mainModel);
    m_tabWidget->addTab(m_mainWidget, "Основная задача");
}

MainWindow::~MainWindow()
{
}
