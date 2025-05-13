#include "mainwindow.h"
#include "DirichletSolverModel.hpp"
#include "DirichletDisplayWidget.hpp"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    m_tabWidget(new QTabWidget(this)),
    m_testModel(new DirichletSolverModel(this)),
    m_mainModel(new DirichletSolverModel(this))
{
    setCentralWidget(m_tabWidget);
    setMinimumSize(800, 600);
    setWindowTitle("Лабораторная работа №1 - Задача Дирихле (вариант 4)");

    setupTestProblem();
    setupMainProblem();
    showMaximized();
}

void MainWindow::setupTestProblem()
{
    auto *display = new DirichletDisplayWidget(m_testModel, true, this);
    m_tabWidget->addTab(display, "Тестовая задача");
}

void MainWindow::setupMainProblem()
{
    auto *display = new DirichletDisplayWidget(m_mainModel, false, this);
    m_tabWidget->addTab(display, "Основная задача");
}

MainWindow::~MainWindow()
{
}
