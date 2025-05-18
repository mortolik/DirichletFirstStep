#include "mainwindow.h"
#include "DirichletSolverModel.hpp"
#include "DirichletDisplayWidget.hpp"

#include <QDockWidget>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    m_tabWidget(new QTabWidget(this)),
    m_testModel(new DirichletSolverModel(this)),
    m_mainModel(new DirichletSolverModel(this)),
    m_tableDock(new QDockWidget(this))
{
    setCentralWidget(m_tabWidget);
    setMinimumSize(800, 600);
    setWindowTitle("Лабораторная работа №1 - Задача Дирихле (вариант 4)");

    m_tableDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    addDockWidget(Qt::RightDockWidgetArea, m_tableDock);

    setupTestProblem();
    setupMainProblem();

    connect(m_tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateDockWidget);
    updateDockWidget(m_tabWidget->currentIndex());

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

void MainWindow::updateDockWidget(int index)
{
    auto *currentDisplay = qobject_cast<DirichletDisplayWidget*>(m_tabWidget->widget(index));
    if (currentDisplay)
    {
        m_tableDock->setWidget(currentDisplay->tableWidget());
    }
}
