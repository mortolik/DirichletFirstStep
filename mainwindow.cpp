#include "mainwindow.h"
#include "DirichletSolverModel.hpp"
#include "DirichletDisplayWidget.hpp"

#include <QDockWidget>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    m_tabWidget(new QTabWidget(this)),
    m_stepOneTabs(new QTabWidget(this)),
    m_testModel(new DirichletSolverModel(this)),
    m_mainModel(new DirichletSolverModel(this)),
    m_tableDock(new QDockWidget(this))
{
    setCentralWidget(m_tabWidget);
    setMinimumSize(800, 600);
    setWindowTitle("Лабораторная работа №1 - Задача Дирихле (вариант 4)");

    m_tableDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    addDockWidget(Qt::RightDockWidgetArea, m_tableDock);

    setupFirstTestProblem();
    setupFirstMainProblem();

    QWidget *stepOneContainer = new QWidget(this);
    QVBoxLayout *stepOneLayout = new QVBoxLayout(stepOneContainer);
    stepOneLayout->setContentsMargins(0, 0, 0, 0);
    stepOneLayout->addWidget(m_stepOneTabs);

    m_tabWidget->addTab(stepOneContainer, "Первая ступень");

    connect(m_tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateDockWidget);
    updateDockWidget();


    showMaximized();
}

void MainWindow::setupFirstTestProblem()
{
    auto *display = new DirichletDisplayWidget(m_testModel, true, this);
    m_stepOneTabs->addTab(display, "Тестовая задача");
}

void MainWindow::setupFirstMainProblem()
{
    auto *display = new DirichletDisplayWidget(m_mainModel, false, this);
    m_stepOneTabs->addTab(display, "Основная задача");
}

MainWindow::~MainWindow()
{
}

void MainWindow::updateDockWidget()
{
    if (m_stepOneTabs)
    {
        int innerIndex = m_stepOneTabs->currentIndex();
        auto* currentDisplay = qobject_cast<DirichletDisplayWidget*>(m_stepOneTabs->widget(innerIndex));
        if (currentDisplay)
        {
            m_tableDock->setWidget(currentDisplay->tableWidget());
        }

        disconnect(m_stepOneTabs, nullptr, this, nullptr);
        connect(m_stepOneTabs, &QTabWidget::currentChanged, this, &MainWindow::updateDockWidget);
    }
}

