#include "mainwindow.h"
#include "Dirichle3StepModel.hpp"
#include "DirichletSolverModel.hpp"
#include "DirichletDisplayWidget.hpp"

#include <QDockWidget>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    m_tabWidget(new QTabWidget(this)),
    m_stepOneTabs(new QTabWidget(this)),
    m_stepTwo1Tabs(new QTabWidget(this)),
    m_stepTwo2Tabs(new QTabWidget(this)),
    m_stepThreeTabs(new QTabWidget(this)),
    m_testModel(new DirichletSolverModel(this)),
    m_mainModel(new DirichletSolverModel(this)),
    m_test21Model(new DirichletSolverModel2(this)),
    m_main21Model(new DirichletSolverModel2(this)),
    m_test22Model(new Dirichlet2step(this)),
    m_main22Model(new Dirichlet2step(this)),
    m_test3Model(new Dirichle3StepModel(this)),
    m_main3Model(new Dirichle3StepModel(this)),
    m_tableDock(new QDockWidget(this))
{
    setCentralWidget(m_tabWidget);
    setMinimumSize(800, 600);
    setWindowTitle("Лабораторная работа №1 - Задача Дирихле (вариант 4)");

    m_tableDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    addDockWidget(Qt::RightDockWidgetArea, m_tableDock);

    setupFirstTestProblem();
    setupFirstMainProblem();

    setupSecond1TestProblem();
    setupSecond1MainProblem();

    setupSecond2TestProblem();
    setupSecond2MainProblem();

    setupThirdTestProblem();
    setupThirdMainProblem();

    QWidget *stepOneContainer = new QWidget(this);
    QVBoxLayout *stepOneLayout = new QVBoxLayout(stepOneContainer);
    stepOneLayout->setContentsMargins(0, 0, 0, 0);
    stepOneLayout->addWidget(m_stepOneTabs);

    QWidget *stepTwo1Container = new QWidget(this);
    QVBoxLayout *stepTwo1Layout = new QVBoxLayout(stepTwo1Container);
    stepTwo1Layout->setContentsMargins(0, 0, 0, 0);
    stepTwo1Layout->addWidget(m_stepTwo1Tabs);

    QWidget *stepTwo2Container = new QWidget(this);
    QVBoxLayout *stepTwo2Layout = new QVBoxLayout(stepTwo2Container);
    stepTwo2Layout->setContentsMargins(0, 0, 0, 0);
    stepTwo2Layout->addWidget(m_stepTwo2Tabs);

    QWidget *stepThreeContainer = new QWidget(this);
    QVBoxLayout *stepThreeLayout = new QVBoxLayout(stepThreeContainer);
    stepThreeLayout->setContentsMargins(0, 0, 0, 0);
    stepThreeLayout->addWidget(m_stepThreeTabs);

    m_tabWidget->addTab(stepOneContainer, "Первая ступень");
    m_tabWidget->addTab(stepTwo1Container, "Вторая ступень (МПИ)");
    m_tabWidget->addTab(stepTwo2Container, "Вторая ступень (МСГ)");
    m_tabWidget->addTab(stepThreeContainer, "Третья ступень");

    connect(m_tabWidget, &QTabWidget::currentChanged, this, &MainWindow::updateDockWidget);
    updateDockWidget(m_tabWidget->currentIndex());


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

void MainWindow::setupSecond1TestProblem()
{
    auto *display = new DirichletDisplayWidget(m_test21Model, true, this);
    m_stepTwo1Tabs->addTab(display, "Тестовая задача");
}

void MainWindow::setupSecond1MainProblem()
{
    auto *display = new DirichletDisplayWidget(m_main21Model, false, this);
    m_stepTwo1Tabs->addTab(display, "Основная задача");
}

void MainWindow::setupSecond2TestProblem()
{
    auto *display = new DirichletDisplayWidget(m_test22Model, true, this);
    m_stepTwo2Tabs->addTab(display, "Тестовая задача");
}

void MainWindow::setupSecond2MainProblem()
{
    auto *display = new DirichletDisplayWidget(m_main22Model, false, this);
    m_stepTwo2Tabs->addTab(display, "Основная задача");
}

void MainWindow::setupThirdTestProblem()
{
    auto *display = new DirichletDisplayWidget(m_test3Model, true, this);
    m_stepThreeTabs->addTab(display, "Тестова задача");
}

void MainWindow::setupThirdMainProblem()
{
    auto *display = new DirichletDisplayWidget(m_main3Model, false, this);
    m_stepThreeTabs->addTab(display, "Основная задача");
}

MainWindow::~MainWindow()
{
}

void MainWindow::updateDockWidget(int index)
{
    QTabWidget* currentStepTabs = nullptr;

    int outerIndex = m_tabWidget->currentIndex();

    switch (outerIndex)
    {
    case 0: currentStepTabs = m_stepOneTabs; break;
    case 1: currentStepTabs = m_stepTwo1Tabs; break;
    case 2: currentStepTabs = m_stepTwo2Tabs; break;
    case 3: currentStepTabs = m_stepThreeTabs; break;
    }

    if (currentStepTabs)
    {
        int innerIndex = currentStepTabs->currentIndex();
        auto* currentDisplay = qobject_cast<DirichletDisplayWidget*>(currentStepTabs->widget(innerIndex));
        if (currentDisplay)
        {
            m_tableDock->setWidget(currentDisplay->tableWidget());
        }

        // Подключаем обновление на смену внутренней вкладки
        disconnect(currentStepTabs, nullptr, this, nullptr);
        connect(currentStepTabs, &QTabWidget::currentChanged, this, &MainWindow::updateDockWidget);
    }
}

