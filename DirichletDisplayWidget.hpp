#pragma once
#include <QWidget>
#include <QHBoxLayout>
#include <QTableWidget>
#include "DirichletWidget.hpp"

class DirichletDisplayWidget : public QWidget
{
    Q_OBJECT

public:
    DirichletDisplayWidget(DirichletSolverModel *model, bool isTest, QWidget *parent = nullptr);
    QTableWidget* tableWidget();
private slots:
    void onSolutionUpdated();
private:
    DirichletSolverModel *m_model;
    bool m_isTest;

    DirichletWidget *m_chart;
    QTableWidget *m_table;
    QHBoxLayout *m_layout;

    QGroupBox* createTable();
    void fillTable(const QVector<QVector<double>> &data);
    void updateTable();
};

