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

private:
    DirichletSolverModel *m_model;
    bool m_isTest;

    DirichletWidget *m_chart;
    QTableWidget *m_table;
    QHBoxLayout *m_layout;

    void fillTable(const QVector<QVector<double>> &data);
};

