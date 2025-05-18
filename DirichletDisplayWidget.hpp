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
    void updateTable();
    QTabWidget* tableWidget();
private slots:
    void onSolutionUpdated();
private:
    DirichletSolverModel *m_model;
    bool m_isTest;

    DirichletWidget *m_chart;
    QTabWidget *m_tableTabs;
    QHBoxLayout *m_layout;

    QTableWidget *m_tableUStar;
    QTableWidget *m_tableV;
    QTableWidget *m_tableDiff;

    QTableWidget *m_tableV2;
    QTableWidget *m_tableVDiff;

    QGroupBox *createTableBox();
    QTableWidget *createEmptyTable();
    void fillTable(QTableWidget *table, const QVector<QVector<double>> &data);
};

