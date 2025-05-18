#include <QGroupBox>
#include <QHeaderView>
#include "DirichletSolverModel.hpp"
#include "DirichletDisplayWidget.hpp"

DirichletDisplayWidget::DirichletDisplayWidget(DirichletSolverModel *model, bool isTest, QWidget *parent)
    : QWidget(parent),
    m_model(model),
    m_isTest(isTest),
    m_layout(new QHBoxLayout(this))
{
    m_layout->setContentsMargins(10, 10, 10, 10);
    m_chart = new DirichletWidget(m_model, m_isTest, this);
    m_layout->addWidget(m_chart, 1);
    createTableBox();

    connect(m_chart, &DirichletWidget::solutionUpdated, this, &DirichletDisplayWidget::onSolutionUpdated);
    //updateTable();
}

QTabWidget *DirichletDisplayWidget::tableWidget()
{
    return m_tableTabs;
}

void DirichletDisplayWidget::fillTable(QTableWidget *table, const QVector<QVector<double>> &data)
{
    int cols = data.size();
    int rows = data[0].size();

    table->setRowCount(rows);
    table->setColumnCount(cols + 1);
    table->setHorizontalHeaderItem(0, new QTableWidgetItem("Y \\ X"));

    for (int i = 0; i < cols; ++i)
    {
        double x = m_model->a() + i * m_model->h();
        table->setHorizontalHeaderItem(i + 1, new QTableWidgetItem(QString::number(x, 'f', 2)));
    }

    for (int j = 0; j < rows; ++j)
    {
        double y = m_model->c() + j * m_model->k();
        table->setItem(j, 0, new QTableWidgetItem(QString::number(y, 'f', 2)));

        for (int i = 0; i < cols; ++i)
        {
            QString val = QString::number(data[i][j], 'f', 4);
            auto *item = new QTableWidgetItem(val);
            item->setTextAlignment(Qt::AlignCenter);
            table->setItem(j, i + 1, item);
        }
    }

    table->horizontalHeader()->setStretchLastSection(true);
    table->resizeColumnsToContents();
    table->resizeRowsToContents();
}

void DirichletDisplayWidget::updateTable()
{
    if (m_isTest)
    {
        fillTable(m_tableUStar, m_model->exactSolution());
        fillTable(m_tableV, m_model->solution());
        fillTable(m_tableDiff, m_model->error());
    }
    else
    {
        fillTable(m_tableV, m_model->solution());
    }
}

void DirichletDisplayWidget::onSolutionUpdated()
{
    updateTable();
}

QGroupBox *DirichletDisplayWidget::createTableBox()
{
    QGroupBox* box = new QGroupBox("Таблица значений");
    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setContentsMargins(5, 5, 5, 5);

    m_tableTabs = new QTabWidget(this);
    m_tableUStar = createEmptyTable();
    m_tableV = createEmptyTable();
    m_tableDiff = createEmptyTable();

    m_tableTabs->addTab(m_tableUStar, "u*(x,y)");
    m_tableTabs->addTab(m_tableV, "v(x,y)");
    m_tableTabs->addTab(m_tableDiff, "|u* - v|");

    layout->addWidget(m_tableTabs);
    return box;
}

QTableWidget *DirichletDisplayWidget::createEmptyTable()
{
    auto *table = new QTableWidget;
    table->setStyleSheet("QTableWidget { border: 1px solid #ccc; }");
    return table;
}
