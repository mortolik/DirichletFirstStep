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
    int rows = data.size();
    int cols = data[0].size();

    table->setRowCount(rows);
    table->setColumnCount(cols + 1);

    // Определяем, нужно ли использовать половинный шаг (только для m_tableV2)
    bool isV2Table = (table == m_tableV2);
    double hStep = isV2Table ? m_model->h() / 2 : m_model->h();
    double kStep = isV2Table ? m_model->k() / 2 : m_model->k();

    table->setHorizontalHeaderItem(0, new QTableWidgetItem("Y \\ X"));
    for (int i = 0; i < cols; ++i)
    {
        double x = m_model->a() + i * hStep; // Используем hStep (обычный или половинный)
        table->setHorizontalHeaderItem(i + 1, new QTableWidgetItem(QString::number(x, 'f', 3)));
    }

    for (int j = 0; j < rows; ++j)
    {
        double y = m_model->c() + j * kStep; // Используем kStep (обычный или половинный)
        table->setItem(j, 0, new QTableWidgetItem(QString::number(y, 'f', 3)));


    for (int i = 0; i < cols; ++i)
    {
        QString val = QString::number(data[j][i], 'f', 5);
        auto *item = new QTableWidgetItem(val);
        item->setTextAlignment(Qt::AlignCenter);
        table->setItem(i, j + 1, item);
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
        auto result = m_model->computeFinerGridComparison();
        fillTable(m_tableV, result.v);
        fillTable(m_tableV2, result.v2);
        fillTable(m_tableVDiff, result.diff);
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

    if (m_isTest)
    {
        m_tableUStar = createEmptyTable();
        m_tableV = createEmptyTable();
        m_tableDiff = createEmptyTable();

        m_tableTabs->addTab(m_tableUStar, "u*(x,y)");
        m_tableTabs->addTab(m_tableV, "v(x,y)");
        m_tableTabs->addTab(m_tableDiff, "|u* - v|");
    }
    else
    {
        m_tableV = createEmptyTable();
        m_tableV2 = createEmptyTable();
        m_tableVDiff = createEmptyTable();

        m_tableTabs->addTab(m_tableV, "v(x,y)");
        m_tableTabs->addTab(m_tableV2, "v₂(x,y)");
        m_tableTabs->addTab(m_tableVDiff, "|v - v₂|");
    }

    layout->addWidget(m_tableTabs);
    return box;
}

QTableWidget *DirichletDisplayWidget::createEmptyTable()
{
    auto *table = new QTableWidget;
    table->setStyleSheet("QTableWidget { border: 1px solid #ccc; }");
    return table;
}
