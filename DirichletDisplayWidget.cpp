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
    createTable();

    connect(m_chart, &DirichletWidget::solutionUpdated, this, &DirichletDisplayWidget::onSolutionUpdated);
    updateTable();

}

QTableWidget *DirichletDisplayWidget::tableWidget()
{
    return m_table;
}

void DirichletDisplayWidget::fillTable(const QVector<QVector<double>> &data)
{
    int cols = data.size();
    int rows = data[0].size();

    m_table->setRowCount(rows);
    m_table->setColumnCount(cols + 1);

    m_table->setHorizontalHeaderItem(0, new QTableWidgetItem("Y \\ X"));
    for (int i = 0; i < cols; ++i)
    {
        double x = m_model->a() + i * m_model->h();
        m_table->setHorizontalHeaderItem(i + 1, new QTableWidgetItem(QString::number(x, 'f', 2)));
    }

    for (int j = 0; j < rows; ++j)
    {
        double y = m_model->c() + j * m_model->k();
        m_table->setItem(j, 0, new QTableWidgetItem(QString::number(y, 'f', 2)));

        for (int i = 0; i < cols; ++i)
        {
            QString val = QString::number(data[i][j], 'f', 4);
            auto *item = new QTableWidgetItem(val);
            item->setTextAlignment(Qt::AlignCenter);
            m_table->setItem(j, i + 1, item);
        }
    }

    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->resizeColumnsToContents();
    m_table->resizeRowsToContents();
}

void DirichletDisplayWidget::updateTable()
{
    QVector<QVector<double>> data;

    if (m_isTest)
        data = m_model->exactSolution(); // Можно использовать solution() или exactSolution()
    else
        data = m_model->solution();

    if (data.isEmpty()) {
        m_table->clear();
        m_table->setRowCount(0);
        m_table->setColumnCount(0);
        m_table->setHorizontalHeaderItem(0, new QTableWidgetItem("Нет данных"));
        m_table->setVerticalHeaderItem(0, new QTableWidgetItem("Нажмите 'Запустить'"));
        return;
    }

    fillTable(data);
}

void DirichletDisplayWidget::onSolutionUpdated()
{
    updateTable();
}

QGroupBox *DirichletDisplayWidget::createTable()
{
    QGroupBox* box = new QGroupBox();
    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setContentsMargins(5, 5, 5, 5);

    m_table = new QTableWidget;
    m_table->setStyleSheet(
        "QTableWidget {"
        "    border-style: solid;"
        "    border-width: 1px;"
        "    border-color: #dcdcdc;"
        "}"
        );
    layout->addWidget(m_table);

    return box;
}
