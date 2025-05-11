#include "DirichletWidget.hpp"
#include "DirichletSolverModel.hpp"
#include <QLabel>
#include <QSurface3DSeries>
#include <QSurfaceDataItem>
#include <QDebug>

DirichletWidget::DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent)
    : QWidget(parent),
    m_model(model),
    m_isTest(isTest),
    m_surface(new Q3DSurface),
    m_container(QWidget::createWindowContainer(m_surface, this)),
    m_mainLayout(new QHBoxLayout(this)),
    m_chartLayout(new QVBoxLayout),
    m_reportLabel(new QLabel)
{
    m_chartLayout->addWidget(m_container);
    m_mainLayout->addLayout(m_chartLayout);

    m_reportLabel->setWordWrap(true);
    m_reportLabel->setFixedWidth(250);
    m_mainLayout->addWidget(m_reportLabel, 0, Qt::AlignTop);

    setupChart();

    if (m_model)
        m_reportLabel->setText(m_model->reportString(m_isTest));
}

DirichletWidget::DirichletWidget(const QVector<QVector<double>> *data, double a, double b, double c, double d, QWidget *parent, const QString &reportText)
    : QWidget(parent),
    m_dataOnly(data),
    m_surface(new Q3DSurface), m_container(QWidget::createWindowContainer(m_surface, this)), m_mainLayout(new QHBoxLayout(this)), m_chartLayout(new QVBoxLayout),
    m_reportLabel(new QLabel),
    m_a(a),
    m_b(b),
    m_c(c),
    m_d(d),
    m_reportText(reportText)
{
    m_chartLayout->addWidget(m_container);
    m_mainLayout->addLayout(m_chartLayout);

    m_reportLabel->setText(m_reportText);
    m_reportLabel->setWordWrap(true);
    m_reportLabel->setMinimumWidth(250);
    m_mainLayout->addWidget(m_reportLabel);

    setupChart();
}

void DirichletWidget::setupChart()
{
    m_surface->setAxisX(new QValue3DAxis);
    m_surface->setAxisY(new QValue3DAxis);
    m_surface->setAxisZ(new QValue3DAxis);

    m_surface->axisX()->setTitle("x");
    m_surface->axisZ()->setTitle("y");
    m_surface->axisY()->setTitle("u");

    m_surface->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
    m_surface->activeTheme()->setType(Q3DTheme::ThemeQt);

    if (m_dataOnly)
    {
        auto *series = createSeries(*m_dataOnly, "Ошибка |u* - v|");
        series->setBaseColor(Qt::blue);
        m_surface->addSeries(series);
        return;
    }

    if (m_isTest)
    {
        m_surface->addSeries(createSeries(m_model->exactSolution(), "u*(x, y)"));
        m_surface->addSeries(createSeries(m_model->solution(), "v(x, y)"));
    }
    else
    {
        m_surface->addSeries(createSeries(m_model->solution(), "Численное решение"));
    }
}


QSurface3DSeries *DirichletWidget::createSeries(const QVector<QVector<double>> &data, const QString &name)
{
    int rows = data.size();
    int cols = data[0].size();

    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(rows);

    double a = 1.0;
    double b = 2.0;
    double c = 2.0;
    double d = 3.0;

    double h = (b - a) / (rows - 1);
    double k = (d - c) / (cols - 1);

    for (int i = 0; i < rows; ++i)
    {
        QSurfaceDataRow *row = new QSurfaceDataRow(cols);
        double x = a + i * h;

        for (int j = 0; j < cols; ++j)
        {
            double y = c + j * k;
            (*row)[j].setPosition(QVector3D(x, data[i][j], y));
        }

        dataArray->append(row);
    }

    QSurfaceDataProxy *proxy = new QSurfaceDataProxy;
    proxy->resetArray(dataArray);

    QSurface3DSeries *series = new QSurface3DSeries(proxy);
    series->setItemLabelFormat(name + ": @xLabel, @zLabel = @yLabel");
    series->setDrawMode(QSurface3DSeries::DrawSurface);
    series->setName(name);
    return series;
}
