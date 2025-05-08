#include "DirichletWidget.hpp"
#include "DirichletSolverModel.hpp"
#include <QSurface3DSeries>
#include <QSurfaceDataItem>
#include <QDebug>

DirichletWidget::DirichletWidget(DirichletSolverModel *model, QWidget *parent)
    : QWidget(parent),
    m_model(model),
    m_surface(new Q3DSurface),
    m_container(QWidget::createWindowContainer(m_surface, this)),
    m_layout(new QVBoxLayout(this))
{
    m_layout->addWidget(m_container);
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

    m_surface->axisX()->setLabelFormat("%.2f");
    m_surface->axisY()->setLabelFormat("%.2f");
    m_surface->axisZ()->setLabelFormat("%.2f");

    m_surface->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
    m_surface->activeTheme()->setType(Q3DTheme::ThemeQt);

    auto *series = createSeries(m_model->solution(), "Численное решение");
    m_surface->addSeries(series);
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
