#include "DirichletWidget.hpp"
#include "DirichletSolverModel.hpp"

#include <QLabel>
#include <QDebug>
#include <QSpinBox>
#include <QGroupBox>
#include <QBoxLayout>
#include <QVBoxLayout>
#include <QSurface3DSeries>
#include <QSurfaceDataItem>

DirichletWidget::DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent)
    : QWidget(parent),
    m_model(model),
    m_isTest(isTest),
    m_surface(new Q3DSurface),
    m_container(QWidget::createWindowContainer(m_surface, this)),
    m_mainLayout(new QVBoxLayout(this)),
    m_chartLayout(new QVBoxLayout),
    m_reportLabel(new QLabel)
{
    m_chartLayout->addWidget(m_container);

    m_reportLabel->setWordWrap(true);
    m_reportLabel->setFixedWidth(250);
    m_mainLayout->addWidget(createSettingsGroup());
    m_mainLayout->addWidget(m_reportLabel, 0, Qt::AlignTop);

    m_mainLayout->addLayout(m_chartLayout);

    setupChart();

    if (m_model)
        m_reportLabel->setText(m_model->reportString(m_isTest));
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

QGroupBox* DirichletWidget::createSettingsGroup()
{
    QGroupBox* box = new QGroupBox("Параметры", this);
    box->setMaximumWidth(175);

    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setSpacing(5);
    layout->setContentsMargins(5, 5, 5, 5);

    m_nParam = new QSpinBox();
    m_mParam = new QSpinBox();
    m_epsParam = new QDoubleSpinBox();
    m_stepsParam = new QSpinBox();
    m_omegaParam = new QDoubleSpinBox();

    m_nParam->setRange(1, 1000);
    m_nParam->setValue(5);

    m_mParam->setRange(1, 1000);
    m_mParam->setValue(5);

    m_epsParam->setDecimals(7);
    m_epsParam->setValue(0.0000005);
    m_epsParam->setSingleStep(0.001);
    m_epsParam->setRange(0.0000001, 1.0);

    m_omegaParam->setRange(0.1, 2.0);
    m_omegaParam->setSingleStep(0.1);
    m_omegaParam->setValue(1.5);

    auto createLabeledWidget = [](const QString& labelText, QWidget* widget, int width) -> QWidget* {
        QLabel* label = new QLabel(labelText);
        label->setFixedWidth(width);
        QWidget* container = new QWidget();
        QHBoxLayout* hLayout = new QHBoxLayout(container);
        hLayout->addWidget(label, 0, Qt::AlignLeft);
        hLayout->addWidget(widget, 1);
        hLayout->setSpacing(5);
        hLayout->setContentsMargins(0, 0, 0, 0);
        return container;
    };

    QWidget* nmContainer = new QWidget();
    QHBoxLayout* nmLayout = new QHBoxLayout(nmContainer);
    nmLayout->addWidget(createLabeledWidget("n =", m_nParam, 20));
    nmLayout->addWidget(createLabeledWidget("m =", m_mParam, 20));
    nmLayout->setSpacing(10);
    nmLayout->setContentsMargins(0, 0, 0, 0);

    layout->addWidget(nmContainer);
    layout->addWidget(createLabeledWidget("Точность ", m_epsParam, 80));
    layout->addWidget(createLabeledWidget("Кол-во шагов ", m_stepsParam, 80));
    layout->addWidget(createLabeledWidget("Omega ", m_omegaParam, 80));

    return box;
}
