#include "DirichletWidget.hpp"
#include "DirichletSolverModel.hpp"

#include <QLabel>
#include <QDebug>
#include <QSpinBox>
#include <QGroupBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QTextEdit>
#include <QBoxLayout>
#include <QPushButton>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QSurface3DSeries>
#include <QSurfaceDataItem>

DirichletWidget::DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent)
    : QWidget(parent),
    m_model(model),
    m_isTest(isTest),
    m_surface(new Q3DSurface),
    m_container(nullptr),
    m_mainLayout(new QHBoxLayout(this)),
    m_chartLayout(new QVBoxLayout)
{
    m_mainLayout->setContentsMargins(0, 0, 0, 0);
    QVBoxLayout* chartAndReportLayot = new QVBoxLayout();
    chartAndReportLayot->setContentsMargins(0, 0, 0, 0);
    chartAndReportLayot->addWidget(createReportBox(), 0, Qt::AlignTop);
    chartAndReportLayot->addWidget(createChartBox());
    chartAndReportLayot->setStretch(0, 0);
    chartAndReportLayot->setStretch(1, 1);
    m_mainLayout->addWidget(createLeftLayout());
    m_mainLayout->addLayout(chartAndReportLayot);

    setupChart();

    connect(m_solveBtn, &QPushButton::clicked, this, &DirichletWidget::onSolveButtonClicked);
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
    if(!m_model->solution().isEmpty())
    {
        if (m_isTest && !m_model->exactSolution().isEmpty())
        {
            m_surface->addSeries(createSeries(m_model->exactSolution(), "u*(x, y)"));
            m_surface->addSeries(createSeries(m_model->solution(), "v(x, y)"));
        }
        else
        {
            m_surface->addSeries(createSeries(m_model->solution(), "Численное решение"));
        }
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

    double h = (b - a) / (cols - 1);
    double k = (d - c) / (rows - 1);

    for (int i = 0; i < rows; ++i)
    {
        QSurfaceDataRow *row = new QSurfaceDataRow;
        double z = c + i * k;

        for (int j = 0; j < cols; ++j)
        {
            double x = a + j * h;
            double y = data[i][j];
            row->append(QSurfaceDataItem(QVector3D(x, y, z)));
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
    box->setFixedSize(190, 150);

    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setSpacing(5);
    layout->setContentsMargins(5, 5, 5, 5);

    m_nParam = new QSpinBox();
    m_mParam = new QSpinBox();
    m_epsParam = new QDoubleSpinBox();
    m_stepsParam = new QSpinBox();
    m_omegaParam = new QDoubleSpinBox();

    m_nParam->setRange(1, 1000);
    m_nParam->setValue(100);

    m_mParam->setRange(1, 1000);
    m_mParam->setValue(100);

    m_stepsParam->setRange(0, 1000000);
    m_stepsParam->setValue(500);

    m_epsParam->setDecimals(8);
    m_epsParam->setValue(0.00000001);
    m_epsParam->setSingleStep(0.001);
    m_epsParam->setRange(0.00000001, 1.0);

    m_optimalCheckBox = new QCheckBox();
    m_optimalCheckBox->setChecked(true);

    m_omegaParam->setRange(0.1, 2.0);
    m_omegaParam->setSingleStep(0.1);
    m_omegaParam->setValue(1.8);

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

    layout->addWidget(createLabeledWidget("Оптимальный ω ", m_optimalCheckBox, 80));
    layout->addWidget(createLabeledWidget("Omega ", m_omegaParam, 80));

    return box;
}

QGroupBox *DirichletWidget::createLeftLayout()
{
    QGroupBox* box = new QGroupBox(this);
    box->setMaximumWidth(200);
    QVBoxLayout* layout = new QVBoxLayout(box);

    layout->addWidget(createSettingsGroup(), 0, Qt::AlignTop);

    m_solveBtn = new QPushButton("Запустить");
    layout->addWidget(m_solveBtn, 0, Qt::AlignTop);
    layout->addStretch();

    return box;
}

void DirichletWidget::onSolveButtonClicked()
{
    int n = m_nParam->value();
    int m = m_mParam->value();
    double eps = m_epsParam->value();
    int maxIter = m_stepsParam->value();

    double omega = m_omegaParam->value();
    m_model->setup(n, m, omega, eps, maxIter);
    if (m_optimalCheckBox->isChecked())
    {
        double omegaOpt = m_model->computeOptimalOmega();
        m_model->setup(n, m, omegaOpt, eps, maxIter);
    }

    if (m_isTest)
        m_model->solveTestProblem();
    else
        m_model->solveMainProblem();

    updateChart();
    emit solutionUpdated();

    QString report = m_model->reportString(m_isTest);
    setReportText(report);
}

void DirichletWidget::updateChart()
{
    auto seriesList = m_surface->seriesList();
    for (auto series : seriesList)
        m_surface->removeSeries(series);
\
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

QGroupBox* DirichletWidget::createReportBox()
{
    QGroupBox* box = new QGroupBox();
    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setContentsMargins(5, 5, 5, 5);

    m_reportEdit = new QTextEdit;
    m_reportEdit->setStyleSheet(
        "QTextEdit {"
        "    border-style: solid;"
        "    border-width: 1px;"
        "    border-color: #dcdcdc;"
        "}"
        );
    m_reportEdit->setReadOnly(true);
    m_reportEdit->setFont(QFont("Courier"));
    
    m_reportEdit->setHtml(R"(
    <font color="#888">Здесь будет справка...</font><br>
    <font color="#888">Нажмите "Запустить", чтобы начать расчёт.</font>
    )");

    layout->addWidget(m_reportEdit);

    return box;
}

void DirichletWidget::setReportText(const QString& text)
{
    m_reportEdit->setPlainText(text);
}

QGroupBox* DirichletWidget::createChartBox()
{
    QGroupBox* box = new QGroupBox();
    QVBoxLayout* layout = new QVBoxLayout(box);
    layout->setContentsMargins(5, 5, 5, 5);

    m_container = QWidget::createWindowContainer(m_surface, this);
    layout->addWidget(m_container);

    return box;
}
