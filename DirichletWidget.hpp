#pragma once
#include <QWidget>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QVBoxLayout>

QT_FORWARD_DECLARE_CLASS(DirichletSolverModel)

using namespace QtDataVisualization;

class DirichletWidget : public QWidget
{
    Q_OBJECT

public:

    DirichletWidget(const QVector<QVector<double>> *data, double a, double b, double c, double d, QWidget *parent = nullptr);
    DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent = nullptr);

private:
    DirichletSolverModel *m_model = nullptr;
    const QVector<QVector<double>> *m_dataOnly = nullptr;
    bool m_isTest = false;
    Q3DSurface *m_surface;
    QWidget *m_container;
    QVBoxLayout *m_layout;
    double m_a = 0, m_b = 1, m_c = 0, m_d = 1;

    void setupChart();
    QSurface3DSeries *createSeries(const QVector<QVector<double>> &data, const QString &name);
};
