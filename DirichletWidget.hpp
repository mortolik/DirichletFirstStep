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
    explicit DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent = nullptr);

private:
    DirichletSolverModel *m_model;
    Q3DSurface *m_surface;
    QWidget *m_container;
    QVBoxLayout *m_layout;

    bool m_isTest;

    void setupChart();
    QSurface3DSeries *createSeries(const QVector<QVector<double>> &data, const QString &name);
};
