#pragma once
#include "qlabel.h"
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

    DirichletWidget(const QVector<QVector<double>> *data, double a, double b, double c, double d, QWidget *parent = nullptr, const QString &reportText = nullptr);
    DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent = nullptr);

private:
    DirichletSolverModel *m_model = nullptr;
    const QVector<QVector<double>> *m_dataOnly = nullptr;
    bool m_isTest = false;
    Q3DSurface *m_surface;
    QWidget *m_container;
    QHBoxLayout *m_mainLayout;
    QVBoxLayout *m_chartLayout;
    QLabel *m_reportLabel;
    double m_a = 0, m_b = 1, m_c = 0, m_d = 1;
    QString m_reportText;

    void setupChart();
    QSurface3DSeries *createSeries(const QVector<QVector<double>> &data, const QString &name);
};
