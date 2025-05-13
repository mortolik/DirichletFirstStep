#pragma once
#include "DirichletSolverModel.hpp"
#include "qtextedit.h"
#include <QWidget>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtDataVisualization/QSurfaceDataProxy>

QT_FORWARD_DECLARE_CLASS(QLabel)
QT_FORWARD_DECLARE_CLASS(QSpinBox)
QT_FORWARD_DECLARE_CLASS(QTextEdit)
QT_FORWARD_DECLARE_CLASS(QLineEdit)
QT_FORWARD_DECLARE_CLASS(QGroupBox)
QT_FORWARD_DECLARE_CLASS(QVBoxLayout)
QT_FORWARD_DECLARE_CLASS(QHBoxLayout)
QT_FORWARD_DECLARE_CLASS(QPushButton)
QT_FORWARD_DECLARE_CLASS(QDoubleSpinBox)
QT_FORWARD_DECLARE_CLASS(DirichletSolverModel)

using namespace QtDataVisualization;

class DirichletWidget : public QWidget
{
    Q_OBJECT

public:

    DirichletWidget(DirichletSolverModel *model, bool isTest, QWidget *parent = nullptr);

private slots:
    void onSolveButtonClicked();

    void setReportText(const QString &text);
signals:
    void solutionUpdated();

private:
    DirichletSolverModel *m_model = nullptr;
    const QVector<QVector<double>> *m_dataOnly = nullptr;
    bool m_isTest = false;
    Q3DSurface *m_surface;
    QWidget *m_container;
    QHBoxLayout *m_mainLayout;
    QVBoxLayout *m_chartLayout;
    double m_a = 0, m_b = 1, m_c = 0, m_d = 1;
    QString m_reportText;

    QSpinBox* m_nParam {nullptr};
    QSpinBox* m_mParam {nullptr};
    QDoubleSpinBox* m_epsParam {nullptr};
    QSpinBox* m_stepsParam {nullptr};
    QDoubleSpinBox* m_omegaParam {nullptr};

    QPushButton* m_solveBtn {nullptr};

    QTextEdit* m_reportEdit;

    void setupChart();
    QSurface3DSeries *createSeries(const QVector<QVector<double>> &data, const QString &name);

    QGroupBox* createSettingsGroup();
    QGroupBox* createLeftLayout();
    void updateChart();
    QGroupBox *createReportBox();
    QGroupBox* createChartBox();
};
