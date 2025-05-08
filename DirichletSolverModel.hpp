#pragma once
#include <QObject>
#include <QVector>

QT_FORWARD_DECLARE_CLASS(QString)

class DirichletSolverModel : public QObject
{
    Q_OBJECT

public:
    explicit DirichletSolverModel(QObject *parent = nullptr);
    void setup(int n, int m, double omega, double eps, int maxIter);
    void solveTestProblem();
    void solveMainProblem();

    double maxError() const;

    const QVector<QVector<double>> &solution() const;
    const QVector<QVector<double>> &exactSolution() const;

    QVector<QVector<double>> error() const;

private:
    int m_n;
    int m_m;
    int m_maxIter;

    double m_a;
    double m_b;
    double m_c;
    double m_d;
    double m_h;
    double m_k;
    double m_omega;
    double m_eps;

    QVector<QVector<double>> m_u;
    QVector<QVector<double>> m_uExact;


    void initialize();

    double f(double x, double y) const;
    double uStar(double x, double y) const;

    double mu1(double y) const;
    double mu2(double y) const;
    double mu3(double x) const;
    double mu4(double x) const;
};
