#pragma once
#include <QObject>
#include <QVector>

QT_FORWARD_DECLARE_CLASS(QString)

class DirichletSolverModel : public QObject
{
    Q_OBJECT

public:
    explicit DirichletSolverModel(QObject *parent = nullptr);
    void setup(int n, int m, double eps, int maxIter);
    void solveTestProblem();
    void solveMainProblem();

    double maxError() const;

    const QVector<QVector<double>> &solution() const;
    const QVector<QVector<double>> &exactSolution() const;

    QVector<QVector<double>> error() const;
    QVector<QVector<double>> compareWithFinerGrid(int finerN, int finerM, double &eps2Out) const;

    struct ReportData
    {
        int n = 0;
        int m = 0;
        double omega = 1.0;
        double eps = 1e-6;
        int maxIter = 10000;
        double maxError = -1.0;
        double accuracy = -1.0; // ε2
        bool isTest = true;
    };

    QString reportString(bool isTestTask = true) const;
    ReportData generateReportData(bool isTestTask = true) const;
    double a() const { return m_a; }
    double b() const { return m_b; }
    double c() const { return m_c; }
    double d() const { return m_d; }
    double h() const { return m_h; }
    double k() const { return m_k; }

    void setOmega(double omega) { m_omega = omega; }

    double computeOptimalOmega() const;

    struct FinerGridResult
    {
        QVector<QVector<double>> v;
        QVector<QVector<double>> v2;
        QVector<QVector<double>> diff;
        int iterations = 0;
        QPair<double,double> maxPt{0,0};
    };

    DirichletSolverModel::FinerGridResult computeFinerGridComparison() const;

    double lastResidual() const;
    int lastIter() const;
    const FinerGridResult& lastFinerGridResult() const { return m_result; }
    bool hasFinerCached() const { return m_resultCached; }


private:
    int m_n;
    int m_m;
    int m_maxIter;
    int m_lastIter;
    mutable int m_last2Iter;
    double m_lastResidual;
    mutable bool m_skipInitialization = false;

    double m_a;
    double m_b;
    double m_c;
    double m_d;
    double m_h;
    double m_k;
    double m_omega;
    double m_eps;

    mutable FinerGridResult m_result;
    mutable bool m_resultCached = false;

    QVector<QVector<double>> m_u;
    QVector<QVector<double>> m_uExact;


    double m_invH2;    // добавить
    double m_invK2;    // добавить
    double m_denom;
    std::vector<double> m_uPrevFlat;  // добавить

    // Удобный доступ к плоскому индексу
    inline int flatIdx(int i, int j) const { return i*(m_m+1) + j; }

    void initialize();

    double f(double x, double y) const;
    double uStar(double x, double y) const;

    double mu1(double y) const;
    double mu2(double y) const;
    double mu3(double x) const;
    double mu4(double x) const;
    QPair<double, double> maxErrorPoint() const;

    double fTest(double x, double y) const;

    QPair<double, double> maxErrorPointCompare() const;

    double computeInitialResidual() const;
    void initializeInterior();
    void applyInterpolatedInitialGuess(const QVector<QVector<double> > &coarseU, int coarseN, int coarseM);
};
