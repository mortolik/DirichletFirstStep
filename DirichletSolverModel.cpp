#include "DirichletSolverModel.hpp"

#include <QtMath>
#include <QTextStream>

DirichletSolverModel::DirichletSolverModel(QObject *parent)
    : QObject(parent),
    m_n(0),
    m_m(0),
    m_maxIter(10000),
    m_a(1.0),
    m_b(2.0),
    m_c(2.0),
    m_d(3.0),
    m_omega(1.0),
    m_eps(1e-6)
{
}

void DirichletSolverModel::setup(int n, int m, double omega, double eps, int maxIter)
{
    m_n = n;
    m_m = m;
    m_omega = omega;
    m_eps = eps;
    m_maxIter = maxIter;

    m_h = (m_b - m_a) / m_n;
    m_k = (m_d - m_c) / m_m;

    m_u = QVector<QVector<double>>(m_n + 1, QVector<double>(m_m + 1, 0.0));
    m_uExact = QVector<QVector<double>>(m_n + 1, QVector<double>(m_m + 1, 0.0));

    initialize();
}

double DirichletSolverModel::f(double x, double y) const
{
    return -qExp(-x * y * y);
}

double DirichletSolverModel::uStar(double x, double y) const
{
    return qSin(M_PI * x * y);
}

double DirichletSolverModel::mu1(double y) const
{
    return (y - 2.0) * (y - 3.0);
}

double DirichletSolverModel::mu2(double y) const
{
    return y * (y - 2.0) * (y - 3.0);
}

double DirichletSolverModel::mu3(double x) const
{
    return (x - 1.0) * (x - 2.0);
}

double DirichletSolverModel::mu4(double x) const
{
    return x * (x - 1.0) * (x - 2.0);
}

void DirichletSolverModel::initialize()
{
    for (int i = 0; i <= m_n; ++i)
    {
        double x = m_a + i * m_h;
        m_u[i][0] = mu3(x);
        m_u[i][m_m] = mu4(x);
        m_uExact[i][0] = uStar(x, m_c);
        m_uExact[i][m_m] = uStar(x, m_d);
    }

    for (int j = 0; j <= m_m; ++j)
    {
        double y = m_c + j * m_k;
        m_u[0][j] = mu1(y);
        m_u[m_n][j] = mu2(y);
        m_uExact[0][j] = uStar(m_a, y);
        m_uExact[m_n][j] = uStar(m_b, y);
    }
}

void DirichletSolverModel::solveTestProblem()
{
    int iter = 0;
    bool stop = false;

    while (!stop && iter < m_maxIter)
    {
        stop = true;

        for (int i = 1; i < m_n; ++i)
        {
            for (int j = 1; j < m_m; ++j)
            {
                double x = m_a + i * m_h;
                double y = m_c + j * m_k;

                double rhs = (m_u[i + 1][j] + m_u[i - 1][j]) / (m_h * m_h)
                             + (m_u[i][j + 1] + m_u[i][j - 1]) / (m_k * m_k)
                             - f(x, y);

                double denom = 2.0 * (1.0 / (m_h * m_h) + 1.0 / (m_k * m_k));
                double uNew = (1.0 - m_omega) * m_u[i][j] + m_omega * rhs / denom;

                if (qAbs(uNew - m_u[i][j]) > m_eps)
                {
                    stop = false;
                }

                m_u[i][j] = uNew;
                m_uExact[i][j] = uStar(x, y);
            }
        }

        ++iter;
    }
}

void DirichletSolverModel::solveMainProblem()
{
    int iter = 0;
    bool stop = false;

    while (!stop && iter < m_maxIter)
    {
        stop = true;

        for (int i = 1; i < m_n; ++i)
        {
            for (int j = 1; j < m_m; ++j)
            {
                double x = m_a + i * m_h;
                double y = m_c + j * m_k;

                double rhs = (m_u[i + 1][j] + m_u[i - 1][j]) / (m_h * m_h)
                             + (m_u[i][j + 1] + m_u[i][j - 1]) / (m_k * m_k)
                             - f(x, y);

                double denom = 2.0 * (1.0 / (m_h * m_h) + 1.0 / (m_k * m_k));
                double uNew = (1.0 - m_omega) * m_u[i][j] + m_omega * rhs / denom;

                if (qAbs(uNew - m_u[i][j]) > m_eps)
                {
                    stop = false;
                }

                m_u[i][j] = uNew;
            }
        }

        ++iter;
    }
}

QVector<QVector<double>> DirichletSolverModel::error() const
{
    QVector<QVector<double>> diff(m_n + 1, QVector<double>(m_m + 1));
    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            diff[i][j] = qAbs(m_uExact[i][j] - m_u[i][j]);
        }
    }
    return diff;
}

double DirichletSolverModel::maxError() const
{
    double maxErr = 0.0;

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            maxErr = qMax(maxErr, qAbs(m_u[i][j] - m_uExact[i][j]));
        }
    }

    return maxErr;
}

const QVector<QVector<double>> &DirichletSolverModel::solution() const
{
    return m_u;
}

const QVector<QVector<double>> &DirichletSolverModel::exactSolution() const
{
    return m_uExact;
}

QVector<QVector<double>> DirichletSolverModel::compareWithFinerGrid(int finerN, int finerM, double &eps2Out) const
{
    DirichletSolverModel fineModel;
    fineModel.setup(finerN, finerM, m_omega, m_eps, m_maxIter);
    fineModel.solveMainProblem(); // основная задача без точного решения

    double hRatio = double(finerN) / m_n;
    double kRatio = double(finerM) / m_m;

    QVector<QVector<double>> fineU = fineModel.solution();
    QVector<QVector<double>> coarseVsFine(m_n + 1, QVector<double>(m_m + 1, 0.0));

    eps2Out = 0.0;

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            int fi = i * hRatio;
            int fj = j * kRatio;

            double diff = qAbs(m_u[i][j] - fineU[fi][fj]);
            coarseVsFine[i][j] = diff;
            eps2Out = qMax(eps2Out, diff);
        }
    }

    return coarseVsFine;
}

DirichletSolverModel::ReportData DirichletSolverModel::generateReportData(bool isTestTask, double extraError) const
{
    ReportData data;
    data.n = m_n;
    data.m = m_m;
    data.omega = m_omega;
    data.eps = m_eps;
    data.maxIter = m_maxIter;
    data.isTest = isTestTask;

    if (isTestTask)
        data.maxError = maxError();

    if (!isTestTask)
        data.accuracy = extraError;

    return data;
}

QString DirichletSolverModel::reportString(bool isTestTask, double extraError) const
{
    ReportData data = generateReportData(isTestTask, extraError);

    QStringList lines;
    lines << (data.isTest ? "Справка по тестовой задаче:" : "Справка по основной задаче:");
    lines << QString("Сетка: n = %1, m = %2").arg(data.n).arg(data.m);
    lines << QString("Метод верхней релаксации: ω = %1").arg(data.omega);
    lines << QString("Точность метода: εмет = %1, максимум итераций: %2").arg(data.eps).arg(data.maxIter);
    lines << "Начальное приближение: интерполяция граничных условий";

    if (data.isTest && data.maxError >= 0)
    {
        lines << QString("Максимальная погрешность ε₁ = %1").arg(data.maxError, 0, 'e', 3);
        auto [xmax, ymax] = maxErrorPoint();
        lines << QString("Максимальное отклонение в точке: x = %1, y = %2").arg(xmax).arg(ymax);
    }

    if (!data.isTest && data.accuracy >= 0)
    {
        lines << QString("Оценка точности ε₂ = %1").arg(data.accuracy, 0, 'e', 3);
        lines << "Оценка: ε₂ = max |v(x, y) − v₂(x, y)| между сетками (n, m) и (2n, 2m)";
    }

    return lines.join("\n");
}

QPair<double, double> DirichletSolverModel::maxErrorPoint() const
{
    double maxErr = 0.0;
    int maxI = 0, maxJ = 0;

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            double err = qAbs(m_uExact[i][j] - m_u[i][j]);
            if (err > maxErr)
            {
                maxErr = err;
                maxI = i;
                maxJ = j;
            }
        }
    }

    double x = m_a + maxI * m_h;
    double y = m_c + maxJ * m_k;
    return qMakePair(x, y);
}


