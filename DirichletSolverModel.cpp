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
    m_u = QVector<QVector<double>>();
    m_uExact = QVector<QVector<double>>();
}

void DirichletSolverModel::setup(int n, int m, double eps, int maxIter)
{
    m_n = n;
    m_m = m;
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

double DirichletSolverModel::fTest(double x, double y) const
{
    return M_PI*M_PI*(x*x + y*y) * qSin(M_PI * x * y);
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
    // граничные условия — именно так, как было изначально
    for (int i = 0; i <= m_n; ++i)
    {
        double x = m_a + i * m_h;
        m_u[i][0]       = mu3(x);
        m_u[i][m_m]     = mu4(x);
        m_uExact[i][0]  = uStar(x, m_c);
        m_uExact[i][m_m]= uStar(x, m_d);
    }
    for (int j = 0; j <= m_m; ++j)
    {
        double y = m_c + j * m_k;
        m_u[0][j]       = mu1(y);
        m_u[m_n][j]     = mu2(y);
        m_uExact[0][j]  = uStar(m_a, y);
        m_uExact[m_n][j]= uStar(m_b, y);
    }
}

void DirichletSolverModel::solveTestProblem()
{
    double maxResidual = 0.0;
    for (int i = 0; i <= m_n; ++i)
    {
        double x = m_a + i * m_h;
        m_u[i][0]      = uStar(x, m_c);
        m_u[i][m_m]    = uStar(x, m_d);
    }
    for (int j = 0; j <= m_m; ++j)
    {
        double y = m_c + j * m_k;
        m_u[0][j]      = uStar(m_a, y);
        m_u[m_n][j]    = uStar(m_b, y);
    }
    for (int i = 1; i < m_n; ++i)
        for (int j = 1; j < m_m; ++j)
            m_u[i][j] = 0.0;

    for (int i = 0; i <= m_n; ++i)
        for (int j = 0; j <= m_m; ++j)
        {
            double x = m_a + i*m_h;
            double y = m_c + j*m_k;
            m_uExact[i][j] = uStar(x, y);
        }

    double h2 = m_h*m_h, k2 = m_k*m_k;
    double denom = 2.0*(1.0/h2 + 1.0/k2);

    int iter = 0;
    bool stop = false;
    while (!stop && iter < m_maxIter)
    {
        double maxResidual = 0.0;
        stop = true;
        for (int i = 1; i < m_n; ++i)
        {
            for (int j = 1; j < m_m; ++j)
            {
                double x = m_a + i*m_h;
                double y = m_c + j*m_k;

                double fstar = -fTest(x, y);

                double rhs = (m_u[i+1][j] + m_u[i-1][j])/h2
                             + (m_u[i][j+1] + m_u[i][j-1])/k2
                             - fstar;

                double uNew = (1.0 - m_omega)*m_u[i][j]
                              + m_omega*(rhs/denom);

                double diff = qAbs(uNew - m_u[i][j]);
                if (diff > m_eps)
                {
                    stop = false;

                }
                maxResidual = qMax(maxResidual, diff);

                m_u[i][j] = uNew;
            }
        }
        ++iter;
    }
    m_lastIter = iter;
    m_lastResidual = maxResidual;
}

double DirichletSolverModel::computeInitialResidual() const
{
    double maxR = 0.0;
    for (int i = 1; i < m_n; ++i)
    {
        for (int j = 1; j < m_m; ++j)
        {
            double x = m_a + i * m_h;
            double y = m_c + j * m_k;

            double laplace =
                (m_u[i+1][j] - 2*m_u[i][j] + m_u[i-1][j]) / (m_h * m_h) +
                (m_u[i][j+1] - 2*m_u[i][j] + m_u[i][j-1]) / (m_k * m_k);

            double R = laplace - f(x, y);
            maxR = qMax(maxR, qAbs(R));
        }
    }
    return maxR;
}

void DirichletSolverModel::solveMainProblem()
{
    initializeInterior();

    int iter = 0;
    bool stop = false;
    double maxResidual = 0.0;

    while (!stop && iter < m_maxIter)
    {
        stop = true;
        maxResidual = 0.0;

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

                double diff = qAbs(uNew - m_u[i][j]);
                if (diff > m_eps)
                {
                    stop = false;

                }
                maxResidual = qMax(maxResidual, diff);

                m_u[i][j] = uNew;
            }
        }

        ++iter;
    }
    m_lastIter = iter;
    m_lastResidual = maxResidual;
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

    if (!m_u.isEmpty() && !m_uExact.isEmpty())
    {
        for (int i = 0; i <= m_n; ++i)
        {
            for (int j = 0; j <= m_m; ++j)
            {
                maxErr = qMax(maxErr, qAbs(m_u[i][j] - m_uExact[i][j]));
            }
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

void DirichletSolverModel::initializeInterior()
{
    for (int i = 1; i < m_n; ++i)
    {
        double x = m_a + i * m_h;
        for (int j = 1; j < m_m; ++j)
        {
            double y = m_c + j * m_k;

            double left = mu1(y);
            double right = mu2(y);
            double bottom = mu3(x);
            double top = mu4(x);

            double tx = (x - m_a) / (m_b - m_a);
            double ty = (y - m_c) / (m_d - m_c);

            m_u[i][j] = (1 - tx) * left + tx * right + (1 - ty) * bottom + ty * top;
            m_u[i][j] /= 2.0;
        }
    }
}

QVector<QVector<double>> DirichletSolverModel::compareWithFinerGrid(int finerN, int finerM, double &eps2Out) const
{
    DirichletSolverModel fineModel;
    fineModel.setup(finerN, finerM, m_eps, m_maxIter);
    fineModel.solveMainProblem();

    QVector<QVector<double>> fineU = fineModel.solution();
    double fineH = (m_b - m_a) / finerN;
    double fineK = (m_d - m_c) / finerM;

    QVector<QVector<double>> coarseVsFine(m_n + 1, QVector<double>(m_m + 1, 0.0));
    eps2Out = 0.0;

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            double x = m_a + i * m_h;
            double y = m_c + j * m_k;

            int fi = static_cast<int>((x - m_a) / fineH);
            int fj = static_cast<int>((y - m_c) / fineK);

            if (fi >= finerN || fj >= finerM) continue;

            double dx = (x - (m_a + fi * fineH)) / fineH;
            double dy = (y - (m_c + fj * fineK)) / fineK;

            double v00 = fineU[fi][fj];
            double v10 = fineU[fi + 1][fj];
            double v01 = fineU[fi][fj + 1];
            double v11 = fineU[fi + 1][fj + 1];

            double interpVal = (1 - dx) * (1 - dy) * v00
                               + dx * (1 - dy) * v10
                               + (1 - dx) * dy * v01
                               + dx * dy * v11;

            double diff = qAbs(m_u[i][j] - interpVal);
            coarseVsFine[i][j] = diff;
            eps2Out = qMax(eps2Out, diff);
        }
    }

    return coarseVsFine;
}

DirichletSolverModel::ReportData DirichletSolverModel::generateReportData(bool isTestTask) const
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
    {
        double eps2 = 0.0;
        compareWithFinerGrid(m_n * 2, m_m * 2, eps2);
        data.accuracy = eps2;
    }

    return data;
}

QString DirichletSolverModel::reportString(bool isTestTask) const
{
    ReportData data = generateReportData(isTestTask);

    QStringList lines;
    lines << (data.isTest ? "Справка по тестовой задаче:" : "Справка по основной задаче:");
    lines << QString("Сетка: n = %1, m = %2").arg(data.n).arg(data.m);
    lines << QString("Метод верхней релаксации: ω = %1").arg(data.omega);
    lines << QString("Количество итераций: %1").arg(m_lastIter);
    lines << QString("Точность метода: εмет = %1, максимум итераций: %2").arg(data.eps).arg(data.maxIter);
    lines << QString("Невязка СЛАУ на начальном приближении R(0): %1").arg(computeInitialResidual(), 0, 'e', 3);
    lines << QString("Достигнутая точность: %1").arg(static_cast<double>(lastResidual()), 0, 'e', 5);
    if (data.isTest && data.maxError >= 0)
    {
        lines << "Начальное приближение: нулевое";
        lines << QString("Максимальная погрешность ε₁ = %1").arg(data.maxError, 0, 'e', 3);
        auto [xmax, ymax] = maxErrorPoint();
        lines << QString("Максимальное отклонение в точке: x = %1, y = %2").arg(xmax).arg(ymax);
    }

    if (!data.isTest && data.accuracy >= 0)
    {
        lines << "Начальное приближение: среднее";
        lines << QString("Оценка точности ε₂ = max |v(x, y) − v₂(x, y)|  = %1").arg(data.accuracy, 0, 'e', 3);
        auto [xmax, ymax] = maxErrorPointCompare();
        lines << QString("Максимальное отклонение в точке: x = %1, y = %2").arg(xmax).arg(ymax);
    }

    return lines.join("\n");
}

double DirichletSolverModel::lastResidual() const
{
    return m_lastResidual;
}

QPair<double, double> DirichletSolverModel::maxErrorPointCompare() const
{
    DirichletSolverModel fineModel;
    fineModel.setup(m_n * 2, m_m * 2, m_eps, m_maxIter);
    fineModel.solveMainProblem();

    QVector<QVector<double>> fineU = fineModel.solution();

    double maxErr = 0.0;
    int maxI = 0, maxJ = 0;

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            double diff = qAbs(m_u[i][j] - fineU[i * 2][j * 2]);
            if (diff > maxErr)
            {
                maxErr = diff;
                maxI = i;
                maxJ = j;
            }
        }
    }

    double x = m_a + maxI * m_h;
    double y = m_c + maxJ * m_k;
    return qMakePair(x, y);
}

QPair<double, double> DirichletSolverModel::maxErrorPoint() const
{
    double maxErr = 0.0;
    int maxI = 0, maxJ = 0;

    if(!m_u.isEmpty() && !m_uExact.isEmpty())
    {
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
    }

    double x = m_a + maxI * m_h;
    double y = m_c + maxJ * m_k;
    return qMakePair(x, y);
}

double DirichletSolverModel::computeOptimalOmega()
{
    double rho = std::cos(M_PI / m_n) + std::cos(M_PI / m_m);
    double value = std::sqrt(1 - std::pow(rho / 2.0, 2));
    m_omega = 2.0 / (1.0 + value);
    return m_omega;
}

DirichletSolverModel::FinerGridResult DirichletSolverModel::computeFinerGridComparison() const
{
    DirichletSolverModel fineModel;
    fineModel.setup(m_n * 2, m_m * 2, m_eps, m_maxIter);
    fineModel.solveMainProblem();

    QVector<QVector<double>> v2Full = fineModel.solution();

    QVector<QVector<double>> v(m_n + 1, QVector<double>(m_m + 1));
    QVector<QVector<double>> diff(m_n + 1, QVector<double>(m_m + 1));

    for (int i = 0; i <= m_n; ++i)
    {
        for (int j = 0; j <= m_m; ++j)
        {
            v[i][j] = m_u[i][j];
            double v2value = v2Full[i * 2][j * 2];
            diff[i][j] = qAbs(v[i][j] - v2value);
        }
    }

    return { v, v2Full, diff };
}

