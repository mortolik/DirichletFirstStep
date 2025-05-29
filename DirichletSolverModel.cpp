#include "DirichletSolverModel.hpp"
#include "qdebug.h"


#include <QtMath>
#include <QTextStream>
#include <QElapsedTimer>

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
    // Сохраняем основные параметры:
    m_n       = n;
    m_m       = m;
    m_eps     = eps;
    m_maxIter = maxIter;

    // Считаем шаги по x/y:
    m_h = (m_b - m_a) / double(m_n);
    m_k = (m_d - m_c) / double(m_m);

    // Предвычисляем константы:
    m_invH2 = 1.0 / (m_h * m_h);
    m_invK2 = 1.0 / (m_k * m_k);
    m_denom = 2.0 * (m_invH2 + m_invK2);

    // Инициализируем контейнер m_u и m_uExact нужного размера и нулями:
    m_u.resize(m_n + 1);
    m_uExact.resize(m_n + 1);
    for (int i = 0; i <= m_n; ++i) {
        m_u[i].fill(0.0, m_m + 1);
        m_uExact[i].fill(0.0, m_m + 1);
    }

    // Плоский буфер:
    m_uPrevFlat.resize((m_n + 1) * (m_m + 1));

    // Заполняем границы и точное решение:
    initialize();  // ваша функция, которая делает applyBoundary() и заполняет m_uExact

    // Копируем gраничные и начальные значения в плоский буфер:
    for (int i = 0; i <= m_n; ++i)
        for (int j = 0; j <= m_m; ++j)
            m_uPrevFlat[flatIdx(i, j)] = m_u[i][j];

    // Сбрасываем отчётные поля:
    m_lastIter     = 0;
    m_lastResidual = 0.0;
    m_resultCached = false;
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
    qDebug() << maxR;
    return maxR;
}

void DirichletSolverModel::solveMainProblem()
{
    QElapsedTimer timer;
    timer.start();
    if (!m_skipInitialization)
        initializeInterior();
    // 1) Вычисляем шаги и константы только один раз
    m_h     = (m_b - m_a) / m_n;
    m_k     = (m_d - m_c) / m_m;
    m_invH2 = 1.0 / (m_h * m_h);
    m_invK2 = 1.0 / (m_k * m_k);
    m_denom = 2.0 * (m_invH2 + m_invK2);

    int rows = m_n + 1, cols = m_m + 1;
    // 2) Выделяем/заполняем плоский буфер
    m_uPrevFlat.resize(rows * cols);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            m_uPrevFlat[flatIdx(i,j)] = m_u[i][j];

    int iter = 0;
    double maxDiff;
    do {
        maxDiff = 0.0;
#pragma omp parallel for reduction(max:maxDiff)
        for(int i = 1; i < m_n; ++i) {
            double x = m_a + i * m_h;
            for(int j = 1; j < m_m; ++j) {
                int idx = flatIdx(i,j);
                double y   = m_c + j * m_k;
                double rhs = (m_uPrevFlat[idx + cols] + m_uPrevFlat[idx - cols]) * m_invH2
                             + (m_uPrevFlat[idx + 1    ] + m_uPrevFlat[idx - 1    ]) * m_invK2
                             - f(x, y);

                double uNew = (1.0 - m_omega)*m_uPrevFlat[idx]
                              + m_omega * (rhs / m_denom);
                double diff = std::abs(uNew - m_uPrevFlat[idx]);
                m_uPrevFlat[idx] = uNew;
                if (diff > maxDiff) maxDiff = diff;
            }
        }
        ++iter;
    } while (maxDiff > m_eps && iter < m_maxIter);

    // 3) Копируем обратно в m_u
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            m_u[i][j] = m_uPrevFlat[flatIdx(i,j)];

    m_lastIter     = iter;
    qDebug() << m_lastIter;
    qDebug() <<"Omega = " << m_omega;
    qDebug() <<"Вычислено за" << timer.elapsed() << "ms или"
             << timer.elapsed()/ 60000.0 << "min";
    m_lastResidual = maxDiff;
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
    DirichletSolverModel finerModel;
    finerModel.setup(finerN, finerM, m_eps, m_maxIter);
    finerModel.solveMainProblem();

    const auto &fineU = finerModel.solution();
    int scaleI = finerN / m_n;
    int scaleJ = finerM / m_m;

    double maxError = 0.0;
    QVector<QVector<double>> diff(m_n + 1, QVector<double>(m_m + 1, 0.0));
    for (int i = 0; i <= m_n; ++i) {
        for (int j = 0; j <= m_m; ++j) {
            double delta = std::abs(m_u[i][j] - fineU[i * scaleI][j * scaleJ]);
            diff[i][j] = delta;
            if (delta > maxError)
                maxError = delta;
        }
    }

    eps2Out = maxError;
    return diff;
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
        m_result = computeFinerGridComparison();
        m_resultCached = true;
        for (int i = 0; i <= m_n; ++i)
            for (int j = 0; j <= m_m; ++j)
                eps2 = std::max(eps2, m_result.diff[i][j]);

        data.accuracy = eps2;
        qDebug() << eps2;
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
        lines << QString("Количество итераций на удвоенной сетке: %1").arg(m_last2Iter);
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

int DirichletSolverModel::lastIter() const
{
    return m_lastIter;
}

QPair<double, double> DirichletSolverModel::maxErrorPointCompare() const
{
    const auto &R = computeFinerGridComparison();
    return R.maxPt;
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

double DirichletSolverModel::computeOptimalOmega() const
{
    double rho = std::cos(M_PI / m_n) + std::cos(M_PI / m_m);
    double value = std::sqrt(1 - std::pow(rho / 2.0, 2));
    return 2.0 / (1.0 + value);
}

void DirichletSolverModel::applyInterpolatedInitialGuess(
    const QVector<QVector<double>>& coarseU, int coarseN, int coarseM)
{
    // 1) билинейная интерполяция для внутренних узлов
    for (int i = 1; i < m_n; ++i) {
        double x = m_a + i * m_h;
        double u = (x - m_a) / (m_b - m_a) * coarseN;
        int i0 = std::clamp(int(floor(u)), 0, coarseN-1);
        double dx = u - i0;

        for (int j = 1; j < m_m; ++j) {
            double y = m_c + j * m_k;
            double v = (y - m_c) / (m_d - m_c) * coarseM;
            int j0 = std::clamp(int(floor(v)), 0, coarseM-1);
            double dy = v - j0;

            double f00 = coarseU[i0  ][j0  ];
            double f10 = coarseU[i0+1][j0  ];
            double f01 = coarseU[i0  ][j0+1];
            double f11 = coarseU[i0+1][j0+1];

            double w0 = (1 - dx)*(1 - dy);
            double w1 = dx*(1 - dy);
            double w2 = (1 - dx)*dy;
            double w3 = dx*dy;

            m_u[i][j] = f00*w0 + f10*w1 + f01*w2 + f11*w3;
        }
    }

    // 2) Снова жёстко задаём граничные условия
    for (int i = 0; i <= m_n; ++i) {
        double x = m_a + i * m_h;
        m_u[i][0]   = mu3(x);
        m_u[i][m_m] = mu4(x);
    }
    for (int j = 0; j <= m_m; ++j) {
        double y = m_c + j * m_k;
        m_u[0][j]   = mu1(y);
        m_u[m_n][j] = mu2(y);
    }

    // 3) Плоский буфер
    for (int i = 0; i <= m_n; ++i)
        for (int j = 0; j <= m_m; ++j)
            m_uPrevFlat[flatIdx(i,j)] = m_u[i][j];
}

DirichletSolverModel::FinerGridResult DirichletSolverModel::computeFinerGridComparison() const
{
    if (m_resultCached)
    {
        return m_result;
    }

    int finerN = m_n * 2;
    int finerM = m_m * 2;

    // 1. Создаём и решаем на утончённой сетке
    DirichletSolverModel finer;
    finer.setup(finerN, finerM, m_eps, m_maxIter);
    double omega = computeOptimalOmega();
    finer.setOmega(omega);
    finer.applyInterpolatedInitialGuess(m_u, m_n, m_m);
    finer.m_skipInitialization = true;
    finer.solveMainProblem();
    m_last2Iter = finer.lastIter();

    auto const &fineU = finer.solution();
    int scaleI = finerN / m_n;
    int scaleJ = finerM / m_m;

    // 2) Готовим результат
    FinerGridResult R;
    R.iterations = finer.lastIter();
    R.v    = m_u;
    R.v2   .resize(m_n + 1);
    R.diff .resize(m_n + 1);

    double maxErr = 0.0;
    int maxI = 0, maxJ = 0;

// 3) Считаем разности и запоминаем точку максимума
#pragma omp parallel for reduction(max:maxErr)
    for (int i = 0; i <= m_n; ++i) {
        R.v2[i].resize(m_m + 1);
        R.diff[i].resize(m_m + 1);
        for (int j = 0; j <= m_m; ++j) {
            double fineVal = fineU[i * scaleI][j * scaleJ];
            double delta   = std::abs(m_u[i][j] - fineVal);
            R.v2  [i][j] = fineVal;
            R.diff[i][j] = delta;
            if (delta > maxErr) {
                maxErr = delta;
                maxI = i;
                maxJ = j;
            }
        }
    }

    // 4) Сохраняем точку максимального отклонения
    R.maxPt.first  = m_a + maxI * m_h;
    R.maxPt.second = m_c + maxJ * m_k;

    m_result       = R;
    m_resultCached = true;
    return m_result;
}
