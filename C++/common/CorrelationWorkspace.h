#ifndef POSSUM_CORRELATION_WORKSPACE_H
#define POSSUM_CORRELATION_WORKSPACE_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace possum {
// Sparse inputs hold one triangle, with diagonal values halved by the caller.
// Drop empty/constant rows before forming Pearson correlations. Recompute
// moments after removing rows so the two solvers use the same observations.
template<class Index, class Value>
int correlationStatistics(long m, Index *i, Index *j, Value *x, int n,
                          std::vector<double> &good, std::vector<double> &mean,
                          std::vector<double> &inverseSd) {
    good.assign(n, 1.0);
    mean.resize(n);
    inverseSd.resize(n);
    std::vector<double> sum(n), squares(n);
    for (long p=0; p<m; ++p) {
        if (i[p] >= static_cast<Index>(n) || j[p] >= static_cast<Index>(n) ||
            static_cast<long long>(i[p]) < 0 || static_cast<long long>(j[p]) < 0)
            throw std::runtime_error("contact coordinate outside chromosome");
        if (!std::isfinite(x[p])) x[p] = 0;
    }
    int active = n;
    for (;;) {
        std::fill(sum.begin(), sum.end(), 0.0);
        std::fill(squares.begin(), squares.end(), 0.0);
        for (long p=0; p<m; ++p) {
            if (!good[i[p]] || !good[j[p]]) { x[p] = 0; continue; }
            double value = x[p];
            sum[i[p]] += value;
            sum[j[p]] += value;
            squares[i[p]] += value*value;
            squares[j[p]] += value*value;
            if (i[p] == j[p]) squares[i[p]] += 2*value*value;
        }
        int retained = 0;
        for (int p=0; p<n; ++p) {
            mean[p] = active > 0 ? sum[p]/active : 0;
            double varianceSum = squares[p] - sum[p]*mean[p];
            if (good[p] && varianceSum > 1e-12*squares[p]) {
                inverseSd[p] = 1/std::sqrt(varianceSum);
                ++retained;
            } else {
                good[p] = mean[p] = inverseSd[p] = 0;
            }
        }
        if (retained == active || retained == 0) return retained;
        active = retained;
    }
}

inline double vectorDot(const double *x, const double *y, int n) {
    double sum = 0;
    for (int p=0; p<n; ++p) sum += x[p]*y[p];
    return sum;
}
inline double vectorNorm(const double *x, int n) {
    return std::sqrt(vectorDot(x, x, n));
}
inline bool unitVector(double *x, int n) {
    double length = vectorNorm(x, n);
    if (!(length > 0) || !std::isfinite(length)) return false;
    for (int p=0; p<n; ++p) x[p] /= length;
    return true;
}
} // namespace possum
#endif
