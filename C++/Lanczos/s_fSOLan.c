#include <lapacke.h>
#include "../common/CorrelationWorkspace.h"

void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads) {
    if (m < 0 || nv <= 0 || nv > maxiter-3 || threads <= 0 || !std::isfinite(tol) || tol <= 0 ||
        !std::isfinite(eps) || eps <= 0 || *N > unsigned(std::numeric_limits<int>::max()))
        throw std::runtime_error("invalid Lanczos solver arguments");
    int n = *N;
    if (!n) {
        for (long p=0; p<m; ++p) n = std::max(n, int(std::max(i[p],j[p]))+1);
        *N = n;
    }
    if (n <= 0) return 300001;
    std::vector<double> good, mean, d;
    int active = possum::correlationStatistics(m,i,j,x,n,good,mean,d);
    if (active < nv) return 300001;
    std::vector<std::vector<double> > buffers(threads, std::vector<double>(n));
    std::vector<double *> space(threads);
    for (int t=0; t<threads; ++t) space[t] = buffers[t].data();
    std::vector<double> temp1(n), temp2(n), u(n), y(n);
    auto multiply = [&](const double *a, double *result) {
        for (int p=0; p<n; ++p) temp1[p] = d[p]*a[p];
        utmvMul(i,j,x,m,temp1.data(),n,temp2.data(),threads,space.data());
        double z = possum::vectorDot(temp1.data(),mean.data(),n);
        for (int p=0; p<n; ++p) temp2[p] -= good[p]*z;
        utmvMul(i,j,x,m,temp2.data(),n,temp1.data(),threads,space.data());
        z = possum::vectorDot(good.data(),temp2.data(),n);
        for (int p=0; p<n; ++p) result[p] = d[p]*(temp1[p]-mean[p]*z);
    };
    // A deterministic, nonconstant start avoids an all-ones vector lying in
    // the nullspace of a centered correlation matrix.
    unsigned int seed = 1234563;
    for (int p=0; p<n; ++p) {
        seed = 1664525u*seed+1013904223u;
        r[p] = good[p]*(double(seed)/4294967296.0-0.5);
    }
    possum::unitVector(r,n);
    const int limit = std::min(maxiter,active);
    std::vector<std::vector<double> > q;
    std::vector<double> alpha, beta, diagonal, offDiagonal, z;
    double residual = std::numeric_limits<double>::infinity();
    int steps = 0;
    auto diagonalize = [&]() {
        diagonal = alpha;
        offDiagonal = beta;
        z.resize(size_t(steps)*steps);
        // A symmetric tridiagonal routine accepts semidefinite and rank-
        // deficient matrices; dpteqr incorrectly required positive definiteness.
        return LAPACKE_dstev(LAPACK_COL_MAJOR,'V',steps,diagonal.data(),offDiagonal.data(),z.data(),steps);
    };
    for (int step=0; step<limit; ++step) {
        q.push_back(std::vector<double>(r,r+n));
        multiply(q.back().data(),u.data());
        for (int p=0; p<n; ++p) r[p] = u[p]-(step ? beta[step-1]*q[step-1][p] : 0);
        double a = possum::vectorDot(q.back().data(),r,n);
        alpha.push_back(a);
        for (int p=0; p<n; ++p) r[p] -= a*q.back()[p];
        // Local reorthogonalization controls cancellation in the recurrence.
        for (int previous=std::max(0,step-1); previous<=step; ++previous) {
            double projection = possum::vectorDot(q[previous].data(),r,n);
            for (int p=0; p<n; ++p) r[p] -= projection*q[previous][p];
        }
        double b = possum::vectorNorm(r,n);
        beta.push_back(b);
        steps = step+1;
        int info = diagonalize();
        if (info) return 100000+info;
        if (steps >= nv) {
            residual = 0;
            for (int v=0; v<nv; ++v)
                residual = std::max(residual,b*std::fabs(z[size_t(steps-1-v)*steps+steps-1]));
            double gap = steps > nv ? diagonal[steps-nv]-diagonal[steps-nv-1] : 0;
            if (gap > 0 && residual <= tol*gap) break;
        }
        if (b <= 64*std::numeric_limits<double>::epsilon()*std::max(1.0,std::fabs(a))) break;
        // Selectively remove converged Ritz vectors from the next residual.
        for (int v=0; v<std::min(steps,nv+2); ++v) {
            int column = steps-1-v;
            if (b*std::fabs(z[size_t(column)*steps+steps-1]) > eps*std::max(1.0,std::fabs(diagonal.back()))) continue;
            std::fill(y.begin(),y.end(),0.0);
            for (int k=0; k<steps; ++k)
                for (int p=0; p<n; ++p) y[p] += q[k][p]*z[size_t(column)*steps+k];
            if (!possum::unitVector(y.data(),n)) continue;
            double projection = possum::vectorDot(r,y.data(),n);
            for (int p=0; p<n; ++p) r[p] -= projection*y[p];
        }
        b = possum::vectorNorm(r,n);
        beta.back() = b;
        if (b <= 64*std::numeric_limits<double>::epsilon()*std::max(1.0,std::fabs(a))) break;
        for (int p=0; p<n; ++p) r[p] /= b;
    }
    if (steps < nv) return 300001;
    for (int v=0; v<nv+2; ++v) {
        int column = steps-1-v;
        lam[v] = column >= 0 ? diagonal[column] : 0.0;
        for (int p=0; p<n; ++p) {
            ev[v][p] = good[p] ? 0.0 : std::numeric_limits<double>::quiet_NaN();
            if (good[p] && column >= 0)
                for (int k=0; k<steps; ++k) ev[v][p] += q[k][p]*z[size_t(column)*steps+k];
        }
    }
    *er = residual;
    return steps;
}
