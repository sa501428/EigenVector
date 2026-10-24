#include "../common/CorrelationWorkspace.h"

void utmvMul(int *i,int *j,double *x,long m,double *v,int k,double *res, int threads, double **space);

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads) {
    if (m < 0 || threads <= 0 || maxiter <= 0 || !std::isfinite(tol) || tol <= 0)
        throw std::runtime_error("invalid Power Method solver arguments");
    int n = *N;
    if (n < 0) {
        n = 0;
        for (long p=0; p<m; ++p) n = std::max(n, std::max(i[p],j[p])+1);
    }
    if (n <= 0) throw std::runtime_error("empty matrix dimension");
    *N = n;
    std::vector<double> good, mean, d;
    int active = possum::correlationStatistics(m,i,j,x,n,good,mean,d);
    *l1 = *l2 = *er = std::numeric_limits<double>::quiet_NaN();
    std::fill(ev, ev+n, std::numeric_limits<double>::quiet_NaN());
    if (!active) return 0;
    std::vector<std::vector<double> > buffers(threads, std::vector<double>(n));
    std::vector<double *> space(threads);
    for (int t=0; t<threads; ++t) space[t] = buffers[t].data();
    std::vector<double> temp1(n), temp2(n), image(n);
    auto multiply = [&](const double *a, double *result) {
        for (int p=0; p<n; ++p) temp1[p] = d[p]*a[p];
        utmvMul(i,j,x,m,temp1.data(),n,temp2.data(),threads,space.data());
        double z = possum::vectorDot(temp1.data(),mean.data(),n);
        for (int p=0; p<n; ++p) temp2[p] -= good[p]*z;
        utmvMul(i,j,x,m,temp2.data(),n,temp1.data(),threads,space.data());
        z = possum::vectorDot(good.data(),temp2.data(),n);
        for (int p=0; p<n; ++p) result[p] = d[p]*(temp1[p]-mean[p]*z);
    };
    for (int p=0; p<n; ++p) { a1[p] *= good[p]; a2[p] *= good[p]; }
    if (!possum::unitVector(a1,n))
        for (int p=0; p<n; ++p) a1[p] = good[p]*std::sin(p+1.0);
    multiply(a1,image.data());
    if (possum::vectorNorm(image.data(),n) < 1e-14) {
        for (int p=0; p<n; ++p) a1[p] = good[p]*std::sin(p+1.0);
        possum::unitVector(a1,n);
    }
    if (!possum::unitVector(a2,n)) {
        for (int p=0; p<n; ++p) a2[p] = good[p]*std::cos(p+1.0);
        possum::unitVector(a2,n);
    }
    int iterations = 0;
    double lambda1 = 0, lambda2 = 0, residual = 0;
    for (; iterations < maxiter; ++iterations) {
        multiply(a1,image.data());
        if (!possum::unitVector(image.data(),n)) break;
        std::copy(image.begin(),image.end(),a1);
        double projection = possum::vectorDot(a1,a2,n);
        for (int p=0; p<n; ++p) a2[p] -= projection*a1[p];
        multiply(a2,image.data());
        projection = possum::vectorDot(a1,image.data(),n);
        for (int p=0; p<n; ++p) image[p] -= projection*a1[p];
        if (possum::unitVector(image.data(),n)) std::copy(image.begin(),image.end(),a2);
        else std::fill(a2,a2+n,0.0); // A rank-one matrix has no second direction.
        multiply(a1,image.data());
        lambda1 = possum::vectorDot(a1,image.data(),n);
        residual = 0;
        for (int p=0; p<n; ++p) residual += std::pow(image[p]-lambda1*a1[p],2);
        residual = std::sqrt(residual);
        multiply(a2,image.data());
        lambda2 = possum::vectorDot(a2,image.data(),n);
        if (lambda2 > lambda1) {
            for (int p=0; p<n; ++p) std::swap(a1[p],a2[p]);
            std::swap(lambda1,lambda2);
            continue; // Recompute the residual for the new leading vector.
        }
        if (residual <= tol*std::max(lambda1-lambda2, std::numeric_limits<double>::epsilon()*lambda1)) {
            ++iterations;
            break;
        }
    }
    // Always report the residual and vector from the same final iterate.
    multiply(a1,image.data());
    lambda1 = possum::vectorDot(a1,image.data(),n);
    residual = 0;
    for (int p=0; p<n; ++p) {
        residual += std::pow(image[p]-lambda1*a1[p],2);
        if (good[p]) ev[p] = a1[p];
    }
    *l1 = lambda1;
    *l2 = lambda2;
    *er = std::sqrt(residual);
    return iterations;
}
