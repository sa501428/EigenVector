#include <lapacke.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef POWER_SOLVER
int bestEigen(long,int*,int*,double*,int*,double*,double*,double*,double*,double*,double*,double,int,int);
#else
int SOLan(long,unsigned int*,unsigned int*,float*,unsigned int*,double*,int,double*,double**,double*,double,double,int,int);
#endif
void check(bool ok, const char *message) { if (!ok) throw std::runtime_error(message); }

void compare(std::vector<double> matrix, int n, int threads, bool rankOne=false, int vectorCount=3) {
    // Independently remove empty/constant rows and compute dense Pearson
    // correlations, without using any of the sparse implementation helpers.
    std::vector<int> retained;
    for (int i=0; i<n; ++i) retained.push_back(i);
    for (;;) {
        std::vector<int> next;
        for (int i : retained) {
            double mean=0,variance=0;
            for (int j : retained) mean += matrix[i*n+j]/retained.size();
            for (int j : retained) variance += std::pow(matrix[i*n+j]-mean,2);
            if (variance > 1e-10) next.push_back(i);
        }
        if (next==retained) break;
        retained=next;
    }
    std::vector<double> correlation(n*n,0),mean(n),sd(n);
    for (int i : retained) {
        for (int j : retained) mean[i] += matrix[i*n+j]/retained.size();
        for (int j : retained) sd[i] += std::pow(matrix[i*n+j]-mean[i],2);
        sd[i]=std::sqrt(sd[i]);
    }
    for (int i : retained) for (int j : retained)
        for (int k : retained) correlation[i*n+j] +=
            (matrix[i*n+k]-mean[i])*(matrix[j*n+k]-mean[j])/(sd[i]*sd[j]);
    std::vector<double> reference=correlation,eigenvalues(n);
    check(LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,reference.data(),n,eigenvalues.data())==0,"dense LAPACK oracle failed");
#ifdef POWER_SOLVER
    std::vector<int> rows,cols;
    std::vector<double> values;
#else
    std::vector<unsigned int> rows,cols;
    std::vector<float> values;
#endif
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
        rows.push_back(i); cols.push_back(j);
        values.push_back(matrix[i*n+j]*(i==j ? .5 : 1));
    }
    const int requested=rankOne ? 1 : vectorCount;
    std::vector<std::vector<double> > vectors(requested+2,std::vector<double>(n));
    std::vector<double *> output;
    for (auto &v : vectors) output.push_back(v.data());
    std::vector<double> lambda(requested+2);
    double error;
#ifdef POWER_SOLVER
    int dimension=n;
    std::vector<double> a1(n),a2(n);
    for (int i=0; i<n; ++i) { a1[i]=std::sin(i+1.0); a2[i]=std::cos(i+2.0); }
    int iterations=bestEigen(values.size(),rows.data(),cols.data(),values.data(),&dimension,
        &lambda[0],&lambda[1],a1.data(),a2.data(),output[0],&error,1e-8,1000,threads);
    const int count=1;
#else
    unsigned int dimension=n;
    std::vector<double> start(n,1.0);
    int iterations=SOLan(values.size(),rows.data(),cols.data(),values.data(),&dimension,
        start.data(),requested,lambda.data(),output.data(),&error,1e-8,1e-10,100,threads);
    const int count=requested;
#endif
    check(iterations>=0 && iterations<90000,"sparse eigensolver failed");
    for (int ev=0; ev<count; ++ev) {
        double expected=eigenvalues[n-1-ev];
        if (std::fabs(lambda[ev]-expected) > 2e-5*std::max(1.,expected)) {
            std::cerr << "EV" << ev+1 << " got " << lambda[ev] << " expected " << expected << '\n';
            throw std::runtime_error("eigenvalue disagrees with dense Pearson matrix");
        }
        double alignment=0,norm=0,residual=0;
        for (int i : retained) {
            alignment += vectors[ev][i]*reference[i*n+n-1-ev];
            norm += vectors[ev][i]*vectors[ev][i];
            double image=0;
            for (int j : retained) image += correlation[i*n+j]*vectors[ev][j];
            residual += std::pow(image-lambda[ev]*vectors[ev][i],2);
        }
        check(std::fabs(std::fabs(alignment)-1)<1e-4,"eigenvector disagrees with dense reference");
        check(std::fabs(norm-1)<1e-5,"eigenvector norm");
        check(std::sqrt(residual)<1e-4,"eigenvector residual");
        for (int i=0; i<n; ++i) if (std::find(retained.begin(),retained.end(),i)==retained.end())
            check(std::isnan(vectors[ev][i]),"masked bins must be NaN internally");
    }
}
int main() try {
    const int n=28;
    std::vector<double> matrix(n*n);
    unsigned int seed=42;
    for (int i=0; i<n; ++i) for (int j=i; j<n; ++j) {
        seed=1664525u*seed+1013904223u;
        double value=1+double(seed%1000)/300.;
        if ((i<n/2)==(j<n/2)) value += 5;
        matrix[i*n+j]=matrix[j*n+i]=value;
    }
    for (int threads : {1,3}) {
        compare(matrix,n,threads);
        compare(matrix,n,threads,false,10);
        auto holes=matrix;
        for (int i=0;i<n;++i) holes[i*n+4]=holes[4*n+i]=0;
        compare(holes,n,threads);
        auto constant=matrix;
        for (int i=0;i<n;++i) constant[i*n+4]=constant[4*n+i]=2;
        compare(constant,n,threads);
        std::vector<double> rankOne={3,3,1,1, 3,3,1,1, 1,1,3,3, 1,1,3,3};
        compare(rankOne,4,threads,true);
    }
    std::cout << "Sparse solver agrees with dense LAPACK Pearson eigenpairs (dense, missing bins, constant bins, rank one; 1/3 threads)\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
}
