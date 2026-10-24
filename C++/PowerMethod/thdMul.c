#include "../common/SparseMultiply.h"
void utmvMul(int *i,int *j,double *x,long m,double *v,int n,double *res,int threads,double **space) {
    possum::sparseMultiply(i,j,x,m,v,n,res,threads,space);
}
