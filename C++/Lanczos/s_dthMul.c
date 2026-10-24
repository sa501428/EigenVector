#include "../common/SparseMultiply.h"
void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int n,double *res,int threads,double **space) {
    possum::sparseMultiply(i,j,x,m,v,n,res,threads,space);
}
