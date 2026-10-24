#ifndef POSSUM_SPARSE_MULTIPLY_H
#define POSSUM_SPARSE_MULTIPLY_H
#include <pthread.h>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace possum {
template<class Index, class Value> struct MultiplyTask {
    Index *i, *j;
    Value *x;
    double *v, *result;
    long count;
    static void *run(void *argument) {
        MultiplyTask &a = *static_cast<MultiplyTask *>(argument);
        for (long p=0; p<a.count; ++p) {
            a.result[a.i[p]] += a.x[p]*a.v[a.j[p]];
            a.result[a.j[p]] += a.x[p]*a.v[a.i[p]];
        }
        return NULL;
    }
};

template<class Index, class Value>
void sparseMultiply(Index *i, Index *j, Value *x, long m, double *v, int n,
                    double *result, int threads, double **space) {
    if (threads <= 0) throw std::runtime_error("thread count must be positive");
    // Use scratch even for one thread: callers may alias v and result.
    if (threads == 1) {
        std::fill(space[0],space[0]+n,0.0);
        MultiplyTask<Index,Value> task = {i,j,x,v,space[0],m};
        MultiplyTask<Index,Value>::run(&task);
        std::copy(space[0],space[0]+n,result);
        return;
    }
    std::vector<MultiplyTask<Index,Value> > tasks(threads);
    std::vector<pthread_t> handles(threads);
    std::vector<bool> started(threads,false);
    const long slice = m/threads;
    for (int t=0; t<threads; ++t) {
        std::fill(space[t],space[t]+n,0.0);
        tasks[t] = {i+t*slice,j+t*slice,x+t*slice,v,space[t],t==threads-1 ? m-t*slice : slice};
        started[t] = pthread_create(&handles[t],NULL,MultiplyTask<Index,Value>::run,&tasks[t]) == 0;
        if (!started[t]) MultiplyTask<Index,Value>::run(&tasks[t]);
    }
    for (int t=0; t<threads; ++t) if (started[t]) pthread_join(handles[t],NULL);
    std::fill(result,result+n,0.0);
    for (int t=0; t<threads; ++t)
        for (int p=0; p<n; ++p) result[p] += space[t][p];
}
}
#endif
