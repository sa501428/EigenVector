#include "../common/CompartmentOrientation.h"

// Compatibility entry point for existing users; both solvers share this logic.
int flipSign(const char *genome, double *x, int n, char *chr, int binsize) {
    if (!genome || !chr) return 0;
    return possum::CompartmentOrientation(genome).orient(x, n, chr, binsize);
}
