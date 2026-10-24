#include <straw.h>
#include "../common/CompartmentOrientation.h"

// Straw returns axes in the requested chromosome order. Inferring a transpose
// from chromosome lengths and the largest observed contact misassigned sparse
// trans contacts, sometimes placing them outside a chromosome.
unsigned int getMatrix(std::string fname, int binsize, std::string norm, std::string ob,
                       bool interOnly, unsigned int **i0, unsigned int **j0, float **x0,
                       long *m, std::vector<std::string> &chroms, std::vector<int64_t> &lengths,
                       std::string &genomeID) {
    std::ifstream in(fname.c_str(),std::ios::binary);
    if (!in) throw std::runtime_error("cannot open Hi-C file: " + fname);
    int64_t master=0,nviPosition=0,nviLength=0;
    int32_t count=0,version=0;
    auto chromosomes = readHeader(in,master,genomeID,count,version,nviPosition,nviLength);
    std::vector<unsigned int> offsets;
    unsigned int dimension=0;
    for (const auto &entry : chromosomes) {
        const chromosome &c = entry.second;
        if (possum::isAllChromosome(c.name)) continue;
        int bins = possum::binCount(c.length,binsize);
        if (dimension > unsigned(std::numeric_limits<int>::max()-bins))
            throw std::runtime_error("genome matrix has too many bins");
        offsets.push_back(dimension);
        dimension += bins;
        chroms.push_back(c.name);
        lengths.push_back(c.length);
    }
    std::vector<unsigned int> rows,columns;
    std::vector<float> values;
    for (size_t a=0; a<chroms.size(); ++a) {
        for (size_t b=a+(interOnly ? 1 : 0); b<chroms.size(); ++b) {
            auto records = straw(ob,norm,fname,chroms[a],chroms[b],"BP",binsize);
            for (const auto &record : records) {
                if (!std::isfinite(record.counts) || record.counts == 0) continue;
                if (record.binX < 0 || record.binY < 0 || record.binX >= lengths[a] || record.binY >= lengths[b])
                    throw std::runtime_error("trans contact coordinate outside chromosome");
                rows.push_back(offsets[a]+record.binX/binsize);
                columns.push_back(offsets[b]+record.binY/binsize);
                values.push_back(record.counts);
            }
        }
    }
    // Size from actual returned records, rather than an unrelated pre-count
    // that can underestimate derived/normalized matrices and overrun allocations.
    *m = static_cast<long>(values.size());
    *i0 = static_cast<unsigned int *>(std::malloc(values.size()*sizeof(unsigned int)));
    *j0 = static_cast<unsigned int *>(std::malloc(values.size()*sizeof(unsigned int)));
    *x0 = static_cast<float *>(std::malloc(values.size()*sizeof(float)));
    if (!values.empty() && (!*i0 || !*j0 || !*x0)) throw std::bad_alloc();
    std::copy(rows.begin(),rows.end(),*i0);
    std::copy(columns.begin(),columns.end(),*j0);
    std::copy(values.begin(),values.end(),*x0);
    return dimension;
}
