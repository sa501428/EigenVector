#ifndef POSSUM_BUNDLED_PHASING_H
#define POSSUM_BUNDLED_PHASING_H
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

namespace possum {
class PhasingGzip {
    gzFile file_;
public:
    explicit PhasingGzip(const std::string &path) : file_(gzopen(path.c_str(), "rb")) {
        if (!file_) throw std::runtime_error("cannot open bundled reference: " + path);
    }
    ~PhasingGzip() { gzclose(file_); }
    PhasingGzip(const PhasingGzip &) = delete;
    PhasingGzip &operator=(const PhasingGzip &) = delete;
    int read(char *buffer, unsigned int size) {
        int got = gzread(file_, buffer, size);
        if (got < 0) throw std::runtime_error("corrupt gzip phasing reference");
        return got;
    }
    void exact(char *buffer, size_t size) {
        while (size) {
            int got = read(buffer, static_cast<unsigned int>(std::min<size_t>(size, 1 << 20)));
            if (!got) throw std::runtime_error("truncated bundled phasing reference");
            buffer += got; size -= got;
        }
    }
    bool line(std::string &text) {
        text.clear();
        char buffer[4096];
        while (gzgets(file_, buffer, sizeof(buffer))) {
            text += buffer;
            if (!text.empty() && text.back() == '\n') return true;
        }
        if (!gzeof(file_)) throw std::runtime_error("corrupt gzip phasing reference");
        return !text.empty();
    }
};

struct PhasingBins {
    std::vector<double> gc, genes, custom;
    bool gcExact = true;
    bool lengthMismatch = false;
};

class BundledPhasing {
    std::string folder_;
    std::vector<std::pair<std::string, int64_t> > chromosomes_;
    std::map<std::string, int64_t> sizes_;
    mutable std::map<std::string, std::vector<unsigned char> > counts_;
    mutable std::map<std::string, std::vector<std::pair<int64_t, int64_t> > > genes_;
    mutable bool genesLoaded_ = false;

    std::string expectedHeader() const {
        // Version 1 is written canonically by scripts/phasing_reference.py.
        std::ostringstream out;
        out << "{\"version\":1,\"bin_size\":100,\"columns\":[\"gc_bases\",\"acgt_bases\"],"
               "\"dtype\":\"uint8\",\"layout\":\"interleaved\",\"chromosomes\":[";
        for (size_t i = 0; i < chromosomes_.size(); ++i) {
            if (i) out << ',';
            out << "[\"" << chromosomes_[i].first << "\"," << chromosomes_[i].second << ']';
        }
        out << "]}";
        return out.str();
    }

    const std::vector<unsigned char> &counts(const std::string &chrom) const {
        auto found = counts_.find(chrom);
        if (found != counts_.end()) return found->second;
        PhasingGzip in(folder_ + "/gc.100bp.counts.bin.gz");
        char prefix[12]; in.exact(prefix, sizeof(prefix));
        if (std::string(prefix, 8) != "GCCOUNT1") throw std::runtime_error("invalid bundled GC magic");
        uint32_t headerSize = 0;
        for (int i = 0; i < 4; ++i) headerSize |= uint32_t(static_cast<unsigned char>(prefix[8+i])) << (8*i);
        if (!headerSize || headerSize > 1000000) throw std::runtime_error("invalid bundled GC header size");
        std::string header(headerSize, '\0'); in.exact(&header[0], headerSize);
        if (header != expectedHeader()) throw std::runtime_error("bundled GC header/chrom.sizes mismatch");
        for (const auto &item : chromosomes_) {
            size_t n = static_cast<size_t>((item.second + 99) / 100);
            std::vector<unsigned char> data(2*n);
            in.exact(reinterpret_cast<char *>(data.data()), data.size());
            if (item.first == chrom) {
                for (size_t i = 0; i < n; ++i) {
                    int width = static_cast<int>(std::min<int64_t>(100, item.second-int64_t(i)*100));
                    if (data[2*i] > data[2*i+1] || data[2*i+1] > width)
                        throw std::runtime_error("invalid bundled GC counts on " + chrom);
                }
                counts_[chrom] = std::move(data);
            }
        }
        char tail;
        if (in.read(&tail, 1)) throw std::runtime_error("unexpected trailing bundled GC data");
        return counts_.at(chrom);
    }

    void loadGenes() const {
        if (genesLoaded_) return;
        PhasingGzip in(folder_ + "/refseq.coverage.bedGraph.gz");
        std::map<std::string, int64_t> cursors;
        std::string line;
        while (in.line(line)) {
            std::istringstream row(line);
            std::string chrom, extra;
            int64_t start, end;
            int value;
            if (!(row >> chrom >> start >> end >> value) || (row >> extra) ||
                !sizes_.count(chrom) || start != cursors[chrom] || end <= start || end > sizes_.at(chrom) ||
                (value != 0 && value != 1)) throw std::runtime_error("invalid bundled RefSeq coverage");
            cursors[chrom] = end;
            if (value) genes_[chrom].push_back({start, end});
        }
        if (cursors != sizes_) throw std::runtime_error("incomplete bundled RefSeq coverage");
        genesLoaded_ = true;
    }

public:
    explicit BundledPhasing(const std::string &folder) : folder_(folder) {
        std::ifstream in(folder + "/chrom.sizes");
        if (!in) throw std::runtime_error("cannot find bundled references at " + folder +
            "; set --phasing-dir or POSSUM_PHASING_DIR, or explicitly use --phasing-method none");
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string chrom, extra;
            int64_t size;
            if (!(row >> chrom >> size) || (row >> extra) || size <= 0 || sizes_.count(chrom))
                throw std::runtime_error("invalid bundled chrom.sizes");
            sizes_[chrom] = size;
            chromosomes_.push_back({chrom, size});
        }
        if (!in.eof() || sizes_.empty()) throw std::runtime_error("empty/unreadable bundled chrom.sizes");
    }

    int64_t length(const std::string &chrom) const {
        auto it = sizes_.find(chrom);
        return it == sizes_.end() ? 0 : it->second;
    }

    PhasingBins bin(const std::string &chrom, int n, int resolution, int64_t inputLength,
                    bool wantGC, bool wantGenes, bool approximate) const {
        PhasingBins result;
        result.gcExact = resolution % 100 == 0;
        double missing = std::numeric_limits<double>::quiet_NaN();
        result.gc.assign(n, missing); result.genes.assign(n, missing);
        int64_t size = length(chrom);
        if (!size) return result;
        if ((inputLength > 0 && inputLength != size) || (size-1)/resolution+1 != n) {
            if (chrom == "chrM") { result.lengthMismatch = true; return result; }
            throw std::runtime_error("Hi-C/reference chromosome length mismatch for " + chrom +
                "; verify --assembly or supply a matching --phasing-track");
        }
        if (wantGC && (result.gcExact || approximate)) {
            const auto &data = counts(chrom);
            for (int b = 0; b < n; ++b) {
                int64_t start = int64_t(b)*resolution, end = std::min(size, start+resolution);
                long double gc = 0, valid = 0;
                for (int64_t tile = start/100; tile <= (end-1)/100; ++tile) {
                    int64_t a = tile*100, z = std::min(size, a+100);
                    long double weight = static_cast<long double>(std::min(end,z)-std::max(start,a)) / (z-a);
                    gc += weight * data[2*tile]; valid += weight * data[2*tile+1];
                }
                if (valid > 0 && valid/(end-start) >= 0.5L) result.gc[b] = static_cast<double>(gc/valid);
            }
        }
        if (wantGenes) {
            loadGenes();
            const auto &spans = genes_[chrom];
            size_t cursor = 0;
            for (int b = 0; b < n; ++b) {
                int64_t start = int64_t(b)*resolution, end = std::min(size,start+resolution), covered = 0;
                while (cursor < spans.size() && spans[cursor].second <= start) ++cursor;
                for (size_t j = cursor; j < spans.size() && spans[j].first < end; ++j)
                    covered += std::min(end,spans[j].second)-std::max(start,spans[j].first);
                result.genes[b] = static_cast<double>(covered)/(end-start);
            }
        }
        return result;
    }
};
} // namespace possum
#endif
