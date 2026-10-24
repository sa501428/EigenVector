#ifndef POSSUM_PHASING_OPTIONS_H
#define POSSUM_PHASING_OPTIONS_H
#include "Arguments.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// Shared long options and validation for every Hi-C entry point.
#define POSSUM_PHASING_LONG_OPTIONS \
    {"assembly", required_argument, NULL, 2100}, \
    {"phasing-dir", required_argument, NULL, 2101}, \
    {"phasing-method", required_argument, NULL, 2102}, \
    {"phasing-min-correlation", required_argument, NULL, 2103}, \
    {"phasing-min-bins", required_argument, NULL, 2104}, \
    {"phasing-allow-approximate-gc", no_argument, NULL, 2105},

namespace possum {
struct PhasingOptions {
    std::string assembly, directory, method = "auto";
    double minCorrelation = 0.1;
    int minBins = 20;
    bool approximateGC = false;

    bool parse(int option, const char *value) {
        switch (option) {
        case 2100: assembly = value; break;
        case 2101: directory = value; break;
        case 2102:
            method = value;
            if (method != "auto" && method != "gc" && method != "refseq" && method != "none")
                throw std::runtime_error("--phasing-method must be auto, gc, refseq, or none");
            break;
        case 2103: {
            char *end = NULL;
            minCorrelation = std::strtod(value, &end);
            if (end == value || *end || !std::isfinite(minCorrelation) || minCorrelation < 0 || minCorrelation > 1)
                throw std::runtime_error("--phasing-min-correlation must be in [0,1]");
            break;
        }
        case 2104:
            minBins = positiveInt(value);
            if (minBins < 2) throw std::runtime_error("--phasing-min-bins must be at least 2");
            break;
        case 2105: approximateGC = true; break;
        default: return false;
        }
        return true;
    }

    std::vector<std::string> arguments() const {
        std::ostringstream threshold;
        threshold << std::setprecision(17) << minCorrelation;
        std::vector<std::string> args = {"--phasing-method", method,
            "--phasing-min-correlation", threshold.str(),
            "--phasing-min-bins", std::to_string(minBins)};
        if (!assembly.empty()) { args.push_back("--assembly"); args.push_back(assembly); }
        if (!directory.empty()) { args.push_back("--phasing-dir"); args.push_back(directory); }
        if (approximateGC) args.push_back("--phasing-allow-approximate-gc");
        return args;
    }
};

inline void phasingUsage() {
    std::fprintf(stderr,
        "\nPhasing options (each Hi-C output also writes an orientation TSV):\n"
        "  --assembly NAME     Override assembly detected from genome ID/chromosome sizes (hg19/hg38/mm9/mm10)\n"
        "  --phasing-dir DIR   Bundled references (or POSSUM_PHASING_DIR)\n"
        "  --phasing-method M  auto (GC + RefSeq check), gc, refseq, or none\n"
        "  --phasing-min-correlation X\n"
        "                      Minimum absolute Spearman correlation (default: 0.1)\n"
        "  --phasing-min-bins N Minimum paired valid bins (default: 20)\n"
        "  --phasing-allow-approximate-gc\n"
        "                      Allow boundary-tile GC estimates off the 100 bp grid\n");
}
} // namespace possum
#endif
