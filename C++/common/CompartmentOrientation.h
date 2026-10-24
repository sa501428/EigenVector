#ifndef POSSUM_COMPARTMENT_ORIENTATION_H
#define POSSUM_COMPARTMENT_ORIENTATION_H
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "PhasingOptions.h"
#include "BundledPhasing.h"

namespace possum {
inline std::string chromosomeKey(std::string name) {
    if (name.size() >= 3 && (name.substr(0,3) == "chr" || name.substr(0,3) == "CHR")) name.erase(0,3);
    if (name == "MT" || name == "mt" || name == "m") name = "M";
    if (name == "x") name = "X";
    if (name == "y") name = "Y";
    return name;
}
inline std::string displayChromosome(const std::string &name) { return "chr" + chromosomeKey(name); }
inline bool isAllChromosome(const std::string &name) { return name == "ALL" || name == "All" || name == "all"; }
inline int binCount(int64_t length, int resolution) {
    if (length <= 0 || resolution <= 0 || (length-1)/resolution >= std::numeric_limits<int>::max())
        throw std::runtime_error("invalid chromosome length/resolution or too many bins");
    return static_cast<int>((length-1)/resolution+1);
}
inline std::vector<std::string> assemblyCandidates(std::string name) {
    std::transform(name.begin(),name.end(),name.begin(),[](unsigned char c) { return std::tolower(c); });
    std::vector<std::string> result;
    if (name.find("hg19") != std::string::npos || name.find("grch37") != std::string::npos || name.find("b37") != std::string::npos) result.push_back("hg19");
    if (name.find("hg38") != std::string::npos || name.find("grch38") != std::string::npos) result.push_back("hg38");
    if (name.find("mm9") != std::string::npos || name.find("ncbi37") != std::string::npos || name.find("ncbim37") != std::string::npos) result.push_back("mm9");
    if (name.find("mm10") != std::string::npos || name.find("grcm38") != std::string::npos) result.push_back("mm10");
    return result;
}
inline std::string assemblyKey(const std::string &name) {
    auto candidates = assemblyCandidates(name);
    return candidates.size() == 1 ? candidates.front() : name;
}
inline bool bundledAssembly(const std::string &name) {
    return name == "hg19" || name == "hg38" || name == "mm9" || name == "mm10";
}
using PhasingChromosomes = std::map<std::string,int64_t>;
template<class ChromosomeMap>
inline PhasingChromosomes phasingChromosomes(const ChromosomeMap &chromosomes) {
    PhasingChromosomes sizes;
    for (const auto &item : chromosomes)
        if (!isAllChromosome(item.second.name)) sizes[item.second.name] = item.second.length;
    return sizes;
}
inline PhasingChromosomes phasingChromosomes(const std::vector<std::string> &names,const std::vector<int64_t> &lengths) {
    if (names.size() != lengths.size()) throw std::runtime_error("chromosome names/sizes count mismatch");
    PhasingChromosomes sizes;
    for (size_t i=0;i<names.size();++i) if (!isAllChromosome(names[i])) sizes[names[i]] = lengths[i];
    return sizes;
}
struct PhaseInterval { int64_t start, end; double value; };
struct PhaseSlice { const double *values; int n; std::string chrom; int64_t length; };
struct PhaseCorrelation {
    double spearman = std::numeric_limits<double>::quiet_NaN();
    double pearson = std::numeric_limits<double>::quiet_NaN();
    size_t bins = 0;
};
struct PhaseDecision {
    int sign = 0;
    std::string source = "none", status = "weak_evidence";
    PhaseCorrelation gc, genes, custom;
    bool gcExact = true;
};
inline double phasePearson(const std::vector<double> &x, const std::vector<double> &y) {
    long double mx = 0, my = 0, xx = 0, yy = 0, xy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        long double dx = x[i]-mx, dy = y[i]-my;
        mx += dx/(i+1); my += dy/(i+1);
        xx += dx*(x[i]-mx); yy += dy*(y[i]-my); xy += dx*(y[i]-my);
    }
    if (x.size() < 2 || xx <= 0 || yy <= 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(std::max(-1.0L,std::min(1.0L,xy/std::sqrt(xx*yy))));
}
inline std::vector<double> phaseRanks(const std::vector<double> &values) {
    std::vector<size_t> order(values.size()); std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b) { return values[a] < values[b]; });
    std::vector<double> ranks(values.size());
    for (size_t i = 0; i < order.size();) {
        size_t j = i+1;
        while (j < order.size() && values[order[j]] == values[order[i]]) ++j;
        for (size_t k = i; k < j; ++k) ranks[order[k]] = (i+j-1)/2.0;
        i = j;
    }
    return ranks;
}
inline PhaseCorrelation phaseCorrelation(const std::vector<double> &x,const std::vector<double> &y) {
    PhaseCorrelation result;
    result.bins = x.size(); result.pearson = phasePearson(x,y);
    result.spearman = phasePearson(phaseRanks(x),phaseRanks(y));
    return result;
}

class CompartmentOrientation {
    std::string genome_;
    std::string assemblySource_, referencePath_;
    PhasingOptions options_;
    bool custom_;
    std::map<std::string,std::vector<PhaseInterval> > intervals_;
    std::unique_ptr<BundledPhasing> bundle_;
    mutable std::map<std::string,PhasingBins> binned_;
    mutable std::ofstream report_;

    struct SizeMatch { bool compatible = true; size_t matched = 0; };
    static SizeMatch matchSizes(const BundledPhasing &reference,const PhasingChromosomes &chromosomes) {
        SizeMatch result;
        for (const auto &item : chromosomes) {
            std::string key = chromosomeKey(item.first);
            // Mitochondrial references vary independently of the nuclear build
            // (notably b37 vs UCSC hg19). Infer assembly from nuclear chromosomes.
            bool canonical = key == "X" || key == "Y";
            for (int i=1;i<=22 && !canonical;++i) canonical = key == std::to_string(i);
            // Random, alternate and unplaced scaffolds are not in these bundles.
            if (!canonical) continue;
            if (reference.length("chr"+key) != item.second) { result.compatible = false; return result; }
            ++result.matched;
        }
        return result;
    }

    void loadTrack(const std::string &path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot read phasing track: " + path);
        std::string line;
        size_t number = 0;
        while (std::getline(in,line)) {
            ++number;
            std::istringstream row(line);
            std::string chrom,text;
            PhaseInterval v;
            if (!(row >> chrom) || chrom[0] == '#' || chrom == "track" || chrom == "browser") continue;
            if (!(row >> v.start >> v.end >> text) || v.start < 0 || v.end <= v.start)
                throw std::runtime_error("invalid bedGraph interval at " + path + ":" + std::to_string(number));
            char *end = NULL; v.value = std::strtod(text.c_str(),&end);
            if (end == text.c_str() || *end) throw std::runtime_error("invalid bedGraph value at " + path + ":" + std::to_string(number));
            // Missing intervals still participate in overlap validation.
            intervals_[chromosomeKey(chrom)].push_back(v);
        }
        if (!in.eof()) throw std::runtime_error("failed reading phasing track: " + path);
        size_t finite = 0;
        for (auto &item : intervals_) {
            auto &spans = item.second;
            std::sort(spans.begin(),spans.end(),[](const PhaseInterval &a,const PhaseInterval &b) { return a.start < b.start; });
            for (size_t i = 0; i < spans.size(); ++i) {
                if (i && spans[i].start < spans[i-1].end)
                    throw std::runtime_error("overlapping phasing intervals on " + item.first);
                finite += std::isfinite(spans[i].value);
            }
        }
        if (!finite) throw std::runtime_error("phasing track contains no finite intervals: " + path);
    }
    const PhasingBins &bins(const PhaseSlice &slice,int resolution) const {
        std::string key = chromosomeKey(slice.chrom) + ":" + std::to_string(resolution) + ":" +
                          std::to_string(slice.n) + ":" + std::to_string(slice.length);
        auto found = binned_.find(key);
        if (found != binned_.end()) return found->second;
        PhasingBins result;
        if (custom_) {
            result.custom.assign(slice.n,std::numeric_limits<double>::quiet_NaN());
            auto item = intervals_.find(chromosomeKey(slice.chrom));
            if (item != intervals_.end()) {
                int64_t limit = slice.length > 0 ? std::min(slice.length,int64_t(slice.n)*resolution) : int64_t(slice.n)*resolution;
                std::vector<long double> sums(slice.n,0),covered(slice.n,0);
                for (const auto &v : item->second) {
                    if (!std::isfinite(v.value)) continue;
                    int64_t start = v.start,end = std::min(v.end,limit);
                    while (start < end) {
                        int b = static_cast<int>(start/resolution);
                        int64_t stop = std::min(end,(int64_t(b)+1)*resolution);
                        sums[b] += (stop-start)*static_cast<long double>(v.value);
                        covered[b] += stop-start; start = stop;
                    }
                }
                for (int b = 0; b < slice.n; ++b) if (covered[b] > 0) result.custom[b] = static_cast<double>(sums[b]/covered[b]);
            }
        } else if (bundle_) {
            result = bundle_->bin(displayChromosome(slice.chrom),slice.n,resolution,slice.length,
                options_.method != "refseq",options_.method != "gc",options_.approximateGC);
        }
        return binned_.emplace(key,std::move(result)).first->second;
    }
    bool strong(const PhaseCorrelation &c) const {
        return c.bins >= static_cast<size_t>(options_.minBins) && std::isfinite(c.spearman) &&
               c.spearman != 0 && std::fabs(c.spearman) >= options_.minCorrelation;
    }
public:
    explicit CompartmentOrientation(std::string genome,const std::string &track = "",
                                    const PhasingOptions &options = PhasingOptions(),
                                    const PhasingChromosomes &chromosomes = PhasingChromosomes())
        : genome_(assemblyKey(options.assembly.empty() ? (genome.empty() ? "hg38" : genome) : options.assembly)),
          options_(options),custom_(!track.empty()) {
        assemblySource_ = !options.assembly.empty() ? "override" : genome.empty() ? "default_hg38" : "map_header";
        if (custom_ && options.method != "auto")
            throw std::runtime_error("--phasing-track requires --phasing-method auto (custom track overrides bundled references)");
        if (!options.assembly.empty() && !bundledAssembly(genome_) && !custom_)
            throw std::runtime_error("unsupported --assembly: " + options.assembly);
        if (options.method == "none") return;
        if (custom_) { referencePath_ = track; loadTrack(track); return; }
        std::string root = options.directory;
        if (root.empty()) { const char *env = std::getenv("POSSUM_PHASING_DIR"); if (env) root = env; }
#ifdef POSSUM_PHASING_DIR
        if (root.empty()) root = POSSUM_PHASING_DIR;
#endif
        if (root.empty()) root = "resources/compartment_phasing";
        if (!options.assembly.empty()) {
            referencePath_ = root + "/" + genome_;
            bundle_.reset(new BundledPhasing(referencePath_));
            if (!matchSizes(*bundle_,chromosomes).compatible)
                throw std::runtime_error("--assembly " + genome_ + " conflicts with Hi-C chromosome names/sizes");
            return;
        }
        const auto labels = assemblyCandidates(genome);
        if (bundledAssembly(genome_) && (labels.size() == 1 || chromosomes.empty())) {
            referencePath_ = root + "/" + genome_;
            bundle_.reset(new BundledPhasing(referencePath_));
            if (matchSizes(*bundle_,chromosomes).compatible) return;
            bundle_.reset();
        }
        if (chromosomes.empty()) return;
        std::vector<std::string> matches;
        std::unique_ptr<BundledPhasing> matching;
        for (const std::string name : {"hg19","hg38","mm9","mm10"}) {
            std::unique_ptr<BundledPhasing> candidate(new BundledPhasing(root+"/"+name));
            SizeMatch signature = matchSizes(*candidate,chromosomes);
            if (signature.compatible && signature.matched) {
                matches.push_back(name); matching = std::move(candidate);
            }
        }
        if (matches.size() == 1) {
            if (labels.size() == 1 && labels.front() != matches.front())
                std::fprintf(stderr,"Genome ID %s suggests %s, but Hi-C chromosome names/sizes identify %s; using %s\n",
                    genome.c_str(),labels.front().c_str(),matches.front().c_str(),matches.front().c_str());
            genome_ = matches.front(); assemblySource_ = "chromosome_sizes";
            referencePath_ = root+"/"+genome_; bundle_ = std::move(matching);
        } else if (labels.size() == 1) {
            throw std::runtime_error("genome ID " + genome + " conflicts with Hi-C chromosome names/sizes; no unique bundled assembly matches");
        } else {
            genome_ = genome.empty() ? "unknown" : genome;
            referencePath_.clear();
            assemblySource_ = matches.empty() ? "unresolved_chromosome_sizes" : "ambiguous_chromosome_sizes";
        }
    }
    const std::string &assembly() const { return genome_; }
    void setReport(const std::string &path) {
        report_.open(path);
        if (!report_) throw std::runtime_error("cannot create orientation report: " + path);
        report_ << "assembly\tchromosome\tresolution\teigenvector\tstatus\tsource\tapplied_sign"
                   "\tgc_spearman\tgc_pearson\tgc_bins\trefseq_spearman\trefseq_pearson\trefseq_bins"
                   "\ttrack_spearman\ttrack_pearson\ttrack_bins\tgc_exact\tmin_bins\tmin_abs_spearman"
                   "\tdecision_sign\tassembly_source\treference_path\n" << std::setprecision(10);
    }
    PhaseDecision evaluate(const std::vector<PhaseSlice> &slices,int resolution) const {
        PhaseDecision result;
        if (options_.method == "none") { result.status = "disabled"; return result; }
        if (!custom_ && !bundle_) { result.status = "unsupported_assembly"; return result; }
        if (resolution <= 0) { result.status = "invalid_resolution"; return result; }
        if (!custom_ && options_.method == "gc" && resolution % 100 && !options_.approximateGC)
            throw std::runtime_error("exact bundled GC requires a resolution divisible by 100 bp; use auto/refseq, a custom exact track, or --phasing-allow-approximate-gc");
        std::vector<double> gx,gy,rx,ry,tx,ty;
        bool lengthMismatch = false;
        for (const auto &slice : slices) {
            if (!slice.values || slice.n <= 0) continue;
            const auto &reference = bins(slice,resolution);
            lengthMismatch = lengthMismatch || reference.lengthMismatch;
            result.gcExact = result.gcExact && reference.gcExact;
            for (int i = 0; i < slice.n; ++i) {
                // These solvers represent masked/constant contact rows as zero.
                double x = slice.values[i];
                if (!std::isfinite(x) || x == 0) continue;
                if (!reference.gc.empty() && std::isfinite(reference.gc[i])) { gx.push_back(x); gy.push_back(reference.gc[i]); }
                if (!reference.genes.empty() && std::isfinite(reference.genes[i])) { rx.push_back(x); ry.push_back(reference.genes[i]); }
                if (!reference.custom.empty() && std::isfinite(reference.custom[i])) { tx.push_back(x); ty.push_back(reference.custom[i]); }
            }
        }
        result.gc = phaseCorrelation(gx,gy); result.genes = phaseCorrelation(rx,ry); result.custom = phaseCorrelation(tx,ty);
        double correlation = 0;
        if (custom_) {
            if (strong(result.custom)) { correlation = result.custom.spearman; result.source = "custom_track"; }
        } else if (options_.method == "gc") {
            if (strong(result.gc)) { correlation = result.gc.spearman; result.source = "bundled_gc"; }
        } else if (options_.method == "refseq") {
            if (strong(result.genes)) { correlation = result.genes.spearman; result.source = "bundled_refseq"; }
        } else {
            if (strong(result.gc) && strong(result.genes) && result.gc.spearman*result.genes.spearman < 0) {
                result.status = "conflicting_references"; return result;
            }
            if (strong(result.gc)) { correlation = result.gc.spearman; result.source = "bundled_gc"; }
            else if (strong(result.genes)) { correlation = result.genes.spearman; result.source = "bundled_refseq"; }
        }
        if (correlation) { result.sign = correlation > 0 ? 1 : -1; result.status = "oriented"; }
        else if (!result.gc.bins && !result.genes.bins && !result.custom.bins)
            result.status = lengthMismatch ? "chromosome_length_mismatch" : "no_reference_overlap";
        return result;
    }
    void record(const PhaseDecision &d,const std::string &chrom,int resolution,int eigenvector,int verbosity = 0) const {
        if (report_.is_open()) {
            report_ << genome_ << '\t' << chrom << '\t' << resolution << '\t' << eigenvector << '\t'
                    << d.status << '\t' << d.source << '\t' << (d.sign < 0 ? -1 : 1) << '\t'
                    << d.gc.spearman << '\t' << d.gc.pearson << '\t' << d.gc.bins << '\t'
                    << d.genes.spearman << '\t' << d.genes.pearson << '\t' << d.genes.bins << '\t'
                    << d.custom.spearman << '\t' << d.custom.pearson << '\t' << d.custom.bins << '\t'
                    << (d.gcExact ? 1 : 0) << '\t' << options_.minBins << '\t' << options_.minCorrelation
                    << '\t' << d.sign << '\t' << assemblySource_ << '\t' << referencePath_ << '\n';
            report_.flush();
            if (!report_) throw std::runtime_error("failed writing orientation report");
        }
        if (verbosity || (!d.sign && d.status != "disabled"))
            std::fprintf(stderr,"A/B orientation %s EV%d: %s; assembly=%s source=%s status=%s GC_rho=%.4g RefSeq_rho=%.4g track_rho=%.4g\n",
                chrom.c_str(),eigenvector,d.sign > 0 ? "kept" : d.sign < 0 ? "flipped" : "unresolved; sign unchanged",
                genome_.c_str(),d.source.c_str(),d.status.c_str(),d.gc.spearman,d.genes.spearman,d.custom.spearman);
    }
    long double evidence(const double *x,int n,const std::string &chrom,int resolution,int64_t length = 0) const {
        return evaluate({{x,n,chrom,length}},resolution).sign;
    }
    int orient(double *x,int n,const std::string &chrom,int resolution,
               int verbosity = 0,int eigenvector = 1,int64_t length = 0) const {
        PhaseDecision d = evaluate({{x,n,chrom,length}},resolution);
        if (d.sign < 0) for (int i = 0; i < n; ++i) x[i] = -x[i];
        record(d,chrom,resolution,eigenvector,verbosity);
        return d.sign;
    }
    int orientGlobal(const std::vector<PhaseSlice> &slices,double *x,size_t n,int resolution,
                     int verbosity,int eigenvector) const {
        PhaseDecision d = evaluate(slices,resolution);
        if (d.sign < 0) for (size_t i = 0; i < n; ++i) x[i] = -x[i];
        record(d,"genome",resolution,eigenvector,verbosity);
        return d.sign;
    }
};
} // namespace possum
#endif
