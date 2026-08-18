#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#ifndef _WIN32
#include <sys/wait.h>
#include <unistd.h>
#else
#include <direct.h>
#include <io.h>
#include <process.h>
#endif
#include <vector>

#include <straw.h>

using namespace std;

namespace {

const int kDefaultResolution = 5000;
const int kLowResolution = 50000;
const int kRescueEigenvectors = 10;

struct Options {
    string normalization = "SCALE";
    bool observed = false;
    bool rescue = false;
    bool keepTemp = false;
    int resolution = kDefaultResolution;
    int jobs = 0;
    int solverThreads = 1;
    int maxIterations = 200;
    int verbosity = 1;
    double tolerance = 1.0e-7;
    double epsilon = 1.0e-8;
    double rescueThreshold = 0.8;
    bool coverageFilter = true;
    double coverageZCutoff = -2.5;
    double maxTopOnePercentEnergy = 0.5;
};

struct ChromosomeInfo {
    string name;
    long length;
    string base;
};

struct Job {
    string label;
    string log;
    vector<string> arguments;
};

struct EigenDiagnostic {
    int firstK;
    double lambdaK;
    double lambdaKPlus1;
    double lambdaKPlus2;
    double gapK;
    double gapKPlus1;
    double ratioK;
    double ratioKPlus1;
    double minRatio;
    double tolerance;
    double estimatedRelativeError;
    bool coverageFilterEnabled;
    double coverageZCutoff;
    double coverageMeanLog;
    double coverageSdLog;
    double coverageMinimumEntries;
    long coverageRemovedBins;
    long coverageRetainedBins;
};

void usage(const char *program) {
    fprintf(stderr,
        "Usage: %s [options] <hicfile> <outbase> [resolution]\n\n"
        "Runs chromosome-specific Lanczos calculations in parallel and merges them.\n\n"
        "Options:\n"
        "  --rescue            Enable 50 kb-guided EV1-EV10 rescue\n"
        "  -r, --resolution N  Output resolution (default: 5000)\n"
        "  -n, --normalization NAME\n"
        "                      Hi-C normalization (default: SCALE)\n"
        "  -o, --observed      Use observed instead of observed/expected\n"
        "  -j, --jobs N        Concurrent chromosomes (default: min(CPUs, 4))\n"
        "  -T, --threads N     Threads inside each chromosome solver (default: 1)\n"
        "  -t, --tolerance X   Lanczos tolerance (default: 1.0e-7)\n"
        "  -e, --epsilon X     Lanczos epsilon (default: 1.0e-8)\n"
        "  -I, --max-iter N    Maximum Lanczos iterations (default: 200)\n"
        "      --threshold X   Rescue absolute-correlation threshold (default: 0.8)\n"
        "      --coverage-z X  Filter high-resolution bins below this log-count z-score (default: -2.5)\n"
        "      --no-coverage-filter\n"
        "                      Disable high-resolution low-coverage filtering\n"
        "      --max-top1-energy X\n"
        "                      Reject localized EVs above this fraction (default: 0.5)\n"
        "      --keep-temp     Retain per-chromosome files after success\n"
        "  -v, --verbose N     Verbosity (default: 1)\n"
        "  -h, --help          Show this help\n\n"
        "Normal output: <outbase>.wig\n"
        "Rescue outputs: <outbase>.wig, .lowres_50kb.wig, .Ev1.wig-.Ev10.wig,\n"
        "                .eigenvalues.tsv, and .rescue.tsv\n",
        program);
}

string toString(int value) {
    ostringstream out;
    out << value;
    return out.str();
}

string toString(double value) {
    ostringstream out;
    out << setprecision(17) << value;
    return out.str();
}

string directoryName(const string &path) {
    size_t slash = path.find_last_of('/');
    if (slash == string::npos) return ".";
    if (slash == 0) return "/";
    return path.substr(0, slash);
}

bool directoryExists(const string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) return false;
#ifdef _WIN32
    return (info.st_mode & S_IFDIR) != 0;
#else
    return S_ISDIR(info.st_mode);
#endif
}

void createDirectory(const string &path) {
    if (path.empty() || path == "." || directoryExists(path)) return;
#ifdef _WIN32
    if (_mkdir(path.c_str()) != 0 && errno != EEXIST)
#else
    if (mkdir(path.c_str(), 0777) != 0 && errno != EEXIST)
#endif
        throw runtime_error("cannot create output directory " + path + ": " +
                            string(strerror(errno)));
    if (!directoryExists(path))
        throw runtime_error("output parent exists but is not a directory: " + path);
}

void ensureOutputDirectory(const string &outbase) {
    string parent = directoryName(outbase);
    if (parent == "." || directoryExists(parent)) return;

    string current;
    size_t position = 0;
    if (!parent.empty() && parent[0] == '/') {
        current = "/";
        position = 1;
    }
#ifdef _WIN32
    if (parent.size() >= 2 && parent[1] == ':') {
        current = parent.substr(0, 2);
        position = 2;
        if (position < parent.size() && (parent[position] == '/' || parent[position] == '\\')) {
            current += parent[position++];
        }
    }
#endif
    while (position <= parent.size()) {
        size_t slash = parent.find_first_of("/\\", position);
        string component = parent.substr(position, slash == string::npos
                                                       ? string::npos
                                                       : slash - position);
        if (!component.empty()) {
            if (!current.empty() && current[current.size() - 1] != '/' &&
                current[current.size() - 1] != '\\') current += '/';
            current += component;
            createDirectory(current);
        }
        if (slash == string::npos) break;
        position = slash + 1;
    }
}

bool regularExecutable(const string &path) {
#ifdef _WIN32
    return _access(path.c_str(), 0) == 0;
#else
    return access(path.c_str(), X_OK) == 0;
#endif
}

bool readable(const string &path) {
#ifdef _WIN32
    return _access(path.c_str(), 4) == 0;
#else
    return access(path.c_str(), R_OK) == 0;
#endif
}

int processId() {
#ifdef _WIN32
    return _getpid();
#else
    return static_cast<int>(getpid());
#endif
}

string naturalKey(const string &name) {
    string key = name;
    if (key.size() >= 3 && (key.substr(0, 3) == "chr" || key.substr(0, 3) == "CHR"))
        key = key.substr(3);
    bool numeric = !key.empty();
    for (size_t i = 0; i < key.size(); ++i)
        if (!isdigit(static_cast<unsigned char>(key[i]))) numeric = false;
    if (numeric) {
        ostringstream out;
        out << '0' << setw(9) << setfill('0') << atoi(key.c_str());
        return out.str();
    }
    return "1" + key;
}

bool chromosomeLess(const ChromosomeInfo &a, const ChromosomeInfo &b) {
    string ak = naturalKey(a.name);
    string bk = naturalKey(b.name);
    if (ak != bk) return ak < bk;
    return a.name < b.name;
}

string displayChromosome(const string &name) {
    string result = name;
    if (result.compare(0, 3, "chr") != 0) result = "chr" + result;
    if (result == "chrMT") result = "chrM";
    return result;
}

vector<ChromosomeInfo> readChromosomes(const string &hicfile) {
    ifstream fin(hicfile.c_str(), ios::in | ios::binary);
    if (!fin) throw runtime_error("cannot open Hi-C file: " + hicfile);

    int64_t master = 0;
    string genomeID;
    int numChromosomes = 0;
    int version = 0;
    int64_t nviPosition = 0;
    int64_t nviLength = 0;
    map<string, chromosome> chromosomeMap =
        readHeader(fin, master, genomeID, numChromosomes, version,
                   nviPosition, nviLength);

    vector<ChromosomeInfo> result;
    for (map<string, chromosome>::const_iterator it = chromosomeMap.begin();
         it != chromosomeMap.end(); ++it) {
        if (it->second.name == "ALL" || it->second.name == "All") continue;
        ChromosomeInfo chr;
        chr.name = it->second.name;
        chr.length = it->second.length;
        result.push_back(chr);
    }
    sort(result.begin(), result.end(), chromosomeLess);
    return result;
}

string makeTempDirectory(const string &outbase) {
#ifdef _WIN32
    for (int attempt = 0; attempt < 1000; ++attempt) {
        string candidate = outbase + ".lan_tmp." + toString(processId()) + "." +
                           toString(attempt);
        if (_mkdir(candidate.c_str()) == 0) return candidate;
        if (errno != EEXIST) break;
    }
    throw runtime_error("cannot create temporary directory near output: " +
                        string(strerror(errno)));
#else
    string pattern = outbase + ".lan_tmp.XXXXXX";
    vector<char> writable(pattern.begin(), pattern.end());
    writable.push_back('\0');
    char *created = mkdtemp(&writable[0]);
    if (created == NULL)
        throw runtime_error("cannot create temporary directory near output: " +
                            string(strerror(errno)));
    return string(created);
#endif
}

#ifndef _WIN32
pid_t launch(const Job &job) {
    pid_t pid = fork();
    if (pid < 0) return pid;
    if (pid == 0) {
        FILE *log = fopen(job.log.c_str(), "w");
        if (log != NULL) {
            dup2(fileno(log), STDOUT_FILENO);
            dup2(fileno(log), STDERR_FILENO);
            fclose(log);
        }
        vector<char *> argv;
        for (size_t i = 0; i < job.arguments.size(); ++i)
            argv.push_back(const_cast<char *>(job.arguments[i].c_str()));
        argv.push_back(NULL);
        execv(argv[0], &argv[0]);
        _exit(127);
    }
    return pid;
}

bool runJobs(const vector<Job> &jobs, int parallelism, int verbosity,
             bool tolerateFailures = false) {
    map<pid_t, size_t> active;
    size_t next = 0;
    bool ok = true;
    bool anyFailure = false;

    while (next < jobs.size() || !active.empty()) {
        while (ok && next < jobs.size() &&
               static_cast<int>(active.size()) < parallelism) {
            pid_t pid = launch(jobs[next]);
            if (pid < 0) {
                cerr << "Unable to start " << jobs[next].label << ": "
                     << strerror(errno) << endl;
                ok = false;
                break;
            }
            active[pid] = next++;
        }

        if (active.empty()) break;
        int status = 0;
        pid_t finished = waitpid(-1, &status, 0);
        if (finished < 0) {
            if (errno == EINTR) continue;
            cerr << "waitpid failed: " << strerror(errno) << endl;
            return false;
        }
        size_t index = active[finished];
        active.erase(finished);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            cerr << "Chromosome job failed: " << jobs[index].label
                 << " (log: " << jobs[index].log << ")" << endl;
            anyFailure = true;
            if (!tolerateFailures) ok = false;
        } else if (verbosity > 0) {
            cerr << "Completed " << jobs[index].label << " ("
                 << (index + 1) << "/" << jobs.size() << ")" << endl;
        }
    }
    if (tolerateFailures && anyFailure && verbosity > 0)
        cerr << "Continuing: sparse chromosome failures will be zero-filled or explicitly reported" << endl;
    return ok && next == jobs.size();
}
#else
string windowsQuote(const string &value) {
    string result = "\"";
    for (size_t i = 0; i < value.size(); ++i) {
        if (value[i] == '\"') result += "\\\"";
        else result += value[i];
    }
    return result + "\"";
}

bool runJobs(const vector<Job> &jobs, int parallelism, int verbosity,
             bool tolerateFailures = false) {
    atomic<size_t> next(0);
    atomic<bool> failed(false);
    mutex outputMutex;
    int workerCount = min(parallelism, static_cast<int>(jobs.size()));
    vector<thread> workers;
    for (int workerIndex = 0; workerIndex < workerCount; ++workerIndex) {
        workers.push_back(thread([&]() {
            while (true) {
                if (failed.load() && !tolerateFailures) return;
                size_t index = next.fetch_add(1);
                if (index >= jobs.size()) return;
                string command;
                for (size_t i = 0; i < jobs[index].arguments.size(); ++i) {
                    if (!command.empty()) command += ' ';
                    command += windowsQuote(jobs[index].arguments[i]);
                }
                command += " > " + windowsQuote(jobs[index].log) + " 2>&1";
                int status = system(command.c_str());
                lock_guard<mutex> lock(outputMutex);
                if (status != 0) {
                    cerr << "Chromosome job failed: " << jobs[index].label
                         << " (log: " << jobs[index].log << ")" << endl;
                    failed.store(true);
                } else if (verbosity > 0) {
                    cerr << "Completed " << jobs[index].label << " ("
                         << (index + 1) << "/" << jobs.size() << ")" << endl;
                }
            }
        }));
    }
    for (size_t i = 0; i < workers.size(); ++i) workers[i].join();
    if (tolerateFailures && failed.load() && verbosity > 0)
        cerr << "Continuing: sparse chromosome failures will be zero-filled or explicitly reported" << endl;
    return tolerateFailures || !failed.load();
}
#endif

bool onlyExpectedSparseFailures(const vector<Job> &jobs) {
    for (size_t i = 0; i < jobs.size(); ++i) {
        if (jobs[i].arguments.size() < 3) return false;
        string base = jobs[i].arguments[jobs[i].arguments.size() - 3];
        if (readable(base + ".Ev1.wig")) continue;
        ifstream log(jobs[i].log.c_str());
        string contents((istreambuf_iterator<char>(log)), istreambuf_iterator<char>());
        // LAPACK info=5 is the observed failure for a matrix with too few
        // independent entries to produce the requested eigensystem.
        if (contents.find("return code is 100005") == string::npos) return false;
    }
    return true;
}

vector<string> workerArguments(const string &worker, const Options &options,
                               const string &hicfile, const string &chromosome,
                               const string &base, int resolution, int nv,
                               bool applyCoverageFilter) {
    vector<string> args;
    args.push_back(worker);
    if (options.observed) args.push_back("-o");
    args.push_back("-n"); args.push_back(options.normalization);
    args.push_back("-t"); args.push_back(toString(options.tolerance));
    args.push_back("-e"); args.push_back(toString(options.epsilon));
    args.push_back("-I"); args.push_back(toString(options.maxIterations));
    args.push_back("-T"); args.push_back(toString(options.solverThreads));
    if (applyCoverageFilter) {
        args.push_back("-c"); args.push_back(toString(options.coverageZCutoff));
    } else {
        args.push_back("-C");
    }
    args.push_back("-v"); args.push_back("0");
    args.push_back(hicfile);
    args.push_back(chromosome);
    args.push_back(base);
    args.push_back(toString(resolution));
    args.push_back(toString(nv));
    return args;
}

vector<Job> makeJobs(vector<ChromosomeInfo> &chromosomes, const string &worker,
                     const Options &options, const string &hicfile,
                     const string &tempDirectory, const string &phase,
                     int resolution, int nv) {
    vector<Job> jobs;
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        ostringstream stem;
        stem << tempDirectory << '/' << phase << "_chr_"
             << setw(3) << setfill('0') << i;
        string base = stem.str();
        if (phase == "high" || phase == "normal") chromosomes[i].base = base;
        Job job;
        job.label = phase + " " + chromosomes[i].name;
        job.log = base + ".log";
        job.arguments = workerArguments(worker, options, hicfile,
                                        chromosomes[i].name, base,
                                        resolution, nv,
                                        options.coverageFilter && phase != "low");
        jobs.push_back(job);
    }
    return jobs;
}

vector<double> readWig(const string &path) {
    ifstream in(path.c_str());
    if (!in) throw runtime_error("cannot read worker WIG: " + path);
    vector<double> values;
    string line;
    while (getline(in, line)) {
        if (line.empty() || line.compare(0, 5, "track") == 0 ||
            line.compare(0, 9, "fixedStep") == 0) continue;
        char *end = NULL;
        double value = strtod(line.c_str(), &end);
        if (end == line.c_str())
            throw runtime_error("invalid value in worker WIG: " + path);
        values.push_back(value);
    }
    return values;
}

vector<string> splitTabs(const string &line) {
    vector<string> fields;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        fields.push_back(line.substr(start, tab == string::npos ? string::npos
                                                               : tab - start));
        if (tab == string::npos) break;
        start = tab + 1;
    }
    return fields;
}

vector<EigenDiagnostic> readEigenDiagnostics(const string &path) {
    ifstream in(path.c_str());
    if (!in) throw runtime_error("cannot read eigenvalues: " + path);
    vector<EigenDiagnostic> values;
    string line;
    getline(in, line);
    while (getline(in, line)) {
        vector<string> fields = splitTabs(line);
        if (fields.size() < 18)
            throw runtime_error("invalid eigenvalue diagnostics: " + path);
        EigenDiagnostic value;
        value.firstK = atoi(fields[0].c_str());
        value.lambdaK = strtod(fields[1].c_str(), NULL);
        value.lambdaKPlus1 = strtod(fields[2].c_str(), NULL);
        value.lambdaKPlus2 = strtod(fields[3].c_str(), NULL);
        value.gapK = strtod(fields[4].c_str(), NULL);
        value.gapKPlus1 = strtod(fields[5].c_str(), NULL);
        value.ratioK = strtod(fields[6].c_str(), NULL);
        value.ratioKPlus1 = strtod(fields[7].c_str(), NULL);
        value.minRatio = strtod(fields[8].c_str(), NULL);
        value.tolerance = strtod(fields[9].c_str(), NULL);
        value.estimatedRelativeError = strtod(fields[10].c_str(), NULL);
        value.coverageFilterEnabled = atoi(fields[11].c_str()) != 0;
        value.coverageZCutoff = strtod(fields[12].c_str(), NULL);
        value.coverageMeanLog = strtod(fields[13].c_str(), NULL);
        value.coverageSdLog = strtod(fields[14].c_str(), NULL);
        value.coverageMinimumEntries = strtod(fields[15].c_str(), NULL);
        value.coverageRemovedBins = atol(fields[16].c_str());
        value.coverageRetainedBins = atol(fields[17].c_str());
        values.push_back(value);
    }
    return values;
}

void replaceFile(const string &partial, const string &destination) {
    if (rename(partial.c_str(), destination.c_str()) != 0) {
        remove(partial.c_str());
        throw runtime_error("cannot install output " + destination + ": " +
                            string(strerror(errno)));
    }
}

void mergeWigs(const vector<ChromosomeInfo> &chromosomes,
               const string &baseKind, int eigenvector, int resolution,
               const string &destination, const string &description,
               bool allowMissing = false) {
    string partial = destination + ".partial." + toString(processId());
    ofstream out(partial.c_str());
    if (!out) throw runtime_error("cannot create output: " + partial);
    out << "track type=wiggle_0 name=\"" << description
        << "\" description=\"" << description << "\"\n";
    out << setprecision(17);
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        string base = chromosomes[i].base;
        if (baseKind == "low") {
            size_t slash = base.find_last_of('/');
            string dir = slash == string::npos ? "." : base.substr(0, slash);
            ostringstream low;
            low << dir << "/low_chr_" << setw(3) << setfill('0') << i;
            base = low.str();
        }
        string wigPath = base + ".Ev" + toString(eigenvector) + ".wig";
        vector<double> values;
        if (readable(wigPath)) {
            values = readWig(wigPath);
        } else if (allowMissing) {
            size_t bins = static_cast<size_t>(ceil(chromosomes[i].length /
                                                   static_cast<double>(resolution)));
            values.assign(bins, 0.0);
        } else {
            throw runtime_error("cannot read worker WIG: " + wigPath);
        }
        out << "fixedStep chrom=" << displayChromosome(chromosomes[i].name)
            << " start=1 step=" << resolution << " span=" << resolution << "\n";
        for (size_t j = 0; j < values.size(); ++j)
            out << (isfinite(values[j]) ? values[j] : 0.0) << '\n';
    }
    out.close();
    if (!out) throw runtime_error("failed while writing output: " + partial);
    replaceFile(partial, destination);
}

double correlationAtHighResolution(const vector<double> &high,
                                   const vector<double> &low,
                                   int highResolution) {
    // Expand the 50 kb reference onto high-resolution bins. A high-resolution
    // bin receives the value of the 50 kb bin containing its genomic center.
    // This also handles resolutions that do not divide 50 kb exactly.
    double meanX = 0.0, meanY = 0.0, sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
    size_t count = 0;
    for (size_t i = 0; i < high.size(); ++i) {
        long long center = static_cast<long long>(i) * highResolution + highResolution / 2;
        size_t lowIndex = static_cast<size_t>(center / kLowResolution);
        if (lowIndex >= low.size() || !isfinite(high[i]) || !isfinite(low[lowIndex]) ||
            high[i] == 0.0 || low[lowIndex] == 0.0) continue;
        ++count;
        double sampleCount = static_cast<double>(count);
        double dx = high[i] - meanX;
        meanX += dx / sampleCount;
        double dy = low[lowIndex] - meanY;
        meanY += dy / sampleCount;
        sumXX += dx * (high[i] - meanX);
        sumYY += dy * (low[lowIndex] - meanY);
        sumXY += dx * (low[lowIndex] - meanY);
    }
    if (count < 3 || sumXX <= 0.0 || sumYY <= 0.0)
        return numeric_limits<double>::quiet_NaN();
    return sumXY / sqrt(sumXX * sumYY);
}

struct VectorMetrics {
    double inverseParticipationRatio;
    double topOnePercentEnergy;
    double autocorrelation50kb;
    double autocorrelation1mb;
    bool localized;
};

double lagCorrelation(const vector<double> &values, size_t lag) {
    if (lag == 0 || lag >= values.size())
        return numeric_limits<double>::quiet_NaN();
    double meanX = 0.0, meanY = 0.0, sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
    size_t count = 0;
    for (size_t i = 0; i + lag < values.size(); ++i) {
        double x = values[i];
        double y = values[i + lag];
        if (!isfinite(x) || !isfinite(y) || x == 0.0 || y == 0.0) continue;
        ++count;
        double sampleCount = static_cast<double>(count);
        double dx = x - meanX;
        meanX += dx / sampleCount;
        double dy = y - meanY;
        meanY += dy / sampleCount;
        sumXX += dx * (x - meanX);
        sumYY += dy * (y - meanY);
        sumXY += dx * (y - meanY);
    }
    if (count < 3 || sumXX <= 0.0 || sumYY <= 0.0)
        return numeric_limits<double>::quiet_NaN();
    return sumXY / sqrt(sumXX * sumYY);
}

VectorMetrics calculateVectorMetrics(const vector<double> &values, int resolution,
                                     double maxTopOnePercentEnergy) {
    vector<double> squared;
    double sumSquares = 0.0;
    double sumFourth = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        if (!isfinite(values[i]) || values[i] == 0.0) continue;
        double square = values[i] * values[i];
        squared.push_back(square);
        sumSquares += square;
        sumFourth += square * square;
    }
    sort(squared.begin(), squared.end(), greater<double>());
    size_t topCount = squared.empty() ? 0 :
        max<size_t>(1, static_cast<size_t>(ceil(0.01 * squared.size())));
    double topEnergy = 0.0;
    for (size_t i = 0; i < topCount; ++i) topEnergy += squared[i];

    VectorMetrics metrics;
    metrics.inverseParticipationRatio = sumSquares > 0.0
        ? sumFourth / (sumSquares * sumSquares)
        : numeric_limits<double>::quiet_NaN();
    metrics.topOnePercentEnergy = sumSquares > 0.0
        ? topEnergy / sumSquares
        : numeric_limits<double>::quiet_NaN();
    size_t lag50kb = max<size_t>(1, static_cast<size_t>(llround(50000.0 / resolution)));
    size_t lag1mb = max<size_t>(1, static_cast<size_t>(llround(1000000.0 / resolution)));
    metrics.autocorrelation50kb = lagCorrelation(values, lag50kb);
    metrics.autocorrelation1mb = lagCorrelation(values, lag1mb);
    metrics.localized = !isfinite(metrics.topOnePercentEnergy) ||
                        metrics.topOnePercentEnergy > maxTopOnePercentEnergy;
    return metrics;
}

void writeRescue(const vector<ChromosomeInfo> &chromosomes,
                 const string &tempDirectory, const Options &options,
                 const string &outbase) {
    string wigDestination = outbase + ".wig";
    string wigPartial = wigDestination + ".partial." + toString(processId());
    string reportDestination = outbase + ".rescue.tsv";
    string reportPartial = reportDestination + ".partial." + toString(processId());
    string eigenDestination = outbase + ".eigenvalues.tsv";
    string eigenPartial = eigenDestination + ".partial." + toString(processId());

    ofstream wig(wigPartial.c_str());
    ofstream report(reportPartial.c_str());
    ofstream eigen(eigenPartial.c_str());
    if (!wig || !report || !eigen) throw runtime_error("cannot create rescue outputs");

    wig << "track type=wiggle_0 name=\"Lanczos rescued EV\" "
        << "description=\"50kb-guided rescued eigenvector\"\n" << setprecision(17);
    report << "chromosome\tselected_eigenvector\tcorrelation\tabs_correlation\tstatus"
           << "\tlambda_k\tlambda_k_plus_1\tlambda_k_plus_2\tgap_k\tgap_k_plus_1"
           << "\tlambda_k_over_gap_k\tlambda_k_plus_1_over_gap_k_plus_1\tmin_ratio"
           << "\ttolerance\testimated_relative_error_first_k"
           << "\tinverse_participation_ratio\ttop_1pct_energy\tautocorrelation_50kb"
           << "\tautocorrelation_1mb\tlocalized"
           << "\tcoverage_filter_enabled\tcoverage_z_cutoff\tcoverage_mean_log1p_nonzero"
           << "\tcoverage_sd_log1p_nonzero\tcoverage_min_nonzero_entries"
           << "\tcoverage_removed_bins\tcoverage_retained_bins";
    for (int ev = 1; ev <= kRescueEigenvectors; ++ev) report << "\tcorr_EV" << ev;
    for (int ev = 1; ev <= kRescueEigenvectors; ++ev) report << "\tlocalized_EV" << ev;
    report << '\n' << setprecision(17);
    eigen << "chromosome\tresolution\tfirst_k\tlambda_k\tlambda_k_plus_1\tlambda_k_plus_2"
          << "\tgap_k\tgap_k_plus_1\tlambda_k_over_gap_k"
          << "\tlambda_k_plus_1_over_gap_k_plus_1\tmin_ratio\ttolerance"
          << "\testimated_relative_error_first_k\tinverse_participation_ratio"
          << "\ttop_1pct_energy\tautocorrelation_50kb\tautocorrelation_1mb\tlocalized"
          << "\tcoverage_filter_enabled\tcoverage_z_cutoff\tcoverage_mean_log1p_nonzero"
          << "\tcoverage_sd_log1p_nonzero\tcoverage_min_nonzero_entries"
          << "\tcoverage_removed_bins\tcoverage_retained_bins\n" << setprecision(17);

    for (size_t c = 0; c < chromosomes.size(); ++c) {
        ostringstream lowBase;
        lowBase << tempDirectory << "/low_chr_" << setw(3) << setfill('0') << c;
        string lowPath = lowBase.str() + ".Ev1.wig";
        bool hasLowReference = readable(lowPath);
        vector<double> low;
        if (hasLowReference) low = readWig(lowPath);
        vector<double> correlations(kRescueEigenvectors,
                                    numeric_limits<double>::quiet_NaN());
        vector<vector<double> > high(kRescueEigenvectors);
        vector<VectorMetrics> metrics(kRescueEigenvectors);
        int selected = -1;
        int best = -1;
        int bestOverall = -1;
        double bestAbsolute = -1.0;
        double bestOverallAbsolute = -1.0;
        for (int ev = 0; ev < kRescueEigenvectors; ++ev) {
            high[ev] = readWig(chromosomes[c].base + ".Ev" + toString(ev + 1) + ".wig");
            correlations[ev] = correlationAtHighResolution(high[ev], low, options.resolution);
            metrics[ev] = calculateVectorMetrics(high[ev], options.resolution,
                                                 options.maxTopOnePercentEnergy);
            if (isfinite(correlations[ev])) {
                double absolute = fabs(correlations[ev]);
                if (absolute > bestOverallAbsolute) {
                    bestOverallAbsolute = absolute;
                    bestOverall = ev;
                }
                if (!metrics[ev].localized && absolute > bestAbsolute) {
                    bestAbsolute = absolute;
                    best = ev;
                }
                if (!metrics[ev].localized && selected < 0 &&
                    absolute > options.rescueThreshold) selected = ev;
            }
        }
        string status = "passed";
        if (!hasLowReference) {
            selected = 0;
            status = "no_low_resolution_reference";
        } else if (selected < 0) {
            if (best >= 0) {
                selected = best;
                status = "fallback_best";
            } else if (bestOverall >= 0) {
                selected = bestOverall;
                status = "fallback_localized_no_eligible_vector";
            } else {
                selected = 0;
                status = "fallback_EV1_no_finite_correlation";
            }
        }
        double selectedCorrelation = correlations[selected];
        double sign = isfinite(selectedCorrelation) && selectedCorrelation < 0.0 ? -1.0 : 1.0;
        vector<EigenDiagnostic> diagnostics =
            readEigenDiagnostics(chromosomes[c].base + ".eigenvalues.tsv");
        if (diagnostics.size() < kRescueEigenvectors)
            throw runtime_error("incomplete eigenvalue diagnostics for " + displayChromosome(chromosomes[c].name));

        string displayChrom = displayChromosome(chromosomes[c].name);
        wig << "fixedStep chrom=" << displayChrom << " start=1 step="
            << options.resolution << " span=" << options.resolution << '\n';
        for (size_t i = 0; i < high[selected].size(); ++i)
            wig << (isfinite(high[selected][i]) ? sign * high[selected][i] : 0.0) << '\n';

        const EigenDiagnostic &selectedDiagnostic = diagnostics[selected];
        report << displayChrom << '\t' << "EV" << (selected + 1) << '\t'
               << selectedCorrelation << '\t'
               << (isfinite(selectedCorrelation) ? fabs(selectedCorrelation)
                                                  : numeric_limits<double>::quiet_NaN())
               << '\t' << status
               << '\t' << selectedDiagnostic.lambdaK
               << '\t' << selectedDiagnostic.lambdaKPlus1
               << '\t' << selectedDiagnostic.lambdaKPlus2
               << '\t' << selectedDiagnostic.gapK
               << '\t' << selectedDiagnostic.gapKPlus1
               << '\t' << selectedDiagnostic.ratioK
               << '\t' << selectedDiagnostic.ratioKPlus1
               << '\t' << selectedDiagnostic.minRatio
               << '\t' << selectedDiagnostic.tolerance
               << '\t' << selectedDiagnostic.estimatedRelativeError
               << '\t' << metrics[selected].inverseParticipationRatio
               << '\t' << metrics[selected].topOnePercentEnergy
               << '\t' << metrics[selected].autocorrelation50kb
               << '\t' << metrics[selected].autocorrelation1mb
               << '\t' << (metrics[selected].localized ? 1 : 0)
               << '\t' << (selectedDiagnostic.coverageFilterEnabled ? 1 : 0)
               << '\t' << selectedDiagnostic.coverageZCutoff
               << '\t' << selectedDiagnostic.coverageMeanLog
               << '\t' << selectedDiagnostic.coverageSdLog
               << '\t' << selectedDiagnostic.coverageMinimumEntries
               << '\t' << selectedDiagnostic.coverageRemovedBins
               << '\t' << selectedDiagnostic.coverageRetainedBins;
        for (int ev = 0; ev < kRescueEigenvectors; ++ev)
            report << '\t' << correlations[ev];
        for (int ev = 0; ev < kRescueEigenvectors; ++ev)
            report << '\t' << (metrics[ev].localized ? 1 : 0);
        report << '\n';
        for (int ev = 0; ev < kRescueEigenvectors; ++ev) {
            const EigenDiagnostic &d = diagnostics[ev];
            eigen << displayChrom << '\t' << options.resolution << '\t' << d.firstK
                  << '\t' << d.lambdaK << '\t' << d.lambdaKPlus1 << '\t' << d.lambdaKPlus2
                  << '\t' << d.gapK << '\t' << d.gapKPlus1
                  << '\t' << d.ratioK << '\t' << d.ratioKPlus1 << '\t' << d.minRatio
                  << '\t' << d.tolerance << '\t' << d.estimatedRelativeError
                  << '\t' << metrics[ev].inverseParticipationRatio
                  << '\t' << metrics[ev].topOnePercentEnergy
                  << '\t' << metrics[ev].autocorrelation50kb
                  << '\t' << metrics[ev].autocorrelation1mb
                  << '\t' << (metrics[ev].localized ? 1 : 0)
                  << '\t' << (d.coverageFilterEnabled ? 1 : 0)
                  << '\t' << d.coverageZCutoff << '\t' << d.coverageMeanLog
                  << '\t' << d.coverageSdLog << '\t' << d.coverageMinimumEntries
                  << '\t' << d.coverageRemovedBins << '\t' << d.coverageRetainedBins << '\n';
        }
    }
    wig.close(); report.close(); eigen.close();
    if (!wig || !report || !eigen) throw runtime_error("failed while writing rescue outputs");
    replaceFile(wigPartial, wigDestination);
    replaceFile(reportPartial, reportDestination);
    replaceFile(eigenPartial, eigenDestination);
}

void cleanup(const vector<ChromosomeInfo> &chromosomes,
             const string &tempDirectory, bool rescue) {
    for (size_t c = 0; c < chromosomes.size(); ++c) {
        vector<string> phases;
        phases.push_back(rescue ? "high" : "normal");
        if (rescue) phases.push_back("low");
        for (size_t p = 0; p < phases.size(); ++p) {
            ostringstream base;
            base << tempDirectory << '/' << phases[p] << "_chr_"
                 << setw(3) << setfill('0') << c;
            int count = phases[p] == "high" ? kRescueEigenvectors : 1;
            for (int ev = 1; ev <= count; ++ev)
                remove((base.str() + ".Ev" + toString(ev) + ".wig").c_str());
            remove((base.str() + ".eigenvalues.tsv").c_str());
            remove((base.str() + ".log").c_str());
        }
    }
#ifdef _WIN32
    _rmdir(tempDirectory.c_str());
#else
    rmdir(tempDirectory.c_str());
#endif
}

} // namespace

int main(int argc, char **argv) {
    Options options;
    unsigned int cpus = thread::hardware_concurrency();
    if (cpus == 0) cpus = 1;
    options.jobs = static_cast<int>(min(4u, cpus));

    enum { OPT_RESCUE = 1000, OPT_THRESHOLD, OPT_KEEP_TEMP, OPT_COVERAGE_Z,
           OPT_NO_COVERAGE_FILTER, OPT_MAX_TOP1_ENERGY };
    static struct option longOptions[] = {
        {"rescue", no_argument, NULL, OPT_RESCUE},
        {"resolution", required_argument, NULL, 'r'},
        {"normalization", required_argument, NULL, 'n'},
        {"observed", no_argument, NULL, 'o'},
        {"jobs", required_argument, NULL, 'j'},
        {"threads", required_argument, NULL, 'T'},
        {"tolerance", required_argument, NULL, 't'},
        {"epsilon", required_argument, NULL, 'e'},
        {"max-iter", required_argument, NULL, 'I'},
        {"threshold", required_argument, NULL, OPT_THRESHOLD},
        {"coverage-z", required_argument, NULL, OPT_COVERAGE_Z},
        {"no-coverage-filter", no_argument, NULL, OPT_NO_COVERAGE_FILTER},
        {"max-top1-energy", required_argument, NULL, OPT_MAX_TOP1_ENERGY},
        {"keep-temp", no_argument, NULL, OPT_KEEP_TEMP},
        {"verbose", required_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int option = 0;
    while ((option = getopt_long(argc, argv, "r:n:oj:T:t:e:I:v:h", longOptions, NULL)) != -1) {
        switch (option) {
            case OPT_RESCUE: options.rescue = true; break;
            case OPT_THRESHOLD: options.rescueThreshold = atof(optarg); break;
            case OPT_KEEP_TEMP: options.keepTemp = true; break;
            case OPT_COVERAGE_Z:
                options.coverageZCutoff = atof(optarg);
                options.coverageFilter = true;
                break;
            case OPT_NO_COVERAGE_FILTER: options.coverageFilter = false; break;
            case OPT_MAX_TOP1_ENERGY: options.maxTopOnePercentEnergy = atof(optarg); break;
            case 'r': options.resolution = atoi(optarg); break;
            case 'n': options.normalization = optarg; break;
            case 'o': options.observed = true; break;
            case 'j': options.jobs = atoi(optarg); break;
            case 'T': options.solverThreads = atoi(optarg); break;
            case 't': options.tolerance = atof(optarg); break;
            case 'e': options.epsilon = atof(optarg); break;
            case 'I': options.maxIterations = atoi(optarg); break;
            case 'v': options.verbosity = atoi(optarg); break;
            case 'h': usage(argv[0]); return EXIT_SUCCESS;
            default: usage(argv[0]); return EXIT_FAILURE;
        }
    }

    if (argc - optind < 2 || argc - optind > 3) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    string hicfile = argv[optind++];
    string outbase = argv[optind++];
    if (optind < argc) options.resolution = atoi(argv[optind++]);
    if (options.resolution <= 0 || options.jobs <= 0 || options.solverThreads <= 0 ||
        options.maxIterations <= (options.rescue ? kRescueEigenvectors + 2 : 3) ||
        options.rescueThreshold < 0.0 || options.rescueThreshold > 1.0 ||
        options.coverageZCutoff >= 0.0 || options.maxTopOnePercentEnergy <= 0.0 ||
        options.maxTopOnePercentEnergy > 1.0) {
        cerr << "Invalid numeric option. Resolution/jobs/threads must be positive, "
             << "max-iter must exceed the requested eigenvector count plus two, "
             << "threshold must be in [0,1], coverage-z must be negative, and "
             << "max-top1-energy must be in (0,1]." << endl;
        return EXIT_FAILURE;
    }
    if (options.rescue && options.resolution > kLowResolution) {
        cerr << "--rescue requires a requested resolution of 50 kb or finer." << endl;
        return EXIT_FAILURE;
    }

    string worker = directoryName(argv[0]) + "/LanChr.exe";
    if (!regularExecutable(worker)) {
        cerr << "Cannot find chromosome worker executable: " << worker << endl;
        return EXIT_FAILURE;
    }

    string tempDirectory;
    try {
        vector<ChromosomeInfo> chromosomes = readChromosomes(hicfile);
        if (chromosomes.empty()) throw runtime_error("Hi-C file has no chromosomes");
        ensureOutputDirectory(outbase);
        tempDirectory = makeTempDirectory(outbase);
        if (options.verbosity > 0) {
            cerr << "Running " << chromosomes.size() << " chromosomes at "
                 << options.resolution << " bp with " << options.normalization
                 << (options.observed ? " observed" : " observed/expected")
                 << "; up to " << options.jobs << " chromosome jobs in parallel"
                 << (options.coverageFilter ? "; high-resolution coverage z cutoff " +
                     toString(options.coverageZCutoff) : "; coverage filter disabled")
                 << endl;
        }

        if (!options.rescue) {
            vector<Job> jobs = makeJobs(chromosomes, worker, options, hicfile,
                                        tempDirectory, "normal",
                                        options.resolution, 1);
            if (!runJobs(jobs, options.jobs, options.verbosity, true) ||
                !onlyExpectedSparseFailures(jobs))
                throw runtime_error("one or more chromosome calculations failed unexpectedly; temporary files retained at " + tempDirectory);
            mergeWigs(chromosomes, "normal", 1, options.resolution,
                      outbase + ".wig", "Lanczos genome-wide EV1", true);
        } else {
            vector<Job> lowJobs = makeJobs(chromosomes, worker, options, hicfile,
                                           tempDirectory, "low",
                                           kLowResolution, 1);
            if (!runJobs(lowJobs, options.jobs, options.verbosity, true) ||
                !onlyExpectedSparseFailures(lowJobs))
                throw runtime_error("one or more 50 kb calculations failed unexpectedly; temporary files retained at " + tempDirectory);
            vector<Job> highJobs = makeJobs(chromosomes, worker, options, hicfile,
                                            tempDirectory, "high",
                                            options.resolution,
                                            kRescueEigenvectors);
            if (!runJobs(highJobs, options.jobs, options.verbosity))
                throw runtime_error("one or more high-resolution calculations failed; temporary files retained at " + tempDirectory);

            mergeWigs(chromosomes, "low", 1, kLowResolution,
                      outbase + ".lowres_50kb.wig", "Lanczos 50kb reference EV1", true);
            for (int ev = 1; ev <= kRescueEigenvectors; ++ev)
                mergeWigs(chromosomes, "high", ev, options.resolution,
                          outbase + ".Ev" + toString(ev) + ".wig",
                          "Lanczos high-resolution EV" + toString(ev));
            writeRescue(chromosomes, tempDirectory, options, outbase);
        }

        if (!options.keepTemp) cleanup(chromosomes, tempDirectory, options.rescue);
        else if (options.verbosity > 0) cerr << "Temporary files retained at " << tempDirectory << endl;
        cout << "Wrote " << outbase << ".wig" << endl;
        return EXIT_SUCCESS;
    } catch (const exception &error) {
        cerr << "Lan.exe: " << error.what() << endl;
        return EXIT_FAILURE;
    }
}
