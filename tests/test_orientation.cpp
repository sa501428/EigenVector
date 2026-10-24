#include "../C++/common/CompartmentOrientation.h"
#include <filesystem>
#include <iostream>

void check(bool ok,const std::string &message) { if (!ok) throw std::runtime_error(message); }
void zipped(const std::string &path,const std::string &text) {
    gzFile out = gzopen(path.c_str(),"wb");
    check(out != NULL,"create fixture gzip");
    check(gzwrite(out,text.data(),static_cast<unsigned int>(text.size())) == static_cast<int>(text.size()),"write gzip");
    check(gzclose(out) == Z_OK,"close gzip");
}
void fixture(const std::string &root,bool conflict = false,bool gaps = false) {
    std::filesystem::create_directories(root + "/hg38");
    std::ofstream(root + "/hg38/chrom.sizes") << "chr1\t1003\n";
    std::string header = "{\"version\":1,\"bin_size\":100,\"columns\":[\"gc_bases\",\"acgt_bases\"],\"dtype\":\"uint8\",\"layout\":\"interleaved\",\"chromosomes\":[[\"chr1\",1003]]}";
    std::string data = "GCCOUNT1";
    for (int i = 0; i < 4; ++i) data += static_cast<char>((header.size() >> (8*i)) & 255);
    data += header;
    std::ostringstream genes;
    for (int i = 0; i < 10; ++i) {
        data += static_cast<char>(gaps && i == 0 ? 0 : gaps && i == 1 ? 10 : (i+1)*10);
        data += static_cast<char>(gaps && i == 0 ? 0 : gaps && i == 1 ? 20 : 100);
        int covered = conflict ? 100-(i+1)*10 : (i+1)*10;
        if (covered) genes << "chr1\t" << i*100 << '\t' << i*100+covered << "\t1\n";
        if (covered < 100) genes << "chr1\t" << i*100+covered << '\t' << (i+1)*100 << "\t0\n";
    }
    data += static_cast<char>(1); data += static_cast<char>(3);
    genes << (conflict ? "chr1\t1000\t1002\t1\nchr1\t1002\t1003\t0\n" : "chr1\t1000\t1001\t1\nchr1\t1001\t1003\t0\n");
    zipped(root + "/hg38/gc.100bp.counts.bin.gz",data);
    zipped(root + "/hg38/refseq.coverage.bedGraph.gz",genes.str());
}
int main(int argc,char **argv) try {
    check(argc == 2,"temporary directory argument required");
    std::string root = std::string(argv[1])+"/phasing-fixture";
    fixture(root);
    possum::PhasingOptions options;
    options.directory = root; options.minBins = 2;
    possum::BundledPhasing bundle(root + "/hg38");
    auto bins = bundle.bin("chr1",11,100,1003,true,true,false);
    check(std::fabs(bins.gc.back()-1.0/3) < 1e-12,"terminal GC count");
    check(std::fabs(bins.genes.back()-1.0/3) < 1e-12,"terminal gene interval");
    auto coarse = bundle.bin("chr1",3,500,1003,true,true,false);
    check(std::fabs(coarse.gc[0]-0.3) < 1e-12,"GC must sum counts");
    check(std::fabs(coarse.gc[1]-0.8) < 1e-12,"coarse GC");
    fixture(root+"/gaps",false,true);
    possum::BundledPhasing gaps(root+"/gaps/hg38");
    auto missing = gaps.bin("chr1",11,100,1003,true,false,false);
    check(std::isnan(missing.gc[0]) && std::isnan(missing.gc[1]),"N-only and poorly covered tiles are missing");
    auto weighted = gaps.bin("chr1",3,500,1003,true,false,false);
    check(std::fabs(weighted.gc[0]-130.0/320) < 1e-12,"GC weights callable base counts, not tile fractions");
    for (const std::string name : {"hg38","GRCh38",""}) {
        possum::CompartmentOrientation orientation(name,"",options);
        check(orientation.assembly() == "hg38","empty assembly default/alias");
        std::vector<double> x = bins.gc;
        for (auto &v : x) v = -v;
        check(orientation.orient(x.data(),x.size(),"1",100,0,1,1003) == -1,"negative auto correlation flips");
        check(orientation.orient(x.data(),x.size(),"chr1",100,0,1,1003) == 1,"orientation is idempotent");
        std::vector<double> constant(11,1);
        check(orientation.evidence(constant.data(),11,"chr1",100,1003) == 0,"constant vector unresolved");
    }
    auto offgrid = bundle.bin("chr1",11,99,1003,true,true,false);
    std::vector<double> x = offgrid.genes;
    for (auto &v : x) v = -v;
    possum::CompartmentOrientation automatic("hg38","",options);
    auto decision = automatic.evaluate({{x.data(),11,"chr1",1003}},99);
    check(decision.sign == -1 && decision.source == "bundled_refseq" && !decision.gcExact,
          "off-grid auto uses exact RefSeq, not silent GC approximation");
    bool rejected = false;
    try { automatic.evidence(x.data(),10,"chr1",100,1000); }
    catch (const std::runtime_error &) { rejected = true; }
    check(rejected,"assembly/length mismatch must fail");
    auto gcOptions = options; gcOptions.method = "gc";
    rejected = false;
    try { possum::CompartmentOrientation("hg38","",gcOptions).evidence(x.data(),11,"chr1",99,1003); }
    catch (const std::runtime_error &) { rejected = true; }
    check(rejected,"off-grid forced GC must require explicit approximation");
    gcOptions.approximateGC = true;
    auto approximate = bundle.bin("chr1",11,99,1003,true,false,true);
    x = approximate.gc;
    for (auto &v : x) v = -v;
    decision = possum::CompartmentOrientation("hg38","",gcOptions).evaluate({{x.data(),11,"chr1",1003}},99);
    check(decision.sign == -1 && !decision.gcExact,"explicit approximate GC is recorded");
    auto weakOptions = options; weakOptions.minBins = 20;
    check(possum::CompartmentOrientation("hg38","",weakOptions).evidence(x.data(),11,"chr1",100,1003) == 0,
          "insufficient bin count is unresolved");
    auto disabled = options; disabled.method = "none"; disabled.directory = "missing";
    check(possum::CompartmentOrientation("hg38","",disabled).evidence(x.data(),11,"chr1",100,1003) == 0,
          "disabled mode needs no resources");
    fixture(root+"/conflict",true);
    auto conflicting = options; conflicting.directory = root+"/conflict";
    x = bins.gc;
    decision = possum::CompartmentOrientation("hg38","",conflicting).evaluate({{x.data(),11,"chr1",1003}},100);
    check(decision.sign == 0 && decision.status == "conflicting_references","strong disagreement must be unresolved");
    check(possum::CompartmentOrientation("custom").assembly() == "custom","unknown map must not become hg38");
    auto override = options; override.assembly = "GRCh38";
    check(possum::CompartmentOrientation("custom","",override).assembly() == "hg38","explicit override");

    const std::string track = root+"/activity.bedGraph";
    { std::ofstream out(track); out << "chrMT 1000 1500 3\nchrMT 0 500 1\nchrMT 500 1000 2\nchrY 0 1500 0\n"; }
    possum::CompartmentOrientation custom("custom",track,options);
    custom.setReport(root+"/custom.orientation.tsv");
    double mt[3] = {3,2,1};
    check(custom.orient(mt,3,"M",500) == -1,"custom aliases and Spearman flip");
    check(custom.orient(mt,3,"chrMT",500) == 1,"custom idempotence");
    double coarseCustom[2] = {2,1};
    check(custom.orient(coarseCustom,2,"MT",750) == -1,"custom arbitrary interval overlaps");
    double constant[3] = {1,2,3};
    check(custom.orient(constant,3,"Y",500) == 0,"zero reference is constant, not missing/inactive evidence");
    // A trans vector may be kept or flipped only as a whole, even when its
    // individual chromosome slices have opposite reference correlations.
    double global[6] = {-3,-2,-1,1,2,3};
    std::vector<double> before(global,global+6);
    std::vector<possum::PhaseSlice> slices = {{global,3,"MT",1500},{global+3,3,"MT",1500}};
    int sign = custom.orientGlobal(slices,global,6,500,0,2);
    for (int i=0;i<6;++i) check(global[i] == before[i]*(sign < 0 ? -1 : 1),"trans must preserve relative chromosome signs");
    { std::ifstream report(root+"/custom.orientation.tsv"); std::string contents((std::istreambuf_iterator<char>(report)),{});
      check(contents.find("gc_spearman") != std::string::npos && contents.find("genome") != std::string::npos,
            "machine-readable diagnostics include the global decision"); }
    for (const std::string text : {"chr1 -1 5 2\n","chr1 5 5 2\n","chr1 0 5 bad\n",
         "chr1 0 10 1\n1 5 15 2\n","chr1 0 10 NaN\nchr1 5 15 2\n","chr1 0 5 NaN\n"}) {
        { std::ofstream out(track); out << text; }
        rejected = false;
        try { possum::CompartmentOrientation invalid("custom",track,options); }
        catch (const std::runtime_error &) { rejected = true; }
        check(rejected,"malformed/empty/overlapping tracks must fail");
    }
    // Exercise the actual packaged format for all assemblies and aliases at a
    // resolution that was impossible for the historical 100 kb anchor caller.
#ifdef POSSUM_PHASING_DIR
    for (const auto &item : std::vector<std::pair<std::string,std::string>>{
            {"hg19","GRCh37"},{"hg38","GRCh38"},{"mm9","NCBIM37"},{"mm10","GRCm38"}}) {
        possum::PhasingOptions real; real.method = "gc";
        possum::BundledPhasing actual(std::string(POSSUM_PHASING_DIR)+"/"+item.first);
        int64_t length = actual.length("chr1");
        int n = possum::binCount(length,7000);
        auto reference = actual.bin("chr1",n,7000,length,true,false,false);
        for (auto &v : reference.gc) v = -v;
        possum::CompartmentOrientation orientation(item.second,"",real);
        check(orientation.orient(reference.gc.data(),n,"1",7000,0,1,length) == -1,"real bundled assembly must flip");
        possum::PhasingChromosomes metadata = {{"1",length},{"chrUn_scaffold",123}};
        check(possum::CompartmentOrientation("sample_"+item.first+"_suffix","",real,metadata).assembly() == item.first,
              "genome ID substring with size confirmation");
        check(possum::CompartmentOrientation("custom_reference","",real,metadata).assembly() == item.first,
              "unique chromosome names/sizes fallback");
        check(possum::CompartmentOrientation("","",real,metadata).assembly() == item.first,
              "missing genome ID uses chromosome sizes for every assembly");
    }
    possum::PhasingOptions real;
    possum::BundledPhasing hg19(std::string(POSSUM_PHASING_DIR)+"/hg19");
    possum::PhasingChromosomes b37 = {{"chr1",hg19.length("chr1")},{"MT",16569}};
    possum::CompartmentOrientation b37Orientation("b37_custom","",real,b37);
    check(b37Orientation.assembly() == "hg19","mitochondrial variants do not change nuclear assembly identification");
    std::vector<double> mtValues(17,1);
    auto mtDecision = b37Orientation.evaluate({{mtValues.data(),17,"MT",16569}},1000);
    check(mtDecision.sign == 0 && mtDecision.status == "chromosome_length_mismatch","incompatible mitochondrial reference stays unresolved");
    possum::PhasingChromosomes onlyMito = {{"MT",16299}};
    check(!possum::bundledAssembly(possum::CompartmentOrientation("unknown","",real,onlyMito).assembly()),
          "mitochondria alone cannot identify a nuclear assembly");
#endif
    const auto ranks = possum::phaseRanks({1,1,9,3});
    check(ranks == std::vector<double>({0.5,0.5,3,2}),"Spearman ties get average ranks");
    possum::PhasingOptions precision; precision.minCorrelation = 1e-9;
    check(std::stod(precision.arguments()[3]) == 1e-9,"worker threshold forwarding preserves precision");
    std::filesystem::remove_all(root);
    std::cout << "Orientation passed: four real bundles/aliases, exact and off-grid bins, weak/conflicting evidence, custom tracks, global sign and diagnostics\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
}
