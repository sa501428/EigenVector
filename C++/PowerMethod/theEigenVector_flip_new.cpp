#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
#include "../common/CompartmentOrientation.h"
#include "../common/Arguments.h"
using namespace std;

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);



static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o (observed)][-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <chromosome> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chromosome>: chromosomee\n");
  fprintf(stderr, "  <outfile>: eigenvector output  file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n  -p, --phasing-track FILE: activity bedGraph (higher = A)\n");
  possum::phasingUsage();
}

int main(int argc, char *argv[]) try {
	string norm("NONE");
	string unit("BP");
	string ob("oe");
	ifstream fin;
	struct timespec t0,t1;

	double tol=1.0e-6;
	int maxiter=100;
	int threads = 1;
	int verb = 1;
	int N;
	int opt;
	string phasingTrack;
    possum::PhasingOptions phasingOptions;

    static option longOptions[] = {
        POSSUM_PHASING_LONG_OPTIONS
        {"phasing-track", required_argument, NULL, 'p'},
        {NULL, 0, NULL, 0}
    };
	while ((opt = getopt_long(argc, argv, "p:ot:I:T:n:v:h", longOptions, NULL)) != -1) {
		if (phasingOptions.parse(opt,optarg)) continue;
        switch (opt) {
                case 'p': phasingTrack = optarg; break;
				case 'o':
					ob = "observed";
					break;
				case 't':
					tol = atof(optarg);
					break;
				case 'I':
					maxiter=possum::positiveInt(optarg);
					break;
				case 'n':
					norm=optarg;
					break;
				case 'T':
					threads=possum::positiveInt(optarg);
					break;
				case 'v':
					verb=atoi(optarg);
					break;
				case 'h':
					usage(argv[0]);
					exit(EXIT_SUCCESS);
				default:
					usage(argv[0]);
					exit(EXIT_FAILURE);
		}
	}

	if (argc - optind < 4) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	string fname = argv[optind++];
	fin.open(fname, ios::in | ios::binary);
	if (fin.is_open()) {
		if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	}
	else {
		fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	string chrom(argv[optind++]);
	char *out_name = argv[optind++];
	int binsize = possum::positiveInt(argv[optind++]);
    if (binsize <= 0 || threads <= 0 || maxiter <= 0 || !std::isfinite(tol) || tol <= 0)
        throw std::runtime_error("resolution, threads, iterations, and tolerance must be positive");
	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	}
	else {
		if (verb) printf("Writing eigenvector to %s\n", out_name);
	}


	printf("\n");

// chromosome map for finding matrix
    int64_t master = 0;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0;
    int64_t nviLength = 0;

    chromosomeMap = readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
    map<string,chromosome>::iterator itr0 = chromosomeMap.find(chrom);
    if (itr0 == chromosomeMap.end()) {
        for (itr0 = chromosomeMap.begin(); itr0 != chromosomeMap.end(); ++itr0)
            if (possum::chromosomeKey(itr0->first) == possum::chromosomeKey(chrom)) break;
    }
    if (itr0 != chromosomeMap.end()) chrom = itr0->first;
    if (itr0 != chromosomeMap.end()) N = possum::binCount(itr0->second.length, binsize);
    	else {
		cout << "chromosome " << chrom << " is not found" << endl;
		exit(1);
	}
	fin.close();

    possum::CompartmentOrientation orientation(genomeID, phasingTrack, phasingOptions, possum::phasingChromosomes(chromosomeMap));
    orientation.setReport(std::string(out_name) + ".orientation.tsv");

	clock_gettime(CLOCK_REALTIME,&t0);
	vector<contactRecord> records = straw(ob, norm, fname, chrom, chrom, unit, binsize);
	long nonZer = records.size();
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	double *xx = (double *) malloc(nonZer*sizeof(double));
	for (long k=0; k<nonZer; k++) {
        if (records[k].binX < 0 || records[k].binY < 0 ||
            records[k].binX >= itr0->second.length || records[k].binY >= itr0->second.length)
            throw std::runtime_error("contact coordinate outside chromosome");
				ii[k] = records[k].binX/binsize; 
				jj[k] = records[k].binY/binsize; 
				xx[k] = (double) records[k].counts;
				if (std::isnan(xx[k])) xx[k] = 0;
	}
	clock_gettime(CLOCK_REALTIME,&t1);
	if (verb) printf("took %10.3f seconds for %ld records\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9,nonZer);
	records.clear();
	records.shrink_to_fit();

	long p;
	long m = nonZer;
	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;
	double l1,l2,er;
	double *a1 = (double *) malloc(N*sizeof(double));
	double *a2 = (double *) malloc(N*sizeof(double));
	double *ev = (double *) malloc(N*sizeof(double));

	long init_seed = 1234563;
	srand48(init_seed);
	for (p=0;p<N;p++) {
		a1[p] = drand48();
		a2[p] = drand48();
	}

	clock_gettime(CLOCK_REALTIME,&t0);
	int iter;
	iter = bestEigen(m,ii,jj,xx,&N,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
	clock_gettime(CLOCK_REALTIME,&t1);

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %10.3f seconds\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9);
		printf("lam1 = %g; lam2 = %g; lam1/lam2 = %g; er = %g; error in EV = %g\n",l1,l2,l1/l2,er,er/(l1-l2));
	}

    orientation.orient(ev, N, chrom, binsize, verb, 1, itr0->second.length);

	for (int j=0;j<N;j++) {
		if (!std::isnan(ev[j])) fprintf(fout,"%.17g\n",ev[j]);
		else fprintf(fout,"%20s\n","NaN");
	}
	fclose(fout);
	// Convergence iterations are diagnostic output, not a process status.
	// Returning the iteration count made every successful scripted run look
	// like a failure to shells using `set -e`.
	return 0;
} catch (const std::exception &error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return EXIT_FAILURE;
}
