#include <cstring>
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
  fprintf(stderr, "Usage: %s [-o observed][-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <outfile>: eigenvectors output file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n  -p, --phasing-track FILE: activity bedGraph (higher = A)\n");
  possum::phasingUsage();
}

int main(int argc, char *argv[]) try {
        string norm("NONE");
        string unit("BP");
        ifstream fin;
        struct timespec t0,t1,st,en;

        double tol=1.0e-6;
        int maxiter=500;
        int threads = 1;
        int verb = 1;
	string ob("oe");
        int opt;
	string phasingTrack;
    possum::PhasingOptions phasingOptions;

	clock_gettime(CLOCK_REALTIME,&st);

    static option longOptions[] = {
        POSSUM_PHASING_LONG_OPTIONS
        {"phasing-track", required_argument, NULL, 'p'},
        {NULL, 0, NULL, 0}
    };
	while ((opt = getopt_long(argc, argv, "p:o:t:I:T:n:v:h", longOptions, NULL)) != -1) {
        if (phasingOptions.parse(opt,optarg)) continue;
        switch (opt) {
                case 'p': phasingTrack = optarg; break;
                                case 'o':
                                        if (strcmp(optarg,"observed") == 0) ob = "observed";
					else {
						usage(argv[0]);
                                        	exit(EXIT_FAILURE);
					}
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

	if (argc - optind < 3) {
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

	char *out_name = argv[optind++];
	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else if (verb) printf("Writing eigenvectors to %s\n\n", argv[optind-1]);
	int binsize = possum::positiveInt(argv[optind++]);
    if (binsize <= 0 || threads <= 0 || maxiter <= 0 || !std::isfinite(tol) || tol <= 0)
        throw std::runtime_error("resolution, threads, iterations, and tolerance must be positive");
	fprintf(fout,"track type=wiggle_0\n");

	int *ii;
	int *jj;
	double *xx;
	long p,q;
	long m;
	int k;
	int i;

	string hg19("hg19");
	string hg38("hg38");

// chromosome map for finding matrix
    int64_t master = 0;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0;
    int64_t nviLength = 0;

    chromosomeMap = readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
    fin.close();
    if (chromosomeMap.empty()) throw std::runtime_error("Hi-C file has no chromosomes");

    possum::CompartmentOrientation orientation(genomeID, phasingTrack, phasingOptions, possum::phasingChromosomes(chromosomeMap));
    orientation.setReport(std::string(out_name) + ".orientation.tsv");

	vector<contactRecord> records;
	int iter0 = 0;

//	loop accross all chromosomes
	for (map<string,chromosome>::iterator itr =  chromosomeMap.begin(); itr !=  chromosomeMap.end(); ++itr) {
		string chrom = itr->second.name;
		char *chr = const_cast<char*> (chrom.c_str());
		if (strcmp(chr,"ALL")==0 || strcmp(chr,"All")==0 || strcmp(chr,"all")==0) continue;
		clock_gettime(CLOCK_REALTIME,&t0);
        	records = straw(ob, norm, fname, chrom, chrom, unit, binsize);

        	m = records.size();
		ii = (int *) malloc(m*sizeof(int));
		jj = (int *) malloc(m*sizeof(int));
		xx = (double *) malloc(m*sizeof(double));
		p=0;
                for (q=0; q<m; q++) {
                        if (std::isfinite((double) records[q].counts) && records[q].binX >= 0 && records[q].binY >= 0) {
                		ii[p] = records[q].binX/binsize;
                        	jj[p] = records[q].binY/binsize;
                        	xx[p] = (double) records[q].counts;
				p++;
			}
                }
        	records.clear();
        	records.shrink_to_fit();
		clock_gettime(CLOCK_REALTIME,&t1);
		if (verb > 1) printf("chromosome %s:  took %10.3f seconds for %ld records\n",chr,((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9,m);
		m = p;
		for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

		int N = possum::binCount(itr->second.length, binsize);

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
		if (verb > 1) {
			printf("total %d iterations\n",iter);
                	printf("iterations took %10.3f seconds\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9);
			printf("lam1 = %g; lam2 = %g; lam1/lam2 = %g; er = %g; error in EV = %g\n",l1,l2,l1/l2,er,er/(l1-l2));
			printf("                           -------------------- \n");
		}

        std::string displayChr = possum::displayChromosome(chrom);
        orientation.orient(ev, N, chrom, binsize, verb, 1, itr->second.length);
        fprintf(fout,"fixedStep chrom=%s start=1 step=%d span=%d\n", displayChr.c_str(), binsize, binsize);
        for (int j=0; j<N; ++j) {
            // A partial terminal bin needs its actual span; exact multiples
            // retain the original section (never a zero step/span).
            long long remainder = itr->second.length % binsize;
            if (j == N-1 && remainder > 0)
                fprintf(fout,"fixedStep chrom=%s start=%lld step=%d span=%lld\n",
                        displayChr.c_str(), static_cast<long long>(j)*binsize+1, binsize, remainder);
            fprintf(fout,"%.17g\n", std::isfinite(ev[j]) ? ev[j] : 0.0);
        }
        free(a1);
        free(a2);
        free(ev);

		free(ii);
		free(jj);
		free(xx);
	}
	fclose(fout);
	clock_gettime(CLOCK_REALTIME,&en);
	if (verb) printf("\n**************    all together took %10.3f seconds\n",((double) (en.tv_sec - st.tv_sec)) + ((double) (en.tv_nsec - st.tv_nsec))/1e9);
	return(0);
} catch (const std::exception &error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return EXIT_FAILURE;
}
