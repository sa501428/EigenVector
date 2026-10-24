#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
#include "../common/CompartmentOrientation.h"
#include "../common/Arguments.h"
using namespace std;

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads);

unsigned int getMatrix(string fname, int binsize, string norm, string ob, bool interOnly, unsigned int **i, unsigned int **j, float **x, long *m, vector<std::string> &CH,vector<int64_t> &CHL, string &genomeID);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [options] <hicfile> <outbase> <resolution> [nv]\n\n", argv0);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -f           Use full matrix instead of inter-chromosomal only\n");
  fprintf(stderr, "  -t <float>   Set tolerance (default: 1.0e-7)\n");
  fprintf(stderr, "  -e <float>   Set epsilon (default: 1.0e-8)\n");
  fprintf(stderr, "  -I <int>     Set maximum iterations (default: 200)\n");
  fprintf(stderr, "  -T <int>     Set number of threads (default: 1)\n");
  fprintf(stderr, "  -v <int>     Set verbosity level (default: 1)\n");
  fprintf(stderr, "  -p, --phasing-track FILE\n               Activity bedGraph (higher = A)\n  -h           Show this help message\n\n");
  fprintf(stderr, "Required arguments:\n");
  fprintf(stderr, "  <hicfile>    Path to .hic file\n");
  fprintf(stderr, "  <outbase>    Base name for output eigenvector files\n");
  fprintf(stderr, "  <resolution> Resolution in base pairs\n");
  fprintf(stderr, "  [nv]        Number of eigenvectors (optional, default: 2)\n\n");
  fprintf(stderr, "Output:\n");
  fprintf(stderr, "  Creates .wig files for each eigenvector with the naming pattern: <outbase>_Ev<n>.wig\n");
  possum::phasingUsage();
}

int main(int argc, char *argv[]) try {
	string norm("NONE");
	string unit("BP");
	string ob("observed");
	ifstream fin;
	time_t t0,t1;

	double tol=1.0e-7;
	double eps=1.0e-8;
	int maxiter=200;
	int threads = 1;
	int verb = 1;
	bool interOnly = true;
	unsigned int N;
	int opt;
	int nv=2;
    string phasingTrack;
    possum::PhasingOptions phasingOptions;

    static option longOptions[] = {
        POSSUM_PHASING_LONG_OPTIONS
        {"phasing-track",required_argument,NULL,'p'}, {NULL,0,NULL,0}};
	while ((opt = getopt_long(argc, argv, "p:ft:e:I:T:v:h",longOptions,NULL)) != -1) {
		if (phasingOptions.parse(opt,optarg)) continue;
        switch (opt) {
                case 'p': phasingTrack=optarg; break;
				case 'f':
					interOnly = false;
					break;
				case 't':
					tol = atof(optarg);
					break;
				case 'e':
					eps = atof(optarg);
					break;
				case 'I':
					maxiter=possum::positiveInt(optarg);
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
	if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	char *out_name = argv[optind++];
	if (verb) printf("Eigenvectors base name is: %s\n", out_name);
	int binsize = possum::positiveInt(argv[optind++]);
	if (argc > optind) nv=possum::positiveInt(argv[optind++]);

    if (binsize <= 0 || threads <= 0 || nv <= 0 || nv > maxiter-3 ||
        !isfinite(tol) || tol <= 0 || !isfinite(eps) || eps <= 0)
        throw runtime_error("invalid resolution or solver options");
	printf("\n");

	unsigned int *i;
        unsigned int *j;
        float *x;
        long p;
        long m;
        vector<std::string> chroms;
        vector<int64_t> chrLen;
        string genomeID;
        time(&t0);
        N = getMatrix(fname, binsize, norm, ob, interOnly, &i, &j, &x, &m,chroms,chrLen,genomeID);
        time(&t1);
        if (N <= 0) {
                cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
                exit(EXIT_FAILURE);
        }
        possum::CompartmentOrientation orientation(genomeID,phasingTrack,phasingOptions,possum::phasingChromosomes(chroms,chrLen));
        orientation.setReport(std::string(out_name) + ".orientation.tsv");
        if (verb) printf("took %ld seconds for %ld records\n",(long) (t1 - t0),m);

        for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;
        int iter;
	
        double er;
        double *r = (double *) malloc(N*sizeof(double));
        for (p=0;p<N;p++) r[p] = 1.0;
        double *lam = (double *) malloc((nv+2)*sizeof(double));
        double **ev = (double **) malloc((nv+2)*sizeof(double *));
        for (p=0;p<nv+2;p++) ev[p] = (double *) malloc(N*sizeof(double));

	time(&t0);
	iter = SOLan(m,i,j,x,&N,r,nv,lam,ev,&er,tol,eps,maxiter,threads);
	time(&t1);
	if (iter > 90000) {
		printf("return code is %d\n",iter);
		exit(EXIT_FAILURE);
	}

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %ld seconds\n",t1 - t0);
	}

	for (int j0=0;j0<nv+2;j0++) printf("%lg ",lam[j0]);
        printf("\n");
        // Trans EVs form one genome-wide eigensystem: pool reference/EV pairs
        // across chromosomes, then apply one global sign.
        for (int j0=0;j0<nv;j0++) {
                std::vector<possum::PhaseSlice> slices;
                size_t offset=0;
                for (size_t c=0; c<chroms.size(); ++c) {
                    int bins=possum::binCount(chrLen[c],binsize);
                    slices.push_back({ev[j0]+offset,bins,chroms[c],chrLen[c]});
                    offset += bins;
                }
                orientation.orientGlobal(slices,ev[j0],N,binsize,verb,j0+1);
                char *curout = (char *) malloc((32+strlen(out_name))*sizeof(char));
                char *temp = (char *) malloc(50);
                snprintf(temp, 50, "_Ev%d", j0+1);
                strcpy(curout,out_name);
                strcat(curout,temp);
                strcat(curout,".wig");
                FILE *fout = fopen(curout,"w");
                if (fout==NULL) {
                        fprintf(stderr, "Error! File %s cannot be opened for writing\n", curout);
                        exit(EXIT_FAILURE);
                }

                fprintf(fout, "track type=wiggle_0 name=\"Eigenvector %d\" description=\"Eigenvector %d from genome-wide analysis\"\n", j0+1, j0+1);
                int a = 0;
                for (int b=0; b<chroms.size(); b++) {
                        string displayChr = possum::displayChromosome(chroms.at(b));
                        const char *chr1 = displayChr.c_str();

                        // Calculate number of bins for this chromosome
                        int numBins = possum::binCount(chrLen.at(b),binsize);

                        // Write all bins with uniform size
                        fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", chr1, binsize, binsize);
                        for (int c=0; c<numBins; c++) {
                                if (!std::isnan(ev[j0][a])) {
                                        fprintf(fout,"%.17g\n", ev[j0][a]);
                                } else {
                                        fprintf(fout,"0\n");
                                }
                                a++;
                        }

                }
                fclose(fout);
                free(curout);
                free(temp);
        }
        return EXIT_SUCCESS;
} catch (const std::exception &error) {
    cerr << "Error: " << error.what() << endl;
    return EXIT_FAILURE;
}
