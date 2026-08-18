#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads);

int flipSign(const char *genome, double *x, int n, char *chr, int binsize);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [options] <hicfile> <chromosome> <outbase> <resolution> [nv]\n\n", argv0);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -o           Use observed matrix instead of observed/expected (o/e) matrix\n");
  fprintf(stderr, "  -t <float>   Set tolerance (default: 1.0e-7)\n");
  fprintf(stderr, "  -e <float>   Set epsilon (default: 1.0e-8)\n");
  fprintf(stderr, "  -I <int>     Set maximum iterations (default: 200)\n");
  fprintf(stderr, "  -n <string>  Set normalization method (default: SCALE)\n");
  fprintf(stderr, "  -T <int>     Set number of threads (default: 1)\n");
  fprintf(stderr, "  -v <int>     Set verbosity level (default: 1)\n");
  fprintf(stderr, "  -h           Show this help message\n\n");
  fprintf(stderr, "Required arguments:\n");
  fprintf(stderr, "  <hicfile>     Path to .hic file\n");
  fprintf(stderr, "  <chromosome>  Chromosome name (e.g., chr1)\n");
  fprintf(stderr, "  <outbase>     Base name for output eigenvector files\n");
  fprintf(stderr, "  <resolution>  Resolution in base pairs\n");
  fprintf(stderr, "  [nv]         Number of eigenvectors to write (optional, default: 1)\n");
}

int main(int argc, char *argv[]) {
	string norm("SCALE");
	string unit("BP");
	string ob("oe");
	ifstream fin;
	time_t t0,t1;

	double tol=1.0e-7;
	double eps=1.0e-8;
	int maxiter=200;
	int threads = 1;
	int verb = 1;
	unsigned int N;
	int opt;

	while ((opt = getopt(argc, argv, "ot:e:I:T:n:v:h")) != -1) {
		switch (opt) {
				case 'o':
					ob = "observed";
					break;
				case 't':
					tol = atof(optarg);
					break;
				case 'e':
					eps = atof(optarg);
					break;
				case 'I':
					maxiter=atoi(optarg);
					break;
				case 'n':
					norm=optarg;
					break;
				case 'T':
					threads=atoi(optarg);
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
	fin.open(fname, fstream::in);
	if (fin.is_open()) {
		if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	}
	else {
		fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	int nv=1;
	string chrom(argv[optind++]);
	char *out_name = argv[optind++];
	if (verb) printf("Eigenvectors base name is: %s\n", out_name);
	int binsize = atoi(argv[optind++]);
	if (argc > optind) nv=atoi(argv[optind++]);
	if (binsize <= 0 || nv <= 0 || maxiter <= nv+2 || threads <= 0) {
		fprintf(stderr, "Error: resolution, eigenvector count, and threads must be positive; max iterations must exceed nv+2\n");
		exit(EXIT_FAILURE);
	}

	printf("\n");

// chromosome map for finding matrix
    int64_t master = 0LL;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0LL;
    int64_t nviLength = 0LL;
    int64_t totalFileSize;

	chromosomeMap = readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
	map<string,chromosome>::iterator itr0 = chromosomeMap.find(chrom);
	if (itr0 != chromosomeMap.end()) N = (int) ceil(itr0->second.length/((double) binsize));
	else {
		cout << "chromosome " << chrom << " is not found" << endl;
		exit(1);
	}
	fin.close();

        string hg19("hg19");
        string hg38("hg38");
        if (genomeID != hg19 && genomeID != hg38) if (verb) {
                cout << genomeID;
                cout << " is not currently supported; no sign flip will be attempted!" << endl << endl;
        }

	time(&t0);
	vector<contactRecord> records = straw(ob, norm, fname, chrom, chrom, unit, binsize);
	long nonZer = records.size();
	unsigned int *i = (unsigned int *) malloc(nonZer*sizeof(int));
	unsigned int *j = (unsigned int *) malloc(nonZer*sizeof(int));
	float *x = (float *) malloc(nonZer*sizeof(float));
	for (long k=0; k<nonZer; k++) {
				i[k] = records[k].binX/binsize; 
				j[k] = records[k].binY/binsize; 
				x[k] = (float) records[k].counts;
				if (isnan(x[k])) x[k] = 0;
	}
	time(&t1);
	if (verb) printf("took %ld seconds for %ld records\n",t1 - t0,nonZer);
	records.clear();
	records.shrink_to_fit();

	long p;
	long m = nonZer;
	for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;
	double er;
	double *r = (double *) malloc(N*sizeof(double));
	for (p=0;p<N;p++) r[p] = 1.0;
	double *lam = (double *) malloc((nv+2)*sizeof(double));
	double **ev = (double **) malloc((nv+2)*sizeof(double *));
	for (p=0;p<nv+2;p++) ev[p] = (double *) malloc(N*sizeof(double));

	time(&t0);
	int iter;
	iter = SOLan( m,i,j,x,&N,r,nv,lam,ev,&er,tol,eps,maxiter,threads);
	time(&t1);
	if (iter > 90000) {
		printf("return code is %d\n",iter);
		exit(EXIT_FAILURE);
	}

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %ld seconds\n",t1-t0);
	}

        char *chr = const_cast<char*> (chrom.c_str());
        char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
        if (!strstr(chr,"chr")) strcpy(chr1,"chr");
        else strcpy(chr1,"");
        strcat(chr1,chr);
        if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
       	char *genome1 = const_cast<char*> (genomeID.c_str());
        if (strcmp(chr1,"chrY")!=0 && strcmp(chr1,"chrM")!=0 && (100000 % binsize == 0)) {
		int junk;
		for (p=0;p<nv+2;p++) junk = flipSign(genome1,ev[p],N,chr1,binsize);
	}

	for (int j0=0;j0<nv+2;j0++) printf("%lg ",lam[j0]);
	printf("\n");

	char *eigenvalue_name = (char *) malloc((20+strlen(out_name))*sizeof(char));
	strcpy(eigenvalue_name, out_name);
	strcat(eigenvalue_name, ".eigenvalues.tsv");
	FILE *eigenvalue_out = fopen(eigenvalue_name, "w");
	if (eigenvalue_out==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", eigenvalue_name);
		exit(EXIT_FAILURE);
	}
	fprintf(eigenvalue_out,
		"first_k\tlambda_k\tlambda_k_plus_1\tlambda_k_plus_2\tgap_k\tgap_k_plus_1"
		"\tlambda_k_over_gap_k\tlambda_k_plus_1_over_gap_k_plus_1\tmin_ratio"
		"\ttolerance\testimated_relative_error_first_k\n");
	for (int j0=0;j0<nv;j0++) {
		double gap_k = lam[j0] - lam[j0+1];
		double gap_k_plus_1 = lam[j0+1] - lam[j0+2];
		double ratio_k = lam[j0] / gap_k;
		double ratio_k_plus_1 = lam[j0+1] / gap_k_plus_1;
		double min_ratio = ratio_k < ratio_k_plus_1 ? ratio_k : ratio_k_plus_1;
		double relative_error = tol * min_ratio;
		fprintf(eigenvalue_out,
			"%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
			j0+1, lam[j0], lam[j0+1], lam[j0+2], gap_k, gap_k_plus_1,
			ratio_k, ratio_k_plus_1, min_ratio, tol, relative_error);
	}
	fclose(eigenvalue_out);
	free(eigenvalue_name);

	// SOLan uses nv+2 Ritz vectors internally, but only the requested nv vectors
	// are part of the public output.
	for (int j0=0;j0<nv;j0++) {
		char *curout = (char *) malloc((10+strlen(out_name))*sizeof(char));
		char *temp = (char *) malloc(50);
		snprintf(temp, 50, ".Ev%d.wig", j0+1);
		strcpy(curout,out_name);
		strcat(curout,temp);
		FILE *fout = fopen(curout,"w");
		if (fout==NULL) {
			fprintf(stderr, "Error! File %s cannot be opened for writing\n", curout);
			exit(EXIT_FAILURE);
		}

                // Write track header
                fprintf(fout, "track type=wiggle_0 name=\"Eigenvector %d\" description=\"Eigenvector %d for %s\"\n", 
                        j0+1, j0+1, chrom.c_str());

                // Convert chromosome name to standard format
                char *chr = const_cast<char*> (chrom.c_str());
                char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
                if (!strstr(chr,"chr")) strcpy(chr1,"chr");
                else strcpy(chr1,"");
                strcat(chr1,chr);
                if (strcmp(chr1,"chrMT") == 0) strcpy(chr1,"chrM");

                // Write all bins with uniform size
                fprintf(fout, "fixedStep chrom=%s start=1 step=%d span=%d\n", chr1, binsize, binsize);
                for (p=0; p<N; p++) {
                        if (!isnan(ev[j0][p])) {
                                fprintf(fout,"%lg\n", ev[j0][p]);
                        } else {
                                fprintf(fout,"0\n");
                        }
                }

                free(chr1);
		fclose(fout);
		free(curout);
		free(temp);
	}
	return(EXIT_SUCCESS);
}
