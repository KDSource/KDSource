#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<math.h>

#include "kdsource.h"


void display_usage(){
	fprintf(stderr, "Usage: kdtool resample sourcefile [options]\n\n");
	fprintf(stderr, "Resample particles from source defined in XML file sourcefile, and save them in\n");
	fprintf(stderr, "a MCPL file.\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-o outfile: Name of MCPL file with new samples\n");
	fprintf(stderr, "\t            (default: \"resampled.mcpl\").\n");
	fprintf(stderr, "\t-n N:       Number of new samples (default: 1E5).\n");
	fprintf(stderr, "\t-h, --help: Display usage instructions.\n");

}

int resample_parse_args(int argc, char **argv, const char** filename, const char** outfilename, long int* N){
	*filename = 0;
	*outfilename = 0;
	*N = 1E5;
	int i;
	for(i=1; i<argc; i++){
		if(argv[i][0] == '\0')
			continue;
		if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0){
			display_usage();
			exit(0);
		}
		if(strcmp(argv[i],"-o") == 0){
			*outfilename = argv[++i];
			continue;
		}
		if(strcmp(argv[i],"-n") == 0){
			*N = atof(argv[++i]);
			continue;
		}
		if(argv[i][0] == '-'){
			fprintf(stderr, "Error: Invalid argument: %s\n",argv[i]);
			exit(1);
		}
		if(!*filename){
			*filename = argv[i];
			continue;
		}
		fprintf(stderr, "Too many arguments. Use -h or --help for help.\n");
		exit(1);
	}
	if(!*filename){
		fprintf(stderr, "No source file. Use -h or --help for help.\n");
		exit(1);
	}
	if(!*outfilename) *outfilename = "resampled.mcpl";
	return 0;
}

int main(int argc, char *argv[]){
	const char *filename;
	const char *outfilename;
	long int N;

	if(resample_parse_args(argc, argv, &filename, &outfilename, &N)) return 1;

	KDSource* kds = KDS_open(filename);
	mcpl_particle_t part;

	mcpl_outfile_t file = mcpl_create_outfile(outfilename);
	mcpl_hdr_set_srcname(file, "KDSource resample");

	double w_crit = KDS_w_mean(kds, 1000, NULL);

	fprintf(stderr, "Resampling...\n");
	long int i;
	for(i=0; i<N; i++){
		KDS_sample2(kds, &part, 1, w_crit, NULL, 1);
		mcpl_add_particle(file, &part);
	}
	mcpl_closeandgzip_outfile(file);
	fprintf(stderr, "Successfully sampled %ld particles.\n", N);

	return 0;
}
