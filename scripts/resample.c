#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"

#define N_DEFAULT 1E5


int resample_parse_args(int argc, char **argv, char* filename, char* outfilename, int* N){
	if(argc <= 1){
		printf("Usage: resample filename [-o outfilename -n Nresamples]\n");
		return 1;
	}
	strcpy(filename, argv[1]);
	strcpy(outfilename, filename);
	outfilename[strcspn(outfilename, ".")] = 0;
	strcat(outfilename, "_resampled");
	*N = 1E5;
	int i;
	for(i=2; i<argc; i++){
		if(argv[i][0] == '\0') continue;
		if(strcmp(argv[i],"-o") == 0) strcpy(outfilename, argv[i+1]);
		if(strcmp(argv[i],"-n") == 0) *N = atoi(argv[i+1]);
	}
	return 0;
}

int main(int argc, char *argv[]){
	char filename[NAME_MAX_LEN] = "";
	char outfilename[NAME_MAX_LEN];
	int N;

	if(resample_parse_args(argc, argv, filename, outfilename, &N)) return 1;

    KSource* ks = KS_open(filename, 0);
	mcpl_particle_t part;

	mcpl_outfile_t file = mcpl_create_outfile(outfilename);
	mcpl_hdr_set_srcname(file, "KSource resample");

	printf("Resampleando...\n");
	int i;
	for(i=0; i<N; i++){
		KS_sample(ks, &part, 1, NULL);
		mcpl_add_particle(file, &part);
	}
	mcpl_closeandgzip_outfile(file);
	printf("Resampleo exitoso\n");

    return 0;
}