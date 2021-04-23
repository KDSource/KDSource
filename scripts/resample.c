#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"

#define N_DEFAULT 1E5


void display_usage(){
	printf("Uso: ksource resample sourcefile [opciones]\n\n");
	printf("Resamplea muestras de la fuente definida en sourcefile, y las guarda\n");
	printf("en un archivo MCPL.\n\n");
	printf("Opciones:\n");
	printf("\t-o outfile:  nombre del archivo MCPL con las nuevas muestras\n");
	printf("\t             (default: \"resampled.mcpl\")\n");
	printf("\t-n N:        cantidad de nuevas muestras (default: 1E5).\n");

}

int resample_parse_args(int argc, char **argv, const char** filename, const char** outfilename, int* N){
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
			*N = atoi(argv[++i]);
			continue;
		}
    	if(argv[i][0] == '-'){
			printf("Error: Argumento invalido: %s\n",argv[i]);
			exit(1);
    	}
		if(!*filename){
			*filename = argv[i];
			continue;
		}
		printf("Demasiados argumentos. Use -h o --help para ayuda.\n");
		exit(1);
	}
	if(!*filename){
		printf("No se especifico archivo de fuente. Use -h o --help para ayuda.\n");
		exit(1);
	}
	if(!*outfilename) *outfilename = "resampled.mcpl";
	return 0;
}

int main(int argc, char *argv[]){
	const char *filename;
	const char *outfilename;
	int N;

	if(resample_parse_args(argc, argv, &filename, &outfilename, &N)) return 1;

    KSource* ks = KS_open(filename);
	mcpl_particle_t part;

	mcpl_outfile_t file = mcpl_create_outfile(outfilename);
	mcpl_hdr_set_srcname(file, "KSource resample");

	double w_crit = KS_w_mean(ks, 1000, NULL);

	printf("Resampleando...\n");
	int i;
	for(i=0; i<N; i++){
		KS_sample(ks, &part, 1, w_crit, NULL);
		mcpl_add_particle(file, &part);
	}
	mcpl_closeandgzip_outfile(file);
	printf("Resampleo exitoso\n");

    return 0;
}