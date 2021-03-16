#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "../ksource_c/ksource.h"

int main(){
    
    double xylims[4] = {-30,30,-30,30};
	PList* plist = PListSimple_create('n', "tracks.ssv", SSV_read, NULL, NULL, 0);
    Metric* metrics[3] = {Metric_create(0, NULL, Let_perturb, 0, NULL),
                          Metric_create(2, NULL, SurfXY_perturb, 4, xylims),
                          Metric_create(3, NULL, Isotrop_perturb, 0, NULL)};
    MetricSepVar* metric = MetricSepVar_create(3, metrics, "bw.txt", 0, NULL, NULL);
	KSource* ks = KS_create(1, plist, metric);

	int N=1e5;
	Part part;
	double w;
	char pt;
	FILE* file = fopen("resampled.ssv", "w");
	int i;
	for(i=0; i<N; i++){
		KS_sample(ks, &pt, &part, &w, 1, NULL);
		fprintf(file, "%le %le %le %le %le %le %le %le\n",
			part.E, part.pos[0], part.pos[1], part.pos[2], part.dir[0], part.dir[1], part.dir[2], w);
	}
	fclose(file);
	printf("Resampleo exitoso\n");
    
    return 0;
}