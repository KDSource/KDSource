#include<stdio.h>
#include<stdlib.h>

#include "ksource.h"

int main(){
    srand(time(NULL));
    
    PList* plist = PListSimple_create('n', NULL, NULL, "GF12p_inicio", PTRAC_read);

    double bw_E[] = {0.1};
    double bw_pos[] = {0.1,0.1,0.1};
    double bw_dir[] = {0.1};
    double* bw[] = {bw_E, bw_pos, bw_dir};
    double* gp[] = {NULL, NULL, NULL};
    PerturbFun perturb[] = {Let_perturb, Vol_perturb, Dir_perturb};
    Metric* metric = Geom_create(gp, "gaussian", bw, perturb);

	KSource* ks = KS_create(1, plist, metric);

	char pt;
	Part part;
	double w;

    int i;
    for(i=0; i<10; i++){
        KS_sample(ks, &pt, &part, &w, 1);
        printf("%c, %.2e, (%.2lf, %.2lf, %.2lf), (%.2lf, %.2lf, %.2lf), %.2lf\n",
            pt, part.E, part.pos[0], part.pos[1], part.pos[2], part.dir[0], part.dir[1], part.dir[2], w);
    }
    
    return 0;
}