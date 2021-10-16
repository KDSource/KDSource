#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<math.h>
#include<libxml/parser.h>

#include "ksource.h"


void KS_error(const char* msg){
	printf("KSource error: %s\n", msg);
	exit(EXIT_FAILURE);
}
void KS_end(const char* msg){
	printf("KSource terminate: %s\n", msg);
	exit(EXIT_SUCCESS);
}

KSource* KS_create(double J, PList* plist, Geometry* geom){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J;
	ks->plist = plist;
	ks->geom = geom;
	return ks;
}

KSource* KS_open(const char* xmlfilename){
	xmlKeepBlanksDefault(0);

	int i, j, n;
	double J;
	char *buf;

	char pt;    
	char mcplfile[NAME_MAX_LEN];
	double *trasl_plist=NULL, *rot_plist=NULL;

	int order;
	double *trasl_geom=NULL, *rot_geom=NULL;
	int switch_x2z, variable_bw;
	char* bwfilename=NULL;
	double bw=0;

	// Read file
	printf("Reading xmlfile %s...\n", xmlfilename);
	xmlDocPtr doc = xmlReadFile(xmlfilename, NULL, 0);
	if(doc == NULL){
		printf("Could not open file %s\n", xmlfilename);
		KS_error("Error in KS_open");
	}
	xmlNodePtr root = xmlDocGetRootElement(doc);
	if(strcmp((char*)root->name, "KSource") != 0){
		printf("Invalid format in source XML file %s\n", xmlfilename);
		KS_error("Error in KS_open");
	}
	xmlNodePtr node = root->children; // Node: J
	sscanf((char*)xmlNodeGetContent(node), "%lf", &J); // Read J

	// PList
	xmlNodePtr pltree = node->next;
	node = pltree->children; // Node: pt
	sscanf((char*)xmlNodeGetContent(node), "%c", &pt); // Read pt
	node = node->next; // Node: mcplname
	if(strlen((char*)xmlNodeGetContent(node)) > NAME_MAX_LEN){
		printf("mcpl file name %s exceeds NAME_MAX_LEN=%d", (char*)xmlNodeGetContent(node), NAME_MAX_LEN);
		KS_error("Error in KS_open");
	}
	strcpy(mcplfile, (char*)xmlNodeGetContent(node)); // Read mcpl file name
	node = node->next; // Node: trasl
	if(strlen((char*)xmlNodeGetContent(node)) > 1){
		trasl_plist = (double*)malloc(3 * sizeof(double));
		sscanf((char*)xmlNodeGetContent(node), "%lf %lf %lf", &trasl_plist[0], &trasl_plist[1], &trasl_plist[2]);
	}
	node = node->next; // Node: rot
	if(strlen((char*)xmlNodeGetContent(node)) > 1){
		rot_plist = (double*)malloc(3 * sizeof(double));
		sscanf((char*)xmlNodeGetContent(node), "%lf %lf %lf", &rot_plist[0], &rot_plist[1], &rot_plist[2]);
	}
	node = node->next; // Node: x2z
	sscanf((char*)xmlNodeGetContent(node), "%d", &switch_x2z); // Read switch_x2z

	// Geometry
	xmlNodePtr gtree = pltree->next;
	sscanf((char*)xmlGetProp(gtree, (const xmlChar*)"order"), "%d", &order); // Read order
	int dims[order], nps[order];
	char metricnames[order][NAME_MAX_LEN];
	double *params[order];
	double *scalings[order];
	PerturbFun perturbs[order];
	xmlNodePtr mtree = gtree->children;
	for(i=0; i<order; i++){
		strcpy(metricnames[i], (char*)mtree->name); // Read metricname
		node = mtree->children; // Node: dim
		sscanf((char*)xmlNodeGetContent(node), "%d", &dims[i]); // Read dim
		scalings[i] = (double*)malloc(dims[i]*sizeof(double));
		node = node->next; // Node: params
		sscanf((char*)xmlGetProp(node,(const xmlChar*)"nps"), "%d", &nps[i]); // Read ngp
		params[i] = (double*)malloc(nps[i]*sizeof(double));
		buf = (char*)xmlNodeGetContent(node);
		for(j=0; j<nps[i]; j++){ // Read params
			sscanf(buf, "%lf %n", &params[i][j], &n);
			buf += n;
		}
		mtree = mtree->next;
	}
	node = mtree; // Node: trasl
	if(strlen((char*)xmlNodeGetContent(node)) > 1){
		trasl_geom = (double*)malloc(3 * sizeof(double));
		sscanf((char*)xmlNodeGetContent(node), "%lf %lf %lf", &trasl_geom[0], &trasl_geom[1], &trasl_geom[2]);
	}
	node = node->next; // Node: rot
	if(strlen((char*)xmlNodeGetContent(node)) > 1){
		rot_geom = (double*)malloc(3 * sizeof(double));
		sscanf((char*)xmlNodeGetContent(node), "%lf %lf %lf", &rot_geom[0], &rot_geom[1], &rot_geom[2]);
	}
	node = gtree->next; // Node: scaling
	buf = (char*)xmlNodeGetContent(node);
	for(i=0; i<order; i++){ // Read scalings
		for(j=0; j<dims[i]; j++){
			sscanf(buf, "%lf %n", &scalings[i][j], &n);
			buf += n;
		}
	}
	node = node->next; // Node: BW
	sscanf((char*)xmlGetProp(node,(const xmlChar*)"variable"), "%d", &variable_bw); // Read variable_bw
	if(variable_bw){
		bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		if(strlen((char*)xmlNodeGetContent(node)) > NAME_MAX_LEN){
			printf("BW file name %s exceeds NAME_MAX_LEN=%d", (char*)xmlNodeGetContent(node), NAME_MAX_LEN);
			KS_error("Error in KS_open");
		}
		strcpy(bwfilename, (char*)xmlNodeGetContent(node)); // Read BW file name
	}
	else
		sscanf((char*)xmlNodeGetContent(node), "%lf", &bw); // Read BW

	// Create PList
	PList* plist = PList_create(pt, mcplfile, trasl_plist, rot_plist, switch_x2z);
	// Create Metric
	for(i=0; i<order; i++){
		for(j=0; j<_n_metrics; j++){
			if(strcmp(metricnames[i], _metric_names[j]) == 0){
				perturbs[i] = _metric_perturbs[j];
				break;
			}
		}
		if(j == _n_metrics){
			printf("Invalid %s metric.\n", metricnames[i]);
			KS_error("Error in KS_open");
		}
	}
	Metric* metrics[order];
	for(i=0; i<order; i++) metrics[i] = Metric_create(dims[i], scalings[i], perturbs[i], nps[i], params[i]);
	Geometry* geom = Geom_create(order, metrics, bw, bwfilename, trasl_geom, rot_geom);
	// Create KSource
	KSource* s = KS_create(J, plist, geom);

	printf("Done.\n");

	// Free allocated variables
	free(trasl_plist); free(rot_plist);
	free(trasl_geom); free(rot_geom);
	free(bwfilename);
	for(i=0; i<order; i++){ free(params[i]); free(scalings[i]); }
	xmlFreeDoc(doc);
	xmlCleanupParser();

	return s;
}

int KS_sample2(KSource* ks, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop){
	int ret=0;
	if(w_crit <= 0){
		PList_get(ks->plist, part);
		if(perturb) Geom_perturb(ks->geom, part);
		ret += PList_next(ks->plist, loop);
		ret += Geom_next(ks->geom, loop);
	}
	else{ // Normalize w to 1
		double bs;
		int resamples = 0;
		while(1){
			PList_get(ks->plist, part);
			if(bias) bs = bias(part);
			else bs = 1;
			if(part->weight*bs > w_crit){ // If w*bs>w_crit, use w_crit/w*bs as prob of going forward in list
				if(perturb) Geom_perturb(ks->geom, part);
				if(rand() < w_crit/(part->weight*bs)*RAND_MAX){
					ret += PList_next(ks->plist, loop);
					ret += Geom_next(ks->geom, loop);
				}
				break;
			}
			else{ // If w*bs<w_crit, use w*bs/w_crit as prob of taking particle
				int take = 0;
				if(rand() < (part->weight*bs)/w_crit*RAND_MAX){
					take = 1;
					if(perturb) Geom_perturb(ks->geom, part);
				}
				ret += PList_next(ks->plist, loop);
				ret += Geom_next(ks->geom, loop);
				if(take) break;
			}
			if(resamples++ > MAX_RESAMPLES)
				KS_error("Error in KS_sample: MAX_RESAMPLES reached.");
		}
		part->weight = 1/bs;
	}
	if(ret==1 && ks->geom->bwfile)
		printf("Warning: Particle list and bandwidths file have different size.\n");
	return ret;
}

int KS_sample(KSource* ks, mcpl_particle_t* part){
	return KS_sample2(ks, part, 1, 1, NULL, 1);
}

double KS_w_mean(KSource* ks, int N, WeightFun bias){
	int i;
	char pt;
	mcpl_particle_t part;
	double w_mean=0;;
	for(i=0; i<N; i++){
		KS_sample2(ks, &part, 0, -1, NULL, 1);
		if(bias) w_mean += part.weight * bias(&part);
		else w_mean += part.weight;
	}
	return w_mean / N;
}

void KS_destroy(KSource* ks){
	PList_destroy(ks->plist);
	Geom_destroy(ks->geom);
	free(ks);
}


MultiSource* MS_create(int len, KSource** s, const double* ws){
	MultiSource* ms = (MultiSource*)malloc(sizeof(MultiSource));
	ms->len = len;
	ms->s = (KSource**)malloc(len*sizeof(KSource*));
	ms->ws = (double*)malloc(len*sizeof(double));
	ms->J = 0;
	int i;
	for(i=0; i<len; i++){
		ms->s[i] = s[i];
		ms->ws[i] = ws[i];
		ms->J += s[i]->J;
	}
	ms->cdf = (double*)malloc(ms->len*sizeof(double));
	for(i=0; i<ms->len; i++) ms->cdf[i] = ms->ws[i];
	for(i=1; i<ms->len; i++) ms->cdf[i] += ms->cdf[i-1];
	return ms;
}

MultiSource* MS_open(int len, const char** xmlfilenames, const double* ws){
	KSource* s[len];
	int i;
	for(i=0; i<len; i++) s[i] = KS_open(xmlfilenames[i]);
	return MS_create(len, s, ws);
}

int MS_sample2(MultiSource* ms, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop){
	double y = rand() / ((double)RAND_MAX+1);
	int i, ret;
	if(ms->cdf[ms->len-1] <= 0) i = (int)(y*ms->len);
	else for(i=0; y*ms->cdf[ms->len-1]>ms->cdf[i]; i++);
	ret = KS_sample2(ms->s[i], part, perturb, w_crit, bias, loop);
	if(ms->cdf[ms->len-1] > 0) part->weight *= (ms->s[i]->J/ms->J) / (ms->ws[i]/ms->cdf[ms->len-1]);
	else part->weight *= (ms->s[i]->J/ms->J) * ms->len;
	return ret;
}

int MS_sample(MultiSource* ms, mcpl_particle_t* part){
	return MS_sample2(ms, part, 1, 1, NULL, 1);
}

double MS_w_mean(MultiSource* ms, int N, WeightFun bias){
	double w_mean=0;
	int i;
	for(i=0; i<ms->len; i++) w_mean += ms->ws[i] * KS_w_mean(ms->s[i], N, bias);
	return w_mean / ms->cdf[ms->len-1];
}

void MS_destroy(MultiSource* ms){
	int i;
	for(i=0; i<ms->len; i++) KS_destroy(ms->s[i]);
	free(ms->s);
	free(ms->ws);
	free(ms->cdf);
	free(ms);
}
