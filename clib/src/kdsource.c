#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml/parser.h>
#include <math.h>

#include "kdsource/kdsource.h"

void KDS_error(const char *msg) {
  printf("KDSource error: %s\n", msg);
  exit(EXIT_FAILURE);
}

void *KDS_malloc(size_t n_bytes) {
  void *a = malloc(n_bytes ? n_bytes : 1);
  if (!a)
    KDS_error("memory allocation failed");
  return a;
}

void *KDS_calloc(size_t n_elem, size_t elem_size) {
  if (elem_size == 0)
    KDS_error("trying to allocate space for elements of size 0");
  void *a = calloc(n_elem, elem_size);
  if (!a)
    KDS_error("memory allocation failed");
  return a;
}

void *KDS_calloc_i(int n_elem, size_t elem_size) {
  if (n_elem < 0)
    KDS_error("trying to allocate space for negative number of elements");
  return KDS_calloc((size_t)(n_elem), elem_size);
}

double KDS_default_rngfct(void *rngstate) {
  // This is a horrible RNG. Also, simply uses global state!
  (void)rngstate;
  return rand() / (RAND_MAX + 1.0);
}

void KDS_end(const char *msg) {
  printf("KDSource terminate: %s\n", msg);
  exit(EXIT_SUCCESS);
}

void KDS_strcpy(char *dst, size_t dst_buflen, const char *src) {
  // Safe + portable + C99 compatible implementation of strcpy
  if (!dst || !src || !dst_buflen) {
    KDS_error("invalid strcpy (null arg)");
    return;
  }
  size_t srclen = 0;
  while (srclen < dst_buflen && src[srclen])
    ++srclen;
  if (!(srclen < dst_buflen)) {
    KDS_error("invalid strcpy (too small buffer)");
    return;
  }
  if (srclen > 1)
    memcpy(dst, src, srclen);
  dst[srclen] = 0;
}

void KDS_strcpyname(char *dst, const char *src) {
  KDS_strcpy(dst, NAME_MAX_LEN, src);
}

KDSource *KDS_create(double J, char kernel, PList *plist, Geometry *geom) {
  KDSource *kds = (KDSource *)KDS_malloc(sizeof(KDSource));
  kds->J = J;
  kds->kernel = kernel;
  kds->plist = plist;
  kds->geom = geom;
  return kds;
}

KDSource *KDS_open(const char *xmlfilename) {
  xmlKeepBlanksDefault(0);

  int i, j, n;
  double J;
  char kernel = 'g';
  char *buf;

  char pt;
  char mcplfile[NAME_MAX_LEN];
  double *trasl_plist = NULL, *rot_plist = NULL;

  int order;
  double *trasl_geom = NULL, *rot_geom = NULL;
  int switch_x2z, variable_bw;
  char *bwfilename = NULL;
  double bw = 0;

  // Read file
  printf("Reading xmlfile %s...\n", xmlfilename);
  xmlDocPtr doc = xmlReadFile(xmlfilename, NULL, 0);
  if (doc == NULL) {
    printf("Could not open file %s\n", xmlfilename);
    KDS_error("Error in KDS_open");
  }
  xmlNodePtr root = xmlDocGetRootElement(doc);
  if (strcmp((char *)root->name, "KDSource") != 0) {
    printf("Invalid format in source XML file %s\n", xmlfilename);
    KDS_error("Error in KDS_open");
  }
  xmlNodePtr node = root->children;                   // Node: J
  sscanf((char *)xmlNodeGetContent(node), "%lf", &J); // Read J

  node = node->next; // Node: kernel
  if (strcmp((char *)node->name, "kernel") == 0) {
    sscanf((char *)xmlNodeGetContent(node), "%s", &kernel); // Read kernel
    node = node->next;
  } else
    printf("No kernel specified. Using gaussian as default.\n");

  // PList
  xmlNodePtr pltree = node;
  node = pltree->children;                            // Node: pt
  sscanf((char *)xmlNodeGetContent(node), "%c", &pt); // Read pt
  node = node->next;                                  // Node: mcplname
  if (strlen((char *)xmlNodeGetContent(node)) > NAME_MAX_LEN) {
    printf("mcpl file name %s exceeds NAME_MAX_LEN=%d",
           (char *)xmlNodeGetContent(node), NAME_MAX_LEN);
    KDS_error("Error in KDS_open");
  }
  KDS_strcpyname(mcplfile,
                 (char *)xmlNodeGetContent(node)); // Read mcpl file name
  node = node->next;                               // Node: trasl
  if (strlen((char *)xmlNodeGetContent(node)) > 1) {
    trasl_plist = (double *)KDS_malloc(3 * sizeof(double));
    sscanf((char *)xmlNodeGetContent(node), "%lf %lf %lf", &trasl_plist[0],
           &trasl_plist[1], &trasl_plist[2]);
  }
  node = node->next; // Node: rot
  if (strlen((char *)xmlNodeGetContent(node)) > 1) {
    rot_plist = (double *)KDS_malloc(3 * sizeof(double));
    sscanf((char *)xmlNodeGetContent(node), "%lf %lf %lf", &rot_plist[0],
           &rot_plist[1], &rot_plist[2]);
  }
  node = node->next;                                          // Node: x2z
  sscanf((char *)xmlNodeGetContent(node), "%d", &switch_x2z); // Read switch_x2z

  // Geometry
  xmlNodePtr gtree = pltree->next;
  sscanf((char *)xmlGetProp(gtree, (const xmlChar *)"order"), "%d",
         &order); // Read order
  if (order < 1)
    KDS_error("Error in KDS_open: invalid value of \"order\"");
  const size_t uorder = (unsigned)(order);

  int *dims = KDS_calloc(uorder, sizeof(int));
  int *nps = KDS_calloc(uorder, sizeof(int));
  double **params = KDS_calloc(uorder, sizeof(double *));
  double **scalings = KDS_calloc(uorder, sizeof(double *));
  char(*metricnames)[NAME_MAX_LEN] =
      KDS_calloc(uorder, NAME_MAX_LEN * sizeof(char));
  PerturbFun *perturbs = KDS_calloc(uorder, sizeof(PerturbFun));
  Metric **metrics = KDS_calloc(uorder, sizeof(Metric *));

  if (!dims || !nps || !params || !scalings || !perturbs || !metricnames ||
      !metrics) {
    free(dims);
    free(nps);
    free(params);
    free(scalings);
    free(metricnames);
    free(perturbs);
    free(metrics);
    dims = NULL;
    KDS_error("Error in KDS_open: memory allocation failure");
  }
  if (!dims)
    return NULL;

  xmlNodePtr mtree = gtree->children;
  for (i = 0; i < order; i++) {
    KDS_strcpyname(metricnames[i], (char *)mtree->name);     // Read metricname
    node = mtree->children;                                  // Node: dim
    sscanf((char *)xmlNodeGetContent(node), "%d", &dims[i]); // Read dim

    if (dims[i] < 1)
      KDS_error("Error in KDS_open: invalid value of \"dims[i]\"");
    scalings[i] = (double *)KDS_calloc_i(dims[i], sizeof(double));
    node = node->next; // Node: params
    sscanf((char *)xmlGetProp(node, (const xmlChar *)"nps"), "%d",
           &nps[i]); // Read ngp
    params[i] = (double *)KDS_calloc_i(nps[i], sizeof(double));
    buf = (char *)xmlNodeGetContent(node);
    for (j = 0; j < nps[i]; j++) { // Read params
      sscanf(buf, "%lf %n", &params[i][j], &n);
      buf += n;
    }
    mtree = mtree->next;
  }
  node = mtree; // Node: trasl
  if (strlen((char *)xmlNodeGetContent(node)) > 1) {
    trasl_geom = (double *)KDS_malloc(3 * sizeof(double));
    sscanf((char *)xmlNodeGetContent(node), "%lf %lf %lf", &trasl_geom[0],
           &trasl_geom[1], &trasl_geom[2]);
  }
  node = node->next; // Node: rot
  if (strlen((char *)xmlNodeGetContent(node)) > 1) {
    rot_geom = (double *)KDS_malloc(3 * sizeof(double));
    sscanf((char *)xmlNodeGetContent(node), "%lf %lf %lf", &rot_geom[0],
           &rot_geom[1], &rot_geom[2]);
  }
  node = gtree->next; // Node: scaling
  buf = (char *)xmlNodeGetContent(node);
  for (i = 0; i < order; i++) { // Read scalings
    for (j = 0; j < dims[i]; j++) {
      sscanf(buf, "%lf %n", &scalings[i][j], &n);
      buf += n;
    }
  }
  node = node->next; // Node: BW
  sscanf((char *)xmlGetProp(node, (const xmlChar *)"variable"), "%d",
         &variable_bw); // Read variable_bw
  if (variable_bw) {
    bwfilename = (char *)KDS_malloc(NAME_MAX_LEN * sizeof(char));
    if (strlen((char *)xmlNodeGetContent(node)) > NAME_MAX_LEN) {
      printf("BW file name %s exceeds NAME_MAX_LEN=%d",
             (char *)xmlNodeGetContent(node), NAME_MAX_LEN);
      KDS_error("Error in KDS_open");
    }
    KDS_strcpyname(bwfilename,
                   (char *)xmlNodeGetContent(node)); // Read BW file name
  } else
    sscanf((char *)xmlNodeGetContent(node), "%lf", &bw); // Read BW

  // Create PList
  PList *plist = PList_create(pt, mcplfile, trasl_plist, rot_plist, switch_x2z);
  // Create Metric
  for (i = 0; i < order; i++) {
    for (j = 0; j < _n_metrics; j++) {
      if (strcmp(metricnames[i], _metric_names[j]) == 0) {
        perturbs[i] = _metric_perturbs[j];
        break;
      }
    }
    if (j == _n_metrics) {
      printf("Invalid %s metric.\n", metricnames[i]);
      KDS_error("Error in KDS_open");
    }
  }
  for (i = 0; i < order; i++)
    metrics[i] =
        Metric_create(dims[i], scalings[i], perturbs[i], nps[i], params[i]);
  Geometry *geom =
      Geom_create(order, metrics, bw, bwfilename, kernel, trasl_geom, rot_geom);
  // Create KDSource
  KDSource *s = KDS_create(J, kernel, plist, geom);

  printf("Done.\n");

  // Free allocated variables
  free(trasl_plist);
  free(rot_plist);
  free(trasl_geom);
  free(rot_geom);
  free(bwfilename);
  for (i = 0; i < order; i++) {
    free(params[i]);
    free(scalings[i]);
  }
  xmlFreeDoc(doc);
  xmlCleanupParser();

  free(dims);
  free(nps);
  free(params);
  free(scalings);
  free(metricnames);
  free(perturbs);
  free(metrics);

  return s;
}

int KDS_rand_sample2(kds_rng_fct_t rng, void *rngstate, KDSource *kds,
                     mcpl_particle_t *part, int perturb, double w_crit,
                     WeightFun bias, int loop) {
  int ret = 0;
  if (w_crit <= 0) {
    PList_get(kds->plist, part);
    if (perturb)
      Geom_perturb(rng, rngstate, kds->geom, part);
    ret += PList_next(kds->plist, loop);
    ret += Geom_next(kds->geom, loop);
  } else { // Normalize w to 1
    double bs;
    int resamples = 0;
    while (1) {
      PList_get(kds->plist, part);
      if (bias)
        bs = bias(part);
      else
        bs = 1;
      if (part->weight * bs > w_crit) { // If w*bs>w_crit, use w_crit/w*bs as
                                        // prob of going forward in list
        if (perturb)
          Geom_perturb(rng, rngstate, kds->geom, part);
        if (rng(rngstate) < w_crit / (part->weight * bs)) {
          ret += PList_next(kds->plist, loop);
          ret += Geom_next(kds->geom, loop);
        }
        break;
      } else { // If w*bs<w_crit, use w*bs/w_crit as prob of taking particle
        int take = 0;
        if (rng(rngstate) < (part->weight * bs) / w_crit) {
          take = 1;
          if (perturb)
            Geom_perturb(rng, rngstate, kds->geom, part);
        }
        ret += PList_next(kds->plist, loop);
        ret += Geom_next(kds->geom, loop);
        if (take)
          break;
      }
      if (resamples++ > MAX_RESAMPLES)
        KDS_error("Error in KDS_sample: MAX_RESAMPLES reached.");
    }
    part->weight = 1 / bs;
  }
  if (ret == 1 && kds->geom->bwfile)
    printf("Warning: Particle list and bandwidths file have different size.\n");
  return ret;
}

int KDS_rand_sample(kds_rng_fct_t rng, void *rngstate, KDSource *kds,
                    mcpl_particle_t *part) {
  return KDS_rand_sample2(rng, rngstate, kds, part, 1, 1, NULL, 1);
}

int KDS_sample2(KDSource *kds, mcpl_particle_t *part, int perturb,
                double w_crit, WeightFun bias, int loop) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return KDS_rand_sample2(rng, NULL, kds, part, perturb, w_crit, bias, loop);
}

int KDS_sample(KDSource *kds, mcpl_particle_t *part) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return KDS_rand_sample(rng, NULL, kds, part);
}

double KDS_w_mean(KDSource *kds, int N, WeightFun bias) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return KDS_rand_w_mean(rng, NULL, kds, N, bias);
}

double KDS_rand_w_mean(kds_rng_fct_t rng, void *rngstate, KDSource *kds, int N,
                       WeightFun bias) {
  int i;
  // char pt;
  mcpl_particle_t part;
  double w_mean = 0;
  for (i = 0; i < N; i++) {
    KDS_rand_sample2(rng, rngstate, kds, &part, 0, -1, NULL, 1);
    if (bias)
      w_mean += part.weight * bias(&part);
    else
      w_mean += part.weight;
  }
  return w_mean / N;
}

void KDS_destroy(KDSource *kds) {
  PList_destroy(kds->plist);
  Geom_destroy(kds->geom);
  free(kds);
}

MultiSource *MS_create(int len, KDSource **s, const double *ws) {
  MultiSource *ms = (MultiSource *)KDS_malloc(sizeof(MultiSource));

  if (len < 1)
    KDS_error("Error in MS_create: invalid value of \"len\"");
  ms->len = len;
  ms->s = (KDSource **)KDS_calloc_i(len, sizeof(KDSource *));
  ms->ws = (double *)KDS_calloc_i(len, sizeof(double));
  ms->J = 0;
  int i;
  for (i = 0; i < len; i++) {
    ms->s[i] = s[i];
    ms->ws[i] = ws[i];
    ms->J += s[i]->J;
  }
  ms->cdf = (double *)KDS_calloc_i(len, sizeof(double));
  for (i = 0; i < ms->len; i++)
    ms->cdf[i] = ms->ws[i];
  for (i = 1; i < ms->len; i++)
    ms->cdf[i] += ms->cdf[i - 1];
  return ms;
}

MultiSource *MS_open(int len, const char **xmlfilenames, const double *ws) {
  if (len < 1)
    KDS_error("Error in MS_open: invalid value of \"len\"");
  KDSource **s = KDS_calloc_i(len, sizeof(KDSource *));
  if (!s) {
    len = 0;
    KDS_error("Error in MS_open: memory allocation failure");
  }
  if (len == 0)
    return NULL;
  int i;
  for (i = 0; i < len; i++)
    s[i] = KDS_open(xmlfilenames[i]);
  MultiSource *result = MS_create(len, s, ws);
  free(s);
  return result;
}

int MS_rand_sample2(kds_rng_fct_t rng, void *rngstate, MultiSource *ms,
                    mcpl_particle_t *part, int perturb, double w_crit,
                    WeightFun bias, int loop) {
  double y = rng(rngstate);
  int i, ret;
  if (ms->cdf[ms->len - 1] <= 0)
    i = (int)(y * ms->len);
  else
    for (i = 0; y * ms->cdf[ms->len - 1] > ms->cdf[i]; i++)
      ;
  ret = KDS_rand_sample2(rng, rngstate, ms->s[i], part, perturb, w_crit, bias,
                         loop);
  if (ms->cdf[ms->len - 1] > 0)
    part->weight *= (ms->s[i]->J / ms->J) / (ms->ws[i] / ms->cdf[ms->len - 1]);
  else
    part->weight *= (ms->s[i]->J / ms->J) * ms->len;
  return ret;
}

int MS_rand_sample(kds_rng_fct_t rng, void *rngstate, MultiSource *ms,
                   mcpl_particle_t *part) {
  return MS_rand_sample2(rng, rngstate, ms, part, 1, 1, NULL, 1);
}

int MS_sample(MultiSource *ms, mcpl_particle_t *part) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return MS_rand_sample(rng, NULL, ms, part);
}

int MS_sample2(MultiSource *ms, mcpl_particle_t *part, int perturb,
               double w_crit, WeightFun bias, int loop) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return MS_rand_sample2(rng, NULL, ms, part, perturb, w_crit, bias, loop);
}

double MS_w_mean(MultiSource *ms, int N, WeightFun bias) {
  kds_rng_fct_t rng = KDS_default_rngfct;
  return MS_rand_w_mean(rng, NULL, ms, N, bias);
}

double MS_rand_w_mean(kds_rng_fct_t rng, void *rngstate, MultiSource *ms, int N,
                      WeightFun bias) {
  double w_mean = 0;
  int i;
  for (i = 0; i < ms->len; i++)
    w_mean += ms->ws[i] * KDS_rand_w_mean(rng, rngstate, ms->s[i], N, bias);
  return w_mean / ms->cdf[ms->len - 1];
}

void MS_destroy(MultiSource *ms) {
  int i;
  for (i = 0; i < ms->len; i++)
    KDS_destroy(ms->s[i]);
  free(ms->s);
  free(ms->ws);
  free(ms->cdf);
  free(ms);
}
