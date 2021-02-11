#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef PLISTS_H
#define PLISTS_H


#define LINE_MAX_LEN 256

typedef struct Part Part;

typedef int (*ReadFun)(char* line, Part* part, double* w);

typedef struct PList{
	char pt; // Tipo de particula (n, p, ...)
	double* trasl; // Traslacion a aplicar
	double* rot; // Rotacion a aplicar
	int x2z; // Si es true, se aplica permutacion x,y,z -> y,z,x

	int n; // Cantidad de archivos de tracks
	FILE** files; // Archivos de tracks (unico o uno por variable)
	char* line; // Buffer para linea de texto
	Part* part; // Buffer para particula
	double w; // Buffer para peso

	ReadFun* read;
} PList;

PList* PList_create(char pt, double* trasl, double* rot, int switch_x2z, int n, FILE** files, ReadFun* read);
int PList_get(PList* plist, Part* part, double* w);
int PList_next(PList* plist);
void PList_destroy(PList* plist);

PList* PListSimple_create(char pt, double* trasl, double* rot, int switch_x2z, char* filename, ReadFun read);
void PListSimple_destroy(PList* plist);

PList* PListSepVar_create(char pt, double* trasl, double* rot, int switch_x2z, char* filename[3], ReadFun read[3]);
void PListSepVar_destroy(PList* plist);

int PTRAC_read(char* line, Part* part, double* w);

int Tripoli_read_part(char* line, Part* part, double* w);

int Decay_read(char* line, Part* part, double* w);
int Tripoli_read_pos(char* line, Part* part, double* w);
int Isotrop_read(char* line, Part* part, double* w);


#endif
