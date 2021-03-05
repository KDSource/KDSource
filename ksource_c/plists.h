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

	int ord; // Cantidad de archivos de tracks
	FILE** files; // Archivos de tracks (unico o uno por variable)
	ReadFun* read; // Funciones de lectura

	double* trasl; // Traslacion a aplicar
	double* rot; // Rotacion a aplicar
	int x2z; // Si es true, se aplica permutacion x,y,z -> y,z,x

	char line[LINE_MAX_LEN]; // Buffer para linea de texto
	Part* part; // Buffer para particula
	double w; // Buffer para peso
} PList;

PList* PList_create(char pt, int ord, char** filenames, ReadFun* read, double* trasl, double* rot, int switch_x2z);
int PList_get(PList* plist, Part* part, double* w);
int PList_next(PList* plist);
void PList_destroy(PList* plist);

PList* PListSimple_create(char pt, char* filename, ReadFun read, double* trasl, double* rot, int switch_x2z);
PList* PListSepVar_create(char pt, char* filenames[3], ReadFun read[3], double* trasl, double* rot, int switch_x2z);


int PTRAC_read(char* line, Part* part, double* w);

int T4stock_read(char* line, Part* part, double* w);

int SSV_read(char* line, Part* part, double* w);

int Decay_read(char* line, Part* part, double* w);
int T4pos_read(char* line, Part* part, double* w);
int Isotrop_read(char* line, Part* part, double* w);


#endif
