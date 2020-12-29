#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef KSOURCE_H
#define KSOURCE_H


#define LINE_MAX_LEN 256

typedef int (*ReadFun)(char* lines[LINE_MAX_LEN], Part* part, double* w);

struct PList{
	char pt;
	FILE** files;
	char* lines[LINE_MAX_LEN];
	Part part;
	double w;
	double pos[3];
	double rot[3];

	ReadFun read;
};

PList* PList_create(char pt_, FILE** files_, char* lines_[LINE_MAX_LEN], double* pos_, double* rot_, ReadFun read_);
int PList_get(PList* plist, Part* part, double* w);
int PList_next(PList* plist);
int PList_rewind(PList* plist);
void PList_destroy(PList* plist);


#endif