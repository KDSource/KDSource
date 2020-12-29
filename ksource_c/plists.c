#include<stdio.h>
#include<stdlib.h>

#include "ksource.h"


PList* PList_create(char pt_, FILE** files_, char* lines_[LINE_MAX_LEN], double* pos_, double* rot_, ReadFun read_){

	PList* plist = (PList*)malloc(sizeof(PList));
	plist->pt = pt_;
	plist->files = files_;
	plist->lines = lines_;
	plist->pos = pos_;
	plist->rot = rot_;
	plist->read = read_;
	return plist;
}

void PList_destroy(PList* plist){
	free(plist);
}