#!/bin/bash

KS=../ksource_c
gcc resample.c $KS/ksource.c $KS/metrics.c $KS/plists.c $KS/aux.c -lm -o resample.out
./resample.out