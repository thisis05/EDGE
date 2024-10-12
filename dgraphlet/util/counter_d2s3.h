#pragma once

#include "digraph.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

struct ThreeSizeInfo {
    Count t1;
    Count t2;
    Count t3;
    Count t4;
    Count t5;
    Count t6;

    ThreeSizeInfo() 
        : t1(0), t2(0), t3(0), t4(0), t5(0), t6(0) {}

};

ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2, int num_threads);
ThreeSizeInfo get3size_b1(CGraph *g, CGraph *g_2, int num_threads);
ThreeSizeInfo get3size_b2(CGraph *gout, CGraph *gout_2, CGraph *gin, CGraph *gin_2, int num_threads);

void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&gcounts)[6], int num_threads);
void countThree_b1(CGraph *cg, CGraph *cg_2, double (&gcounts)[6], int num_threads);
void countThree_b2(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&gcounts)[6], int num_threads);

void comb_three(double (&gcounts)[6]);
void print3size(double (&gcounts)[6]);