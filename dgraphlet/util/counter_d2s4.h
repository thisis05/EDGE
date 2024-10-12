#pragma once

#include "digraph.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

struct TriangleInfo {
    VertexIdx* triend;
    Count count;
    Count maxTri;

    TriangleInfo(Count degree) {
        maxTri = degree;
        count = 0;
        triend = new VertexIdx[maxTri]();
    }

    ~TriangleInfo() {
        delete[] triend;
    }
};

struct FourSizeInfo {
    Count n;
    Count e1;
    Count e2;
    Count t1, t2, t3, t4;
    Count clique1, clique2, clique3, clique4, clique5, clique6, clique7, clique8, clique9, clique10, clique11;
    Count chord1, chord2, chord3, chord4, chord5, chord6, chord7, chord8; 
    Count tailed1, tailed2, tailed3, tailed4, tailed5, tailed6, tailed7, tailed8;
    Count cycle1, cycle2, cycle3;
    Count star1, star2;
    Count path1, path2, path3, path4;

    FourSizeInfo(VertexIdx nVertices, EdgeIdx nEdges, EdgeIdx nEdges2) 
        : t1(0), t2(0), t3(0), t4(0),
        clique1(0), clique2(0), clique3(0), clique4(0), clique5(0), clique6(0), clique7(0), clique8(0), clique9(0), clique10(0), clique11(0),
        chord1(0), chord2(0), chord3(0), chord4(0), chord5(0), chord6(0), chord7(0), chord8(0),
        tailed1(0), tailed2(0), tailed3(0), tailed4(0), tailed5(0), tailed6(0), tailed7(0), tailed8(0),
        cycle1(0), cycle2(0), cycle3(0),
        star1(0), star2(0),
        path1(0), path2(0), path3(0), path4(0) {
        n = nVertices;
        e1 = nEdges;
        e2 = nEdges2;
    }
};

struct StarInfo {
    VertexIdx* star1;
    VertexIdx* star2_1;
    VertexIdx* star2_2;
    VertexIdx* tri3;

    StarInfo(VertexIdx nv) {
        star1 = new VertexIdx[nv+1]();
        star2_1 = new VertexIdx[nv+1]();
        star2_2 = new VertexIdx[nv+1]();
        tri3 = new VertexIdx[nv+1]();
    }

    ~StarInfo() {
        delete[] star1;
        delete[] star2_1;
        delete[] star2_2;
        delete[] tri3;
    }
};


FourSizeInfo get4size(CGraph *gout, CGraph *gin, CGraph *gin_2, CGraph *gout_2, int num_threads);

void countFour(CDAG *dag, CDAG *dag_2, double (&gcounts)[36], int num_threads);
void comb_Four(double (&gcounts)[36]);
void print4size(double (&gcounts)[36]);