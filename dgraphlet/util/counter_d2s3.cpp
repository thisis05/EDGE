#include "counter_d2s3.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <stdexcept>
#include <cstring>
#ifdef _WIN32
    #include <windows.h>
    #include <psapi.h>
#endif
#include <iostream>

using namespace std;

ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2, int num_threads) {

    omp_set_num_threads(num_threads); 
    printf("# using threads : %d\n", num_threads);
    
    ThreeSizeInfo ret{};
    VertexIdx current = 0;
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret{};
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                // T4 : (di, dj, dk) = (1, 1, 1)
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.t4++;
                    }
                }

                // T2 : (1) (di, dj, dk) = (1, 2, 2) + (2) (di, dj, dk) = (2, 1, 2)
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.t2++;
                    }
                }                 
            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];

                    // T1 : (di, dj, dk) = (2, 2, 2)
                    EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.t1++;
                    } else {
                        
                        // T2 : (di, dj, dk) = (2, 2, 1)
                        EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                        if (loc_221 != -1) {
                            local_ret.t2++;
                        }
                    }
                }
            }    
        }

        #pragma omp critical
        {
            ret.t1 += local_ret.t1;
            ret.t2 += local_ret.t2;
            ret.t4 += local_ret.t4;
        }
    }
    return ret;
}

ThreeSizeInfo get3size_b1(CGraph *g, CGraph *g_2, int num_threads) {

    omp_set_num_threads(num_threads); 
    printf("# using threads : %d\n", num_threads);
    
    ThreeSizeInfo ret{};
    VertexIdx current = 0;
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret{};
        
        #pragma omp for schedule(guided)
        for (VertexIdx i = 0; i < g->nVertices; ++i) {
            const EdgeIdx start = g->offsets[i];
            const EdgeIdx end = g->offsets[i+1];
            const EdgeIdx start_2 = g_2->offsets[i];
            const EdgeIdx end_2 = g_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = g->nbors[j];

                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = g->nbors[k];
                    EdgeIdx loc_111 = g->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.t4++;
                    }else{
                        local_ret.t3++;
                    }
                }
            
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = g_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? g_2->getEdgeBinary(end2, end1) : g_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.t2++;
                    }
                    else{
                        EdgeIdx loc121 = end1 > end2 ? g->getEdgeBinary(end2, end1) : g->getEdgeBinary(end1, end2);
                        if (loc121 != -1) {
                            local_ret.t3++;
                        }else{
                            local_ret.t6++;
                        }
                    }
                }

                for (EdgeIdx k = g_2->offsets[end1];  k < g_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_2->nbors[k];
                    if (g_2->getEdgeBinary(i, end2) == -1){
                        if (g->getEdgeBinary(i, end2) == -1){
                            local_ret.t6++;
                        }
                    }
                }

            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = g_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = g_2->nbors[k];
                    EdgeIdx loc_222 = g_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.t1++;
                    }
                    else{
                        EdgeIdx loc_221 = g->getEdgeBinary(end1, end2);
                        if(loc_221 != -1){
                            local_ret.t2++;
                        }else{
                            local_ret.t5++;
                        }
                    }
                }

                for (EdgeIdx k = g_2->offsets[end1];  k < g_2->offsets[end1+1]; ++k) {
                    const VertexIdx end2 = g_2->nbors[k];
                    if (i == end2) {continue;}  
                    if (g_2->getEdgeBinary(i, end2) == -1){
                        if (g->getEdgeBinary(i, end2) == -1){
                            local_ret.t5++;
                        }
                    }
                }

                for (EdgeIdx k = g->offsets[end1];  k < g->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g->nbors[k];
                    if (g_2->getEdgeBinary(i, end2) == -1){
                        if (g->getEdgeBinary(i, end2) == -1){
                            local_ret.t6++;
                        }
                    }
                }
                
            }
        
        }

        #pragma omp critical
        {
            ret.t1 += local_ret.t1;
            ret.t2 += local_ret.t2;
            ret.t3 += local_ret.t3;
            ret.t4 += local_ret.t4;
            ret.t5 += local_ret.t5;
            ret.t6 += local_ret.t6;
        }
    }
    return ret;
}

ThreeSizeInfo get3size_b2(CGraph *gout, CGraph *gout_2, CGraph *gin, CGraph *gin_2, int num_threads) {

    omp_set_num_threads(num_threads);
    ThreeSizeInfo ret{};
    VertexIdx current = 0;
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret{};
        
        #pragma omp for schedule(guided)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.t4++;
                    }
                    else{
                        local_ret.t3++;
                    }
                }

                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.t2++;
                    }
                    else{
                        EdgeIdx loc121 = end1 > end2 ? gout->getEdgeBinary(end2, end1) : gout->getEdgeBinary(end1, end2);
                        if (loc121 != -1) {
                            local_ret.t3++;
                        }
                        else{
                            local_ret.t6++;
                        }
                    }
                }

                for (EdgeIdx k = gout_2->offsets[end1];  k < gout_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_2->nbors[k];
                    if (gout_2->getEdgeBinary(i, end2) == -1){
                        if (gout->getEdgeBinary(i, end2) == -1){
                            local_ret.t6++;
                        }
                    }
                }

                for (EdgeIdx k = gin_2->offsets[end1];  k < gin_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_2->nbors[k];
                    if (i < end2){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.t6++;
                            }
                        }
                    }
                }        
            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];

                    EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.t1++;
                    } else {
                        EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                        if (loc_221 != -1) {
                            local_ret.t2++;
                        }
                        else{
                            local_ret.t5++;
                        }
                    }
                }

                for (EdgeIdx k = gout_2->offsets[end1];  k < gout_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_2->nbors[k];
                    if (gout_2->getEdgeBinary(i, end2) == -1){
                        if (gout->getEdgeBinary(i, end2) == -1){
                            local_ret.t5++;
                        }
                    }
                } 

                for (EdgeIdx k = gin_2->offsets[end1];  k < gin_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_2->nbors[k];
                    if (i < end2){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.t5++;
                            }
                        }
                    }
                } 

                for (EdgeIdx k = gout->offsets[end1];  k < gout->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout->nbors[k];
                    if (gout_2->getEdgeBinary(i, end2) == -1){
                        if (gout->getEdgeBinary(i, end2) == -1){
                            local_ret.t6++;
                        }
                    }
                }

                for (EdgeIdx k = gin->offsets[end1];  k < gin->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin->nbors[k];
                    if (i < end2){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.t6++;
                            }
                        }
                    }
                } 
            }      
        
        }

        #pragma omp critical
        {
            ret.t1 += local_ret.t1;
            ret.t2 += local_ret.t2;
            ret.t3 += local_ret.t3;
            ret.t4 += local_ret.t4;
            ret.t5 += local_ret.t5;
            ret.t6 += local_ret.t6;
        }
    }
    return ret;
}

void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&gcounts)[6], int num_threads){

    double t3 = 0, t5 = 0, t6 = 0;
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); 
        VertexIdx deg_2 = cg_2->degree(i); 

        t3 += deg * (deg-1) / 2; 
        t5 += deg_2 * (deg_2-1) / 2;
        t6 += deg * deg_2; 
    }

    gcounts[2] = t3;
    gcounts[4] = t5;
    gcounts[5] = t6;
    
    ThreeSizeInfo tricount = get3size(&(dag->outlist), &(dag_2->outlist), num_threads);
    gcounts[0] = tricount.t1;
    gcounts[1] = tricount.t2;
    gcounts[3] = tricount.t4;
}

void countThree_b1(CGraph *cg, CGraph *cg_2, double (&gcounts)[6], int num_threads){

    ThreeSizeInfo tcount = get3size_b1(cg, cg_2, num_threads);

    gcounts[0] = tcount.t1/3;
    gcounts[1] = tcount.t2/3;
    gcounts[2] = tcount.t3/3;
    gcounts[3] = tcount.t4/3;
    gcounts[4] = tcount.t5/3;
    gcounts[5] = tcount.t6/3;
}

void countThree_b2(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&gcounts)[6], int num_threads){

    ThreeSizeInfo tcount = get3size_b2(&(dag->outlist), &(dag_2->outlist), &(dag->inlist), &(dag_2->inlist), num_threads);

    gcounts[0] = tcount.t1;
    gcounts[1] = tcount.t2;
    gcounts[2] = tcount.t3;
    gcounts[3] = tcount.t4;
    gcounts[4] = tcount.t5;
    gcounts[5] = tcount.t6;
}


void comb_three(double (&gcounts)[6]){
    double t[6];

    t[0] = gcounts[0];
    t[1] = gcounts[1];
    t[2] = gcounts[2] - 3*gcounts[3];
    t[3] = gcounts[3];
    t[4] = gcounts[4] - 3 * t[0] - t[1];
    t[5] = gcounts[5] - 2 * t[1] - 2 * t[2];

    for (int m = 0; m < 6; ++m){
        gcounts[m] = t[m];
    }
}

void print3size(double (&gcounts)[6]){
    for (int i = 0; i < 6; ++i){
        printf("\"T%d\" : %.1f,\n", i+1, gcounts[i]);
    }
}
