#include "counter_d3s3.h"
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


ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2, CGraph* gout_3, int num_threads) {

    omp_set_num_threads(num_threads);
    ThreeSizeInfo ret = {};
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret = {};
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            const EdgeIdx start_3 = gout_3->offsets[i];
            const EdgeIdx end_3 = gout_3->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                // T9 : (di, dj, dk) = (1, 1, 1)
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.tri9++;
                    }
                }

                // T7 : (1) (di, dj, dk) = (1, 2, 2) + (2) (di, dj, dk) = (2, 1, 2)
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.tri7++;
                    }
                }       

                // T4: (1) (di, dj, dk) = (1, 3, 3) + (2) (di, dj, dk) = (3, 1, 3)
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc133 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc133 != -1) {
                        local_ret.tri4++;
                    }
                }    

            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    // T3 : (di, dj, dk) = (2, 2, 3)
                    EdgeIdx loc_223 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_223 != -1){
                        local_ret.tri3++;
                    }
                    else{
                        // T6 : (di, dj, dk) = (2, 2, 2)
                        EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_222 != -1) {
                            local_ret.tri6++;
                        } else {
                            // T7 : (di, dj, dk) = (2, 2, 1)
                            EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                            if (loc_221 != -1) {
                                local_ret.tri7++;
                            }
                        }
                    }
                }
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    // T2: (1) (di, dj, dk) = (2, 3, 3) + (2) (di, dj, dk) = (3, 2, 3)
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc233 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc233 != -1) {
                        local_ret.tri2++;
                    }
                    else{
                        // T3: (1) (di, dj, dk) = (2, 3, 2) + (2) (di, dj, dk) = (3, 2, 2)
                        EdgeIdx loc232 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                        if (loc232 != -1) {
                            local_ret.tri3++;
                        }
                    }
                }
            }

            for (EdgeIdx j = start_3; j < end_3; ++j) {
                const VertexIdx end1 = gout_3->nbors[j];
                for (EdgeIdx k = j+1; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    // T1 : (di, dj, dk) = (3, 3, 3)
                    EdgeIdx loc_333 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_333 != -1){
                        local_ret.tri1++;
                    }
                    else{
                        // T2 : (di, dj, dk) = (3, 3, 2)
                        EdgeIdx loc_332 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_332 != -1) {
                            local_ret.tri2++;
                        } else {
                            // T4 : (di, dj, dk) = (3, 3, 1)
                            EdgeIdx loc_331 = gout->getEdgeBinary(end1, end2);
                            if (loc_331 != -1) {
                                local_ret.tri4++;
                            }
                        }
                    }
                }
            }     
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri3 += local_ret.tri3;
            ret.tri4 += local_ret.tri4;
            ret.tri6 += local_ret.tri6;
            ret.tri7 += local_ret.tri7;
            ret.tri9 += local_ret.tri9;
        }
    }

    return ret;
}

ThreeSizeInfo get3size_b1(CGraph *g, CGraph *g_2, CGraph* g_3, int num_threads) {

    omp_set_num_threads(num_threads);
    ThreeSizeInfo ret = {};
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret = {};
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < g->nVertices; ++i) {
            const EdgeIdx start = g->offsets[i];
            const EdgeIdx end = g->offsets[i+1];
            const EdgeIdx start_2 = g_2->offsets[i];
            const EdgeIdx end_2 = g_2->offsets[i+1];
            const EdgeIdx start_3 = g_3->offsets[i];
            const EdgeIdx end_3 = g_3->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = g->nbors[j];
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = g->nbors[k];
                    EdgeIdx loc_112 = g_2->getEdgeBinary(end1, end2);
                    if(loc_112 != -1){
                        local_ret.tri8++;
                    }
                    else{
                        EdgeIdx loc_111 = g->getEdgeBinary(end1, end2);
                        if (loc_111 != -1) {
                            local_ret.tri9++;
                        }
                    }    
                }

                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = g_2->nbors[k];
                    EdgeIdx loc123 = end1 > end2 ? g_3->getEdgeBinary(end2, end1) : g_3->getEdgeBinary(end1, end2);
                    if (loc123 != -1){
                        local_ret.tri5++;
                    }
                    else{
                        EdgeIdx loc122 = end1 > end2 ? g_2->getEdgeBinary(end2, end1) : g_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.tri7++;
                        }
                        else{
                            EdgeIdx loc121 = end1 > end2 ? g->getEdgeBinary(end2, end1) : g->getEdgeBinary(end1, end2);
                            if (loc121 != -1){
                                local_ret.tri8++;
                            }
                        }
                    }
                    
                }       

                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = g_3->nbors[k];
                    EdgeIdx loc133 = end1 > end2 ? g_3->getEdgeBinary(end2, end1) : g_3->getEdgeBinary(end1, end2);
                    if (loc133 != -1) {
                        local_ret.tri4++;
                    }
                    else{
                        EdgeIdx loc132 = end1 > end2 ? g_2->getEdgeBinary(end2, end1) : g_2->getEdgeBinary(end1, end2);
                        if (loc132 != -1){
                            local_ret.tri5++;
                        }
                        else{
                            local_ret.tri13++;
                        }
                    }
                } 

                for (EdgeIdx k = g_3->offsets[end1];  k < g_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_3->nbors[k];
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri13++;
                            }
                        }
                    } 
                }   

            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = g_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = g_2->nbors[k];
                    EdgeIdx loc_223 = g_3->getEdgeBinary(end1, end2);
                    if (loc_223 != -1){
                        local_ret.tri3++;
                    }
                    else{
                        EdgeIdx loc_222 = g_2->getEdgeBinary(end1, end2);
                        if (loc_222 != -1) {
                            local_ret.tri6++;
                        } else {
                            EdgeIdx loc_221 = g->getEdgeBinary(end1, end2);
                            if (loc_221 != -1) {
                                local_ret.tri7++;
                            }
                            else{
                                local_ret.tri12++;
                            }
                        }
                    }
                }
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = g_3->nbors[k];
                    EdgeIdx loc233 = end1 > end2 ? g_3->getEdgeBinary(end2, end1) : g_3->getEdgeBinary(end1, end2);
                    if (loc233 != -1) {
                        local_ret.tri2++;
                    }
                    else{
                        EdgeIdx loc232 = end1 > end2 ? g_2->getEdgeBinary(end2, end1) : g_2->getEdgeBinary(end1, end2);
                        if (loc232 != -1) {
                            local_ret.tri3++;
                        }
                        else{
                            EdgeIdx loc231 = end1 > end2 ? g->getEdgeBinary(end2, end1) : g->getEdgeBinary(end1, end2);
                            if (loc231 != -1){
                                local_ret.tri5++;
                            }
                            else{
                                local_ret.tri11++;
                            }
                        }
                    }
                }


                for (EdgeIdx k = g_2->offsets[end1];  k < g_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_2->nbors[k];
                    if (i == end2) {continue;}
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri12++;
                            }
                        }
                    }  
                } 

                for (EdgeIdx k = g_3->offsets[end1];  k < g_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_3->nbors[k];
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri11++;
                            }
                        }
                    } 
                }
                
            }

            for (EdgeIdx j = start_3; j < end_3; ++j) {
                const VertexIdx end1 = g_3->nbors[j];
                for (EdgeIdx k = j+1; k < end_3; ++k) {
                    const VertexIdx end2 = g_3->nbors[k];
                    EdgeIdx loc_333 = g_3->getEdgeBinary(end1, end2);
                    if (loc_333 != -1){
                        local_ret.tri1++;
                    }
                    else{
                        EdgeIdx loc_332 = g_2->getEdgeBinary(end1, end2);
                        if (loc_332 != -1) {
                            local_ret.tri2++;
                        } else {
                            EdgeIdx loc_331 = g->getEdgeBinary(end1, end2);
                            if (loc_331 != -1) {
                                local_ret.tri4++;
                            }
                            else{
                                local_ret.tri10++;
                            }
                        }
                    }
                }

                for (EdgeIdx k = g->offsets[end1];  k < g->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g->nbors[k];
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri13++;
                            }
                        }
                    }  
                } 


                for (EdgeIdx k = g_2->offsets[end1];  k < g_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_2->nbors[k];
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri11++;
                            }
                        }
                    }  
                } 

                for (EdgeIdx k = g_3->offsets[end1];  k < g_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = g_3->nbors[k];
                    if (i == end2) {continue;}  
                    if (g_3->getEdgeBinary(i, end2) == -1){
                        if (g_2->getEdgeBinary(i, end2) == -1){
                            if (g->getEdgeBinary(i, end2) == -1){
                                local_ret.tri10++;
                            }
                        }
                    } 
                }
            }     
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri3 += local_ret.tri3;
            ret.tri4 += local_ret.tri4;
            ret.tri5 += local_ret.tri5;
            ret.tri6 += local_ret.tri6;
            ret.tri7 += local_ret.tri7;
            ret.tri8 += local_ret.tri8;
            ret.tri9 += local_ret.tri9;
            ret.tri10 += local_ret.tri10;
            ret.tri11 += local_ret.tri11;
            ret.tri12 += local_ret.tri12;
            ret.tri13 += local_ret.tri13;
        }
    }

    return ret;
}


ThreeSizeInfo get3size_b2(CGraph *gout, CGraph *gout_2, CGraph *gout_3, CGraph *gin, CGraph *gin_2, CGraph *gin_3, int num_threads) {

    omp_set_num_threads(num_threads);
    ThreeSizeInfo ret = {};
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret = {};
        
        #pragma omp for schedule(dynamic, 64)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            const EdgeIdx start_3 = gout_3->offsets[i];
            const EdgeIdx end_3 = gout_3->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    EdgeIdx loc_112 = gout_2->getEdgeBinary(end1, end2);
                    if(loc_112 != -1){
                        local_ret.tri8++;
                    }
                    else{
                        EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                        if (loc_111 != -1) {
                            local_ret.tri9++;
                        }
                    }    
                }

                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc123 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc123 != -1){
                        local_ret.tri5++;
                    }
                    else{
                        EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.tri7++;
                        }
                        else{
                            EdgeIdx loc121 = end1 > end2 ? gout->getEdgeBinary(end2, end1) : gout->getEdgeBinary(end1, end2);
                            if (loc121 != -1){
                                local_ret.tri8++;
                            }
                        }
                    }
                    
                }       

                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc133 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc133 != -1) {
                        local_ret.tri4++;
                    }
                    else{
                        EdgeIdx loc132 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                        if (loc132 != -1){
                            local_ret.tri5++;
                        }
                        else{
                            local_ret.tri13++;
                        }
                    }
                } 

                for (EdgeIdx k = gout_3->offsets[end1];  k < gout_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_3->nbors[k];
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri13++;
                            }
                        }
                    } 
                }   

                for (EdgeIdx k = gin_3->offsets[end1];  k < gin_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_3->nbors[k];
                    if (i < end2){
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri13++;
                                }
                            }
                        } 
                    }
                }

            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc_223 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_223 != -1){
                        local_ret.tri3++;
                    }
                    else{
                        EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_222 != -1) {
                            local_ret.tri6++;
                        } else {
                            EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                            if (loc_221 != -1) {
                                local_ret.tri7++;
                            }
                            else{
                                local_ret.tri12++;
                            }
                        }
                    }
                }
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc233 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc233 != -1) {
                        local_ret.tri2++;
                    }
                    else{
                        EdgeIdx loc232 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                        if (loc232 != -1) {
                            local_ret.tri3++;
                        }
                        else{
                            EdgeIdx loc231 = end1 > end2 ? gout->getEdgeBinary(end2, end1) : gout->getEdgeBinary(end1, end2);
                            if (loc231 != -1){
                                local_ret.tri5++;
                            }
                            else{
                                local_ret.tri11++;
                            }
                        }
                    }
                }


                for (EdgeIdx k = gout_2->offsets[end1];  k < gout_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_2->nbors[k];
                    if (i == end2) {continue;}
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri12++;
                            }
                        }
                    }  
                } 

                for (EdgeIdx k = gin_2->offsets[end1];  k < gin_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_2->nbors[k];
                    if (i < end2) {
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri12++;
                                }
                            }
                        }
                    }  
                } 

                for (EdgeIdx k = gout_3->offsets[end1];  k < gout_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_3->nbors[k];
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri11++;
                            }
                        }
                    } 
                }

                for (EdgeIdx k = gin_3->offsets[end1];  k < gin_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_3->nbors[k];
                    if (i < end2){
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri11++;
                                }
                            }
                        } 
                    }

                }
                
            }

            for (EdgeIdx j = start_3; j < end_3; ++j) {
                const VertexIdx end1 = gout_3->nbors[j];
                for (EdgeIdx k = j+1; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc_333 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_333 != -1){
                        local_ret.tri1++;
                    }
                    else{
                        EdgeIdx loc_332 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_332 != -1) {
                            local_ret.tri2++;
                        } else {
                            EdgeIdx loc_331 = gout->getEdgeBinary(end1, end2);
                            if (loc_331 != -1) {
                                local_ret.tri4++;
                            }
                            else{
                                local_ret.tri10++;
                            }
                        }
                    }
                }

                for (EdgeIdx k = gout->offsets[end1];  k < gout->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout->nbors[k];
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri13++;
                            }
                        }
                    }  
                } 

                for (EdgeIdx k = gin->offsets[end1];  k < gin->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin->nbors[k];
                    if (i < end2){
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri13++;
                                }
                            }
                        } 
                    }
                     
                } 

                for (EdgeIdx k = gout_2->offsets[end1];  k < gout_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_2->nbors[k];
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri11++;
                            }
                        }
                    }
                } 

                for (EdgeIdx k = gin_2->offsets[end1];  k < gin_2->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_2->nbors[k];
                    if (i < end2){
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri11++;
                                }
                            }
                        }
                    }
                      
                } 

                for (EdgeIdx k = gout_3->offsets[end1];  k < gout_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gout_3->nbors[k];
                    if (gout_3->getEdgeBinary(i, end2) == -1){
                        if (gout_2->getEdgeBinary(i, end2) == -1){
                            if (gout->getEdgeBinary(i, end2) == -1){
                                local_ret.tri10++;
                            }
                        }
                    }     
                }

                for (EdgeIdx k = gin_3->offsets[end1];  k < gin_3->offsets[end1+1]; ++k) {  
                    const VertexIdx end2 = gin_3->nbors[k];
                    if (i < end2){
                        if (gout_3->getEdgeBinary(i, end2) == -1){
                            if (gout_2->getEdgeBinary(i, end2) == -1){
                                if (gout->getEdgeBinary(i, end2) == -1){
                                    local_ret.tri10++;
                                }
                            }
                        } 
                    }     
                }
            }     
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri3 += local_ret.tri3;
            ret.tri4 += local_ret.tri4;
            ret.tri5 += local_ret.tri5;
            ret.tri6 += local_ret.tri6;
            ret.tri7 += local_ret.tri7;
            ret.tri8 += local_ret.tri8;
            ret.tri9 += local_ret.tri9;
            ret.tri10 += local_ret.tri10;
            ret.tri11 += local_ret.tri11;
            ret.tri12 += local_ret.tri12;
            ret.tri13 += local_ret.tri13;
        }
    }

    return ret;
}


void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, CGraph *cg_3, CDAG *dag_3, double (&gcounts)[13], int num_threads){

    double t5 = 0, t8 = 0, t10 = 0, t11 = 0, t12 = 0, t13 = 0;
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); 
        VertexIdx deg_2 = cg_2->degree(i); 
        VertexIdx deg_3 = cg_3->degree(i);

        t5 += deg * deg_2;
        t8 += deg * (deg-1) / 2; 
        t10 += deg_3 * (deg_3 - 1) / 2; 
        t11 += deg_2 * deg_3;
        t12 += deg_2 * (deg_2-1) / 2;
        t13 += deg * deg_3;
    }

    gcounts[4] = t5;
    gcounts[7] = t8;
    gcounts[9] = t10;
    gcounts[10] = t11;
    gcounts[11] = t12;
    gcounts[12] = t13;
    ThreeSizeInfo tricount = get3size(&(dag->outlist), &(dag_2->outlist), &(dag_3->outlist), num_threads);
    gcounts[0] = tricount.tri1;
    gcounts[1] = tricount.tri2;
    gcounts[2] = tricount.tri3;
    gcounts[3] = tricount.tri4;
    gcounts[5] = tricount.tri6;
    gcounts[6] = tricount.tri7;
    gcounts[8] = tricount.tri9;
}

void countThree_b1(CGraph *cg, CGraph *cg_2, CGraph *cg_3, double (&gcounts)[13], int num_threads){

    ThreeSizeInfo tricount = get3size_b1(cg, cg_2, cg_3, num_threads);

    gcounts[0] = tricount.tri1/3;
    gcounts[1] = tricount.tri2/3;
    gcounts[2] = tricount.tri3/3;
    gcounts[3] = tricount.tri4/3;
    gcounts[4] = tricount.tri5/3;
    gcounts[5] = tricount.tri6/3;
    gcounts[6] = tricount.tri7/3;
    gcounts[7] = tricount.tri8/3;
    gcounts[8] = tricount.tri9/3;
    gcounts[9] = tricount.tri10/3;
    gcounts[10] = tricount.tri11/3;
    gcounts[11] = tricount.tri12/3;
    gcounts[12] = tricount.tri13/3;
}

void countThree_b2(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, CGraph *cg_3, CDAG *dag_3, double (&gcounts)[13], int num_threads){

    ThreeSizeInfo tricount = get3size_b2(&(dag->outlist), &(dag_2->outlist), &(dag_3->outlist), &(dag->inlist), &(dag_2->inlist), &(dag_3->inlist), num_threads);

    gcounts[0] = tricount.tri1;
    gcounts[1] = tricount.tri2;
    gcounts[2] = tricount.tri3;
    gcounts[3] = tricount.tri4;
    gcounts[4] = tricount.tri5;
    gcounts[5] = tricount.tri6;
    gcounts[6] = tricount.tri7;
    gcounts[7] = tricount.tri8;
    gcounts[8] = tricount.tri9;
    gcounts[9] = tricount.tri10;
    gcounts[10] = tricount.tri11;
    gcounts[11] = tricount.tri12;
    gcounts[12] = tricount.tri13;
}

void comb_three(double (&gcounts)[13]){
    double t[13];

    t[0] = gcounts[0];
    t[1] = gcounts[1];
    t[2] = gcounts[2];
    t[3] = gcounts[3];

    t[8] = gcounts[8];
    t[7] = gcounts[7] - 3 * t[8];
    t[6] = gcounts[6];

    t[4] = gcounts[4] - 2 * t[6] - 2 * t[7];
    t[5] = gcounts[5];
    
    t[9] = gcounts[9] - t[3] - t[1] - 3 * t[0];
    t[10] = gcounts[10] - 2 * t[1] - 2 * t[2] - t[4];
    t[11] = gcounts[11] - t[2] - 3 * t[5] - t[6];
    t[12] = gcounts[12] - 2 * t[3] - t[4];


    for (int m = 0; m < 13; ++m){
        gcounts[m] = t[m];
    }
}

void print3size(double (&gcounts)[13]){
    for (int i = 0; i < 13; ++i){
        printf("\"T%d\" : %.1f,\n", i+1, gcounts[i]);
    }
}
