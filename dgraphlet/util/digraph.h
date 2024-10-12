#pragma once

#include "graph.h"
#include <algorithm>

struct CDAG
{
    CGraph outlist;
    CGraph inlist;
};

// Structure for comparing nodes according to their degree.
// So u < v if degree of u less than that of v in graph g.

struct DegreeComp
{
    CGraph *g;
    DegreeComp(CGraph *g) { this->g = g;}

    bool operator () (VertexIdx u, VertexIdx v)
    {
        if (u<0 || v<0) // something wrong, so print error message, and exit
        {
            printf("Something wrong in DegreeComp, negative vertices. u is %d, v is %d\n",u,v);
            exit(EXIT_FAILURE);
        }
        VertexIdx degu = g->offsets[u+1] - g->offsets[u];  // Degree of u
        VertexIdx degv = g->offsets[v+1] - g->offsets[v];  // Degree of v
    
        if (degu < degv || (degu == degv && u < v))    // Comparing degrees and breaking ties by id
            return true;
        else
            return false;
    }
};

CDAG degreeOrdered(CGraph *g);
CDAG degreeOrdered2(CGraph *g, CGraph* g1);