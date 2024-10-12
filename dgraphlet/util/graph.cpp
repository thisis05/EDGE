#include "graph.h"
#include "JointSort.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>  
#include <unordered_map>
#include <vector>
#include <omp.h>


Graph Graph::copy() const
{
  Graph ret = newGraph(nVertices, nEdges);
  std::copy(srcs, srcs + nEdges, ret.srcs);
  std::copy(dsts, dsts + nEdges, ret.dsts);
  return ret;
}

Graph newGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new VertexIdx[nEdges], new VertexIdx[nEdges]};
}

CGraph newCGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new EdgeIdx[nVertices + 1], new VertexIdx[nEdges]};
}

void CGraph::sortById() const
{
    for (VertexIdx i=0; i < nVertices; i++){
        std::sort(nbors+offsets[i],nbors+offsets[i+1]);
    }
}

CGraph makeCSR(Graph& g){

    auto tmpG = g.copy();
    auto begin = JSIterator<VertexIdx, VertexIdx>{tmpG.srcs, tmpG.dsts};
    auto end   = begin + tmpG.nEdges;
    std::sort(begin, end);

    //We steal the dsts array from tmpG.
    CGraph ret = {g.nVertices, g.nEdges, new EdgeIdx[g.nVertices + 1], tmpG.dsts};
    //Now we have everything sorted by src, compress:
    VertexIdx cv = 0;
    for (EdgeIdx i = 0; i < tmpG.nEdges; ++i)
    {
        auto src = tmpG.srcs[i];
        while (cv <= src)
            ret.offsets[cv++] = i;
    }
    while (cv <= g.nVertices)
        ret.offsets[cv++] = g.nEdges;

    delete[] tmpG.srcs; //we retain tmpG.dsts in the output

    ret.sortById();

    VertexIdx maxDegree = 0;
    for (VertexIdx v = 0; v < ret.nVertices; ++v) {
        VertexIdx degree = ret.offsets[v + 1] - ret.offsets[v];
        if (degree > maxDegree) {
            maxDegree = degree;
        }
    }
    ret.maxDegree = maxDegree;

    return ret;
}

bool compareDegree(Pair firstPair, Pair nextPair) {
    if (firstPair.second != nextPair.second)
        return firstPair.second < nextPair.second;
    return firstPair.first < nextPair.first;
}

VertexIdx binarySearch(EdgeIdx* array, VertexIdx end, EdgeIdx val)
{
    VertexIdx low = 0;
    VertexIdx high = end-1;
    VertexIdx mid;

    while (low <= high)
    {
        mid = (low+high)/2;
        if (array[mid] == val)
            return mid;
        if (array[mid] > val)
            high = mid-1;
        if (array[mid] < val)
            low = mid+1;
    }
    return -1;
}


CGraph CGraph::renameByDegreeOrder() const
{
    CGraph ret = newCGraph(nVertices, nEdges);
    // std::unordered_map<long long, std::unordered_map<long long, long long>> adjMat;
    Pair *deg_info = new Pair[nVertices];

    VertexIdx *mapping = new VertexIdx[nVertices];
    VertexIdx *inverse = new VertexIdx[nVertices];


    // Construct array of pairs, storing old vertex label and degree
    for (VertexIdx i=0; i < nVertices; i++)
    {
        deg_info[i].first = i;
        deg_info[i].second = offsets[i+1]-offsets[i];
    }

    // sort the pairs by degree (if degree is same, sort by old vertex label)
    std::sort(deg_info,deg_info+nVertices, compareDegree);

    // Construct the mapping of old vertex label to new vertex label
    // So mapping[i] is what i is mapped to
    // And inverse[i] is what maps to i
    for (VertexIdx i=0; i < nVertices; i++)
    {
        mapping[deg_info[i].first] = i;
        inverse[i] = deg_info[i].first;
    }

    // Initialize offsets of output CGraph
    ret.offsets[0] = 0;
    EdgeIdx current = 0;

    Count cur_maxDegree = 0;
    Count new_maxDegree = 0;

    // Loop over new vertices
    for (VertexIdx new_label=0; new_label < nVertices; new_label++)
    {
        VertexIdx old_label = inverse[new_label]; // corresponding old label for new vertices
        for (EdgeIdx pos = offsets[old_label]; pos < offsets[old_label+1]; pos++) // loop over neighbors of old label
        {
            VertexIdx old_nbr = nbors[pos];
            VertexIdx new_nbr = mapping[old_nbr]; //corresponding new neighbor
            ret.nbors[current] = new_nbr; // insert new neighbor in nbors of output
            current++;
        }
        ret.offsets[new_label+1] = current; // all neighbors of new_label have been added, so we set offset for new_label+1
        new_maxDegree = ret.degree(new_label);
        if (cur_maxDegree < new_maxDegree){
            cur_maxDegree = new_maxDegree;
        }
    }
    //ret.adj_matrix = adjMat;
    ret.maxDegree = cur_maxDegree;
    ret.mapping = mapping;
    ret.inverse = inverse;
    return ret;
}

CGraph CGraph::reMapping(VertexIdx *mapping, VertexIdx *inverse) const{

    CGraph ret = newCGraph(nVertices, nEdges);

    ret.offsets[0] = 0;
    EdgeIdx current = 0;

    Count cur_maxDegree = 0;
    Count new_maxDegree = 0;

    for (VertexIdx new_label=0; new_label < nVertices; new_label++)
    {
        VertexIdx old_label = inverse[new_label]; // corresponding old label for new vertices
        for (EdgeIdx pos = offsets[old_label]; pos < offsets[old_label+1]; pos++) // loop over neighbors of old label
        {
            VertexIdx old_nbr = nbors[pos];
            VertexIdx new_nbr = mapping[old_nbr]; //corresponding new neighbor
            ret.nbors[current] = new_nbr; // insert new neighbor in nbors of output

            // adjMat[new_label][new_nbr] = 1;
            // adjMat[new_nbr][new_label] = 1;
            current++;
        }
        ret.offsets[new_label+1] = current; // all neighbors of new_label have been added, so we set offset for new_label+1
        new_maxDegree = ret.degree(new_label);
        if (cur_maxDegree < new_maxDegree){
            cur_maxDegree = new_maxDegree;
        }
    }
    //ret.adj_matrix = adjMat;
    ret.maxDegree = cur_maxDegree;

    return ret;
}

CGraph CGraph::getE2() const {
    
    CGraph ret = newCGraph(nVertices, (nEdges * 100));
    EdgeIdx current = 0;
    Count cur_maxDegree = 0;
    Count new_maxDegree = 0;
    
    ret.offsets[0] = 0;
    for (VertexIdx start = 0; start < nVertices; start++) {
        std::vector<bool> check(nVertices, false);

        for (EdgeIdx idx = offsets[start]; idx < offsets[start + 1]; idx++) {
            VertexIdx nbr = nbors[idx];
            check[nbr] = true;
        }

        for (EdgeIdx idx = offsets[start]; idx < offsets[start + 1]; idx++) {
            VertexIdx nbr = nbors[idx];
            for (EdgeIdx idx2 = offsets[nbr]; idx2 < offsets[nbr + 1]; idx2++) {
                VertexIdx final = nbors[idx2];
                if (start != final && !check[final]) {
                    check[final] = true;
                    ret.nbors[current++] = final;
                }
            }
        }
        ret.offsets[start + 1] = current;
        new_maxDegree = ret.degree(start);
        if (cur_maxDegree < new_maxDegree){
            cur_maxDegree = new_maxDegree;
        }

    }

    ret.nEdges = current;
    ret.nbors = (VertexIdx*)realloc(ret.nbors, current * sizeof(VertexIdx));
    ret.maxDegree = cur_maxDegree;
    
    return ret;
}

CGraph CGraph::getE3(CGraph G1) const {
    CGraph ret = newCGraph(nVertices, (nEdges * 50));
    EdgeIdx current = 0;
    Count cur_maxDegree = 0;
    Count new_maxDegree = 0;

    ret.offsets[0] = 0;
    for (VertexIdx start = 0; start < nVertices; start++) {
        std::vector<bool> check(nVertices, false);

        for (EdgeIdx idx = G1.offsets[start]; idx < G1.offsets[start + 1]; idx++) {
            VertexIdx nbr = G1.nbors[idx];
            check[nbr] = true;
        }

        for (EdgeIdx idx = offsets[start]; idx < offsets[start + 1]; idx++) {
            VertexIdx nbr = nbors[idx];
            check[nbr] = true;
        }

        for (EdgeIdx idx = offsets[start]; idx < offsets[start + 1]; idx++) {
            VertexIdx nbr = nbors[idx];
            for (EdgeIdx idx2 = G1.offsets[nbr]; idx2 < G1.offsets[nbr + 1]; idx2++) {
                VertexIdx final = G1.nbors[idx2];
                if (start != final && !check[final]) {
                    check[final] = true;
                    ret.nbors[current++] = final;
                }
            }
        }

        ret.offsets[start + 1] = current;
        new_maxDegree = ret.degree(start);
        if (cur_maxDegree < new_maxDegree){
            cur_maxDegree = new_maxDegree;
        }
    }

    ret.nEdges = current;
    ret.nbors = (VertexIdx*)realloc(ret.nbors, current * sizeof(VertexIdx));
    ret.maxDegree = cur_maxDegree;

    return ret;
}


EdgeIdx CGraph::getEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;
        if (nbors[mid] == v2)
            return mid;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return -1;
}

bool CGraph::isEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;

        if (nbors[mid] == v2)
            return true;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return false;
}


void read_mtx(const std::string& filename, Graph& graph) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }

    std::istringstream iss(line);
    VertexIdx n;
    EdgeIdx m;
    iss >> n >> n >> m;
    graph.nVertices = n;
    graph.nEdges = 2*m;

    VertexIdx u, v;
    EdgeIdx iEdge = 0;

    graph.srcs = new VertexIdx[graph.nEdges];
    graph.dsts = new VertexIdx[graph.nEdges];

    while (file >> u >> v) {
        u--, v--;
        graph.srcs[iEdge] = v;
        graph.dsts[iEdge] = u;
        ++iEdge;
        

        graph.srcs[iEdge] = u;
        graph.dsts[iEdge] = v;
        ++iEdge;

    }
    file.close();
}

