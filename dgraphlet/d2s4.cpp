#include "./util/counter_d2s4.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstdio> 

using namespace std;
using namespace std::chrono;

//#define PRINT_CSR

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string filename = argv[1];
    int num_threads = atoi(argv[2]);

    Graph g;
    printf("Read Graph\n");
    read_mtx(filename, g);

    printf("Loaded graph\n");
    CGraph pre_cg = makeCSR(g);
    pre_cg.sortById();
    printf("Converted to CSR\n");

    // 2-edge Construction
    auto start_2edge = high_resolution_clock::now();
    printf("Construct 2-edge...\n");
    CGraph pre_cg_2 = pre_cg.getE2();
    auto end_2edge = high_resolution_clock::now();
    auto dur_2edge = duration_cast<milliseconds>(end_2edge - start_2edge);
    double sec_2edge = dur_2edge.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for 2-edge construction: %.3f\n", sec_2edge);

    // DAG construction
    printf("Creating DAG...\n");
    auto start_dag = high_resolution_clock::now();
    CGraph cg_2 = pre_cg_2.renameByDegreeOrder();
    cg_2.sortById();
    CGraph cg = pre_cg.reMapping(cg_2.mapping, cg_2.inverse);
    cg.sortById();
    CDAG dag_2 = degreeOrdered(&cg_2);
    (dag_2.outlist).sortById();
    (dag_2.inlist).sortById();
    CDAG dag = degreeOrdered2(&cg, &cg_2);
    (dag.outlist).sortById();
    (dag.inlist).sortById();
    auto end_dag = high_resolution_clock::now();
    auto dur_dag = duration_cast<milliseconds>(end_dag - start_dag);
    double sec_dag = dur_dag.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for DAG construction : %.3f\n", sec_dag);
    
    //Count (2,3)-graphlet
    double gcounts4[36];
    printf("Count (2,4)-graphlet...\n");
    auto start_count4 = high_resolution_clock::now();
    countFour(&dag, &dag_2, gcounts4, num_threads);
    auto end_count3 = high_resolution_clock::now();
    auto dur_count3 = duration_cast<milliseconds>(end_count3 - start_count4);
    double sec_count3 = dur_count3.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting (2,3)-graphlet : %.3f\n", sec_count3);
    comb_Four(gcounts4);
    print4size(gcounts4);

    printf("Total Execution time for 3-size: %.3f\n", sec_2edge + sec_dag + sec_count3);

    printf("# of Edge (1): %d", cg.nEdges / 2);
    printf("\n# of Edge (2): %d", cg_2.nEdges / 2);

    printf("\nMax Degree (1): %ld", cg.maxDegree);
    printf("\nMax Degree (2): %ld\n", cg_2.maxDegree);
}