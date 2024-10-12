#include "./util/counter_d3s3.h"
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

    // 3-edge Construction
    auto start_3edge = high_resolution_clock::now();
    printf("Construct 3-edge...\n");
    CGraph pre_cg_3 = pre_cg_2.getE3(pre_cg);
    auto end_3edge = high_resolution_clock::now();
    auto dur_3edge = duration_cast<milliseconds>(end_3edge - start_3edge);
    double sec_3edge = dur_3edge.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for 3-edge construction: %.3f\n", sec_3edge);
    
    printf("Creating DAG...\n");
    auto start_dag = high_resolution_clock::now();

    CGraph cg, cg_2, cg_3;
    CDAG dag, dag_2, dag_3;

    cg_3 = pre_cg_3.renameByDegreeOrder();
    cg_3.sortById();
    cg = pre_cg.reMapping(cg_3.mapping, cg_3.inverse);
    cg_2 = pre_cg_2.reMapping(cg_3.mapping, cg_3.inverse);
    cg.sortById();
    cg_2.sortById();

    dag_3 = degreeOrdered(&cg_3);
    (dag_3.outlist).sortById();
    (dag_3.inlist).sortById();

    dag_2 = degreeOrdered2(&cg_2, &cg_3);
    (dag_2.outlist).sortById();
    (dag_2.inlist).sortById();

    dag = degreeOrdered2(&cg, &cg_3);
    (dag.outlist).sortById();
    (dag.inlist).sortById();
    auto end_dag = high_resolution_clock::now();
    auto dur_dag = duration_cast<milliseconds>(end_dag - start_dag);
    double sec_dag = dur_dag.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for DAG construction : %.3f\n", sec_dag);
    
    
    //Count (3,3)-graphlet
    double gcounts3[13];
    printf("Count (3,3)-graphlet...\n");
    auto start_count3 = high_resolution_clock::now();
    countThree_b2(&cg, &dag, &cg_2, &dag_2, &cg_3, &dag_3, gcounts3, num_threads);
    auto end_count3 = high_resolution_clock::now();
    auto dur_count3 = duration_cast<milliseconds>(end_count3 - start_count3);
    double sec_count3 = dur_count3.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting (2,3)-graphlet : %.3f\n", sec_count3);
    print3size(gcounts3);

    printf("Total Execution time for 3-size: %.3f\n", sec_2edge + sec_3edge + sec_dag + sec_count3);

    printf("# of Edge (1): %d", cg.nEdges / 2);
    printf("\n# of Edge (2): %d", cg_2.nEdges / 2);
    printf("\n# of Edge (3): %d", pre_cg_3.nEdges / 2);

    printf("\nMax Degree (1): %ld", cg.maxDegree);
    printf("\nMax Degree (2): %ld", cg_2.maxDegree);
    printf("\nMax Degree (3): %ld\n", cg_3.maxDegree);
}