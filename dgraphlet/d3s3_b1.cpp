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
    CGraph cg = makeCSR(g);
    cg.sortById();
    printf("Converted to CSR\n");

    // 2-edge Construction
    auto start_2edge = high_resolution_clock::now();
    printf("Construct 2-edge...\n");
    CGraph cg_2 = cg.getE2();
    cg_2.sortById();
    auto end_2edge = high_resolution_clock::now();
    auto dur_2edge = duration_cast<milliseconds>(end_2edge - start_2edge);
    double sec_2edge = dur_2edge.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for 2-edge construction: %.3f\n", sec_2edge);

    // 3-edge Construction
    auto start_3edge = high_resolution_clock::now();
    printf("Construct 3-edge...\n");
    CGraph cg_3 = cg_2.getE3(cg);
    cg_3.sortById();
    auto end_3edge = high_resolution_clock::now();
    auto dur_3edge = duration_cast<milliseconds>(end_3edge - start_3edge);
    double sec_3edge = dur_3edge.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for 3-edge construction: %.3f\n", sec_3edge);

    //Count (3,3)-graphlet
    double gcounts3[13];
    
    printf("Count (3,3)-graphlet...\n");
    auto start_count3 = high_resolution_clock::now();
    countThree_b1(&cg, &cg_2, &cg_3, gcounts3, num_threads);
    auto end_count3 = high_resolution_clock::now();
    auto dur_count3 = duration_cast<milliseconds>(end_count3 - start_count3);
    double sec_count3 = dur_count3.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting (2,3)-graphlet : %.3f\n", sec_count3);
    print3size(gcounts3);

    printf("Total Execution time for 3-size: %.3f\n", sec_2edge + sec_3edge + sec_count3);

    printf("# of Edge (1): %d", cg.nEdges / 2);
    printf("\n# of Edge (2): %d", cg_2.nEdges / 2);
    printf("\n# of Edge (3): %d", cg_3.nEdges / 2);

    printf("\nMax Degree (1): %ld", cg.maxDegree);
    printf("\nMax Degree (2): %ld", cg_2.maxDegree);
    printf("\nMax Degree (3): %ld\n", cg_3.maxDegree);
}