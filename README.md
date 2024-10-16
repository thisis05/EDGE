# Official Code Implementation of EDGE

## Paper Information
- **Title**: Beyond Neighbors: Distance-Generalized Graphlets for Enhanced Graph Characterization
- **Code for**: Algorithm **EDGE** for counting distance-generalized graphlets, i.e., $(d,s)$-graphlets.

## Datasets
We provide 13 real-world graph datasets:
- **Collaboration**: ca-dblp-2012, ca-citeseer, ca-HepTh
- **Web**: web-arabic-2005, web-indochina
- **Social-Facebook**: socfb-UCSC68, socfb-UC64, socfb-Middlebury45
- **Tags**: tags-ask-ubuntu, tags-math-sx
- **Road Networks**: roadNet-PA, roadNet-CA, roadNet-TX

Preprocessed versions are included. If you need the original data, refer to the original sources:
- [Network repository](https://networkrepository.com/)
- [SNAP datsets](https://snap.stanford.edu/data/index.html)
- [ARB](https://www.cs.cornell.edu/~arb/data/)

## Explanation

### Main Codes (dgraphlet)
These files are used to execute the EDGE algorithm:
- `d2s3.cpp`: Counts $(2,3)$-graphlets (distance up to 2, size 3 graphlets).
- `d2s3_b1.cpp`: Baseline version 1 for counting $(2,3)$-graphlets (**EDGE-D2** in the paper).
- `d2s3_b2.cpp`: Baseline version 2 for counting $(2,3)$-graphlets (**EDGE-D** in the paper).

- `d3s3.cpp`: Counts $(3,3)$-graphlets (distance up to 3, size 3 graphlets).
- `d3s3_b1.cpp`: Baseline version 1 for counting $(3,3)$-graphlets (**EDGE-D2** in the paper).
- `d3s3_b2.cpp`: Baseline version 2 for counting $(3,3)$-graphlets (**EDGE-D** in the paper).

- `d2s4.cpp`: Counts $(2,4)$-graphlets (distance up to 2, size 4 graphlets).

### Utility Codes
These are supporting codes for the main execution (constructing $d$-edge, $d$-graph, $d$-DAG...).
Some construction code is based on the existing code available [here](https://bitbucket.org/seshadhri/escape/src/master/) of the paper: [ESCAPE: Efficiently Counting All 5-Vertex Subgraphs](https://dl.acm.org/doi/abs/10.1145/3038912.3052597?casa_token=S-hx23geQAAAAAAA:DaLCMp0j28f_TfQkgZbNlBMN6rJQEsctf8LYN5-Iv0QBou0BKOYL-aFUfiq3SHDsPTSBwFqAIBzDaQ).

## Run Instructions

### Input File Format
Currently, we only support the Matrix Market Coordinate Format (symmetric). We will update soon to support all datatypes.
For more details, see: <http://math.nist.gov/MatrixMarket/formats.html#MMformat>.

### Main Commands

1. Navigate to the dgraphlet folder:
   ```bash
   cd dgraphlet
2. Choose the $(d,s)$-graphlet to count (e.g., for $(2,3)$-graphlets, use `d2s3.cpp`).
3. Use the Makefile to compile:
   ```bash
   make
4. Run the executable:
   ```bash
   ./[filename].cpp ../data/[dataset].mtx -n
Here, `n` represents the number of threads you wish to use for processing.

Example command for running the $(2,3)$-graphlet on the ca-HepTh dataset using 6 threads:
   ```bash
   ./d2s3.cpp ../data/ca-HepTh.mtx -6


