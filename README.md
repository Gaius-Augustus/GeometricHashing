# Geometric Hashing

## Prerequisites
* Currently only works on Linux distributions (tested under Ubuntu 18.04)
* You need Boost version 1.70.0 or higher. Either install via your distributions package manager or manually in `lib/boost`. You need the directories `lib/boost/build/lib` and `lib/boost/build/include`.
* You need a compiler that can handle the `-std=c++17` flag, preferrably GCC 8.3.0 or higher

## Installation
* Adjust the optimization in `Makefile` (we reccomend to use `OPTIMIZATION = 3`)
* make

```$ make```

* If desired, run the test cases
```$ make test && ./bin/test_geometricHashing```

## Run
* Run `geometricHashing` with the `--help` flag to see all available options
```$ ./bin/geometricHashing --help```

## Data
* The file `dataWithArtificial.tar.bz` contains the fasta files we used in our experiments (contains also the random sequences) and a JSON file with the orthology relations
* Random sequences have a tag `artificial` in their fasta headers
