# Geometric Hashing

## Prerequisites
* Currently only works on Linux distributions (tested under Ubuntu 18.04)
* You need Boost version 1.70.0 or higher. Either install via your distributions package manager or manually in `lib/boost_1_70_0`. You need the directories `lib/boost_1_70_0/build/lib` and `lib/boost_1_70_0/build/include`.
* You need a compiler that can handle the `-std=c++17` flag, preferrably GCC 8.3.0 or higher

## Installation
* Adjust the optimization in `Makefile` (we reccomend to use `OPTIMIZATION = 3`)
* make
```$ make```
* If desired, run the test cases
```$ make test && ./bin/test_seedFinding```

## Run
* Run `seedFinding` with the `--help` flag to see all available options
```$ ./bin/seedFinding --help```