# Geometric Hashing

### Prerequisites

* GCC 9 or any compiler that supports `--std=c++17`and has `std::execution::par`

  * It might be necessary to also install the `tbb` library (Intel(R) Threading Building Blocks 2018 or later), i.e. `sudo apt install libtbb-dev`

  * If your distribution's repositories do not have the required packages, it might be easiest to install them with [Homebrew](https://docs.brew.sh/Homebrew-on-Linux)

    * If you chose to install GCC with Homebrew, you can put this into your `~/.bash_profile` (or whichever shell configuration file applies in your setup) to make `make` use the correct version

      ```bash
      export CC=$(which gcc-9)
      export CXX=$(which g++-9)
      ```

* Boost 1.70.0 or higher

* Doxygen if you want to generate the documentation

#### Installing Boost manually
In case the Boost packages from your official Linux distribution are too old and you don't want to use Homebrew, you need to manually install Boost in a certain place

* Download the latest Boost version from https://www.boost.org/users/download/
* Unpack the archive to `lib/boost`
* Install using the following commands:
```
$ cd lib/boost/
$ mkdir build
$ ./bootstrap.sh --prefix=./build
$ ./b2 stage threading=multi link=shared
$ ./b2 install threading=multi link=shared
```

Please refer to the official Boost documentation if errors occur during installation of Boost

*Remark* --  You can unpack the Boost archive to wherever you want, just make sure that the installed files end up in the correct place: `$ ./bootstrap.sh --prefix=(/path/to/)GeometricHashing/lib/boost/build`

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
