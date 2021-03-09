# seedFinding

Find accurate alignment seeds fast with geometric hashing.

### Cloning

Get this repo and the external libraries with

```
$ git clone --recurse-submodules https://github.com/Gaius-Augustus/GeometricHashing
$ git submodule update --init --recursive
```

If you already cloned but forgot the submodules, run

`$ git submodule update --init --recursive`

If you cloned with submodules but it was some time ago, you may want to run

`$ git submodule update --remote`

to  get the most recent versions of the external libraries

### Prerequisites

* GCC 9 or any compiler that supports `--std=c++17`and has `std::execution::par`

  * It might be necessary to also install the `tbb` library (Intel(R) Threading Building Blocks 2018 or later), i.e. `sudo apt install libtbb-dev`

  * If your distribution's repositories do not have the required packages, it might be easiest to install them with [Homebrew](https://docs.brew.sh/Homebrew-on-Linux)

    * If you chose to install GCC with Homebrew, you can put this into your `~/.bash_profile` (or whichever shell configuration file applies in your setup) to make sure the correct version is used

      ```bash
      export CC=$(which gcc-9)
      export CXX=$(which g++-9)
      ```

* CMake 3.6.3 or higher

* Boost 1.70.0 or higher

* Doxygen if you want to generate the documentation

### Build

Run the commands below to build.

```
$ mkdir build && cd build
$ cmake .. && make
```

Release mode is built by default. If you want a different build type, pass the following to CMake:

`$ cmake .. -DCMAKE_BUILD_TYPE=[DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]`

### Usage

Run `build/seedFinding` without parameters or with `--help` to get a list of command line parameters.
You should use `--allvsall` by default unless you get memory issues.

### Data

The file `dataWithArtificial.tar.bz` contains the fasta files we used in our experiments (contains also the random sequences) and a JSON file with the orthology relations  
Random sequences have a tag `gid:artificial` in their fasta headers

### Output

By default, the output file is called `seedFindingOutput.json`

It contains a JSON list of the following structure:

````
[
    [ [<pos>, <"0"|"1">, <genome name>, <sequence name>], [<pos>, <0|1>, <genome name>, <sequence name>] ]
    , ...
]
````

i.e. it is a list of lists of lists. 

The outer list holds all matches found by `seedFinding`. 

The inner lists are pairs of _occurrences_, one for `--genome1` and one for `--genome2`. If the input data contained more than two genomes, the additional genomes were used as hints to find better matches but are not reported here.

The occurrence lists can be seen as 4-tuples where `<pos>` is a positive integer, denoting the  _midpoint_ of the matching seed in both genomes. `<"0"|"1">` is either `0` or `1` (as a string type), denoting the strand on which the matching seed lies in the respective genome (`0` -- positive/forward strand). `<genome name>` is the name of the genome as used in `metagraph` (usually the filename with file extension, e.g. `hg38.fa`) or the filename of the input fasta file (in this case without extension, e.g. `hg38`). `<sequence name>` is the fasta header (without the leading `>`) of the respective sequence in the genome.

#### Geometric Hashing Output

If you ran `seedFinding --geometric-hashing --cube-output [1|2]`, another output file is created, by default called `seedFindingOutput.cubes.json`.

It contains a JSON object of the following structure:

```
{
	<cube identifier>: {
		<link identifier 1>: <link count>,
		<link identifier 2>: <link count>,
		<link identifier 3>: <link count>,
		...,
		"rawCount": <int>,
		"regionTuples": [<tuple>, ...],
		"score": <float>
	},
	...
}
```

`rawCount` is the number of links stored in the cube, `regionTupes` is a list of lists (pairs) of tile IDs, i.e. `[ [<start>, <end>], ...]`, where `<start>` and `<end>` are tile IDs (int) w.r.t. the reference genome (`--genome1`) sequence in that cube. `score` is the score that the respective cube reached in geometric hashing.

A tile ID here denotes an interval of positions on a sequence of length `--tilesize`: `[<tile ID> * tilesize, <tile ID + 1> * tilesize[`

`<cube identifier>` and `<link identifier>` are strings but can be parsed to JSON lists.

Use your favorite JSON parser on a `<cube identifier>` string and get:

```
[ [<genome name>, <sequence name>, <tile>, <"0"|"1">], ... ]
```

i.e. a list of lists (4-tuples) with elements similar to the normal output. The `<tile>` is an integer denoting the tile of the respective "side" of the cube (works slightly different than the tile ID above, see the reminder below on how it is computed).

Parsing `<link identifier>` yields

```
[ [<pos>, <"0"|"1">, <genome name>, <sequence name>], ... ]
```

i.e. again a list of lists (4-tuples)  which are seed occurrences just like above, only that now there can be more than two genomes per link (but each genome only once per link!)

---

Reminder: Cubes are defined by their links. Assume you parsed a `cube = json.parse(<cube identifier>)` (pseudo code) and the accompanying links `links = [json.parse(<link identifier 1>), json.parse(<link identifier 2>), json.parse(<link identifier 3>), ...])` (pseudo code) and stored the `--tilesize` used in the geometric hashing run in a variable `tilesize`.

Then, the following is true (Python):

```python
math.floor( (links[0][0] - links[0][0]) / tilesize ) == cube[0][2]
math.floor( (links[1][0] - links[0][0]) / tilesize ) == cube[1][2]
math.floor( (links[2][0] - links[0][0]) / tilesize ) == cube[2][2]
...
```



## Advanced Installation Notes

If you don't have all the dependencies and don't want to use Homebrew, here are some notes on how to install things manually.

You could install everything into your system paths with root access, however here we assume you have some separate private directory `${DIR}` into which everything will go.

### Installing Boost Manually

* Download the latest Boost version from https://www.boost.org/users/download/
* Unpack the archive and install using the following commands:
```
$ cd /path/to/boost_X_XX_X
$ mkdir build
$ ./bootstrap.sh --prefix=${DIR}
$ ./b2 stage threading=multi link=shared
$ ./b2 install threading=multi link=shared
```

Please refer to the official Boost documentation if errors occur during installation of Boost

You later need to tell CMake where to look for it:

`cmake .. -DBOOST_ROOT=/path/to/boost/build`

