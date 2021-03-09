import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description = "Create all testdata",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--alpha",
                        dest = "alpha", 
                        metavar = "FLOAT", 
                        type=float,
                        default=1,
                        help="alpha parameter for region tuple extraction")
scriptArgs.add_argument("--beta",
                        dest = "beta", 
                        metavar = "FLOAT", 
                        type=float,
                        default=1,
                        help="beta parameter for region tuple extraction")
scriptArgs.add_argument("--binsize",
                        dest = "binsize", 
                        metavar = "INT", 
                        type=int,
                        default=100,
                        help="binsize of metagraph")
scriptArgs.add_argument("--chunksize",
                        dest = "chunksize", 
                        metavar = "INT", 
                        type=int,
                        default=200,
                        help="chunksize/length in geometric hashing")
scriptArgs.add_argument("--gamma",
                        dest = "gamma", 
                        metavar = "FLOAT", 
                        type=float,
                        default=1,
                        help="gamma parameter for region tuple extraction")
scriptArgs.add_argument("--metagraph",
                        dest = "metagraph", 
                        metavar = "FILE", 
                        type=argparse.FileType("r"),
                        help="metagraph binary", 
                        required = True)
scriptArgs.add_argument("--nsequences",
                        dest = "nsequences", 
                        metavar = "INT", 
                        type=int,
                        default=10,
                        help="number of sequences to generate per species")
scriptArgs.add_argument("--nspecies",
                        dest = "nspecies", 
                        metavar = "INT", 
                        type=int,
                        default=10,
                        help="number of species to generate")
scriptArgs.add_argument("--nsubtiles",
                        dest = "nsubtiles", 
                        metavar = "INT", 
                        type=int,
                        default=25,
                        help="number of subcubes in GH (approximate)")
scriptArgs.add_argument("--pnorm",
                        dest = "pnorm", 
                        metavar = "INT", 
                        type=int,
                        default=3,
                        help="p for p-norm (GH)")
scriptArgs.add_argument("--redmask",
                        dest = "redmask", 
                        metavar = "BOOL", 
                        type=bool,
                        default=True,
                        help="Discard low-complexity seeds")
scriptArgs.add_argument("--sequence-lengths",
                        dest = "seqlen", 
                        metavar = "INT", 
                        type=int,
                        default=10000,
                        help="length of each sequence")
scriptArgs.add_argument("--similarity",
                        dest = "similarity", 
                        metavar = "FLOAT", 
                        type=float,
                        default=0.85,
                        help="sequence similarity between neighbouring species")
scriptArgs.add_argument("--stringhasher",
                        dest = "stringhasher", 
                        metavar = "FILE", 
                        type=argparse.FileType("r"),
                        help="stringHasher binary", 
                        required = True)
scriptArgs.add_argument("--tilesize",
                        dest = "tilesize", 
                        metavar = "INT", 
                        type=int,
                        default=500,
                        help="tilesize in geometric hashing")
scriptArgs.add_argument("--thinning",
                        dest = "thinning", 
                        metavar = "INT", 
                        type=int,
                        default=5,
                        help="Discard roughly a fraction of 1/thinning k-mers")
args = parser.parse_args()
metagraphBin = args.metagraph.name
stringHasherBin = args.stringhasher.name

assert args.binsize >= 1
assert args.nsequences >= 1
assert args.nspecies >= 2
assert args.redmask in [True, False]
assert args.seqlen >= 1
assert args.similarity > 0
assert args.similarity < 1
assert args.tilesize >= 1



# k fixed to 39 to match masks in expectedSeeds.py
k = 39



def run(command):
    print("Running", command)
    returncode = subprocess.run([command], shell = True, executable = '/bin/bash')
    return(returncode.returncode == 0)



# run the scripts in order
if run("python3 generateTestdata.py --metagraph " + metagraphBin \
    + " --binsize " + str(args.binsize) \
    + " --k " + str(k) \
    + " --nsequences " + str(args.nsequences) \
    + " --nspecies " + str(args.nspecies) \
    + " --sequence-lengths " + str(args.seqlen) \
    + " --similarity " + str(args.similarity)): 
    if run("python3 expectedSeeds.py --binsize " + str(args.binsize) + " --redmask " + str(args.redmask) + " --stringhasher " + stringHasherBin + " --thinning " + str(args.thinning)):
        if run("python3 expectedLinks.py"):
            if run("python3 expectedCubes.py --tilesize " + str(args.tilesize)):
                run("python3 expectedCubeScores.py --tilesize " + str(args.tilesize) \
                    + " --binsize " + str(args.binsize) \
                    + " --chunksize " + str(args.chunksize) \
                    + " --nsequences " + str(args.nsequences) \
                    + " --nspecies " + str(args.nspecies) \
                    + " --nsubtiles " + str(args.nsubtiles) \
                    + " --pnorm " + str(args.pnorm) \
                    + " --sequence-lengths " + str(args.seqlen))
                run("python3 expectedRegionTuples.py --alpha " + str(args.alpha) \
                    + " --beta " + str(args.beta) + " --gamma " + str(args.gamma) \
                    + " --reference-sequence-lengths " + str(args.seqlen) \
                    + " --tilesize " + str(args.tilesize))

sys.exit()