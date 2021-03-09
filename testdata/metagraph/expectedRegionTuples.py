import argparse
import json
import math
import numpy as np
import os

parser = argparse.ArgumentParser(description = "Create expected region tuples",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--alpha",
                        dest = "alpha", 
                        metavar = "FLOAT", 
                        type=float,
                        required=True,
                        help="alpha parameter for region tuple extraction")
scriptArgs.add_argument("--beta",
                        dest = "beta", 
                        metavar = "FLOAT", 
                        type=float,
                        required=True,
                        help="beta parameter for region tuple extraction")
scriptArgs.add_argument("--gamma",
                        dest = "gamma", 
                        metavar = "FLOAT", 
                        type=float,
                        required=True,
                        help="gamma parameter for region tuple extraction")
scriptArgs.add_argument("--reference-sequence-lengths",
                        dest = "refseqlen", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="lengths of reference sequences")
scriptArgs.add_argument("--tilesize",
                        dest = "tilesize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="tilesize in geometric hashing")
args = parser.parse_args()



path = "."
tilesize = args.tilesize
refseqlen = args.refseqlen

alpha = args.alpha
beta = args.beta
gamma = args.gamma
assert refseqlen >= 1
assert tilesize >= 1


def round(pos):
    return math.floor(pos/tilesize)

# occurrence: [genome, sequence, position, reverseStrand]
def tuplesFromLinks(links):
    tuples = []
    ntiles = round(refseqlen-1) + 1
    fi = np.zeros(ntiles)
    for link in links:
        # link: [[occ (i.e. [gen, seq, pos, strand]), occ, ...], span]
        assert link[0][0][0] == "genome0.fa"
        tile = round(link[0][0][2])
        assert tile < ntiles, "tile " + str(tile) + " from link " + str(link) + " greater than ntiles (" + str(ntiles) + ")"
        fi[tile] = fi[tile] + 1

    pprev = {0: 0}
    for i in range(1,ntiles):
        if fi[i-1] > 0:
            pprev[i] = i-1
        else:
            pprev[i] = pprev[i-1]

    # forward
    A = [float('-inf')] * (ntiles+1)
    B = [0] * (ntiles+1)
    m1choosen = [False] * (ntiles+1)
    for i in range(1,ntiles+1):
        tile = i-1
        assert tile in pprev, str(tile) + " not in " + str(pprev)
        m1 = B[i-1] - beta - gamma
        m2 = A[i-1] - beta*(tile - pprev[tile])
        m1choosen[i] = True if m1 > m2 else False
        A[i] = alpha*fi[tile] + max(m1,m2)
        B[i] = max(A[i], B[i-1])

    # backward
    # i in [nitles+1, ntiles, ..., 1]
    intervalOpen = False
    for i in range(ntiles,0,-1):
        assert i < len(A), str(i) + " greater than A (" + str(len(A)) + ")"
        assert i < len(B), str(i) + " greater than B (" + str(len(B)) + ")"
        assert i < len(m1choosen), str(i) + " greater than m1choosen (" + str(len(m1choosen)) + ")"

        if A[i] >= B[i]:
            if not intervalOpen:
                tuples.append([0,i-1])
                intervalOpen = True
            
            if m1choosen[i]:
                tuples[-1][0] = i-1
                intervalOpen = False
        elif intervalOpen:
            tuples[-1][0] = i-1
            intervalOpen = False

    return(tuples)

    

def getRegionTuples(cubes):
    tuples = []
    for i in cubes:
        cube = cubes[i]["cube"]
        links = cubes[i]["links"]
        tuples.append({"cube": cube,
                       "tuples": tuplesFromLinks(links)})

    return tuples



for infile in ["expectedCubesGrouped.json", "expectedCubesGroupedGraph.json"]:
    # load data
    with open(os.path.join(path,infile), "r") as fh:
        expectedCubes = json.load(fh)

    # create cubes
    expectedRegionTuples = getRegionTuples(expectedCubes)

    # save data
    outfile = infile.replace("CubesGrouped", "RegionTuples")
    with open(os.path.join(path,outfile), "w") as fh:
        json.dump(expectedRegionTuples, fh)