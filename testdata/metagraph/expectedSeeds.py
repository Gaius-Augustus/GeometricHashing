import argparse
import json
import math
import os
import re
import subprocess

parser = argparse.ArgumentParser(description = "Create expected seeds",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--binsize",
                        dest = "binsize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="binsize of metagraph")
scriptArgs.add_argument("--stringhasher",
                        dest = "stringhasher", 
                        metavar = "FILE", 
                        type=argparse.FileType("r"),
                        help="stringHasher binary", 
                        required = True)
scriptArgs.add_argument("--redmask",
                        dest = "redmask", 
                        metavar = "BOOL", 
                        type=bool,
                        required=True,
                        help="Discard low-complexity seeds")
scriptArgs.add_argument("--thinning",
                        dest = "thinning", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="Discard roughly a fraction of 1/thinning k-mers")
args = parser.parse_args()
stringHasherBin = args.stringhasher.name



path = "."
masks = ["11101111010110110111111", "11110110011000110101010011111", "111110010010101001001001100011111", "111110001100100100000101000100001011111"]
span = max([len(mask) for mask in masks]) # spoiler: 39
weight = sum([int(b) for b in masks[0]])  # spoiler: 18
binsize = args.binsize # MUST MATCH generateTestdata.py

assert binsize >= 1



fastas = [fa for fa in os.listdir(path) if os.path.isfile(os.path.join(path,fa)) and fa[-3:] == ".fa"]

sequences = {}
for fa in fastas:
    sequences[fa] = {}
    with open(os.path.join(path,fa), "r") as fh:
        head = ""
        for line in fh:
            line = line.rstrip()
            if len(line) > 0:
                if line[0] == ">":
                    head = line[1:]
                else:
                    sequences[fa][head] = line

def keepKmer(kmer):
    run = subprocess.run([stringHasherBin + " " + kmer], shell = True, executable = '/bin/bash', stdout = subprocess.PIPE)
    kmerHash = int(run.stdout)
    if kmerHash % args.thinning == 1:
        return(False)
    else:
        return(True)

def getKmers(seq, span):
    assert len(seq) >= span 
    kmers = []
    for i in range(len(seq) - span + 1):
        kmer = seq[i:i+span]
        # print(kmer, "has hash", kmerHash)
        if re.search("^[ACGT]{"+str(span)+"}$", kmer) and keepKmer(kmer):
            kmers.append([kmer, i])

    return(kmers)

def roundPosition(pos, binsize):
    return(math.floor(pos/binsize) * binsize)

# each kmer only reported once per bin
def filterKmerBins(kmers, binsize):
    filtered = set()
    for kmerTuple in kmers:
        filtered.add(tuple([kmerTuple[0], roundPosition(kmerTuple[1], binsize)]))

    return [list(t) for t in filtered]
    #filtered = []
    #seen = {}
    #for kmerTuple in kmers:
    #    kmer = kmerTuple[0]
    #    p = roundPosition(kmerTuple[1], binsize)
    #    if p not in seen:
    #        seen[p] = set()
    #
    #    if kmer not in seen[p]:
    #        filtered.append([kmer, p])
    #        seen[p].add(kmer)
    #
    #return(filtered)

def getSeed(mask, seq):
    assert len(mask) <= len(seq)
    seed = ""
    for i in range(len(mask)):
        if int(mask[i]):
            seed = seed + seq[i]
            
    return(seed)

def keepSeed(seed):
    if args.redmask:
        bases = set([b for b in seed])
        return(len(bases) > 2)
    else:
        return(True)

def fillExpectedSeeds(kmers, masks):
    seeds = {}
    for kmerTuple in kmers:
        kmer = kmerTuple[0]
        pos = kmerTuple[1]
        for i in range(len(masks)):
            seed = getSeed(masks[i], kmer)
            if (keepSeed(seed)):
                if seed not in seeds:
                    seeds[seed] = [[] for mask in masks]  # [[],[],[],[]]
                
                # [genome, sequence, position, reverseStrand, kmer, span]
                seeds[seed][i].append([genome, sequence, pos, False, kmer, len(masks[i])])
        
    return(seeds)



expectedSeeds = {}
expectedSeedsGraph = {}
for genome in sequences:
    for sequence in sequences[genome]:
        kmers = getKmers(sequences[genome][sequence], span)
        kmersGraph = filterKmerBins(kmers, binsize)

        seeds = fillExpectedSeeds(kmers, masks)
        seedsGraph = fillExpectedSeeds(kmersGraph, masks)
        for seed in seeds:
            if seed not in expectedSeeds:
                expectedSeeds[seed] = seeds[seed]
            else:
                for i in range(len(masks)):
                    expectedSeeds[seed][i].extend(seeds[seed][i])

        for seed in seedsGraph:
            if seed not in expectedSeedsGraph:
                expectedSeedsGraph[seed] = seedsGraph[seed]
            else:
                for i in range(len(masks)):
                    expectedSeedsGraph[seed][i].extend(seedsGraph[seed][i])
        


# store data
with open(os.path.join(path,"expectedSeeds.json"), "w") as fh:
    json.dump(expectedSeeds, fh)

with open(os.path.join(path,"expectedSeedsGraph.json"), "w") as fh:
    json.dump(expectedSeedsGraph, fh)