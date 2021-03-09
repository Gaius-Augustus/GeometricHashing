import argparse
import itertools
import json
import math
import numpy as np
import os

parser = argparse.ArgumentParser(description = "Create expected cube scores",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--binsize",
                        dest = "binsize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="binsize of metagraph")
scriptArgs.add_argument("--chunksize",
                        dest = "chunksize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="chunksize (length) in geometric hashing")
scriptArgs.add_argument("--nsequences",
                        dest = "nsequences", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="number of sequences to generate per species")
scriptArgs.add_argument("--nspecies",
                        dest = "nspecies", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="number of species to generate")
scriptArgs.add_argument("--nsubtiles",
                        dest = "nsubtiles", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="number of subcubes in GH (approximate)")
scriptArgs.add_argument("--pnorm",
                        dest = "pnorm", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="p for p-norm (GH)")
scriptArgs.add_argument("--sequence-lengths",
                        dest = "seqlen", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="length of each sequence")
scriptArgs.add_argument("--tilesize",
                        dest = "tilesize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="tilesize in geometric hashing")
args = parser.parse_args()



path = "."
binsize = args.binsize # MUST MATCH generateTestdata.py
chunksize = args.chunksize
b = args.nsubtiles
nsequences = args.nsequences
nspecies = args.nspecies
pnorm = args.pnorm
tilesize = args.tilesize
assert binsize >= 1
assert chunksize >= 1
assert b >= 1
assert nsequences >= 1
assert nspecies > 1
assert pnorm >= 1
assert tilesize >= 1



def round(pos, tilesize):
    return math.floor(pos/tilesize)

# occurrence: [genome, sequence, position, reverseStrand]
def cubeFromLink(link, tilesize):
    assert len(link) >= 1
    assert link[0][0] == "genome0.fa", "link " + str(link) + " invalid"

    refPos = link[0][2]
    cube = tuple([tuple([o[0],o[1],round(o[2]-refPos, tilesize), o[3]]) for o in link])
    return cube

def chunk(link, chunksize):
    psum = sum([occ[2] for occ in link])
    chnk = math.floor(psum/chunksize)
    assert chnk >= 0, "Chunk "+str(chnk)+" from link "+str(link)
    return chnk

def minChunk(cube, tilesize, chunksize):
    pos = np.zeros(len(cube))
    u = -tilesize * min([a[2] for a in cube])
    i0 = u
    shift = 1-tilesize if i0 > 0 else 0
    for i in range(len(cube)):
        pos[i] = cube[i][2]*tilesize + i0
        if pos[i] > 0:
            pos[i] += shift

    chnk = math.floor(sum(pos)/chunksize)
    assert chnk >= 0, "Chunk "+str(chnk)+" from cube "+str(cube)+", F "+str(tilesize)+", chunksize "+str(chunksize)
    return chnk

def maxChunk(cube, seqlens, tilesize, chunksize):
    pos = np.zeros(len(cube))
    v = min([seqlens[i] - tilesize*cube[i][2] for i in range(len(cube))])
    u = -tilesize * min([a[2] for a in cube])
    i0 = u
    pos[0] = i0 - 1 + v - u
    for i in range(1,len(cube)):
        ii = tilesize*(cube[i][2]+1) + pos[0] - 1
        pos[i] = min(ii, seqlens[i]-1)

    chnk = math.floor(sum(pos)/chunksize)
    assert chnk >= 0, "Chunk "+str(chnk)+" from cube "+str(cube)+", seqlens "+str(seqlens)+", F "+str(tilesize)+", chunksize "+str(chunksize)
    return chnk

# for cubes from a single sequence tuple
def score(cubes, seqlens, linkfrac, tilesize, chunksize, nsubtiles, pnorm):
    #nlinks = sum([len(cubes[cube]) for cube in cubes])
    #linkfrac = nlinks
    #for seqlen in seqlens:
    #    linkfrac /= seqlen # nlinks / (seqlen**cubedim)

    scores = {}
    for cube in cubes:
        nspecies = len(cube)
        assert nspecies >= 2, "Less than two dimensions in cube "+str(cube)

        nsubtilesPerSide = math.floor(nsubtiles**(1/(nspecies-1)))
        subtilewidth = math.floor(tilesize/nsubtilesPerSide)
        nsubtilesPerSide = math.ceil(tilesize/subtilewidth) # update to true value
        nsubtilesTrue = nsubtilesPerSide**(nspecies-1)
        chunkvol = chunksize * subtilewidth**(nspecies-1)

        assert nsubtilesTrue > 0, "nsubtiles "+str(nsubtiles)+", "+"subtilewidth "+str(subtilewidth)+", "+"nsubtilesPerSide "+str(nsubtilesPerSide)
        assert chunkvol > 0, "chunksize "+str(chunksize)+", nspecies "+str(nspecies)

        maxchunk = maxChunk(cube, seqlens, tilesize, chunksize)
        minchunk = minChunk(cube, tilesize, chunksize)
        assert minchunk <= maxchunk, "min "+str(minchunk)+", max "+str(maxchunk)
        nchunks = maxchunk - minchunk + 1

        psum = 0
        chunkCounts = {}
        for link in cubes[cube]:
            subcube = cubeFromLink(link[0], subtilewidth)
            chnk = chunk(link[0], chunksize)
            assert chnk >= minchunk, "cube "+str(cube)+", link "+str(link[0])+", chnk "+str(chnk)+", min "+str(minchunk)
            assert chnk <= maxchunk, "cube "+str(cube)+", link "+str(link[0])+", seqlens "+str(seqlens)+", chnk "+str(chnk)+", max "+str(maxchunk)
            chunkID = tuple([subcube, chnk])
            if chunkID not in chunkCounts:
                chunkCounts[chunkID] = 0
            
            chunkCounts[chunkID] += 1

        for chnk in chunkCounts:
            p = chunkCounts[chnk]**pnorm
            psum += p

        assert psum > 0, str(psum)
        psumRoot = psum**(1/pnorm)

        n = nchunks * nsubtilesTrue
        assert n > 0, "nchunks "+str(nchunks)+", nsubtiles "+str(nsubtilesTrue)
                
        lamb = chunkvol * linkfrac            
        assert lamb > 0, "vol "+str(chunkvol)+", frac "+str(linkfrac)

        norm = lamb * (n**(1/pnorm))
        assert norm > 0, "lambda "+str(lamb)+", n "+str(n)+", p-th root of n "+str(n**(1/pnorm))

        s = psumRoot / norm
        assert s > 0, "psum "+str(psum)+", psumRoot "+str(psumRoot)+", norm "+str(norm)

        # print(cube, "--", s, "=", psum, "/", norm)
        # print("nsubtilesPerSide", nsubtilesPerSide)
        # print("subtilewidth", subtilewidth)
        # print("nsubtilesTrue", nsubtilesTrue)
        # print("chunkvol", chunkvol)
        # print("maxchunk", maxchunk)
        # print("minchunk", minchunk)
        # print("nchunks", nchunks)
        # print("seqlens", seqlens)
        # print("nlinks", nlinks)
        # print("linkfrac", linkfrac)
        # print("lambda", lamb)
        # print()

        scores[cube] = s

    return scores

# {"i": {"cube": cube, "links": [link, ...]}, ...}
# cube = (tiledist, ...), tiledist = (gen, seq, dist, strand)
# link = [[occ, ...], span], occ = [gen, seq, pos, strand]
def scoreCubes(cubes, seqlen, nspecies, nsequences, tilesize, chunksize, nsubtiles, pnorm):
    nlinksTotal = 0
    seqtupleToCubes = {}
    for i in cubes:
        cube = tuple([tuple(td) for td in cubes[i]['cube']])
        seqtuple = tuple([tuple([td[0], td[1], td[3]]) for td in cube])
        if seqtuple not in seqtupleToCubes:
            seqtupleToCubes[seqtuple] = {}

        seqtupleToCubes[seqtuple][cube] = cubes[i]['links']
        nlinksTotal += len(cubes[i]['links'])

    linkfrac = nlinksTotal
    for _ in range(nspecies):
        lensum = seqlen*nsequences
        linkfrac /= lensum # nLinks / nPossibleLinks

    scores = {}
    for seqtuple in seqtupleToCubes:
        seqlens = [seqlen for _ in range(len(seqtuple))]
        seqscores = score(seqtupleToCubes[seqtuple], seqlens, linkfrac, tilesize, chunksize, nsubtiles, pnorm)
        for cube in seqscores:
            scores[cube] = seqscores[cube]

    outscores = {}
    i = 0
    for cube in scores:
        outscores[str(i)] = {'cube': cube, 'score': scores[cube]}
        i += 1

    return outscores

def roundPosition(pos, binsize):
    return(math.floor(pos/binsize) * binsize)



for infile in ["expectedCubes.json", "expectedCubesGraph.json"]:
    # adjust seqlen
    seqlen = args.seqlen
    assert seqlen >= 1
    if infile == "expectedCubesGraph.json":
        seqlen = roundPosition(seqlen, binsize) + 1

    # load data
    with open(os.path.join(path,infile), "r") as fh:
        expectedCubes = json.load(fh)

    # score cubes
    print(infile)
    expectedScores = scoreCubes(expectedCubes, seqlen, nspecies, nsequences, tilesize, chunksize, b, pnorm)

    # save data
    outfile = infile.replace("Cubes", "CubeScores")
    with open(os.path.join(path,outfile), "w") as fh:
        json.dump(expectedScores, fh)