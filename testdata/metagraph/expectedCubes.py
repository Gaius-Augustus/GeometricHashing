import argparse
import itertools
import json
import math
import numpy as np
import os

parser = argparse.ArgumentParser(description = "Create expected cubes",
                                 formatter_class = argparse.RawTextHelpFormatter)
scriptArgs = parser.add_argument_group("Script Arguments")
scriptArgs.add_argument("--tilesize",
                        dest = "tilesize", 
                        metavar = "INT", 
                        type=int,
                        required=True,
                        help="tilesize in geometric hashing")
args = parser.parse_args()



path = "."
tilesize = args.tilesize
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

def groupLinks(links, tilesize):
    # links = [[occurrences, span], ...] with occurrences = [["genome", "sequence", pos, reverse], ...]
    # group links of same diagonals
    diagonals = {}
    for link in links:
        assert cubeFromLink(links[0][0], tilesize) == cubeFromLink(link[0], tilesize), "Not all links from same cube"
        occurrences = link[0]
        assert len(occurrences) >= 1, "Empty link"
        span = link[1]

        refPos = occurrences[0][2]
        diagonal = tuple([o[2] - refPos for o in occurrences])
        if diagonal not in diagonals:
            diagonals[diagonal] = {}

        # if multiple links on same diagonal at same start position, just take max span
        if str(refPos) in diagonals[diagonal]:
            diagonals[diagonal][str(refPos)][1] = max(diagonals[diagonal][str(refPos)][1], span)
        else:
            diagonals[diagonal][str(refPos)] = link

    # combine overlapping seeds by extending spans
    groupedLinks = []
    for diagonal in diagonals:
        groupedDiag = []
        sortedPos = sorted([int(k) for k in diagonals[diagonal]])
        for p in sortedPos:
            curr = diagonals[diagonal][str(p)]
            if len(groupedDiag) == 0:
                groupedDiag.append(curr)
            else:
                prev = groupedDiag[-1]
                prevEnd = prev[0][0][2] + prev[1] - 1 # pos + span - 1
                currEnd = curr[0][0][2] + curr[1] - 1
                overlap = prevEnd - curr[0][0][2]
                if overlap >= 0:
                    if currEnd > prevEnd:
                        newSpan = currEnd - prev[0][0][2] + 1
                        groupedDiag[-1][1] = newSpan
                else:
                    groupedDiag.append(curr)

        groupedLinks.extend(groupedDiag)

    return groupedLinks

def makeCubes(links, tilesize, group=False):
    cubes = {}
    for i in links:
        link = links[i]["link"]
        cube = cubeFromLink(link, tilesize)
        if cube not in cubes:
            cubes[cube] = []

        cubes[cube].append([link, links[i]['span']])

    cubesOutput = {}
    i = 0
    for cube in cubes:
        if group:
            cubesOutput[str(i)] = {"cube": cube,
                                   "links": groupLinks(cubes[cube], tilesize)}
        else:
            cubesOutput[str(i)] = {"cube": cube,
                                   "links": cubes[cube]}
        
        i = i + 1

    return cubesOutput

def linkFromSubcubes(cubesOutput):
    cubeToLinks = {}
    for idx in cubesOutput:
        cubeToLinks[cubesOutput[idx]['cube']] = cubesOutput[idx]['links']

    subcubesOutput = {}
    i = 0
    for cube in cubeToLinks:
        links = []
        subcubes = []
        for r in range(1, len(cube)+1):
            subcubes.extend([c for c in itertools.combinations(cube, r)])

        for subcube in subcubes:
            if subcube in cubeToLinks:
                links.extend(cubeToLinks[subcube])

        subcubesOutput[str(i)] = {"cube": cube,
                                  "links": links}
        i = i + 1

    return subcubesOutput



for infile in ["expectedLinks.json", "expectedLinksGraph.json"]:
    # load data
    with open(os.path.join(path,infile), "r") as fh:
        expectedLinks = json.load(fh)

    for group in [True, False]:
        # create cubes
        expectedCubes = makeCubes(expectedLinks, tilesize, group)
        expectedCubesWithSubcubeLinks = linkFromSubcubes(expectedCubes)

        # save data
        if group:
            outfile = infile.replace("Links", "CubesGrouped")
        else:
            outfile = infile.replace("Links", "Cubes")

        with open(os.path.join(path,outfile), "w") as fh:
            json.dump(expectedCubes, fh)

        if group:
            outfile = infile.replace("Links", "CubesGroupedWithSubcubeLinks")
        else:
            outfile = infile.replace("Links", "CubesWithSubcubeLinks")

        with open(os.path.join(path,outfile), "w") as fh:
            json.dump(expectedCubesWithSubcubeLinks, fh)