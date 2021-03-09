import json
import os

path = "."



# recursive link creation
def linksFromOccDict(occDict, incompleteLinks):
    assert len(occDict) > 0, "Cannot create links from empty dict"

    genomeToOccs = dict(occDict) # make true copy
    genome = list(genomeToOccs.keys())[0]
    occs = genomeToOccs.pop(genome)
    links = []
    if len(incompleteLinks) > 0:
        for linkl in incompleteLinks:
            for occ in occs:
                link = list(linkl) # true copy
                link.append(occ)
                links.append(link)
    else:
        links = [[o] for o in occs]

    if len(genomeToOccs) == 0:
        return(links)
    else:
        return(linksFromOccDict(genomeToOccs, links))



def stripKmerFromOccs(occs):
    # occurrence: [genome, sequence, position, reverseStrand, kmer, span] -> [genome, sequence, position, reverseStrand, span]
    # return unique set of occurrences (i.e. remove duplicates ["g0","s0",0,False,"ACGT",4] and ["g0","s0",0,False,"GCTA",4])
    uniqueOccs = set()
    for occ in occs:
        uniqueOccs.add(tuple([occ[0], occ[1], occ[2], occ[3], occ[5]]))

    return [list(o) for o in uniqueOccs]



# exclude if
#    no occurrence in reference (genome0.fa)
#    no occurrence besides reference
def makeLinks(seeds):
    links = {}
    for seed in seeds:
        for mask in range(len(seeds[seed])):
            genomeToOccs = {}
            occs = stripKmerFromOccs(seeds[seed][mask])
            for occ in occs:#seeds[seed][mask]:
                # [genome, sequence, position, reverseStrand, span]
                genome = occ[0]
                if genome not in genomeToOccs:
                    genomeToOccs[genome] = []

                genomeToOccs[genome].append(occ)

            # create links and count equivalent links
            if (len(genomeToOccs.keys()) > 1) and ("genome0.fa" in genomeToOccs):
                linkv = linksFromOccDict(genomeToOccs, [])
                for link in linkv:
                    # occurrence: [genome, sequence, position, reverseStrand, span]
                    linkt = tuple(sorted([tuple(occ[0:4]) for occ in link]))
                    if linkt not in links:
                        links[linkt] = {'count': 0, 'span': link[0][4]}

                    links[linkt]['count'] = links[linkt]['count'] + 1
                    links[linkt]['span'] = max(links[linkt]['span'], link[0][4])

                    # quality assessment
                    genomeSeen = set()
                    for occ in linkt:
                        gen = occ[0]
                        assert gen not in genomeSeen, "Invalid link " + str(linkt) + ": " + str(links[linkt]) + " created (genomeToOccs: " + str(genomeToOccs) + ")"
                        genomeSeen.add(gen)

    # json cannot handle tuples (as dict keys), so transform
    linkOutput = {}
    i = 0
    for linkt in links:
        linkOutput[str(i)] = {"link": linkt,
                              "count": links[linkt]['count'],
                              "span": links[linkt]['span']}
        i = i + 1

    return(linkOutput)




for infile in ["expectedSeeds.json", "expectedSeedsGraph.json"]:
    # load data
    with open(os.path.join(path,infile), "r") as fh:
        expectedSeeds = json.load(fh)

    # create links
    expectedLinks = makeLinks(expectedSeeds)

    # save data
    outfile = infile.replace("Seeds", "Links")
    with open(os.path.join(path,outfile), "w") as fh:
        json.dump(expectedLinks, fh)