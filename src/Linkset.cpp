#include "Linkset.h"

//! External helper function to deal with template mismatches
void addLinkExt(Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan> & linkset, LinkPtr link) {
    linkset.addLink(link);
}
//! External helper function to deal with template mismatches
void addLinkExt(Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan> & linkset, LinkPtr link) {
    linkset.addLink(*link);
}

template<typename LinkType, typename LinkTypeHash, typename LinkTypeEqual>
void Linkset<LinkType, LinkTypeHash, LinkTypeEqual>::addLink(LinkType link) {
    auto linkIt = linkset_.find(link);
    if (linkIt == linkset_.end()) {
        linkset_.insert({link, 0});
    } else if (linkIt->first.span() < link.span()) {
        auto count = linkIt->second;
        linkset_.erase(linkIt);
        linkset_.insert({link, count});
    }
    linkset_[link] += 1;
}



template<typename LinkType, typename LinkTypeHash, typename LinkTypeEqual>
void Linkset<LinkType, LinkTypeHash, LinkTypeEqual>::createLinks(std::vector<KmerOccurrence> const & occurrences, size_t span) {
    auto processingFunction = [this, span](std::vector<tsl::hopscotch_map<size_t, // seqID
                                                                          tsl::hopscotch_set<KmerOccurrence,
                                                                                             KmerOccurrencePositionHash,
                                                                                             KmerOccurrencePositionEqual>>> const & occurrenceMap,
                                           size_t nPossible) {
        // flatten occurrenceMap for link creation
        std::vector<std::vector<KmerOccurrence>> allOccs{};
        for (size_t i = 0; i < idMapping_->numGenomes(); ++i) {
            std::vector<KmerOccurrence> occs;
            for (auto&& elem : occurrenceMap.at(i)) {
                if (elem.second.size()) {
                    occs.insert(occs.end(), elem.second.begin(), elem.second.end());
                }
            }
            if (occs.size()) {
                allOccs.emplace_back(occs);
            }
        }
        // create links, sampling if too many
        if (nPossible > config_->matchLimit()) {
            tsl::hopscotch_set<size_t> linkIDs;
            std::default_random_engine rng(rd_());  // create default rng with seed from rd_
            std::uniform_int_distribution<size_t> runif(1, nPossible-1); // uniformly distributed random values in given interval
            while (linkIDs.size() < config_->matchLimit()) {    // create matchLimit distinct IDs
                auto randomID = runif(rng);
                linkIDs.emplace(randomID);   // set -> no duplicate IDs
            }
            for (auto linkID : linkIDs) {
                auto tiles = cartesianProductByID(linkID, allOccs);
                addLink(LinkType(tiles, span));   // add link to linkset and increase link count
            }
        } else {
            for (size_t id = 0; id < nPossible; ++id) {
                auto tiles = cartesianProductByID(id, allOccs);
                addLink(LinkType(tiles, span));   // add link to linkset and increase link count
            }
        }
    };
    // create valid links
    auto valid = processOccurrences(occurrences, processingFunction, config_->hasse());
    if (!valid) { ++numDiscarded_; }
}



template<typename LinkType, typename LinkTypeHash, typename LinkTypeEqual>
void Linkset<LinkType, LinkTypeHash, LinkTypeEqual>::createRelevantLinks(std::vector<KmerOccurrence> const & occurrences, size_t span,
                                                                         tsl::hopscotch_set<std::shared_ptr<Cube const>, CubePtrHash, CubePtrEqual> const & relevantCubes) {
    auto processingFunction = [this,
                               &relevantCubes,
                               span](std::vector<tsl::hopscotch_map<size_t, // seqID
                                                                    tsl::hopscotch_set<KmerOccurrence,
                                                                                       KmerOccurrencePositionHash,
                                                                                       KmerOccurrencePositionEqual>>> const & occurrenceMap,
                                      size_t nPossibleGlobal) {
        (void)nPossibleGlobal;
        for (auto&& cube : relevantCubes) { // do not sample noisy links
            // flatten occurrenceMap for link creation
            std::vector<std::vector<KmerOccurrence>> allOccs{};
            size_t nPossible = 1;
            bool cubeHasLinks = true;
            for (auto&& td : cube->tiledistance()) {
                std::vector<KmerOccurrence> occs;
                auto gid = td.genome();
                auto sid = td.sequence();
                if (occurrenceMap.at(gid).find(sid) != occurrenceMap.at(gid).end()
                        && occurrenceMap.at(gid).at(sid).size()) {
                    occs.insert(occs.end(), occurrenceMap.at(gid).at(sid).begin(), occurrenceMap.at(gid).at(sid).end());
                } else {
                    cubeHasLinks = false;
                    break;
                }
                allOccs.emplace_back(occs);
                nPossible *= occs.size(); // [ATTENTION] matchLimit works differently here: on cube level!
            }
            if (!cubeHasLinks) { continue; } // next cube
            // create links, sampling if too many
            std::vector<size_t> linkIDVector;
            if (nPossible > config_->matchLimit()) {
                tsl::hopscotch_set<size_t> linkIDs;
                std::default_random_engine rng(rd_());  // create default rng with seed from rd_
                std::uniform_int_distribution<size_t> runif(1, nPossible-1); // uniformly distributed random values in given interval
                while (linkIDs.size() < config_->matchLimit()) {    // create matchLimit distinct IDs
                    auto randomID = runif(rng);
                    linkIDs.emplace(randomID);   // set -> no duplicate IDs
                }
                linkIDVector = std::vector<size_t>(linkIDs.begin(), linkIDs.end());
            } else {
                linkIDVector.resize(nPossible);
                std::iota(linkIDVector.begin(), linkIDVector.end(), 0); // fill vector with sequence (0 .. numLinks-1) (i.e. create all links)
            }
            for (auto linkID : linkIDVector) {
                auto tiles = cartesianProductByID(linkID, allOccs);
                LinkPtr link(tiles, span);
                if (relevantCubes.find(std::make_shared<Cube>(*link, config_->tileSize())) != relevantCubes.end()) {
                    addLinkExt(*this, link);
                } else if (link.dimensionality() > 2 && config_->hasse()) {
                    // strip 3+ dimension from link and see if link fits
                    if (relevantCubes.find(
                                std::make_shared<Cube>(*(LinkPtr{std::vector<KmerOccurrence>{link.occurrence(0), link.occurrence(1)}, span}),
                                                       config_->tileSize())
                                ) != relevantCubes.end()) {
                        addLinkExt(*this, link);
                    }
                }
            }
        }
    };

    // create valid links
    auto valid = processOccurrences(occurrences, processingFunction, config_->hasse());
    if (!valid) { ++numDiscarded_; }
}



template<typename LinkType, typename LinkTypeHash, typename LinkTypeEqual>
bool Linkset<LinkType, LinkTypeHash, LinkTypeEqual>::processOccurrences(std::vector<KmerOccurrence> const & occurrences,
                                                                        std::function<void(std::vector<tsl::hopscotch_map<size_t, // seqID
                                                                                                                          tsl::hopscotch_set<KmerOccurrence,
                                                                                                                                             KmerOccurrencePositionHash,
                                                                                                                                             KmerOccurrencePositionEqual>>> const &,
                                                                                           size_t)> processingFunction,
                                                                        bool hasse) const {
    std::vector<size_t> genomeOccCount(idMapping_->numGenomes(), 0); // count occurrences per genome
    std::vector<tsl::hopscotch_map<size_t, // seqID
                                   tsl::hopscotch_set<KmerOccurrence,
                                                      KmerOccurrencePositionHash,
                                                      KmerOccurrencePositionEqual>>
               > occurrenceMap(idMapping_->numGenomes()); // one map for each genome
    for (auto&& occ : occurrences) {
        ++(genomeOccCount.at(occ.genome()));
        occurrenceMap.at(occ.genome())[occ.sequence()].insert(occ);
    }
    if (genomeOccCount.at(0) == 0) { return false; }    // not in ref
    if (genomeOccCount.at(0) == occurrences.size()) { return false; }   // only in ref

    // occurrencePerGenome min/max
    bool brk = false;
    for (size_t i = 0; i < idMapping_->numGenomes(); ++i) {
        auto c = genomeOccCount.at(i);
        if ((c > config_->occurrencePerGenomeMax())
                || (c > 0 && c < config_->occurrencePerGenomeMin())) {
            brk = true;
            break;
        }
    }
    if (brk) { return false; }
    // occurrence per sequence max -> delete cases with too many occs
    size_t nPossible = 1;
    size_t nonRefCount = 0;
    for (size_t i = 0; i < idMapping_->numGenomes(); ++i) {
        size_t occCount = 0;
        for (auto&& elem : occurrenceMap.at(i)) {
            if (elem.second.size() > config_->occurrencePerSequenceMax()) {
                occurrenceMap.at(i).at(elem.first).clear();
            } else {
                occCount += elem.second.size();
            }
        }
        if (i == 0) {
            nPossible *= occCount; // if no occs in ref remain, nPossible is zero
        } else {
            nonRefCount += occCount;
            if (hasse) {
                if (occCount) { nPossible *= occCount; } // need to check whether only ref genome remains
            } else {
                nPossible *= occCount; // if no hasse and no occs (remain) in any genome, nPossible is zero
            }
        }
    }
    if (nPossible == 0) { return false; }
    if (nonRefCount == 0) { return false; }
    // match limit
    if (nPossible > config_->matchLimit() && config_->matchLimitDiscardSeeds()) { return false; }

    // at this point, seed is considered valid -> call processing function
    processingFunction(occurrenceMap, nPossible);
    return true;
}



template class Linkset<LinkPtr, LinkPtrHashIgnoreSpan, LinkPtrEqualIgnoreSpan>;
template class Linkset<Link, LinkHashIgnoreSpan, LinkEqualIgnoreSpan>;
