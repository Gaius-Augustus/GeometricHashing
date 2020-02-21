#include "Linkset.h"

void Linkset::addLink(std::shared_ptr<Link const> link) {
    // throw if genome or sequence name are unknown in idMapping_
    for (auto&& occ : link->occurrence()) {
        idMapping_->queryGenomeName(occ.genome());
        idMapping_->querySequenceName(occ.sequence());
    }
    // insert link / increase link counter
    linkset_[link]++;
}



void Linkset::createLinks(std::vector<KmerOccurrence> const & occurrences) {
    std::unordered_map<uint8_t, size_t> occurrenceCount;
    auto excludeFlag = false;
    for (auto&& occ : occurrences) {
        // track number of occurences per genome
        ++occurrenceCount[occ.genome()];
    }

    // exclude if
    //    no occurrence in reference (always ID 0)
    //    no occurrence besides reference
    //    too many or too few occurences in any genome
    if ((occurrenceCount.find(0) == occurrenceCount.end()) || (occurrenceCount.at(0) == 0)) {
        ++discardedNotInReference_;
        excludeFlag = true;
    } else if (occurrenceCount.size() < 2) {
        ++discardedOnlyInReference_;
        excludeFlag = true;
    } else {
        for (auto&& occs : occurrenceCount) {
            if (occs.second > occurrencePerGenomeMax_) {
                ++discardedTooManyOccurrences_;
                excludeFlag = true;
                break;
            } else if ((occs.second > 0) && (occs.second < occurrencePerGenomeMin_)) {
                ++discardedTooFewOccurrences_;
                excludeFlag = true;
                break;
            }
        }
    }

    if (!excludeFlag) { addLinksFromKmer(occurrences); }
}



void Linkset::addLinksFromKmer(std::vector<KmerOccurrence> const & occurrences) {
    // map genome to vector of occurrences
    std::map<uint8_t, std::vector<KmerOccurrence>> genomeToTile;
    size_t numPossibleLinks = (occurrences.size() > 0) ? 1 : 0;
    for (auto&& occ : occurrences) {
        genomeToTile[occ.genome()].emplace_back(occ);
        if (genomeToTile.at(occ.genome()).size() > 1) { numPossibleLinks *= 2; }
    }

    // limit number of link creations
    std::vector<size_t> linkIDVector;
    if (numPossibleLinks > linkLimit_) {
        if (!discardExceeding_) {
            linkIDVector.resize(linkLimit_);
            std::unordered_set<size_t> linkIDs;
            linkIDs.insert(0);              // always create first (all left occurrences)
            linkIDs.insert(linkLimit_-1);   // and last (all right occurrences) links
            std::default_random_engine rng(rd_());                                      // create default rng with seed from rd_
            std::uniform_int_distribution<size_t> runif(1, numPossibleLinks-1);         // uniformly distributed random values in given interval
            while (linkIDs.size() < linkLimit_) {   // create linkLimit_ distinct IDs
                auto randomID = runif(rng);
                linkIDs.emplace(randomID);   // set -> no duplicate IDs
            }
            linkIDVector.insert(linkIDVector.end(), linkIDs.begin(), linkIDs.end());
        }
    } else {
        linkIDVector.resize(numPossibleLinks);
        std::iota(linkIDVector.begin(), linkIDVector.end(), 0); // fill vector with sequence (0 .. numLinks-1) (i.e. create all links)
    }

    // create links
    for (auto linkID : linkIDVector) {
        std::vector<KmerOccurrence> tiles;
        size_t occurrenceProduct = 1;
        for (auto&& elem : genomeToTile) {
            auto const & occurrenceVector = elem.second;
            auto clock = std::floor(static_cast<double>(linkID) / static_cast<double>(occurrenceProduct));
            auto tileID = static_cast<size_t>(clock) % occurrenceVector.size();
            auto const & occurrence = occurrenceVector.at(tileID);
            tiles.emplace_back(occurrence); // use exact positions, let Cubeset take care of tiling
            occurrenceProduct *= occurrenceVector.size();
        }
        auto link = std::make_shared<Link>(tiles);
        linkset_[link]++;   // add link to linkset and increase link count
    }
}
