#ifndef CONTAINERCHUNKS_H
#define CONTAINERCHUNKS_H

#include <cstddef>
#include <tuple>

#include "mabl3/disambiguateSTLContainer.h"

namespace mabl3 {

/* container:
 *     array
 *     (unordered) map
 *     (unordered) set
 *     vector
 * associative:
 *     (unordered) map
 *     (unordered) set
 * sequential:
 *     array
 *     vector
 * map:
 *     (unordered) map */

template <typename ContainerA, typename ContainerB> // for associative container
void containerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator begin, typename ContainerB::const_iterator end, std::true_type) {
    chunk.insert(begin, end);
}

template <typename ContainerA, typename ContainerB> // for sequential container
void containerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator begin, typename ContainerB::const_iterator end, std::false_type) {
    chunk.insert(chunk.end(), begin, end);
}

template <typename ContainerA, typename ContainerB> // merge container
void containerChunkInsert(ContainerA & chunk, typename ContainerB::const_iterator begin, typename ContainerB::const_iterator end) {
    containerChunkInsert_impl<ContainerA, ContainerB>(chunk, begin, end, isAssociativeContainer(chunk));
}



template <typename ContainerA, typename ContainerB> // for map container
void associativeContainerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator elem, std::true_type) {
    chunk.insert({elem->first, elem->second});
}

template <typename ContainerA, typename ContainerB> // for set container
void associativeContainerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator elem, std::false_type) {
    chunk.emplace(*elem);
}

template <typename ContainerA, typename ContainerB> // for associative container
void containerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator elem, std::true_type) {
    associativeContainerChunkInsert_impl<ContainerA, ContainerB>(chunk, elem, isMapContainer(chunk));
}

template <typename ContainerA, typename ContainerB> // for sequential container
void containerChunkInsert_impl(ContainerA & chunk, typename ContainerB::const_iterator elem, std::false_type) {
    chunk.emplace_back(*elem);
}

template <typename ContainerA, typename ContainerB> // single element insert
void containerChunkInsert(ContainerA & chunk, typename ContainerB::const_iterator elem) {
    containerChunkInsert_impl<ContainerA, ContainerB>(chunk, elem, isAssociativeContainer(chunk));
}



//! Return \c std::pair of two \c size_t that give the amount to shift a container.begin() iterator to begin and past-the-end of a chunk from input container, chunks evenly distributed over a number of threads
/*! \c thread count must be zero based! */
template <typename Container>
auto ContainerChunkShift(Container const & container,
                         size_t thread, size_t nThreads) {
    if (thread >= nThreads) { throw std::runtime_error("[ERROR] -- getTupleChunk -- thread count must be zero based and lower than nThreads"); }
    auto baseChunkSize = container.size() / nThreads;
    auto remainder = container.size() % nThreads;
    size_t shiftBegin = 0;
    size_t shiftEnd = 0;
    if (remainder > 0) {
        // some threads get one element more than others
        if (thread < remainder) {
            shiftBegin += thread;
            shiftEnd = 1;
        } else {
            shiftBegin += remainder;
        }
    }
    return std::pair<size_t, size_t>{
        (baseChunkSize * thread) + shiftBegin,                  // amount needed to shift container.begin() to first element of this threads chunk
        (baseChunkSize * (thread + 1)) + shiftBegin + shiftEnd  // amount needed to shift container.begin() to past-the-last element of this threads chunk
    };
}

//! Return an iterator to the begin of a chunk from input container, chunks evenly distributed over a number of threads
template <typename Container>
auto getContainerChunkBegin(Container const & container,
                            size_t thread, size_t nThreads) {
    auto shiftValue = ContainerChunkShift<Container>(container, thread, nThreads);
    auto itBegin = container.begin();
    std::advance(itBegin, shiftValue.first);
    return itBegin;
}

//! Return a const iterator to past-the-end of a chunk from input container, chunks evenly distributed over a number of threads
template <typename Container>
auto getContainerChunkEnd(Container const & container,
                          size_t thread, size_t nThreads) {
    auto shiftValue = ContainerChunkShift<Container>(container, thread, nThreads);
    auto itEnd = container.cbegin();
    std::advance(itEnd, shiftValue.second);
    return itEnd;
}

//! Return chunks from input container, evenly distributed over a number of threads
/*! \c thread count must be zero based! */
template <typename Container>
auto getContainerChunk(Container const & container,
                       size_t thread, size_t nThreads) {
    Container chunk;
    if (container.size() == 0) { return chunk; }
    auto begin = getContainerChunkBegin<Container>(container, thread, nThreads);
    auto end = getContainerChunkEnd<Container>(container, thread, nThreads);
    containerChunkInsert<Container, Container>(chunk, begin, end);
    return chunk;
}

} // NAMESPACE mabl3

#endif // CONTAINERCHUNKS_H
