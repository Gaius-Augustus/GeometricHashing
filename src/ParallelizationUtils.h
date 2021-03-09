#ifndef PARALLELIZATIONUTILS_H
#define PARALLELIZATIONUTILS_H

#include <iostream>
#include <mutex>
#include <thread>

#include "mabl3/ProgressBar.h"
#include "ContainerChunks.h"

//! Divide \param input container into equal chunks and execute \param function on each chunk in parallel
/*! \c function must accept InputContainer::const_iterator as the last two arguments (begin and end of InputContainer chunk). */
template<typename InputContainer, typename F, typename... Args>
void executeParallel(InputContainer const & input,
                     size_t nThreads,
                     F function,
                     Args... args) {
    if (nThreads == 0) { std::cerr << "[WARNING] -- executeParallel -- Zero threads specified, not running" << std::endl; }
    if (nThreads == 1) {
        // do not start extra threads for single threaded execution
        function(args..., input.begin(), input.end());
    } else {
        // start parallel execution
        std::vector<std::thread> threads;
        for (size_t i = 0; i < nThreads; ++i) {
            auto it = getContainerChunkBegin<InputContainer>(input, i, nThreads);
            auto end = getContainerChunkEnd<InputContainer>(input, i, nThreads);
            if (it != end) {
                threads.push_back(std::thread(function, args..., it, end));
            }
        }
        //   wait until all threads are finished
        std::for_each(threads.begin(), threads.end(), [](std::thread & t) { t.join(); });
    }
}



//! Progressbar that can be called from parallel threads
class ParallelProgressBar {
public:
    ParallelProgressBar(size_t total, bool quiet = false, size_t barwidth = 0)
        : mutex_{}, pb_{total, quiet, barwidth} {}

    //! Increase progress count
    void increase(size_t steps = 1) {
        std::unique_lock<std::mutex> lock(mutex_);
        pb_ += steps;
        lock.unlock();
    }
    //! Be careful when using this as there is no thread protection
    auto & unprotectedProgressBar() { return pb_; }


private:
    std::mutex mutex_;
    mabl3::ProgressBar pb_;
};

#endif // PARALLELIZATIONUTILS_H
