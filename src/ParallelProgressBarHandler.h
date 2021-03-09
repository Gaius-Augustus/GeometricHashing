#ifndef PARALLELPROGRESSBARHANDLER_H
#define PARALLELPROGRESSBARHANDLER_H

#include <mutex>
#include "mabl3/ProgressBar.h"

using namespace mabl3;

//! Take care of concurrent ProgressBar updates
/*! Make sure the instance does not live longer than the
 *    ProgressBar and the lock! */
class ParallelProgressBarHandler {
public:
    ParallelProgressBarHandler(ProgressBar & pb,
                               std::unique_lock<std::mutex> & lock)
        : lock_{lock}, pb_{pb}, processed_{0}, step_{pb.step()} {}
    ~ParallelProgressBarHandler() {
        if (!pb_.quiet()) {
            update(); // final update before destruction
        }
    }
    ParallelProgressBarHandler & operator++() {
        if (!pb_.quiet()) {
            ++processed_;
            update();
        }
        return *this;
    }

private:
    void update() {
        if (processed_>= step_) {
            lock_.lock();
            pb_ += processed_;
            lock_.unlock();
            processed_ = 0;
        }
    }

    std::unique_lock<std::mutex> & lock_;
    ProgressBar & pb_;
    size_t processed_;
    size_t const step_;
};

#endif // PARALLELPROGRESSBARHANDLER_H
