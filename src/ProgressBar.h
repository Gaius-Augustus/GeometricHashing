#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

// inspired from here: https://stackoverflow.com/a/14539953/5586698

#include <cstdlib>
#include <cmath>
#include <sys/ioctl.h> // For ioctl, TIOCGWINSZ
#include <unistd.h> // For STDOUT_FILENO

class ProgressBar {
public:
    ProgressBar(size_t total, bool quiet = false, size_t barwidth = 0) :
        barw_{barwidth}, current_{0},
        step_{static_cast<size_t>(std::ceil(static_cast<double>(total) / 100.))},
        next_{step_}, quiet_{quiet},  total_{total} {
        // if barw_ is 0 (default), determine fitting barwidth
        if (!quiet_ && barw_ == 0) {
            struct winsize wsize;
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);
            if (wsize.ws_col < 8) { // not enough space for a bar (min: `[=] 100%`)
                quiet_ = true;
            } else {
                // 7 chars for `[] 100%`
                barw_ = wsize.ws_col - 7;
            }
        }
        if (!quiet_) {
            std::cout << "Bar -- total: " << total_ << ", step: " << step_ << std::endl;
            std::cout << "[";
            for (size_t i = 0; i < barw_; ++i) {
                std::cout << " ";
            }
            std::cout << "] 0%\r";
            std::cout.flush();
        }
    }
    ProgressBar & operator++() {    // prefix increment ( https://docs.microsoft.com/de-de/cpp/cpp/increment-and-decrement-operator-overloading-cpp?view=vs-2017 )
        current_++;
        update();
        return *this;
    }
    ProgressBar operator++(int) {   // postfix increment ( https://docs.microsoft.com/de-de/cpp/cpp/increment-and-decrement-operator-overloading-cpp?view=vs-2017 )
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    ProgressBar & operator+=(int progress) {    // with negative argument also possible to rewind the bar
        current_ += progress;
        update();
        return *this;
    }
    size_t step() const { return step_; }
    void update(size_t idx) {
        current_ = idx;
        update();
    }

private:
    size_t barw_;
    size_t current_;
    size_t step_;
    size_t next_;
    bool quiet_;
    size_t total_;

    void printBar(size_t percent) {
        if (!quiet_) {
            std::cout << "[";
            auto bar = static_cast<size_t>(std::floor(static_cast<double>(barw_ * percent) / 100.));
            for (size_t i = 0; i < barw_; ++i) {
                if (i <= bar) std::cout << "=";
                else std::cout << " ";
            }
            std::cout << "] " << percent << "%\r";
            std::cout.flush();
        }
    }
    void update() {
        if (current_ >= next_) {
            auto progress = static_cast<double>(current_) / static_cast<double>(total_);
            auto percent = static_cast<size_t>(std::floor(progress * 100.));
            next_ += step_;
            if (next_ > total_) { next_ = total_; }
            printBar(percent);
        }
    }
};

#endif // PROGRESSBAR_H
