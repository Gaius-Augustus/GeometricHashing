#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <chrono>
#include <iostream>

//! helper class to easily introduce time measurement into code
class Timestep {
public:
    enum unit {nanoseconds, microseconds, milliseconds, seconds, minutes, hours};

    Timestep(std::string const & description)
        : description_{description},
          end_{std::chrono::system_clock::now()},
          ended_{false},
          start_{std::chrono::system_clock::now()} { std::cout << "Starting '" << description << "'" << std::endl; }
    auto elapsed(Timestep::unit unit) const {
        auto now = (ended_) ? end_ : std::chrono::system_clock::now();
        switch(unit) {
            case nanoseconds : return std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_).count();
            case microseconds : return std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
            case milliseconds : return std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count();
            case seconds : return std::chrono::duration_cast<std::chrono::seconds>(now - start_).count();
            case minutes : return std::chrono::duration_cast<std::chrono::minutes>(now - start_).count();
            case hours : return std::chrono::duration_cast<std::chrono::hours>(now - start_).count();
            default : return std::chrono::duration_cast<std::chrono::seconds>(now - start_).count();    // suppress compiler warning
        }
    }
    void end() {
        end_ = std::chrono::system_clock::now();
        ended_ = true;
    }
    //! automatically chose the timestep
    void endAndPrint() {
        end();
        print();
    }
    void endAndPrint(Timestep::unit unit) {
        end();
        print(unit);
    }
    bool ended() const { return ended_; }
    //! automatically chose the timestep
    void print() const {
        auto now = (ended_) ? end_ : std::chrono::system_clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_).count();
        auto mus = std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count();
        auto s = std::chrono::duration_cast<std::chrono::seconds>(now - start_).count();
        auto m = std::chrono::duration_cast<std::chrono::minutes>(now - start_).count();
        if (ns < 10000) { print(nanoseconds); }
        else if (mus < 10000) { print(microseconds); }
        else if (ms < 10000) { print(milliseconds); }
        else if (s <= 180) { print(seconds); }
        else if (m <= 180) { print(minutes); }
        else { print(hours); }
    }
    void print(Timestep::unit unit) const {
        auto now = (ended_) ? end_ : std::chrono::system_clock::now();
        std::cout << description_ << " (took ";
        switch(unit) {
            case nanoseconds : std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_).count() << " ns)";
            break;
            case microseconds : std::cout << std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count() << " us)";
            break;
            case milliseconds : std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count() << " ms)";
            break;
            case seconds : std::cout << std::chrono::duration_cast<std::chrono::seconds>(now - start_).count() << " s)";
            break;
            case minutes : std::cout << std::chrono::duration_cast<std::chrono::minutes>(now - start_).count() << " min)";
            break;
            case hours : std::cout << std::chrono::duration_cast<std::chrono::hours>(now - start_).count() << " h)";
            break;
        }
        std::cout << std::endl;
    }
    bool running() const { return (!ended_); }

private:
    std::string const description_;
    std::chrono::time_point<std::chrono::system_clock> end_;
    bool ended_;
    std::chrono::time_point<std::chrono::system_clock> start_;
};

#endif // TIMESTEP_H
