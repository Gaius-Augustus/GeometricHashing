#ifndef MEMORYMONITOR_H
#define MEMORYMONITOR_H

#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/sysinfo.h>

namespace MM {

inline int byteToMB(unsigned long byte) {
    auto dbyte = static_cast<double>(byte);
    auto factor = static_cast<double>(1) / static_cast<double>(1024 * 1024);
    return static_cast<int>(std::ceil(dbyte * factor));
}
unsigned long parseLine(char* line);
unsigned long getProcessUsedMemory(bool physical = false);

class MemoryMonitor {
public:
    MemoryMonitor();

    friend std::ostream& operator<<(std::ostream& out, MemoryMonitor const & m) {
        // compute total used virtual mem
        unsigned long virtualMemUsed = m.memInfo_.totalram - m.memInfo_.freeram;
        virtualMemUsed += m.memInfo_.totalswap - m.memInfo_.freeswap;
        virtualMemUsed *= m.memInfo_.mem_unit;

        // compute total used physical mem
        unsigned long physMemUsed = m.memInfo_.totalram - m.memInfo_.freeram;
        physMemUsed *= m.memInfo_.mem_unit;

        // print memory usage
        out << "~~~ MEMORY INFORMATION ~~~" << std::endl;
        out << "Currently Used by Process (virtual/phys) [MB] \t" << byteToMB(getProcessUsedMemory()) << " / " << byteToMB(getProcessUsedMemory(true)) << std::endl;
        out << "Currently Used (virtual/phys) [MB] \t" << byteToMB(virtualMemUsed) << " / " << byteToMB(physMemUsed) << std::endl;
        out << "Total (virtual/phys) [MB] \t" << byteToMB(m.totalVirtualMem_) << " / " << byteToMB(m.totalPhysMem_) << std::endl;
        out << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

        return out;
    }
private:
    struct sysinfo memInfo_;
    unsigned long totalVirtualMem_;
    unsigned long totalPhysMem_;
};

} // namespace MM

#endif // MEMORYMONITOR_H
