#include "MemoryMonitor.h"

MM::MemoryMonitor::MemoryMonitor()
    : memInfo_{}, totalVirtualMem_{}, totalPhysMem_{} {
    sysinfo(&memInfo_);
    totalVirtualMem_ = memInfo_.totalram;
    totalVirtualMem_ += memInfo_.totalswap;
    totalVirtualMem_ *= memInfo_.mem_unit;

    totalPhysMem_ = memInfo_.totalram;
    totalPhysMem_ *= memInfo_.mem_unit;
}

unsigned long MM::parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    unsigned long i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    if (atoi(p) < 0) { throw std::runtime_error("[ERROR] -- MM::parseLine -- integer overflow"); }
    i = static_cast<unsigned long>(atoi(p));
    return i;
}

unsigned long MM::getProcessUsedMemory(bool physical) { //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    unsigned long result = 0;
    char line[128];

    while (fgets(line, 128, file) != NULL) {
        //std::cout << "> " << line;
        if (physical) {
            if (std::strncmp(line, "VmRSS:", 6) == 0) {
                result = parseLine(line);
                break;
            }
        } else {
            if (std::strncmp(line, "VmSize:", 7) == 0) {
                result = parseLine(line);
                break;
            }
        }
    }
    fclose(file);
    return result * 1024;   // convert to byte
}
