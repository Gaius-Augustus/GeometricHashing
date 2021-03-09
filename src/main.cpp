/* Geometric Hashing -- Global, Highly Specic and Fast Filtering of Alignment Seeds
 * Copyright (C) 2020 Matthis Ebel
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/* External Libraries
 * ~~~~~~~~~~~~~~~~~~
 * - Boost 1.70.0
 * - cxx-prettyprint (https://github.com/louisdx/cxx-prettyprint) 
 * - Catch 2 (https://github.com/catchorg/Catch2)
 * - JSON for Modern C++ version 3.1.2 (https://github.com/nlohmann/json)
 * - hopscotch-map v2.2.1 (https://github.com/Tessil/hopscotch-map)
 * ~~~~~~~~~~~~~~~~~~ */

#include <climits>
#include <iostream>
#include <memory>

#include "mabl3/Timestep.h"
#include "Configuration.h"
#include "MemoryMonitor.h"
#include "SeedFinder.h"

using namespace mabl3;



//! Main
int main(int argc, char * argv[]) {
    // parse cmdline options
    auto config = std::make_shared<Configuration>(argc, argv);

    // Print info on parameters
    std::cout << std::endl << "Run Seed Finding" << std::endl
              <<              "~~~~~~~~~~~~~~~~" << std::endl << std::endl;
    std::cout << "Processing files" << std::endl;
    for (auto&& file : config->inputFiles()) { std::cout << "\t" << file << std::endl; }
    std::cout << "With " << *config << std::endl;
    std::cout << std::endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Timestep tsMain("Running Program");

    auto mm = MM::MemoryMonitor();
    std::cout << mm << std::endl << std::endl;
    std::cout << "Starting Program" << std::endl;

    // Run pipeline
    if (config->weight() <= 29) {
        SeedFinder<TwoBitKmerDataShort>(config).run();
    } else if (config->weight() <= 61) {
        SeedFinder<TwoBitKmerDataMedium>(config).run();
    } else {
        SeedFinder<TwoBitKmerDataLong>(config).run();
    }
    std::cout << "Finished Program" << std::endl;
    std::cout << mm << std::endl;
    tsMain.endAndPrint(Timestep::minutes);

    return 0;
}
