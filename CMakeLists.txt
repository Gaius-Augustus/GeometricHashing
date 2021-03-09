cmake_minimum_required(VERSION 3.6.3)
project(seedFinding VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # generate compile_commands.json for e.g. sourcetrail

# Print info about compiler
message("[INFO] -- CXX Compiler version: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message("[INFO] -- BUILD MODE: ${CMAKE_BUILD_TYPE}")

# set verbose make output and extra warnings only for debug builds
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
    if(NOT MSVC)
        add_compile_options(-pedantic)
    endif()
endif()

# set warnings for everything
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Weffc++ -Woverloaded-virtual
        -Wuninitialized -Wmaybe-uninitialized)# -Werror)
endif()

# tell cmake where to put the compiled binaries (usually in build/ dir)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# for telling the unit tests where the testfiles are
set(SEEDFINDING_TESTFILELOCATION "${CMAKE_CURRENT_SOURCE_DIR}/testdata")
message("[INFO] -- SEEDFINDING_TESTFILELOCATION set to: ${SEEDFINDING_TESTFILELOCATION}")

# find libraries
# need at least boost 1.70.0, need the program_options
set(BOOST_ROOT "lib/boost/build" CACHE PATH "Preferred path to boost library")
message("[INFO] -- BOOST_ROOT set to: ${BOOST_ROOT}")

find_package(Boost
             1.70.0
             REQUIRED
             COMPONENTS program_options)

# create interface library targets from other libs
add_library(catch2 INTERFACE)
target_include_directories(catch2 SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/Catch2/single_include")
add_library(cpp-json-outstream INTERFACE)
target_include_directories(cpp-json-outstream SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-json-outstream/include")
add_library(cpp-progress-bar INTERFACE)
target_include_directories(cpp-progress-bar SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-progress-bar/include")
add_library(cpp-timestep INTERFACE)
target_include_directories(cpp-timestep SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-timestep/include")
add_library(cxx-prettyprint INTERFACE)
target_include_directories(cxx-prettyprint SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cxx-prettyprint")
add_library(hopscotch-map INTERFACE)
target_include_directories(hopscotch-map SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/hopscotch-map/include")
add_library(json INTERFACE)
target_include_directories(json SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/json/single_include")



# make targets in subdirs visible
add_subdirectory(lib/metagraphInterface/include/mabl3)
add_subdirectory(src)
add_subdirectory(src/test)
