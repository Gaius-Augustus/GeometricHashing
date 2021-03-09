project(seedFindingProgram)

# sources as library to make them testable
add_library(seedFindingLib STATIC Configuration.cpp Configuration.h computeYassParameters.h
                                  SeedFinder.h
                                  SeedMap.h
                                  ExtractSeeds.h
                                  Linkset.cpp Linkset.h Link.h
                                  Cubeset.cpp Cubeset.h Cube.h
                                  DiagonalMatchesFilter.cpp DiagonalMatchesFilter.h
                                  FastaRepresentation.cpp FastaRepresentation.h FastaCollection.h
                                  MemoryMonitor.cpp MemoryMonitor.h
                                  ContainerChunks.h
                                  CustomHashGeneral.h
                                  IdentifierMapping.h
                                  KmerOccurrence.h
                                  optimalSpacedSeeds.h
                                  ReverseComplement.h
                                  SpacedSeedMask.h SpacedSeedMaskCollection.h
                                  StrongType.h
                                  TwoBitKmer.h
                                  regionTupleExtraction.h)
target_include_directories(seedFindingLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# set C++ standard
target_compile_features(seedFindingLib PUBLIC cxx_std_17)

# link third party libraries
target_include_directories(seedFindingLib SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(seedFindingLib PUBLIC cpp-json-outstream)
target_link_libraries(seedFindingLib PUBLIC cpp-progress-bar)
target_link_libraries(seedFindingLib PUBLIC cpp-timestep)
target_link_libraries(seedFindingLib PUBLIC cxx-prettyprint)
target_link_libraries(seedFindingLib PUBLIC hopscotch-map)
target_link_libraries(seedFindingLib PUBLIC json)
target_link_libraries(seedFindingLib PRIVATE pthread)
target_link_libraries(seedFindingLib PRIVATE tbb)
target_link_libraries(seedFindingLib PRIVATE Boost::program_options)
target_link_libraries(seedFindingLib PRIVATE stdc++fs)

# the main program
add_executable(seedFinding main.cpp)

# link sources
target_link_libraries(seedFinding PRIVATE seedFindingLib)
