project(seedFindingTests)

# small tool for testdata creation (get string hashes)
add_executable(stringHasher stringHash.cpp)

# the main program
add_executable(test_seedFinding EXCLUDE_FROM_ALL testMain.cpp
                                                 ConfigurationGenerator.h
                                                 fillDiagonalMatchesFilterTestdata.h
                                                 fillLinksetTestdata.h
                                                 GeometricHashingTestdata.h
                                                 LoadTestdata.cpp LoadTestdata.h
                                                 testContainerChunks.cpp
                                                 testConfiguration.cpp
                                                 testCubeRelated.cpp
                                                 testDiagonalMatchesFilter.cpp
                                                 testFastaRepresentation.cpp
                                                 testIdentifierMapping.cpp
                                                 testKmerOccurrence.cpp
                                                 testLink.cpp
                                                 testLinkset.cpp
                                                 testSeedFinder.cpp
                                                 testSeedKmerMap.cpp
                                                 testSeedMap.cpp
                                                 testSpacedSeedMask.cpp
                                                 testTwoBitKmer.cpp
                                                 WriteTestdata.h)

target_include_directories(test_seedFinding PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# set testfile location (sets -DSEEDFINDING_TESTFILEPATH=\"/home/.../testfile\" that is then used in LoadTestdata.h)
target_compile_definitions(test_seedFinding PRIVATE SEEDFINDING_TESTFILEPATH="${SEEDFINDING_TESTFILELOCATION}")

# link sources
target_link_libraries(test_seedFinding PRIVATE seedFindingLib)

# link catch
target_link_libraries(test_seedFinding PRIVATE catch2)
