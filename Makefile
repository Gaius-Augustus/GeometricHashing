VERSION = 0.1

# Allowed values: 0, 1, 2, 3
OPTIMIZATION = 0

ifdef CXX
CC      = $(CXX)
else
CC      = g++ -fPIC
endif

rootdir=$(realpath .)

# -isystem instead of -I suppresses compiler warnings for libraries
CFLAGS  = -Wall -Wextra -Weffc++ -Woverloaded-virtual -Wuninitialized -Wmaybe-uninitialized \
          -std=c++17 -g -O$(OPTIMIZATION) -DVERSION=\"$(VERSION)\" \
          -isystem ./lib/ \
          -isystem ./lib/boost/build/include/ \
          -isystem ./lib/hopscotch-map/ \
          -I ./ \
          -I ./src/ \
          -I ./src/test/

LDFLAGS = -L$(rootdir)/lib/boost/build/lib \
          -Wl,-rpath=$(rootdir)/lib/boost/build/lib \
          -lpthread -ltbb \
          -lboost_program_options -lstdc++fs

#TESTFLAGS = --param=max-vartrack-size=0
TESTFLAGS =



SRCdir = src
OBJdir = obj
DEPdir = obj

# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPdir)/$*.Td
POSTCOMPILE = @mv -f $(DEPdir)/$*.Td $(DEPdir)/$*.d && touch $@

FILES = main Configuration SeedMapGeneral \
        SeedMapContiguous SeedMapSpaced \
        ExtractSeeds \
        Linkset Cubeset DiagonalMatchesFilter \
        FastaRepresentation MemoryMonitor
TFILES = Configuration SeedMapGeneral \
         SeedMapContiguous SeedMapSpaced \
         ExtractSeeds \
         Linkset Cubeset DiagonalMatchesFilter \
         FastaRepresentation MemoryMonitor \
         test/testMain test/testSeedMapGeneral test/testSeedMapContiguous \
         test/testSeedMapSpaced test/testDiagonalMatchesFilter \
         test/testContainerChunks test/testDisambiguateSTLContainer \
         test/testFastaRepresentation test/testJsonStream \
         test/testKmerOccurrence test/testTwoBitKmer test/testIdentifierMapping test/testSpacedSeedMask \
         test/testLinkRelated test/testCubeRelated test/testConfiguration

OBJ = $(patsubst %, $(OBJdir)/%.o,  $(FILES))
SRC = $(patsubst %, $(SRCdir)/%.cpp, $(FILES))
DEP = $(patsubst %, $(DEPdir)/%.d, $(FILES))

TOBJ = $(patsubst %, $(OBJdir)/%.o,  $(TFILES))
TSRC = $(patsubst %, $(SRCdir)/%.cpp, $(TFILES))



# first rule (i.e. only running make will do this -> building the program)
geometricHashing: $(OBJ)
	mkdir -p bin
	mkdir -p lib/boost
	$(CC) $(CFLAGS) -o bin/geometricHashing $(OBJ) $(LDFLAGS)



# make test will build the unit test program
test: $(TOBJ)
	mkdir -p bin
	mkdir -p lib/boost
	$(CC) $(CFLAGS) $(TESTFLAGS) -o bin/test_geometricHashing $(TOBJ) $(LDFLAGS)



# rule for any object file including auto-dependency creation (see link above for details)
# $@ is the target name (e.g. main.o), $< is the first prerequisite (main.cpp), $^ would be all prerequisites (does not apply here)
$(OBJdir)/%.o: $(SRCdir)/%.cpp
	mkdir -p $(OBJdir) $(DEPdir) $(OBJdir)/test $(DEPdir)/test
	$(CC) $(DEPFLAGS) $(CFLAGS) $(TESTFLAGS) -c -o $@ $<
	$(POSTCOMPILE)



# auto-dependency stuff
(DEPdir)/%.d: ;
.PRECIOUS: $(DEPdir)/%.d

include $(wildcard $(patsubst %, $(DEPdir)/%.d, $(FILES)))



# make clean will remove all generated object and dependency files and the testlib
clean:
	rm -f $(OBJdir)/*.o $(DEPdir)/*.d $(OBJdir)/test/*.o $(DEPdir)/test/*.d
