#include "Configuration.h"

Configuration::Configuration(int argc, char * argv[])
    : allowOverlap_{false},
      allvsall_{false},
      artificialSequenceSizeFactor_{1},
      batchsize_{1},
      createAllMatches_{false},
      cubeLengthCutoff_{300000000},
      cubeOutput_{0},
      cubeScoreNormalizationParameter_{300000000},
      cubeScoreParameter_{500},
      cubeScoreParameterChunks_{0},
      cubeScoreThreshold_{25.},
      diagonalDelta_{0},
      diagonalRho_{0},
      diagonalThreshold_{2.},
      dynamicArtificialSequences_{false},
      genome1_{},
      genome2_{},
      graphAnnotationFile_{},
      graphFile_{},
      hasse_{false},
      inputFiles_{},
      localSearchArea_{1000},
      maskCollection_{nullptr},
      matchLimit_{10},
      matchLimitDiscardSeeds_{false},
      maxPrefixLength_{15},
      metagraphInterface_{nullptr},
      minMatchDistance_{0},
      nThreads_{1},
      occurrencePerGenomeMax_{ULLONG_MAX},
      occurrencePerGenomeMin_{1},
      occurrencePerSequenceMax_{ULLONG_MAX},
      oldCubeScore_{false},
      optimalSeed_{false},
      output_{},
      outputArtificialSequences_{},
      outputRunInformation_{},
      performDiagonalFiltering_{false},
      performGeometricHashing_{false},
      preAddNeighbouringCubes_{false},
      preHasse_{false},
      preLinkThreshold_{5},
      preMaskCollection_{nullptr},
      preOptimalSeed_{false},
      postSequential_{false},
      redmask_{false},
      thinning_{1},
      tileSize_{0},
      verbose_{2},
      yass_{false},
      yassEpsilon_{0.05},
      yassIndel_{0.08},
      yassMutation_{0.15} {
    po::options_description generalOptions("General Program Options");
    generalOptions.add_options()
            ("allvsall", "Create all Links from all vs. all sequences, if not stated, process one reference sequence at a time (vs. all non-ref sequences)")
            ("artificial-sequence-size-factor", po::value<int>()->default_value(1), "If '--dynamic-artificial-sequences', create artificial sequences of length of this factor times the length of the input sequences")
            ("check-parameters-and-exit", "Evaluate the other command line parameters, output any warnings or errors and exit without actually doing something")
            ("dynamic-artificial-sequences", "For each real input sequence, add an artificial sequence of the same length to the respective genome.")
            ("batchsize", po::value<int>()->default_value(1), "Divide each input fasta (does not work with metagraph) into this number of  batches, run for each possible batch combination (1 for single run, default)")
            ("input,i", po::value<std::vector<std::string>>()->multitoken(), "List of input files (including first and second genome).")
            ("help,h", "Show this message and exit immediately.")
            ("genome1", po::value<std::string>(), "Filename of the first genome. Can be omitted if '--input' and exactly two input genomes.")
            ("genome2", po::value<std::string>(), "Filename of the second genome. Can be omitted if '--input' and exactly two input genomes.")
            ("graph,g", po::value<std::string>(), "Filename of a metagraph graph, assumes '.dbg' etension if no extension is given. Ignored if '--input' is given")
            ("graph-annotation", po::value<std::string>(), "Filename of the graph annotation, defaults to <graph basename>.row.annodbg. Ignored if '--input' is given")
            ("masks", po::value<std::vector<std::string>>()->multitoken(), "Directly define a set of SpacedSeedMasks of equal weight. Space separated strings can only contain `0` and `1`. Overwrites '--optimal-seed' and explicit '--weight'/'--span'.")
            ("match-limit", po::value<int>()->default_value(10), "Create at most this many (randomly chosen) matches from a seed. Corresponds to link limit in geometric hashing setting. Set to 0 for no limit.")
            ("match-limit-discard-exceeding", "If a seed would give more than '--match-limit' matches, discard all matches rather than sampling")
            ("max-prefix-length", po::value<int>()->default_value(15), "The internal seed mapping data structure allocates 8*4^(max-prefix-length) bytes of memory plus additional memory for each seed")
            ("occurrence-per-genome-max", po::value<int>()->default_value(0), "At most this many seed occurrences in any genome. Set to 0 for no threshold. Does not work as expected with --batchsize > 1.")
            ("occurrence-per-genome-min", po::value<int>()->default_value(1), "At least this many seed occurrences in a genome (if any occurrences in the respective genome). Does not work as expected with --batchsize > 1.")
            ("occurrence-per-sequence-max", po::value<int>()->default_value(0), "At most this many seed occurrences in a sequence, otherwise the respective sequence is not considered in seed creation. Set to 0 for no threshold. Does not work as expected with --fast.")
            ("optimal-seed", "Use pre-computed optimal seed of weight '--weight' instead of a randomly generated. Terminates if no such seed is found. Overwrites explicit '--span'.")
            ("output,o", po::value<std::string>()->default_value("seedFindingOutput.json"), "Output goes in this file, using json format.")
            ("output-artificial-sequences", po::value<std::string>(), "Write generated artificial sequences to this file. Does not write if not stated.")
            ("output-run-information", po::value<std::string>(), "Write JSON representation of configuration and run statistics. Does not write if not stated.")
            ("p,p", po::value<int>(), "Number of threads to create when code is executed in parallel (positive integer, default: number of CPU cores available).")
            ("pre-add-neighbouring-cubes", "For pre-filter steop (GH or M1-3). Add all neighbours to found relevant cubes")
            ("pre-link-threshold", po::value<int>()->default_value(5), "For pre-filter step (GH or M1-3). A cube needs at least this many links to pass the pre-filter and be considered in second GH run")
            ("pre-masks", po::value<std::vector<std::string>>()->multitoken(), "For pre-filter step (GH or M1-3). Directly define a set of SpacedSeedMasks of equal weight. Space separated strings can only contain `0` and `1`. Overwrites '--pre-optimal-seed' and explicit '--pre-weight'/'--pre-span'.")
            ("pre-optimal-seed", "For pre-filter (GH or M1-3). Use pre-computed optimal seed of weight '--pre-weight' instead of a randomly generated. Terminates if no such seed is found. Overwrites explicit '--pre-span'.")
            ("pre-seed-set-size", po::value<int>()->default_value(1), "For pre-filter step (GH or M1-3). Number of spaced seeds (if any) to generate. No effect if pre-span equals pre-weight (default).")
            ("pre-span", po::value<int>(), "For pre-filter step (GH or M1-3). Spaced seed length >= weight. Default: same as '--pre-weight', i.e. contiguous seeds. Overwrites '--pre-weight-fraction' if stated.")
            ("pre-weight", po::value<int>(), "For pre-filter step (GH or M1-3). Weight of spaced seed (positive integer).")
            ("pre-weight-fraction", po::value<double>()->default_value(1.), "For pre-filter step (GH or M1-3). Fraction of 'care'-positions in a seed, i.e. pre-span = ceil(pre-weight/pre-weight-fraction). No effect if '--pre-span' is given explicitly.")
            ("redmask", "Apply YASS-like redmask filter, i.e. discard low complexity seeds consisting of only one or two nucleotides.")
            ("seed-set-size", po::value<int>()->default_value(1), "Number of spaced seeds (if any) to generate. No effect if span equals weight (default).")
            ("span", po::value<int>(), "spaced seed length >= weight. Default: same as '--weight', i.e. contiguous seeds. Overwrites '--weight-fraction' if stated.")
            ("thinning", po::value<int>()->default_value(1), "Discard roughly 1/thinning of input k-mers to save memory, set to 1 for not thinning (default)")
            ("verbose", po::value<int>()->default_value(2), "Set verbosity level of terminal output: 0 - no output, 1 - static output messages (use to pipe program output into logfile), 2 - static messages and progress bars (default)")
            ("weight", po::value<int>(), "Weight of spaced seed (positive integer).")
            ("weight-fraction", po::value<double>()->default_value(1.), "Fraction of 'care'-positions in a seed, i.e. span = ceil(weight/weight-fraction). No effect if '--span' is given explicitly.");
    po::options_description diagonalOptions("Diagonal Filtering Options (M4)");
    diagonalOptions.add_options()
            ("yass", "Perform yass filtering on a seedMap before reporting matches")
            ("yass-epsilon", po::value<double>()->default_value(0.05), "Set epsilon (or alpha) parameter of yass")
            ("yass-indel", po::value<double>()->default_value(0.08), "Set indel rate parameter for yass")
            ("yass-mutation", po::value<double>()->default_value(0.15), "Set mutation rate parameter for yass")

            ("allow-overlap", "Allow overlapping matches when searching for neighbouring matches on the same diagonal. No effect if '--diagonal-threshold' is 0 (default).")
            ("diagonal-filtering", "Perform diagonal filtering on a seedMap instead of reporting plain matches")
            ("diagonal-threshold", po::value<double>()->default_value(2.), "At least this many seed-matches must be found on the same diagonal between genome1 and genome2.")
            ("local-search-area", po::value<int>()->default_value(1000), "Search in an area of this many bp for neighbouring matches on the same diagonal")
            ("min-match-distance", po::value<int>()->default_value(0), "Neighbouring matches must be at least this many bp apart. No effect if '--allow-overlap' is stated.");

    po::options_description geometricHashingOptions("Geometric Hashing Options (M5)");
    geometricHashingOptions.add_options()
            ("cube-length-cutoff", po::value<int>()->default_value(30000), "Parameter for Cube score normalization.")
            ("cube-output", po::value<int>()->default_value(0), "Write cubes to file with this verbosity level (0 - don't write, 1 - only link count and score, 2 - link count, score and all the links)")
            ("cube-score-parameter", po::value<int>()->default_value(500), "Parameter for Cube scoring.")
            ("cube-score-parameter-chunks", po::value<int>()->default_value(0), "Parameter 2 for Cube scoring, 0 for no chunking (default).")
            ("cube-score-normalization-parameter", po::value<int>()->default_value(30000), "Parameter for Cube score normalization.")
            ("cube-score-threshold", po::value<double>()->default_value(25.), "Cubes must have at least this score to be considered in match creation.")
            ("geometric-hashing", "Perform geometricHashing on a seedMap instead of reporting plain matches.")
            ("hasse", "Perform geometric hashing with 'Hasse' subcubes (three or more input genomes)")
            ("old-cube-score", "Use the old scoring algorithm (i.e. score diagonals rather than sub-tiles)")
            ("pre-hasse", "For pre-filter step in GH. Perform geometric hashing with 'Hasse' subcubes (three or more input genomes)")
            ("post-sequential", "If running with pre-filter step in GH, perform the second GH steps sequentially rather than in parallel")
            ("tilesize", po::value<int>()->default_value(0), "Cut genome in tiles of this size. Set to zero (default) to skip tiling");

    po::options_description allOptions("Program Options");
    allOptions.add(generalOptions).add(diagonalOptions).add(geometricHashingOptions);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, allOptions), vm);
    boost::program_options::notify(vm);

    // display help if no args (first arg is always program call) or -h (--help) is specified
    if (argc == 1 || vm.count("help")) {
        std::cout << allOptions << std::endl;
        exit(vm.count("help") ? 0 : 1);
    }

    validateOptions(vm);

    if (vm.count("check-parameters-and-exit")) {
        std::cout << "[INFO] -- Parameters checked, not running anything." << std::endl;
        exit(0);
    }
}



void Configuration::validateOptions(po::variables_map & vm) {
    // lambda to really check if a parameter was set by user
    auto userSet = [&vm](std::string const & key) {
        return (vm.count(key) && (!vm[key].defaulted()));
    };
    // lambda to throw if a mandatory option is not there
    auto throwMandatory = [&userSet](std::string const & key) {
        if (!userSet(key)) { throw std::runtime_error("[ERROR] -- '--" + key + "' must be specified"); }
    };
    // lambda to warn if useless parameter combinations
    auto warnUselessIfSet = [&userSet](std::string const & uselessIf, std::string const & isSet) {
        if (userSet(uselessIf) && userSet(isSet)) { std::cerr << "[WARNING] -- '--" << uselessIf << "' is useless if '--" << isSet << "' is set, ignoring '--" << uselessIf << "'!" << std::endl; }
    };
    auto warnUselessIfNotSet = [&userSet](std::string const & uselessIf, std::string const & isNotSet) {
        if (userSet(uselessIf) && !userSet(isNotSet)) { std::cerr << "[WARNING] -- '--" << uselessIf << "' is useless if '--" << isNotSet << "' is not set, ignoring '--" << uselessIf << "'!" << std::endl; }
    };
    // lambda to strip extension from a path string
    auto stripExtension = [&vm](std::string const & pathstr, bool graphInput = false) {
        if (graphInput) { return pathstr; } // metagraph stores complete paths as genome names
        auto path = fs::path(pathstr);
        path.replace_extension("");
        return (path.filename()).string();
    };
    // lambda to deal with mask stuff
    auto setMasks = [this,
                     &vm,
                     &userSet,
                     &throwMandatory,
                     &warnUselessIfSet](std::string const & masksOption,
                                        std::string const & optimalOption,
                                        std::string const & weightOption,
                                        std::string const & weightFractionOption,
                                        std::string const & spanOption,
                                        std::string const & sizeOption,
                                        bool & optimalMember,
                                        std::shared_ptr<SpacedSeedMaskCollection const> & maskMember) {
        size_t weight = 0;
        size_t span = 0;
        size_t seedSetSize = static_cast<size_t>(castWithBoundaryCheck<int, size_t>(vm, sizeOption, 1, INT_MAX));
        if (userSet(masksOption)) {
            warnUselessIfSet(weightOption, masksOption);
            warnUselessIfSet(weightFractionOption, masksOption);
            warnUselessIfSet(sizeOption, masksOption);
            warnUselessIfSet(optimalOption, masksOption);
            warnUselessIfSet(spanOption, masksOption);
            optimalMember = false;
            auto masks = vm[masksOption].as<std::vector<std::string>>();
            maskMember = std::make_shared<SpacedSeedMaskCollection const>(masks); // throws if masks are ill-formed
        } else {
            throwMandatory(weightOption);
            weight = castWithBoundaryCheck<int, size_t>(vm, weightOption, 1, INT_MAX);
            optimalMember = userSet(optimalOption);
            if (optimalMember) {
                warnUselessIfSet(weightFractionOption, optimalOption);
                warnUselessIfSet(spanOption, optimalOption);
                maskMember = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(weight),
                                                                              SpacedSeedMaskCollection::SeedSetSize(seedSetSize));  // throws if not possible to load
            } else {
                auto weightFraction = castWithBoundaryCheck<double, double>(vm, weightFractionOption, 0, 1);
                if (userSet(spanOption)) {
                    warnUselessIfSet(weightFractionOption, spanOption);
                    span = castWithBoundaryCheck<int, size_t>(vm, spanOption, 1, INT_MAX);
                } else if (userSet(weightFractionOption)) {
                    span = static_cast<size_t>(std::ceil(static_cast<double>(weight)/weightFraction));
                } else {
                    span = weight;
                }
                maskMember = std::make_shared<SpacedSeedMaskCollection const>(SpacedSeedMaskCollection::Weight(weight),
                                                                              SpacedSeedMaskCollection::Span(span),
                                                                              SpacedSeedMaskCollection::SeedSetSize(seedSetSize));
            }
        }
    };



    /* set all members according to the program options */
    // general options
    // --allvsall
    allvsall_ = userSet("allvsall");
    // --artificial-sequence-size-factor
    if (userSet("artificial-sequence-size-factor")) { warnUselessIfNotSet("artificial-sequence-size-factor", "dynamic-artificial-sequences"); }
    artificialSequenceSizeFactor_ = castWithBoundaryCheck<int, size_t>(vm, "artificial-sequence-size-factor", 1, INT_MAX);
    // --dynamic-artificial-sequences
    dynamicArtificialSequences_ = userSet("dynamic-artificial-sequences");
    // batchsize
    batchsize_ = castWithBoundaryCheck<int, size_t>(vm, "batchsize", 1, INT_MAX);
    // --input or --graph
    //throwMandatory("input");
    if (userSet("input")) {
        warnUselessIfSet("graph", "input");
        inputFiles_ = vm["input"].as<std::vector<std::string>>();
        auto numberOfInputGenomes = inputFiles_.size();
        if (numberOfInputGenomes < 2) { throw std::runtime_error("[ERROR] -- '--input' or '--graph' need at least two genomes"); }
        // --genome1/2
        if (numberOfInputGenomes > 2) {
            throwMandatory("genome1");
            throwMandatory("genome2");
        }
    } else if (userSet("graph")) {
        warnUselessIfSet("input", "graph");
        throwMandatory("genome1");
        throwMandatory("genome2");
        graphFile_ = fs::path(vm["graph"].as<std::string>());
        if (!graphFile_.has_extension()) { graphFile_.replace_extension(".dbg"); }
        if (userSet("graph-annotation")) {
            graphAnnotationFile_ = fs::path(vm["graph-annotation"].as<std::string>());
            if (!graphAnnotationFile_.has_extension()) { graphAnnotationFile_.replace_extension(".row.annodbg"); }
        } else {
            graphAnnotationFile_ = graphFile_;
            graphAnnotationFile_.replace_extension(".row.annodbg");
        }
        // may take very long, skip if just checking parameteres
        if (!userSet("check-parameters-and-exit")) {
            Timestep ts("Loading graph");
            metagraphInterface_ = std::make_shared<MetagraphInterface const>(graphFile_, graphAnnotationFile_);
            ts.endAndPrint();
        }
    } else {
        throw std::runtime_error("[ERROR] -- '--input' or '--graph' must be specified");
    }
    // if --input, possible that no genome1/2; if --graph, genome1/2 mandatory so else never executed
    if (userSet("genome1")) {
        throwMandatory("genome2");
        genome1_ = stripExtension(vm["genome1"].as<std::string>(), userSet("graph"));
    } else {
        genome1_ = stripExtension(inputFiles_.at(0));
    }
    if (userSet("genome2")) {
        throwMandatory("genome1");
        genome2_ = stripExtension(vm["genome2"].as<std::string>(), userSet("graph"));
    } else {
        genome2_ = stripExtension(inputFiles_.at(1));
    }
    // --match-limit
    matchLimit_ = castWithBoundaryCheck<int, size_t>(vm, "match-limit", 0, INT_MAX);
    matchLimit_ = (matchLimit_ == 0) ? ULLONG_MAX : matchLimit_;
    // --match-limit-discard-exceeding
    warnUselessIfNotSet("match-limit-discard-exceeding", "match-limit");
    matchLimitDiscardSeeds_ = userSet("match-limit-discard-exceeding");
    // --max-prefix-length
    maxPrefixLength_ = castWithBoundaryCheck<int, size_t>(vm, "max-prefix-length", 1, INT_MAX);
    // --occurrence-per-genome-max
    occurrencePerGenomeMax_ = castWithBoundaryCheck<int, size_t>(vm, "occurrence-per-genome-max", 0, INT_MAX);
    occurrencePerGenomeMax_ = (occurrencePerGenomeMax_ == 0) ? ULLONG_MAX : occurrencePerGenomeMax_;
    // --occurrence-per-genome-min
    occurrencePerGenomeMin_ = castWithBoundaryCheck<int, size_t>(vm, "occurrence-per-genome-min", 1, INT_MAX);
    // --occurrence-per-sequence-max
    occurrencePerSequenceMax_ = castWithBoundaryCheck<int, size_t>(vm, "occurrence-per-sequence-max", 0, INT_MAX);
    occurrencePerSequenceMax_ = (occurrencePerSequenceMax_ == 0) ? ULLONG_MAX : occurrencePerSequenceMax_;
    // --output
    output_ = fs::path(vm["output"].as<std::string>());
    if (fs::exists(output_)) {
        std::cerr << "[WARNING] -- File '" << output_.string() << "' already exists and gets overwritten" << std::endl;
    } else if (fs::exists(output_) && !fs::is_regular_file(output_)) {
        throw std::runtime_error("[ERROR] -- Outfile '" + output_.string() + "' is not a regular file");
    } else {
        std::ofstream os(output_);
        if (!os.good()) { throw std::runtime_error("[ERROR] -- Cannot write to outfile '" + output_.string() + "'"); }
    }
    // --output-artificial-sequences
    if (userSet("output-artificial-sequences")) {
        outputArtificialSequences_ = fs::path(vm["output-artificial-sequences"].as<std::string>());
        if (fs::exists(outputArtificialSequences_)) {
            std::cerr << "[WARNING] -- File '" << outputArtificialSequences_.string() << "' already exists and gets overwritten" << std::endl;
        } else if (fs::exists(outputArtificialSequences_) && !fs::is_regular_file(outputArtificialSequences_)) {
            throw std::runtime_error("[ERROR] -- Artificial Sequences Outfile '" + outputArtificialSequences_.string() + "' is not a regular file");
        } else {
            std::ofstream os(outputArtificialSequences_);
            if (!os.good()) { throw std::runtime_error("[ERROR] -- Cannot write to artificial sequences output file '" + outputArtificialSequences_.string() + "'"); }
        }
    }
    // --output-run-information
    if (userSet("output-run-information")) {
        outputRunInformation_ = fs::path(vm["output-run-information"].as<std::string>());
        if (fs::exists(outputRunInformation_)) {
            std::cerr << "[WARNING] -- File '" << outputRunInformation_.string() << "' already exists and gets overwritten" << std::endl;
        } else if (fs::exists(outputRunInformation_) && !fs::is_regular_file(outputRunInformation_)) {
            throw std::runtime_error("[ERROR] -- Run Information Outfile '" + outputRunInformation_.string() + "' is not a regular file");
        } else {
            std::ofstream os(outputRunInformation_);
            if (!os.good()) { throw std::runtime_error("[ERROR] -- Cannot write to run information output file '" + outputRunInformation_.string() + "'"); }
        }
    }
    // --p
    nThreads_ = userSet("p")
            ? castWithBoundaryCheck<int, size_t>(vm, "p", 1, INT_MAX)
            : std::thread::hardware_concurrency();
    if (nThreads_ == 0) {
        std::cerr << "[WARNING] -- Could not determine number of threads automatically, using only one thread" << std::endl;
        nThreads_ = 1;
    }
    // --pre-add-neghbouring-cubes
    preAddNeighbouringCubes_ = userSet("pre-add-neighbouring-cubes");
    // --pre-link-threshold
    preLinkThreshold_ = castWithBoundaryCheck<int, size_t>(vm, "pre-link-threshold", 0, INT_MAX);
    // pre-masks
    if (userSet("pre-masks") || userSet("pre-weight")) {
        setMasks("pre-masks",
                 "pre-optimal-seed",
                 "pre-weight",
                 "pre-weight-fraction",
                 "pre-span",
                 "pre-seed-set-size",
                 preOptimalSeed_, preMaskCollection_);
        if (metagraphInterface_ != nullptr && metagraphInterface_->getK() < preMaskCollection_->maxSpan()) {
            throw std::runtime_error("[ERROR] -- Max. required pre-span (" + std::to_string(preMaskCollection_->maxSpan())
                                     + ") exceeds graph's k (" + std::to_string(metagraphInterface_->getK())
                                     + "). Build graph with higher k or choose lower pre-weight/span");
        }
    }
    // --redmask
    redmask_ = userSet("redmask");
    // masks
    setMasks("masks",
             "optimal-seed",
             "weight",
             "weight-fraction",
             "span",
             "seed-set-size",
             optimalSeed_, maskCollection_);
    if (metagraphInterface_ != nullptr && metagraphInterface_->getK() < maskCollection_->maxSpan()) {
        throw std::runtime_error("[ERROR] -- Max. required span (" + std::to_string(maskCollection_->maxSpan())
                                 + ") exceeds graph's k (" + std::to_string(metagraphInterface_->getK())
                                 + "). Build graph with higher k or choose lower weight/span");
    }
    // --verbose
    verbose_ = castWithBoundaryCheck<int, size_t>(vm, "verbose", 0, 2);

    // diagonal filter options
    // --diagonal-filtering
    performDiagonalFiltering_ = userSet("yass") || userSet("diagonal-filtering");
    yass_ = userSet("yass");
    if (performDiagonalFiltering_ && preMaskCollection_) {
        throw std::runtime_error("[ERROR] -- cannot run pre-filter step with M4!");
    }
    // --yass-epsilon
    warnUselessIfNotSet("yass-epsilon", "yass");
    yassEpsilon_ = castWithBoundaryCheck<double, double>(vm, "yass-epsilon", 0, 1);
    // --yass-indel
    warnUselessIfNotSet("yass-indel", "yass");
    yassIndel_ = castWithBoundaryCheck<double, double>(vm, "yass-indel", 0, 1);
    // --yass-mutation
    warnUselessIfNotSet("yass-mutation", "yass");
    yassMutation_ = castWithBoundaryCheck<double, double>(vm, "yass-mutation", 0, 1);

    diagonalRho_ = computeRho(yassMutation_, maskCollection_->weight(), yassEpsilon_);
    diagonalDelta_ = computeDelta(yassIndel_, diagonalRho_, yassEpsilon_);

    // --diagonal-threshold
    warnUselessIfNotSet("diagonal-threshold", "diagonal-filtering");
    diagonalThreshold_ = castWithBoundaryCheck<double, double>(vm, "diagonal-threshold", 1, DBL_MAX);
    // --local-search-area
    warnUselessIfNotSet("local-search-area", "diagonal-filtering");
    localSearchArea_ = castWithBoundaryCheck<int, size_t>(vm, "local-search-area", 1, INT_MAX);
    // --allow-overlap
    warnUselessIfNotSet("allow-overlap", "diagonal-filtering");
    allowOverlap_ = userSet("allow-overlap");
    if (allowOverlap_ && !performDiagonalFiltering_) { std::cerr << "[WARNING] -- '--allow-overlap' is ignored as diagonal filter is not applied" << std::endl; }
    // --min-match-distance
    warnUselessIfNotSet("min-match-distance", "diagonal-filtering");
    warnUselessIfSet("min-match-distance", "allow-overlap");
    minMatchDistance_ = userSet("allow-overlap")
            ? 0
            : castWithBoundaryCheck<int, size_t>(vm, "min-match-distance", 0, INT_MAX);

    // geometric hashing options
    // --geometric-hashing
    performGeometricHashing_ = userSet("geometric-hashing");
    // --cube-length-cutoff
    warnUselessIfNotSet("cube-length-cutoff", "geometric-hashing");
    cubeLengthCutoff_ = castWithBoundaryCheck<int, size_t>(vm, "cube-length-cutoff", 0, INT_MAX);
    // --cube-output
    cubeOutput_ = castWithBoundaryCheck<int, size_t>(vm, "cube-output", 0, 2);
    // --cube-score-parameter
    warnUselessIfNotSet("cube-score-parameter", "geometric-hashing");
    cubeScoreParameter_ = castWithBoundaryCheck<int, size_t>(vm, "cube-score-parameter", 1, INT_MAX);
    // --cube-score-parameter-chunks
    warnUselessIfNotSet("cube-score-parameter-chunks", "geometric-hashing");
    cubeScoreParameterChunks_ = castWithBoundaryCheck<int, size_t>(vm, "cube-score-parameter-chunks", 0, INT_MAX);
    // --cube-score-normalization-parameter
    warnUselessIfNotSet("cube-score-normalization-parameter", "geometric-hashing");
    cubeScoreNormalizationParameter_ = castWithBoundaryCheck<int, size_t>(vm, "cube-score-normalization-parameter", 1, INT_MAX);
    // --cube-score-threshold
    warnUselessIfNotSet("cube-score-threshold", "geometric-hashing");
    cubeScoreThreshold_ = castWithBoundaryCheck<double, double>(vm, "cube-score-threshold", 0, DBL_MAX);
    warnUselessIfNotSet("hasse", "geometric-hashing");
    hasse_ = userSet("hasse");
    // --old-cube-score
    warnUselessIfNotSet("old-cube-score", "geometric-hashing");
    oldCubeScore_ = userSet("old-cube-score");
    // --pre-hasse
    preHasse_ = userSet("pre-hasse");
    // --post-sequential
    postSequential_ = userSet("post-sequential");
    // --thinning
    thinning_ = castWithBoundaryCheck<int, size_t>(vm, "thinning", 1, INT_MAX);
    // --tilesize
    warnUselessIfNotSet("tilesize", "geometric-hashing");
    tileSize_ = castWithBoundaryCheck<int, size_t>(vm, "tilesize", 0, INT_MAX);

    if ((!performGeometricHashing_) && preMaskCollection_ && tileSize_ > 0) {
        std::cerr << "[WARNING] -- Setting tilesize to zero for pre-filter step as '--geometric-hashing' is not executed!";
        tileSize_ = 0;
    }
}



JsonValue Configuration::configJson() const {
    auto jsonstream = std::stringstream(std::ios_base::out);    // stream json of config into string
    JsonStreamDict map(jsonstream);
    map.addValue("allowOverlap", allowOverlap_);
    map.addValue("allvsall", allvsall_);
    map.addValue("artificialSequenceSizeFactor", artificialSequenceSizeFactor_);
    map.addValue("batchsize", batchsize_);
    map.addValue("createAllMatches", createAllMatches_);
    map.addValue("cubeLengthCutoff", cubeLengthCutoff_);
    map.addValue("cubeOutput", cubeOutput_);
    map.addValue("cubeScoreNormalizationParameter", cubeScoreNormalizationParameter_);
    map.addValue("cubeScoreParameter", cubeScoreParameter_);
    map.addValue("cubeScoreParameterChunks", cubeScoreParameterChunks_);
    map.addValue("cubeScoreThreshold", cubeScoreThreshold_);
    map.addValue("diagonalDelta", diagonalDelta_);
    map.addValue("diagonalRho", diagonalRho_);
    map.addValue("diagonalThreshold", diagonalThreshold_);
    map.addValue("dynamicArtificialSequences", dynamicArtificialSequences_);
    map.addValue("genome1", genome1_);
    map.addValue("genome2", genome2_);
    map.addValue("graph", graphFile_.string());
    map.addValue("graphAnnot", graphAnnotationFile_.string());
    map.addValue("hasse", hasse_);
    map.addValue("inputFiles", inputFiles_);
    map.addValue("localSearchArea", localSearchArea_);
    map.addValue("masks", masks());
    map.addValue("matchLimit", matchLimit_);
    map.addValue("matchLimitDiscardExceeding", matchLimitDiscardSeeds_);
    map.addValue("maxPrefixLength", maxPrefixLength_);
    map.addValue("minMatchDistance", minMatchDistance_);
    map.addValue("nThreads", nThreads_);
    map.addValue("occurrencePerGenomeMax", occurrencePerGenomeMax_);
    map.addValue("occurrencePerGenomeMin", occurrencePerGenomeMin_);
    map.addValue("occurrencePerSequenceMax", occurrencePerSequenceMax_);
    map.addValue("oldCubeScore", oldCubeScore_);
    map.addValue("optimalSeed", optimalSeed_);
    map.addValue("performDiagonalFiltering", performDiagonalFiltering_);
    map.addValue("performGeometricHashing", performGeometricHashing_);
    map.addValue("preAddNeighbouringCubes", preAddNeighbouringCubes_);
    map.addValue("pre-hasse", preHasse_);
    map.addValue("pre-link-threshold", preLinkThreshold_);
    if (preMaskCollection_) {
        map.addValue("pre-masks", preMaskCollection_->masksAsString());
    } else {
        map.addValue("pre-masks", "");
    }
    map.addValue("pre-optimalSeed", preOptimalSeed_);
    if (preMaskCollection_) {
        map.addValue("pre-seedSetSize", preMaskCollection_->size());
        map.addValue("pre-span", preMaskCollection_->maxSpan());
        map.addValue("pre-weight", preMaskCollection_->weight());
    } else {
        map.addValue("pre-seedSetSize", 0);
        map.addValue("pre-span", 0);
        map.addValue("pre-weight", 0);
    }
    map.addValue("post-sequential", postSequential_);
    map.addValue("redmask", redmask_);
    map.addValue("seedSetSize", seedSetSize());
    map.addValue("span", span());
    map.addValue("thinning", thinning_);
    map.addValue("tileSize", tileSize_);
    map.addValue("verbose", verbose_);
    map.addValue("weight", weight());
    map.addValue("yass", yass_);
    map.addValue("yassEpsilon", yassEpsilon_);
    map.addValue("yassIndel", yassIndel_);
    map.addValue("yassMutation", yassMutation_);
    map.close();
    return JsonValue(jsonstream.str(), JsonValue::stringIsValidJson);
}



std::ostream & operator<<(std::ostream & os, Configuration const & conf) {
    os << "Configuration:" << std::endl;
    os << "\t" << "--allow-overlap " << conf.allowOverlap_ << std::endl;
    os << "\t" << "--allvsall " << conf.allvsall_ << std::endl;
    os << "\t" << "--artificial-sequence-size-factor " << conf.artificialSequenceSizeFactor_ << std::endl;
    os << "\t" << "--batchsize " << conf.batchsize_ << std::endl;
    os << "\t" << "createAllMatches_ " << conf.createAllMatches_ << std::endl;
    os << "\t" << "--cube-length-cutoff " << conf.cubeLengthCutoff_ << std::endl;
    os << "\t" << "--cube-output " << conf.cubeOutput_ << std::endl;
    os << "\t" << "--cube-score-parameter " << conf.cubeScoreParameter_ << std::endl;
    os << "\t" << "--cube-score-parameter-chunks " << conf.cubeScoreParameterChunks_ << std::endl;
    os << "\t" << "--cube-score-normalization-parameter " << conf.cubeScoreNormalizationParameter_ << std::endl;
    os << "\t" << "--cube-score-threshold " << conf.cubeScoreThreshold_ << std::endl;
    os << "\t" << "diagonalDelta " << conf.diagonalDelta_ << std::endl;
    os << "\t" << "diagonalRho " << conf.diagonalRho_ << std::endl;
    os << "\t" << "--diagonal-threshold " << conf.diagonalThreshold_ << std::endl;
    os << "\t" << "--dynamic-artificial-sequences " << conf.dynamicArtificialSequences_ << std::endl;
    os << "\t" << "--input " << conf.inputFiles_ << std::endl;
    os << "\t" << "--genome1 " << conf.genome1_ << std::endl;
    os << "\t" << "--genome2 " << conf.genome2_ << std::endl;
    os << "\t" << "--graph " << conf.graphFile_ << std::endl;
    os << "\t" << "--graph-annotation " << conf.graphAnnotationFile_ << std::endl;
    os << "\t" << "--hasse " << conf.hasse_ << std::endl;
    os << "\t" << "--local-search-area " << conf.localSearchArea_ << std::endl;
    os << "\t" << "--masks " << conf.masks() << std::endl;
    os << "\t" << "--match-limit " << conf.matchLimit_ << std::endl;
    os << "\t" << "--match-limit-discard-exceeding " << conf.matchLimitDiscardSeeds_ << std::endl;
    os << "\t" << "--max-prefix-length " << conf.maxPrefixLength_ << std::endl;
    os << "\t" << "--min-match-distance " << conf.minMatchDistance_ << std::endl;
    os << "\t" << "--occurrence-per-genome-max " << conf.occurrencePerGenomeMax_ << std::endl;
    os << "\t" << "--occurrence-per-genome-min " << conf.occurrencePerGenomeMin_ << std::endl;
    os << "\t" << "--occurrence-per-sequence-max " << conf.occurrencePerSequenceMax_ << std::endl;
    os << "\t" << "--old-cube-score " << conf.oldCubeScore_ << std::endl;
    os << "\t" << "--optimal-seed " << conf.optimalSeed_ << std::endl;
    os << "\t" << "--output " << conf.output_ << std::endl;
    os << "\t" << "--output-artificial-sequences " << conf.outputArtificialSequences_ << std::endl;
    os << "\t" << "--p " << conf.nThreads_ << std::endl;
    os << "\t" << "perform diagonal filtering " << conf.performDiagonalFiltering_ << std::endl;
    os << "\t" << "perform geometric-hashing " << conf.performGeometricHashing_ << std::endl;
    os << "\t" << "--pre-add-neighbouring-cubes " << conf.preAddNeighbouringCubes_ << std::endl;
    os << "\t" << "--pre-hasse " << conf.preHasse_ << std::endl;
    os << "\t" << "--pre-link-threshold " << conf.preLinkThreshold_ << std::endl;
    os << "\t" << "--pre-masks ";
    if (conf.preMaskCollection()) { os << conf.preMaskCollection()->masksAsString() << std::endl; } else { os << "" << std::endl; }
    os << "\t" << "--pre-optimal-seed " << conf.preOptimalSeed() << std::endl;
    os << "\t" << "--pre-seed-set-size ";
    if (conf.preMaskCollection()) { os << conf.preMaskCollection()->size() << std::endl; } else { os << "" << std::endl; }
    os << "\t" << "--pre-span ";
    if (conf.preMaskCollection()) { os << conf.preMaskCollection()->maxSpan() << std::endl; } else { os << "" << std::endl; }
    os << "\t" << "--pre-weight ";
    if (conf.preMaskCollection()) { os << conf.preMaskCollection()->weight() << std::endl; } else { os << "" << std::endl; }
    os << "\t" << "--post-sequential " << conf.postSequential_ << std::endl;
    os << "\t" << "--redmask " << conf.redmask_ << std::endl;
    os << "\t" << "--seed-set-size " << conf.seedSetSize() << std::endl;
    os << "\t" << "--span " << conf.span() << std::endl;
    os << "\t" << "--thinning " << conf.thinning_ << std::endl;
    os << "\t" << "--tilesize " << conf.tileSize_ << std::endl;
    os << "\t" << "--verbose " << conf.verbose_ << std::endl;
    os << "\t" << "--weight " << conf.weight() << std::endl;
    os << "\t" << "--yass " << conf.yass_ << std::endl;
    os << "\t" << "--yass-epsilon " << conf.yassEpsilon_ << std::endl;
    os << "\t" << "--yass-indel " << conf.yassIndel_ << std::endl;
    os << "\t" << "--yass-mutation " << conf.yassMutation_ << std::endl;
    return os;
}
