#include "Configuration.h"

Configuration::Configuration(int argc, char * argv[])
    : allowOverlap_{false},
      artificialSequenceSizeFactor_{1},
      createAllMatches_{false},
      cubeScoreMu_{0.},
      cubeScoreThreshold_{0},
      diagonalThreshold_{0.},
      dynamicArtificialSequences_{false},
      genome1_{},
      genome2_{},
      inputFiles_{},
      localSearchArea_{0},
      masks_{},
      matchLimit_{ULLONG_MAX},
      matchLimitDiscardSeeds_{false},
      minMatchDistance_{0},
      noProgressbar_{false},
      nThreads_{1},
      occurrencePerGenomeMax_{ULLONG_MAX},
      occurrencePerGenomeMin_{1},
      optimalSeed_{false},
      output_{},
      outputArtificialSequences_{},
      outputRunInformation_{},
      performDiagonalFiltering_{false},
      performGeometricHashing_{false},
      seedSetSize_{1},
      span_{0},
      tileSize_{1},
      weight_{0} {
    po::options_description generalOptions("General Program Options");
    generalOptions.add_options()
            ("artificial-sequence-size-factor", po::value<int>()->default_value(1), "If '--dynamic-artificial-sequences', create artificial sequences of length of this factor times the length of the input sequences")
            ("check-parameters-and-exit", "Evaluate the other command line parameters, output any warnings or errors and exit without actually doing something")
            ("dynamic-artificial-sequences", "For each real input sequence, add an artificial sequence of the same length to the respective genome.")
            ("input,i", po::value<std::vector<std::string>>()->multitoken(), "List of input files (including first and second genome).")
            ("help,h", "Show this message and exit immediately.")
            ("genome1", po::value<std::string>(), "Filename of the first genome. Can be omitted if exactly two input genomes.")
            ("genome2", po::value<std::string>(), "Filename of the second genome. Can be omitted if exactly two input genomes.")
            ("masks", po::value<std::vector<std::string>>()->multitoken(), "Directly define a set of SpacedSeedMasks of equal weight. Space separated strings can only contain `0` and `1`. Overwrites '--optimal-seed' and explicit '--weight'/'--span'.")
            ("match-limit", po::value<int>()->default_value(0), "Create at most this many (randomly chosen) matches from a seed. Corresponds to link limit in geometric hashing setting. Set to 0 for no limit.")
            ("match-limit-discard-exceeding", "If a seed would give more than '--match-limit' matches, discard all matches rather than sampling")
            ("optimal-seed", "Use pre-computed optimal seed of weight '--l' instead of a randomly generated. Terminates if no such seed is found. Overwrites explicit '--span'.")
            ("output,o", po::value<std::string>()->default_value("seedFindingOutput.json"), "Output goes in this file, using json format.")
            ("output-artificial-sequences", po::value<std::string>(), "Write generated artificial sequences to this file. Does not write if not stated.")
            ("output-run-information", po::value<std::string>(), "Write JSON representation of configuration and run statistics. Does not write if not stated.")
            ("p,p", po::value<int>(), "Number of threads to create when code is executed in parallel (positive integer, default: number of CPU cores available).")
            ("quiet,q", "Run more quietly, i.e. don't create progress bars.")
            ("seed-set-size", po::value<int>()->default_value(1), "Number of spaced seeds (if any) to generate. No effect if span equals weight (default).")
            ("span", po::value<int>(), "spaced seed length >= weight. Default: same as '--weight', i.e. contiguous seeds. Overwrites '--weight-fraction' if stated.")
            ("weight", po::value<int>(), "Weight of spaced seed (positive integer).")
            ("weight-fraction", po::value<double>()->default_value(1), "Fraction of 'care'-positions in a seed, i.e. span = ceil(weight/weight-fraction). No effect if '--span' is given explicitly.");
    po::options_description diagonalOptions("Diagonal Filtering Options (M4)");
    diagonalOptions.add_options()
            ("allow-overlap", "Allow overlapping matches when searching for neighbouring matches on the same diagonal. No effect if '--diagonal-threshold' is 0 (default).")
            ("diagonal-threshold", po::value<double>()->default_value(0.), "If > 0, diagonal filtering is performed: at least this many seed-matches must be found on the same diagonal between genome1 and genome2.")
            ("local-search-area", po::value<int>()->default_value(100), "Search in an area of this many bp for neighbouring matches on the same diagonal. No effect if '--diagonal-threshold' is 0 (default).")
            ("min-match-distance", po::value<int>()->default_value(0), "Neighbouring matches must be at least this many bp apart. No effect if '--diagonal-threshold' is 0 (default) or if --allow-overlap is stated.");
    po::options_description geometricHashingOptions("Geometric Hashing Options (M6)");
    geometricHashingOptions.add_options()
            ("cube-score-parameter", po::value<double>()->default_value(3), "Parameter mu for Cube scoring.")
            ("cube-score-threshold", po::value<double>()->default_value(10), "Cubes must have at least this score to be considered in match creation.")
            ("geometric-hashing", "Perform geometricHashing on a seedMap instead of reporting plain matches.")
            ("occurrence-per-genome-max", po::value<int>()->default_value(0), "At most this many seed occurrences in any genome. Set to 0 for no threshold.")
            ("occurrence-per-genome-min", po::value<int>()->default_value(1), "At least this many seed occurrences in a genome (if any occurrences in the respective genome).")
            ("tilesize", po::value<int>()->default_value(100), "Cut genome in tiles of this size.");

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
    auto stripExtension = [&vm](std::string const & pathstr) {
        auto path = fs::path(pathstr);
        path.replace_extension("");
        return (path.filename()).string();
    };



    /* set all members according to the program options */
    // general options
    // --artificial-sequence-size-factor
    if (userSet("artificial-sequence-size-factor")) { warnUselessIfNotSet("artificial-sequence-size-factor", "dynamic-artificial-sequences"); }
    artificialSequenceSizeFactor_ = castWithBoundaryCheck<int, size_t>(vm, "artificial-sequence-size-factor", 1, INT_MAX);
    // --dynamic-artificial-sequences
    dynamicArtificialSequences_ = userSet("dynamic-artificial-sequences");
    // --input
    throwMandatory("input");
    inputFiles_ = vm["input"].as<std::vector<std::string>>();
    if (inputFiles_.size() < 2) { throw std::runtime_error("[ERROR] -- '--input' needs at least two genomes"); }
    // --genome1/2
    if (inputFiles_.size() > 2) {
        throwMandatory("genome1");
        throwMandatory("genome2");
    }
    if (userSet("genome1")) {
        throwMandatory("genome2");
        genome1_ = stripExtension(vm["genome1"].as<std::string>());
    } else {
        genome1_ = stripExtension(inputFiles_.at(0));
    }
    if (userSet("genome2")) {
        throwMandatory("genome1");
        genome2_ = stripExtension(vm["genome2"].as<std::string>());
    } else {
        genome2_ = stripExtension(inputFiles_.at(1));
    }
    // --masks
    if (userSet("masks")) {
        masks_ = vm["masks"].as<std::vector<std::string>>();
        SpacedSeedMaskCollection{masks_}; // throws if masks are ill-formed
    }
    // --match-limit
    matchLimit_ = castWithBoundaryCheck<int, size_t>(vm, "match-limit", 0, INT_MAX);
    matchLimit_ = (matchLimit_ == 0) ? ULLONG_MAX : matchLimit_;
    // --match-limit-discard-exceeding
    warnUselessIfNotSet("match-limit-discard-exceeding", "match-limit");
    matchLimitDiscardSeeds_ = userSet("match-limit-discard-exceeding");
    // --weight
    if (masks_.size()) {
        warnUselessIfSet("weight", "masks");
        weight_ = SpacedSeedMaskCollection{masks_}.weight();
    } else {
        throwMandatory("weight");
        weight_ = castWithBoundaryCheck<int, size_t>(vm, "weight", 1, INT_MAX);
    }
    // --weight-fraction
    warnUselessIfSet("weight-fraction", "masks");
    warnUselessIfSet("weight-fraction", "optimal-seed");
    auto weightFraction = castWithBoundaryCheck<double, double>(vm, "weight-fraction", 0, 1);
    // --seed-set-size
    if (masks_.size()) {
        warnUselessIfSet("seed-set-size", "masks");
        seedSetSize_ = masks_.size();
    } else {
        seedSetSize_ = castWithBoundaryCheck<int, size_t>(vm, "seed-set-size", 1, INT_MAX);
    }
    // --optimal-seed
    if (masks_.size()) {
        warnUselessIfSet("optimal-seed", "masks");
        optimalSeed_ = false;
    } else {
        optimalSeed_ = userSet("optimal-seed");
        if (optimalSeed_) {
            SpacedSeedMaskCollection{SpacedSeedMaskCollection::Weight(weight_),
                                     SpacedSeedMaskCollection::SeedSetSize(seedSetSize_)};  // throws if not possible to load
        }
    }
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
    // --quiet
    noProgressbar_ = userSet("quiet");
    // --span
    if (masks_.size()) {
        warnUselessIfSet("span", "masks");
        span_ = SpacedSeedMaskCollection{masks_}.maxSpan();
    } else if (optimalSeed_) {
        warnUselessIfSet("span", "optimal-seed");
        span_ = SpacedSeedMaskCollection{SpacedSeedMaskCollection::Weight(weight_),
                                         SpacedSeedMaskCollection::SeedSetSize(seedSetSize_)}.maxSpan();
    } else if (userSet("span")) {
        warnUselessIfSet("weight-fraction", "span");
        span_ = castWithBoundaryCheck<int, size_t>(vm, "span", 1, INT_MAX);
    } else if (userSet("weight-fraction")) {
        span_ = static_cast<size_t>(std::ceil(static_cast<double>(weight_)/weightFraction));
    } else {
        span_ = weight_;
    }
    SpacedSeedMaskCollection{SpacedSeedMaskCollection::Weight(weight_),
                             SpacedSeedMaskCollection::Span(span_),
                             SpacedSeedMaskCollection::SeedSetSize(seedSetSize_)}; // check if parameter combination works

    // diagonal filter options
    // --diagonal-threshold
    diagonalThreshold_ = castWithBoundaryCheck<double, double>(vm, "diagonal-threshold", 0, DBL_MAX);
    performDiagonalFiltering_ = (diagonalThreshold_ > 0);
    // --local-search-area
    localSearchArea_ = castWithBoundaryCheck<int, size_t>(vm, "local-search-area", 1, INT_MAX);
    // --allow-overlap
    allowOverlap_ = userSet("allow-overlap");
    if (allowOverlap_ && !performDiagonalFiltering_) { std::cerr << "[WARNING] -- '--allow-overlap' is ignored as diagonal filter is not applied" << std::endl; }
    // --min-match-distance
    if (userSet("min-match-distance") && !performDiagonalFiltering_) {
        std::cerr << "[WARNING] -- '--min-match-distance' is ignored as diagonal filter is not applied" << std::endl;
    }
    warnUselessIfSet("min-match-distance", "allow-overlap");
    minMatchDistance_ = userSet("allow-overlap")
            ? 0
            : castWithBoundaryCheck<int, size_t>(vm, "min-match-distance", 0, INT_MAX);

    // geometric hashing options
    // --geometric-hashing
    performGeometricHashing_ = userSet("geometric-hashing");
    // --cube-score-parameter
    warnUselessIfNotSet("cube-score-parameter", "geometric-hashing");
    cubeScoreMu_ = castWithBoundaryCheck<double, double>(vm, "cube-score-parameter", 0, DBL_MAX);
    // --cube-score-threshold
    warnUselessIfNotSet("cube-score-threshold", "geometric-hashing");
    cubeScoreThreshold_ = castWithBoundaryCheck<double, double>(vm, "cube-score-threshold", 0, DBL_MAX);
    // --occurrence-per-genome-max
    warnUselessIfNotSet("occurrence-per-genome-max", "geometric-hashing");
    occurrencePerGenomeMax_ = castWithBoundaryCheck<int, size_t>(vm, "occurrence-per-genome-max", 0, INT_MAX);
    occurrencePerGenomeMax_ = (occurrencePerGenomeMax_ == 0) ? ULLONG_MAX : occurrencePerGenomeMax_;
    // --occurrence-per-genome-min
    warnUselessIfNotSet("occurrence-per-genome-min", "geometric-hashing");
    occurrencePerGenomeMin_ = castWithBoundaryCheck<int, size_t>(vm, "occurrence-per-genome-min", 1, INT_MAX);
    // --tilesize
    warnUselessIfNotSet("tilesize", "geometric-hashing");
    tileSize_ = castWithBoundaryCheck<int, size_t>(vm, "tilesize", 1, INT_MAX);
}



JsonValue Configuration::configJson() const {
    auto jsonstream = std::stringstream(std::ios_base::out);    // stream json of config into string
    JsonStreamDict map(jsonstream);
    map.addValue("allowOverlap", allowOverlap_);
    map.addValue("artificialSequenceSizeFactor", artificialSequenceSizeFactor_);
    map.addValue("createAllMatches", createAllMatches_);
    map.addValue("cubeScoreParameter", cubeScoreMu_);
    map.addValue("cubeScoreThreshold", cubeScoreThreshold_);
    map.addValue("diagonalThreshold", diagonalThreshold_);
    map.addValue("dynamicArtificialSequences", dynamicArtificialSequences_);
    map.addValue("genome1", genome1_);
    map.addValue("genome2", genome2_);
    map.addValue("inputFiles", inputFiles_);
    map.addValue("localSearchArea", localSearchArea_);
    map.addValue("masks", masks_);
    map.addValue("matchLimit", matchLimit_);
    map.addValue("matchLimitDiscardExceeding", matchLimitDiscardSeeds_);
    map.addValue("minMatchDistance", minMatchDistance_);
    map.addValue("quiet", noProgressbar_);
    map.addValue("nThreads", nThreads_);
    map.addValue("occurrencePerGenomeMax", occurrencePerGenomeMax_);
    map.addValue("occurrencePerGenomeMin", occurrencePerGenomeMin_);
    map.addValue("optimalSeed", optimalSeed_);
    map.addValue("performDiagonalFiltering", performDiagonalFiltering_);
    map.addValue("performGeometricHashing", performGeometricHashing_);
    map.addValue("seedSetSize", seedSetSize_);
    map.addValue("span", span_);
    map.addValue("tileSize", tileSize_);
    map.addValue("weight", weight_);
    map.close();
    return JsonValue(jsonstream.str(), JsonValue::stringIsValidJson);
}



std::ostream & operator<<(std::ostream & os, Configuration const & conf) {
    os << "Configuration:" << std::endl;
    os << "\t" << "--allow-overlap " << conf.allowOverlap_ << std::endl;
    os << "\t" << "--artificial-sequence-size-factor " << conf.artificialSequenceSizeFactor_ << std::endl;
    os << "\t" << "createAllMatches_ " << conf.createAllMatches_ << std::endl;
    os << "\t" << "--cube-score-parameter " << conf.cubeScoreMu_ << std::endl;
    os << "\t" << "--cube-score-threshold " << conf.cubeScoreThreshold_ << std::endl;
    os << "\t" << "--diagonal-threshold " << conf.diagonalThreshold_ << std::endl;
    os << "\t" << "--dynamic-artificial-sequences " << conf.dynamicArtificialSequences_ << std::endl;
    os << "\t" << "--input " << conf.inputFiles_ << std::endl;
    os << "\t" << "--genome1 " << conf.genome1_ << std::endl;
    os << "\t" << "--genome2 " << conf.genome2_ << std::endl;
    os << "\t" << "--local-search-area " << conf.localSearchArea_ << std::endl;
    os << "\t" << "--masks " << conf.masks_ << std::endl;
    os << "\t" << "--match-limit " << conf.matchLimit_ << std::endl;
    os << "\t" << "--match-limit-discard-exceeding " << conf.matchLimitDiscardSeeds_ << std::endl;
    os << "\t" << "--min-match-distance " << conf.minMatchDistance_ << std::endl;
    os << "\t" << "--occurrence-per-genome-max " << conf.occurrencePerGenomeMax_ << std::endl;
    os << "\t" << "--occurrence-per-genome-min " << conf.occurrencePerGenomeMin_ << std::endl;
    os << "\t" << "--optimal-seed " << conf.optimalSeed_ << std::endl;
    os << "\t" << "--output " << conf.output_ << std::endl;
    os << "\t" << "--output-artificial-sequences " << conf.outputArtificialSequences_ << std::endl;
    os << "\t" << "--p " << conf.nThreads_ << std::endl;
    os << "\t" << "perform diagonal filtering " << conf.performDiagonalFiltering_ << std::endl;
    os << "\t" << "perform geometric-hashing " << conf.performGeometricHashing_ << std::endl;
    os << "\t" << "--quiet " << conf.noProgressbar_ << std::endl;
    os << "\t" << "--seed-set-size " << conf.seedSetSize_ << std::endl;
    os << "\t" << "--span " << conf.span_ << std::endl;
    os << "\t" << "--tilesize " << conf.tileSize_ << std::endl;
    os << "\t" << "--weight " << conf.weight_ << std::endl;
    return os;
}
