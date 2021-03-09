#! /usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use POSIX qw(ceil floor);

# use singleFasta.fa
die "No input file specified" unless ($ARGV[0]);
die "No output file specified" unless ($ARGV[1]);


# Spaced Seeds Parameter
my $valid = 0;
my $spaced = 1;
while (not $valid) {
	print "Spaced seeds? [Y/n]: ";
	my $in = <STDIN>;
	chomp $in;
	$valid = 1 if ($in =~ /^$/ or $in =~ /^[yn]$/i);
	$spaced = 0 if ($valid and $in =~ /^n$/i);
}

# k Parameter
my $k = 20;
if (not $spaced) {
	$valid = 0;
	while (not $valid) {
		print "k? [20]: ";
		my $in = <STDIN>;
		chomp $in;
		$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
		$k = $in if ($valid and $in =~ /^\d+$/);
	}
}

# Link parameters: max occurrences per genome
my $occurrencePerGenomeMax = 1;
$valid = 0;
while (not $valid) {
	print "Linkset: Max. Per Genome Occurrence Threshold? [1]: ";
	my $in = <STDIN>;
	chomp $in;
	$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
	$occurrencePerGenomeMax = $in if ($valid and $in =~ /^\d+$/);
}

# Link parameters: max occurrences per genome
my $occurrencePerGenomeMin = 0;
$valid = 0;
while (not $valid) {
	print "Linkset: Min. Per Genome Occurrence Threshold? [0]: ";
	my $in = <STDIN>;
	chomp $in;
	$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
	$occurrencePerGenomeMin = $in if ($valid and $in =~ /^\d+$/);
}

# Triangulation parameters: distance
my $triangulationDistance = 0;
$valid = 0;
while (not $valid) {
	print "Triangulation: Distance? [0]: ";
	my $in = <STDIN>;
	chomp $in;
	$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
	$triangulationDistance = $in if ($valid and $in =~ /^\d+$/);
}

# diagonal filtering parameters
my $diagonal = 0;
$valid = 0;
while (not $valid) {
	print "Diagonal? [y/N]: ";
	my $in = <STDIN>;
	chomp $in;
	$valid = 1 if ($in =~ /^$/ or $in =~ /^[yn]$/i);
	$diagonal = 1 if ($valid and $in =~ /^y$/i);
}

my $t = 2;
my $searchArea = 1000;
my $minMatchDistance = 0;
my $allowOverlap = 0;
if ($diagonal) {
	$valid = 0;
	while (not $valid) {
		print "t? [2]: ";
		my $in = <STDIN>;
		chomp $in;
		$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
		$t = $in if ($valid and $in =~ /^\d+$/);
	}

	$valid = 0;
	while (not $valid) {
		print "Search Area? [1000]: ";
		my $in = <STDIN>;
		chomp $in;
		$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
		$searchArea = $in if ($valid and $in =~ /^\d+$/);
	}

	$valid = 0;
	while (not $valid) {
		print "Min Match Distance? [0]: ";
		my $in = <STDIN>;
		chomp $in;
		$valid = 1 if ($in =~ /^$/ or $in =~ /^\d+$/);
		$minMatchDistance = $in if ($valid and $in =~ /^\d+$/);
	}

	$valid = 0;
	while (not $valid) {
		print "Allow overlap? [y/N]: ";
		my $in = <STDIN>;
		chomp $in;
		$valid = 1 if ($in =~ /^$/ or $in =~ /^[yn]$/i);
		$allowOverlap = 1 if ($valid and $in =~ /^y$/i);
	}
}
my $halfSearchArea = ceil($searchArea/2);


my @masks = ($spaced) ? ("10011010011010010011", "01101101001101101000") : ("11111111111111111111");
if ($k != 20) { 
	@masks = ("");
	foreach (1..$k) { $masks[0] = $masks[0]."1"; }
}
my $span = length($masks[0]);

#===============#
# Read in Fasta #
#===============#
my %fasta;
my $head = "";

open(my $fh, "<", $ARGV[0]) or die "Cannot open file $ARGV[0]: $!\n";
while (<$fh>) {
	my $line = $_;
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$head = $line;
		warn "Sequence $head already exists, overwrite with new sequence" if ($fasta{$head});
		$fasta{$head} = "";
	} elsif ($head ne "") {
		$fasta{$head} = $fasta{$head}.$line;
	} else {
		die "No head for line '$line'\n";
	}
}
close($fh) or warn "Cannot close file $ARGV[0]: $!\n";

open(my $out, ">", $ARGV[1]) or die "Cannot open output file $ARGV[1]: $!\n";
print "Writing Output to $ARGV[1]\n";

#=========================================#
# Extract all Seeds and their Occurrences #
#=========================================#
# occurrence = [genome, sequence, position, seed-string, kmer-string]
my %seedToOccurrence;	# seed-string => {maskIndex => [occurrence, ...], ...}
my $numKeys = keys %fasta;

foreach my $seqName (keys %fasta) {
	my $seq = $fasta{$seqName};
	print {$out} "\n// processing $seqName\n// $seq\n";
	for (my $i = 0; $i < (length($seq) - $k + 1); $i++) {
		my $kmer = substr($seq, $i, $k);
		if ($kmer =~ /^[ACGT]+$/) {
			$seqName =~ /^([^\|]+)/;
			my $genome = $1."_orthologs";
			
			for (my $m = 0; $m < scalar @masks; $m++) {
				my $mask = $masks[$m];
			#foreach my $mask (@masks) {
				my $seed = getSeed($kmer, $mask);
				print {$out} "// keeping $kmer ($seed)\n";
				
				$seedToOccurrence{$seed} = {} unless ($seedToOccurrence{$seed});
				$seedToOccurrence{$seed}{$m} = [] unless ($seedToOccurrence{$seed}{$m});
				push(@{$seedToOccurrence{$seed}{$m}}, [$genome, $seqName, $i, $seed, $kmer]);
			}
		} else {
			print {$out} "// discarding $kmer\n";
		}
	}
}

#===========================================================#
# Filter Seeds that are too Frequent and Create a Blacklist #
#===========================================================#
my %filteredSeedToOccurrence;	# seed-string => {maskIndex => [occurrence, ...], ...}
foreach my $seed (keys %seedToOccurrence) {
	for (my $m = 0; $m < scalar @masks; $m++) {
		if ($seedToOccurrence{$seed}{$m}) {	# not neccessarily every seed genereated by both masks, if the only genereated gets discarded, the loose key in filtered hash disturbs the valid seeds code!
			my %seqToCount;
			my $discard = 0;
			foreach my $occ (@{$seedToOccurrence{$seed}{$m}}) {
				#my $seed = $$occ[3];
				my $seqName = $$occ[1];
				$seqToCount{$seqName} = 0 unless ($seqToCount{$seqName});
				$seqToCount{$seqName}++;
			}
			if (not $discard) {
				$filteredSeedToOccurrence{$seed} = {} unless ($filteredSeedToOccurrence{$seed});
				$filteredSeedToOccurrence{$seed}{$m} = $seedToOccurrence{$seed}{$m};
			}
		}
	}
}

#===============================#
# Generate Code for Valid Seeds #
#===============================#
sub generateOccurrenceCode {
	my $occ = shift;
	my $genome = $$occ[0];
	my $seqName = $$occ[1];
	my $pos = $$occ[2];
	my $kmer = $$occ[4];
	return 'KmerOccurrence{static_cast<uint8_t>(idMap.queryGenomeIDConst("'.$genome.'")), static_cast<uint32_t>(idMap.querySequenceID("'.$seqName.'")), '.$pos.', false, "'.$kmer.'"}'
}
my @allMatches; # [[occ1, occ2], ...]
my %diagonalsToMatches;	# "$seqName1|$seqName2|dist" => [[occ1, occ2], ...]
my $expectedSeedsCode = "";
my $expectedMatchesCode = "";
my $seedCodeCounter = 0;
my $matchCodeCounter = 0;
foreach my $seed (keys %filteredSeedToOccurrence) {
	for (my $m = 0; $m < scalar @masks; $m++) {
		$expectedSeedsCode = $expectedSeedsCode.'kmerMapExp[TwoBitKmer<TwoBitKmerDataShort>{"'.$seed.'"}].emplace_back();'."\n"; $seedCodeCounter++;
		# expected seeds
		my @codes;
		my $in0 = 0;
		my $in1 = 0;
		foreach my $occ (@{$filteredSeedToOccurrence{$seed}{$m}}) {
			my $genome = $$occ[0];
			my $occCode = generateOccurrenceCode($occ);
			$expectedSeedsCode = $expectedSeedsCode.'kmerMapExp[TwoBitKmer<TwoBitKmerDataShort>{"'.$seed.'"}].at('.$m.').emplace_back('.$occCode.');'."\n"; $seedCodeCounter++;
			
			push(@codes, $occCode);
			$in0 = 1 if $genome =~ /hg38/;
			$in1 = 1 if $genome =~ /mm10/;
		}
		
		# all matches
		for (my $i = 0; $i < scalar @{$filteredSeedToOccurrence{$seed}{$m}}; $i++) {
			for (my $j = ($i + 1); $j < scalar @{$filteredSeedToOccurrence{$seed}{$m}}; $j++) {
				my $occ1 = ${$filteredSeedToOccurrence{$seed}{$m}}[$i];
				my $occ2 = ${$filteredSeedToOccurrence{$seed}{$m}}[$j];
				if ($$occ1[0] eq $$occ2[0]) { next; } # same genomes
				my $match = [];
				if ($$occ1[1] lt $$occ2[1]) {
					$match = [$occ1, $occ2];
				} else {
					$match = [$occ2, $occ1];
				}
				push(@allMatches, $match);

				# generate code
				my $gen1 = $$match[0][0];
				my $gen2 = $$match[1][0];
				my $occCode1 = generateOccurrenceCode($$match[0]);
				my $occCode2 = generateOccurrenceCode($$match[1]);
				#if ($gen1 ne $gen2 && ($gen1 =~ /hg38|mm10/) && ($gen2 =~ /hg38|mm10/)) { # EQUIVALENT TO createAllMatches_ = false
					if ($spaced) {
						$expectedMatchesCode = $expectedMatchesCode.'matchesExp.emplace(KmerOccurrencePair{'.$occCode1.', '.$occCode2.', '.$span.'});'."\n"; $matchCodeCounter++;
					} elsif ($diagonal) {
						$expectedMatchesCode = $expectedMatchesCode.'matchesExp.emplace('.$occCode1.', '.$occCode2.', '.$span.');'."\n"; $matchCodeCounter++;
					} else {
						$expectedMatchesCode = $expectedMatchesCode.'matchesExp.emplace_back('.$occCode1.', '.$occCode2.', '.$span.');'."\n"; $matchCodeCounter++;
					}
				#}
			}
		}

		# diagonals
		for (my $i = 0; $i < scalar @{$filteredSeedToOccurrence{$seed}{$m}}; $i++) {
			for (my $j = ($i + 1); $j < scalar @{$filteredSeedToOccurrence{$seed}{$m}}; $j++) {
				my $occ1 = ${$filteredSeedToOccurrence{$seed}{$m}}[$i];
				my $occ2 = ${$filteredSeedToOccurrence{$seed}{$m}}[$j];
				if ($$occ1[1] eq $$occ2[1]) { next; }	# same sequences
				
				my $match = [];
				if ($$occ1[1] lt $$occ2[1]) {
					$match = [$occ1, $occ2];
				} else {
					$match = [$occ2, $occ1];
				}
				my $distance = ${$$match[1]}[2] - ${$$match[0]}[2];
				my $diagonal = ${$$match[0]}[1]."|".${$$match[1]}[1]."|".$distance;
				$diagonalsToMatches{$diagonal} = [] unless ($diagonalsToMatches{$diagonal});
				push(@{$diagonalsToMatches{$diagonal}}, $match);
				# print "Diagonal $diagonal has matches [".join(",", @{$$match[0]})."], [".join(",", @{$$match[1]})."]\n";
			}
		}
	}
	
	if ($seedCodeCounter >= 100) {
		$expectedSeedsCode = $expectedSeedsCode."\n\n\n";
		$seedCodeCounter = 0;
	}
	if ($matchCodeCounter >= 100) {
		$expectedMatchesCode = $expectedMatchesCode."\n\n\n";
		$matchCodeCounter = 0;
	}
}
print {$out} $expectedSeedsCode;
print {$out} $expectedMatchesCode;

#=============================#
# Generate Triangulation Code #
#=============================#
sub keyFromMatch {
	my $match = shift;
	my $occ0 = $$match[0];
	my $occ1 = $$match[1];
	my $key = join("/", $$occ0[0], $$occ0[1], $$occ0[2])."-".join("/", $$occ1[0], $$occ1[1], $$occ1[2]);	# using also seeds does not work in triangulation as inferred occs don't have one
	return $key;
}
if ($triangulationDistance > 0) {
	my $triangulationCode = "";
	my $triangulationMatchesCode = "";
	my $triangulationCodeCounter = 0;
	my $triangulationMatchCounter = 0;
	my %seenMatches;
	my %seenInferred;
	foreach (@allMatches) {
		$seenMatches{keyFromMatch($_)} = 1;
		$triangulationMatchesCode = $triangulationMatchesCode."triangulationMatches.emplace(".generateOccurrenceCode($$_[0]).", ".generateOccurrenceCode($$_[1]).", $span);\n"; $triangulationMatchCounter++;
		if ($triangulationMatchCounter >= 100) { $triangulationMatchesCode = $triangulationMatchesCode."\n\n\n"; $triangulationMatchCounter = 0; }
	}
	# create all inferrable matches
	foreach my $match1 (@allMatches) {
		my $match1genome1 = ${$$match1[0]}[0];
		my $match1genome2 = ${$$match1[1]}[0];
		
		my $match1het_hg = ($match1genome1 =~ /^hetGla2/ and $match1genome2 =~ /^hg38/) ? 1 : 0;
		my $match1het_mm = ($match1genome1 =~ /^hetGla2/ and $match1genome2 =~ /^mm10/) ? 1 : 0; 
		my $match1hg_mac = ($match1genome1 =~ /^hg38/ and $match1genome2 =~ /^macFas5/) ? 1 : 0;
		my $match1mac_mm = ($match1genome1 =~ /^macFas5/ and $match1genome2 =~ /^mm10/) ? 1 : 0;
		
		if (not ($match1het_hg or $match1het_mm or $match1hg_mac or $match1mac_mm)) { next; }
			
		foreach my $match2 (@allMatches) {
			if (keyFromMatch($match1) eq keyFromMatch($match2)) { next; }
			
			my $match2genome1 = ${$$match2[0]}[0];
			my $match2genome2 = ${$$match2[1]}[0];
			
			my $match2het_hg = ($match2genome1 =~ /^hetGla2/ and $match2genome2 =~ /^hg38/) ? 1 : 0;
			my $match2het_mm = ($match2genome1 =~ /^hetGla2/ and $match2genome2 =~ /^mm10/) ? 1 : 0; 
			my $match2hg_mac = ($match2genome1 =~ /^hg38/ and $match2genome2 =~ /^macFas5/) ? 1 : 0;
			my $match2mac_mm = ($match2genome1 =~ /^macFas5/ and $match2genome2 =~ /^mm10/) ? 1 : 0;
			
			if (($match1het_hg and $match2het_mm) 
			    or ($match1het_mm and $match2het_hg)
			    or ($match1hg_mac and $match2mac_mm)
			    or ($match1mac_mm and $match2hg_mac)) {
				#valid genome constellation at this point
				my $hgOcc;
				my $mmOcc;
				my $thirdOcc1;
				my $thirdOcc2;
				my $hgLeft;
				if ($match1het_hg) {
					$hgOcc = $$match1[1];
					$mmOcc = $$match2[1];
					$thirdOcc1 = $$match1[0];
					$thirdOcc2 = $$match2[0];
					$hgLeft = ($$thirdOcc1[2] <= $$thirdOcc2[2]) ? 1 : 0;
				} elsif ($match1het_mm) {
					$hgOcc = $$match2[1];
					$mmOcc = $$match1[1];
					$thirdOcc1 = $$match1[0];
					$thirdOcc2 = $$match2[0];
					$hgLeft = ($$thirdOcc2[2] <= $$thirdOcc1[2]) ? 1 : 0;
				} elsif ($match1hg_mac) {
					$hgOcc = $$match1[0];
					$mmOcc = $$match2[1];
					$thirdOcc1 = $$match1[1];
					$thirdOcc2 = $$match2[0];
					$hgLeft = ($$thirdOcc1[2] <= $$thirdOcc2[2]) ? 1 : 0;
				} elsif ($match1mac_mm) {
					$hgOcc = $$match2[0];
					$mmOcc = $$match1[1];
					$thirdOcc1 = $$match1[0];
					$thirdOcc2 = $$match2[1];
					$hgLeft = ($$thirdOcc2[2] <= $$thirdOcc1[2]) ? 1 : 0;
				}
				if ($$thirdOcc1[1] eq $$thirdOcc2[1]) {
					# valid sequences, check distance
					my $pos1 = $$thirdOcc1[2];
					my $pos2 = $$thirdOcc2[2];
					my $absDistance = ($pos1 >= $pos2) ? $pos1 - $pos2 : $pos2 - $pos1;
					if ($absDistance <= $triangulationDistance) {
						my $newOcc = ($hgLeft)
							? [$$hgOcc[0],$$hgOcc[1],($$hgOcc[2]+$absDistance),"",""]
							: [$$hgOcc[0],$$hgOcc[1],($$hgOcc[2]-$absDistance),"",""];
						my $newMatch = [$newOcc, $mmOcc];
						if (not $seenMatches{keyFromMatch($newMatch)} and not $seenInferred{keyFromMatch($newMatch)}) {
							print {$out} "// Inferring (hg,".$$newOcc[2].") -- (mm,".$$mmOcc[2].") (distance: $absDistance) from ".join("/",@{$$match1[0]})."--".join("/",@{$$match1[1]})." and ".join("/",@{$$match2[0]})."--".join("/",@{$$match2[1]})."\n";
							$triangulationCode = $triangulationCode."triangulationResult.emplace(".generateOccurrenceCode($newOcc).", ".generateOccurrenceCode($mmOcc).");\n"; $triangulationCodeCounter++;
							if ($triangulationCodeCounter >= 100) { $triangulationCode = $triangulationCode."\n\n\n"; $triangulationCodeCounter = 0; }
							$seenInferred{keyFromMatch($newMatch)} = 1;
						}
					}
				}
			}
		}
	}
	print {$out} $triangulationMatchesCode;
	print {$out} $triangulationCode;
}

#========================#
# Generate Diagonal Code #
#========================#
print ${out} "\\ ".scalar(keys %diagonalsToMatches)." diagonals\n";
my $diagonalMatchesCode = "";
my $notInHumanMouse = 0;
my $tooCloseOrOverlap = 0;
my $tooFewOnDiagonal = 0;
my $tooFewNeighbours = 0;
for my $diagonal (keys %diagonalsToMatches) {
	# only human/mouse
	my $firstDiagMatch = ${$diagonalsToMatches{$diagonal}}[0];
	if (not (${$$firstDiagMatch[0]}[0] =~ /^hg38|mm10/
	         and ${$$firstDiagMatch[1]}[0] =~ /^hg38|mm10/) ) {
		print {$out} "// discarding diagonal $diagonal as first match (match [".join(",", @{$$firstDiagMatch[0]})."], [".join(",", @{$$firstDiagMatch[1]})."] (${$$firstDiagMatch[0]}[4])) not in human and mouse\n";
		$notInHumanMouse += scalar @{$diagonalsToMatches{$diagonal}};
		next;
	}
	
	my @filteredDiagonal;
	my %posToMatch;
	foreach my $match (@{$diagonalsToMatches{$diagonal}}) {
		my $pos = ${$$match[0]}[2];
		die "duplicate position in diagonal $diagonal\n" if $posToMatch{$pos};
		$posToMatch{$pos} = $match;
	}
	my @sortedPositions = sort {$a <=> $b} keys %posToMatch;
	my $nextValidPosition = -1;
	foreach my $pos (@sortedPositions) {
		if ($pos >= $nextValidPosition) {
			print {$out} "// match at position $pos on diagonal $diagonal lies behind 'next valid position' on this diagonal: $nextValidPosition)\n";
			push(@filteredDiagonal, $posToMatch{$pos});
			
			my $currentMatch = $posToMatch{$pos};
			my $matchPos = ${$$currentMatch[0]}[2];
			if ($allowOverlap) {
				$nextValidPosition = ($matchPos + 1);
			} else {
				$nextValidPosition = ($matchPos + $k + $minMatchDistance);
			}
		} else {
			print {$out} "// match at position $pos on diagonal $diagonal overlaps with / too close to previous (next valid position on this diagonal: $nextValidPosition)\n";
			$tooCloseOrOverlap += 1;
		}
	}
	
	if (scalar @filteredDiagonal < $t) {
		print {$out} "// discarding diagonal $diagonal as there are too few matches\n";
		$tooFewOnDiagonal += 1;
	} else {
		foreach my $match (@filteredDiagonal) {
			my $matchMid = ${$$match[0]}[2] + ceil($k / 2) - 1;
			my $leftBorder = $matchMid - $halfSearchArea;
			my $rightBorder = $matchMid + $halfSearchArea;
			
			my %positions;
			foreach my $neighbour (@filteredDiagonal) {
				my @positions = ${$$neighbour[0]}[2]..(${$$neighbour[0]}[2] + $k - 1);
				foreach (@positions) { 
					if ($leftBorder <= $_ and $_ <= $rightBorder) { $positions{$_} = 1; }
				} 
			}
			my $seedWeight = ($spaced) ? 10 : $k;
			my $seedCount = (scalar keys %positions) / $seedWeight;
			
			if ($seedCount >= $t) {
				my $genome1 = ${$$match[0]}[0];
				my $genome2 = ${$$match[1]}[0];
				my $seqName1 = ${$$match[0]}[1];
				my $seqName2 = ${$$match[1]}[1];
				my $pos1 = ${$$match[0]}[2];
				my $pos2 = ${$$match[1]}[2];
				my $kmer1 = ${$$match[0]}[4];
				my $kmer2 = ${$$match[1]}[4];
				if (($genome1 =~ /^hg38/ and $genome2 =~ /^mm10/)
                    or ($genome1 =~ /^mm10/ and $genome2 =~ /^hg38/)) {
					$diagonalMatchesCode = $diagonalMatchesCode."expectedMatches.emplace_back(KmerOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst(\"$genome1\")), static_cast<uint32_t>(idMap.querySequenceID(\"$seqName1\")), $pos1, false, \"$kmer1\"), KmerOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst(\"$genome2\")), static_cast<uint32_t>(idMap.querySequenceID(\"$seqName2\")), $pos2, false, \"$kmer2\"), $span);\n";
					#my $midPos1 = $pos1 + floor($k / 2);
					#my $midPos2 = $pos2 + floor($k / 2);
					#print "($midPos1, $midPos2) in ($seqName1, $seqName2)\n";
				} else {
					print {$out} "// match [".join(",", @{$$match[0]})."], [".join(",", @{$$match[1]})."] (${$$match[0]}[4]) not in human or mouse\n";
					$notInHumanMouse += 1;
				}
				print {$out} "// match [".join(",", @{$$match[0]})."], [".join(",", @{$$match[1]})."] (${$$match[0]}[4]) valid\n";
			} else {
				print {$out} "// match [".join(",", @{$$match[0]})."], [".join(",", @{$$match[1]})."] (${$$match[0]}[4]) too few neighbours\n";
				$tooFewNeighbours += 1;
			}
		}
	}
}
print {$out} $diagonalMatchesCode;
print {$out} "not in human/mouse: $notInHumanMouse\n";
print {$out} "too close/overlap: $tooCloseOrOverlap\n";
print {$out} "too few on diagonal: $tooFewOnDiagonal\n";
print {$out} "too few neighbours: $tooFewNeighbours\n";

#==============#
# Linkset Data #
#==============#
sub addGenomeTile {	# needs a vector-ref of vectors [[occ0, occ1, ...], ...],
	                # link (vector-ref of occurrences to fill) [occ0, ...]
	                # linkVector (vetor-ref of links to fill)
	my $genomeOccurrences = shift;
	my $genomeOccurrencesCopy = [@$genomeOccurrences];
	my $link = shift;
	my $linkCopy = [@$link];
	my $linkVector = shift;
	
	my $thisGenome = shift @$genomeOccurrencesCopy;
	if (not defined $thisGenome) {
		push(@$linkVector, $link);
	} else {
		foreach my $occ (@$thisGenome) {
			push(@$linkCopy, $occ);
			addGenomeTile($genomeOccurrencesCopy, $linkCopy, $linkVector);
		}
	}
}



my $linksetCode = "";
my $linksetCodeCounter = 0;
my $count = 0;
my $discarded = 0;
foreach my $seed (keys %filteredSeedToOccurrence) {
	for (my $m = 0; $m < scalar @masks; $m++) {
		my %genomeCount;
		my $tooMany = 0;
		foreach my $occ (@{$filteredSeedToOccurrence{$seed}{$m}}) {
			my $genome = $$occ[0];
			$genomeCount{$genome} = 0 unless ($genomeCount{$genome});
			$genomeCount{$genome}++;
			$tooMany = 1 if ($genomeCount{$genome} > $occurrencePerGenomeMax);
		}
		
	#	print Dumper(%genomeCount);
		
		$discarded++ if ($tooMany);
		next if ($tooMany);
		
		my $tooFew = 0;
		foreach (keys %genomeCount) {
			$tooFew = 1 if ($genomeCount{$_} < $occurrencePerGenomeMin);
		}
		$discarded++ if ($tooFew);
		next if ($tooFew);
		
		my $notInRef = (not defined $genomeCount{"hg38_orthologs"}) ? 1 : 0;
		$discarded++ if ($notInRef);
		next if ($notInRef);
		
		my $sum = 0;
		foreach (keys %genomeCount) { $sum += $genomeCount{$_}; }
		$sum -= $genomeCount{"hg38_orthologs"};
		my $onlyInRef = ($sum == 0) ? 1 : 0;
		$discarded++ if ($onlyInRef);
		next if ($onlyInRef);
		
	#	print "...valid\n\n";
		
		# create Links from valid seeds
		my %genomeIdx;
		my @genomeOccurrences;
		foreach my $occ (@{$filteredSeedToOccurrence{$seed}{$m}}) {
			my $genome = $$occ[0];
			if (not $genomeIdx{$genome}) {
				push(@genomeOccurrences, [$occ]);
				$genomeIdx{$genome} = (scalar @genomeOccurrences) - 1;
			} else {
				push(@{$genomeOccurrences[$genomeIdx{$genome}]}, $occ);
			}
		}
		
		my $linkVector = [];
		addGenomeTile(\@genomeOccurrences, [], $linkVector);
		
		# LinksetType = std::unordered_map<std::shared_ptr<Link const>, size_t>;
		# Link: insertTile(uint8_t genomeID,  uint32_t sequenceID, size_t tileID, bool reverse, std::string const & kmer)
		# occurrence = [genome, sequence, position, seed-string, kmer-string]
		foreach my $link (@$linkVector) {
			$linksetCode = $linksetCode."auto link$count = std::make_shared<Link>();\n"; $linksetCodeCounter++;
			foreach my $occ (@$link) {
				my $genome = $$occ[0];
				my $seqName = $$occ[1];
				my $pos = $$occ[2];
				my $kmer = $$occ[3];
				$linksetCode = $linksetCode."link$count".'->insertOccurrence(static_cast<uint8_t>(idMap.queryGenomeIDConst("'.$genome.'")), static_cast<uint32_t>(idMap.querySequenceID("'.$seqName.'")), '.$pos.', false, "'.$kmer.'");'."\n"; $linksetCodeCounter++;
			}
			# $linksetCode = $linksetCode."link$count->sortTiles();\n";
			$linksetCode = $linksetCode."if (linkset.find(link$count) != linkset.end()) { ++linkset.at(link$count); } else { linkset[link$count] = 1; }\n"; $linksetCodeCounter++;
			$count++;
			if ($linksetCodeCounter >= 100) { $linksetCode = $linksetCode."\n\n\n";  $linksetCodeCounter = 0; }
		}
	}
}
print {$out} $linksetCode;
print {$out} "discarded kmers in link creation: $discarded\n";
print {$out} "links creatied: $count\n";

close($out) or warn "Cannot close output file $ARGV[1]: $!\n";

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



sub getSeed {
	my $kmer = shift;
	my $mask = shift;
	die "kmer ($kmer, ".length($kmer).") and mask ($mask, ".length($mask).") don't match" if (length($kmer) != length($mask));
	my $seed = "";
	for (my $i = 0; $i < length($mask); $i++) {
		if (substr($mask, $i, 1) eq "1") {
			$seed = $seed.substr($kmer, $i, 1);
		}
	}
	return $seed;
}
