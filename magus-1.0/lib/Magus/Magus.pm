package Magus::Magus;

################################################################################
# *
# * Copyright Carole Dossat  / Institut de Genomique / DSV / CEA
# *                            <cdossat@genoscope.cns.fr>
# *           Jean-Marc Aury / Institut de Genomique / DSV / CEA
# *                            <jmaury@genoscope.cns.fr>
# *			  Amin Madoui    / Institut de Genomique / DSV / CEA
# *                            <amadoui@genoscope.cns.fr>
# *
# * This software, called MAGUS is a computer program whose purpose is to
# * improve assemblies by using maps and NGS sequencing technologies
# *
# * This software is governed by the CeCILL license under French law and
# * abiding by the rules of distribution of free software.  You can  use,
# * modify and/ or redistribute the software under the terms of the CeCILL
# * license as circulated by CEA, CNRS and INRIA at the following URL
# * "http://www.cecill.info".
# *
# * As a counterpart to the access to the source code and  rights to copy,
# * modify and redistribute granted by the license, users are provided only
# * with a limited warranty  and the software's author,  the holder of the
# * economic rights,  and the successive licensors  have only  limited
# * liability.
# *
# * In this respect, the user's attention is drawn to the risks associated
# * with loading,  using,  modifying and/or developing or reproducing the
# * software by the user in light of its specific status of free software,
# * that may mean  that it is complicated to manipulate,  and  that  also
# * therefore means  that it is reserved for developers  and  experienced
# * professionals having in-depth computer knowledge. Users are therefore
# * encouraged to load and test the software's suitability as regards their
# * requirements in conditions enabling the security of their systems and/or
# * data to be ensured and,  more generally, to use and operate it in the
# * same conditions as regards security.
# *
# * The fact that you are presently reading this means that you have had
# * knowledge of the CeCILL license and that you accept its terms.
################################################################################

use strict;
use Magus::AlignReadsMP;
use Magus::AlignTags;
use Magus::Wgp;
use Magus::CloneLink;
use Magus::Connection;

### TO DO perldoc


### FUNCTIONS FOR MAGUS ###

sub getScaffoldslengths ($$) {
	my ($file,$path) = @_;
	my %scaffoldsLengths = ();
	
	my $fastalength = $path."fastalength";
	
	open (LENGTH, "$fastalength -f $file |") 
		or die "[ERROR] :: fill :: Cannot use fastalength : $!\n";

	while (<LENGTH>) {
		if ($_ =~ /^(\d+)\s+(\S+)\s*$/ ) {
			my $scaffoldId = $2;
			my $length = $1; 
			$scaffoldsLengths{$scaffoldId} = $length;
		}
		else {
			die "[ERROR] problem with parsing fastalength line $_\n";
		}		
	}
	close LENGTH;
	return \%scaffoldsLengths;
}

sub checkWgp2map ($$$) {
	my ($help,$mapFile,$tagsFile) = @_;
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus wg2map -w wgpFile -t tags.bam [options] [-h]\n";
	$synopsis .= "\t-w <string>\twgpFile: WGP data\n";
	$synopsis .= "\t-t <string>\ttags.bam: tags alignment on the assembly (BAM)\n\n";
	$synopsis .= "OPTIONAL PARAMETERS:\n";
	$synopsis .= "\t-p <string>\tprefix for output files";
	$synopsis .= "\t(default : magus)\n";
	$synopsis .= "\t-h\t\tthis help\n\n";
	$synopsis .= "EXAMPLE:\n\tmagus wg2map -w wgpFile -t tags.bam -p Arabido\n\n";
	
	
	if ($help) {
		die $synopsis;
	}
	if ($mapFile eq "") {
		die "$synopsis";
	}
	if ($tagsFile eq "") {
		die "$synopsis";
	}
}

sub wgp2map ($$$$$) {
	my ($mapFile,$path,$tagsFile,$checkedTagsFile,$orderedTagsFile) = @_;
	### tags vs assembly
	my $tags = new Magus::AlignTags();
	$tags->binPath($path);
	$tags->fill($tagsFile);
	### tags wgp
	my $map = new Magus::Wgp();
	$map->fill($mapFile,$tags->matches());
	### add ranks and clean bad mapped tags
	$tags->addTagsData($map->tags(),$orderedTagsFile);	
	### create file
	$tags->printAnchors($checkedTagsFile);
}

sub checkMap2links ($$) {
	my ($help,$checkedTagsFile) = @_;
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus map2links -a anchoring_file.txt [options] [-h]\n";
	$synopsis .= "\t-a <string>\tanchoring_file.txt: anchored assembly (wgp2map output)\n\n";
	$synopsis .= "OPTIONAL PARAMETERS:\n";	
	$synopsis .= "\t-l <string>\toutput file containing links between contigs/scaffolds"; 
	$synopsis .= "\t(default: prefix_map_links.txt)\n";
	$synopsis .= "\t-m <string>\tBin path\t(default: \$PATH)\n";
	$synopsis .= "\t-p <string>\tprefix for output files";
	$synopsis .= "\t(default: magus)\n";
	$synopsis .= "\t-h\t\tthis help\n\n";
	$synopsis .= "EXAMPLE:\n\tmagus map2links -a anchoring_file.txt -l links_file.txt -m /env/bin/";
	$synopsis .= "\n\tmagus map2links -a anchoring_file.txt -m /env/bin/ -p Arabido\n\n";
	if ($help) {
		die $synopsis;
	}
	if (! -f $checkedTagsFile) {
		die "$synopsis";
	}
}

sub map2links ($$) {
	my ($checkedTagsFile,$linksFile) = @_;
	open (my $handleLinks, ">$linksFile") # clean old file if exists
		or die "Cannot open $linksFile : $!\n";
	my $tags = new Magus::AlignTags();
	$tags->createLinks($checkedTagsFile);
	$tags->createLinksNotIncluded($checkedTagsFile);
	$tags->createLinksNoNoise($checkedTagsFile);
	#$tags->createLinksMoreTags($checkedTagsFile);
	$tags->printLinks($handleLinks);
	close $handleLinks;
}

sub checkPairs2links ($$$) {
	my ($help,$scaffoldsFile,$readsData) = @_;
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus pairs2links -f assembly.fa -l links_file.txt -b file.bam,m,sd,s [options] [-h]\n";
	$synopsis .= "\t-f <string>\tassembly.fa: assembly file (FASTA)\n";
	$synopsis .= "\t-l <string>\tlinks_file.txt: file containing links between contigs/scaffolds\n";
	$synopsis .= "\t-b <string>\tfile.bam,m,sd,s: paired reads alignment (BAM), library median size (bp), library standart deviation (bp), reads size (bp)\n\n";
	$synopsis .= "OPTIONAL PARAMETERS:\n";	
	$synopsis .= "\t-v <string>\tpath to samtools\t(default: \$PATH)\n";
	$synopsis .= "\t-p <string>\tprefix for output files";
	$synopsis .= "\t(default: magus)\n\n";
	$synopsis .= "EXAMPLE:\n\tmagus pairs2links -f Arabidopsis.fa -l links_file.txt -b mapping_library1.bam,3500,600,101 -b mapping_library2.bam,6000,1000,151 -v /env/bin/PathToSamtools -p Arabido\n\n";
	
#	$synopsis .= "\n\n This module looks for proofs between putative links by using mate pair sequences.\n";
#	$synopsis .= " It needs reads file and assembly file.\n";
#	$synopsis .= "\n Outputs are 2 files : one contains the links with de format (prefix_new_connections.de) and the other contains statistical logs (prefix_new_connections.log).\n\n";
	
	if ($help) {
		die $synopsis;
	}
	if ($scaffoldsFile eq "") {
		die "$synopsis";
	}
	if (scalar(@$readsData) == 0 ) {
		die "$synopsis";
	} 
}

sub pairs2links ($$$$$$$$) {
	my ($samtoolsPath,$fastalengthPath,$readsData,$scaffoldsFile,$linksFile,$noLinkFile,$connectionsFile,$connectionsLog) = @_;
		
	### initialize
	open (CONNECTIONS, ">$connectionsFile")
		or die "Cannot open connectionsFile : $connectionsFile : $!\n";
	close CONNECTIONS;		
	
	## add reads mapping data 
	my $readsAlignment = new Magus::AlignReadsMP();
	$readsAlignment->binPath($samtoolsPath);
	
	foreach my $data (@$readsData) {
		my ($readsFile,$bankMean,$bankDev,$readsSize) = split (/,/, $data);
		$readsAlignment->fill($readsFile,$bankMean,$bankDev,$readsSize);	
	}

	### get scaffolds lengths
	my $scaffoldsLengths = getScaffoldslengths($scaffoldsFile,$fastalengthPath);

	### for stats
	my ($goodLinks,$badLinks,$allLinks) = (0,0,0);

	### create checked connections
	open (LINKS, "<$linksFile")
		or die "Cannot open linksFile : $linksFile : $!\n";
	
	open (NOLINK, ">$noLinkFile")
		or die "Cannot open noLinkFile : $noLinkFile : $!\n";
	
	while (<LINKS>) {
		$allLinks++;
		chomp $_;
		my ($scaff1,$scaff2) = split (/_/, $_);
		my $length1 = $$scaffoldsLengths{$scaff1};
		my $length2 = $$scaffoldsLengths{$scaff2};
	
		if ( (!defined $length1) || (!defined $length2) ) {
			next;
		}
	
		my %connectionsByOrientation = ();
		$connectionsByOrientation{"++"} = new Magus::Connection(orientation => "++", scaffoldId1 => $scaff1, scaffoldId2 => $scaff2);
		$connectionsByOrientation{"+-"} = new Magus::Connection(orientation => "+-", scaffoldId1 => $scaff1, scaffoldId2 => $scaff2);
		$connectionsByOrientation{"-+"} = new Magus::Connection(orientation => "-+", scaffoldId1 => $scaff1, scaffoldId2 => $scaff2);
		$connectionsByOrientation{"--"} = new Magus::Connection(orientation => "--", scaffoldId1 => $scaff1, scaffoldId2 => $scaff2);
		
		my $links = $readsAlignment->createCloneLinks($scaff1,$scaff2,$length1,$length2);				
			
		# if no links
		#if (scalar(@$links) == 0) {
		#	$badLinks++;
		#	next;
		#}
			
		# add to connection
		foreach my $link (@$links) {
			my $or = $link->orientation();
			push @{$connectionsByOrientation{$or}->links()}, $link;		
		} 		
																		
		# calculate gap and choose good orientation and connection														
		my %hash = ();
		foreach my $connect (values %connectionsByOrientation) {
			$connect->sizeGap();
			my $nbLinks = $connect->validLinksNb();
			my $or = $connect->orientation();
			$hash{$or} = $nbLinks;
		}
				
		my $selectedOrientation = "*";
		my @sort = sort { $b <=> $a } (values %hash);
		#if ( ($sort[0] != 0) && ($sort[0] != $sort[1]) ) { # no result or undefined
		if ($sort[0] != 0) { # less stringent
			foreach my $key (keys %hash) {
				if ($hash{$key} == $sort[0]) {
					$selectedOrientation = $key;	
				}
			}
		}
						
		# print results
		if ($selectedOrientation ne "*") {
			$connectionsByOrientation{$selectedOrientation}->printDeFormat($connectionsFile) ;
			$goodLinks++;			
		}	
		else {
			$badLinks++;
			print NOLINK "${scaff1}_$scaff2\n";
		}					
	}
	close LINKS;
	close NOLINK;
	
	open (LOG, ">$connectionsLog")
		or die "Cannot open connection log : $connectionsLog : $!\n";
		
	print LOG "number of links : $allLinks\n";
	print LOG "number of validates links : $goodLinks\n";
	print LOG "number of not validated links : $badLinks\n";
	
	close LOG;
}

sub execCmd ($) {
	my ($cmd) = @_;
	system($cmd) == 0 
		or die "[ERROR] : Cannot execute $cmd : $!\n";
} 

sub checkLinks2scaf ($$) {
	my ($help,$scaffoldsFile) = @_;
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus links2scaf -f assembly.fa -c links.de [options] [-h]\n";
	$synopsis .= "\t-f <string>\tassembly.fa: assembly file (FASTA)\n";
	$synopsis .= "\t-c <string>\tlinks.de: file containing links in DE format\n\n";
	$synopsis .= "OPTIONAL PARAMETERS:\n";	
	$synopsis .= "\t-z <string>\tpath to sga\t(default: \$PATH)\n";
	$synopsis .= "\t-g <string>\tpath to getseq\t(default: \$PATH)\n";
	$synopsis .= "\t-p <string>\tprefix for output files\t(default: magus)\n";
	$synopsis .= "EXAMPLE:\n\tmagus links2scaf -f Arabidopsis.fa -c links.de -z /env/bin/PathToSga -g /env/bin/PathToGetseq -p Arabido\n\n";
#	$synopsis .= "\n This module outputs the final assembly by using validated links.\n";
#	$synopsis .= " It needs the links file and the first assembly.\n";
#	$synopsis .= "\n Outputs are sga output file, sga log file and final assembly file prefix_all_scaffolds.fa\n\n";
	
	if ($help) {
		die $synopsis;
	}
	if ($scaffoldsFile eq "") {
		die "$synopsis";
	}
}

sub getScafNameFromSgaScafFile ($) {
	my ($scafFile) = @_;
	my %hash = ();
	
	open (SGA, "<$scafFile")
		or die "Cannot open sga scaffold file : $scafFile : $!\n";
	
	while (<SGA>) {
		if ($_ =~ /^(\S+)\s+$/ ) {
			$hash{$1} = 1;
		}
		elsif ($_ =~ /^(\S+)([\s+\S+]+)$/) {
			$hash{$1} = 1;
			my @split1 = split(/\s+/,$2);	
			foreach my $s (@split1) {
				if ($s eq "") {
					next;
				}
				my @split2 = split(/,/,$s);
				$hash{$split2[0]} = 1;
			}
		}
		else {
			die "pb getScafNameFromSgaScafFile ! $_ \n";
		}
	}	
	close SGA;
	return \%hash;
}

sub getLostScafNameFromFasta ($$$) {
	my ($scafFile,$sgaScaffFile,$sgaLostFile) = @_;
	
	my $hash = getScafNameFromSgaScafFile($sgaScaffFile);
	
	open (FASTA, "<$scafFile")
		or die "Cannot open sga scaffold file : $scafFile : $!\n";
	
	open (LOST, ">$sgaLostFile")
		or die "Cannot open sga scaffold file : $sgaLostFile : $!\n";
	
	while (<FASTA>) {
		if ($_ =~ /^>(\S+)\s+/ ) {
			if (!defined $$hash{$1}) {
				print LOST "$1\n";
			}
			else {
				next;
			}
		}
	}	
	close FASTA;
	close LOST;
}

sub links2scaf ($$$$$$$$$$$) {
	my ($connectionsFile,$sgaScaffFile,$scaffoldsFile,$sgaScaffLog,$newScaffFile,$newScaffLog,$sgaLostFile,$scaffLostFile,$allScaffFile,$sgaPath,$getseqPath) = @_;
	my $sga = $sgaPath."sga";
	my $getseq = $getseqPath."getseq";
	
	execCmd("$sga scaffold -v --mate-pair $connectionsFile -o $sgaScaffFile $scaffoldsFile > $sgaScaffLog");
	execCmd("$sga scaffold2fasta -o $newScaffFile -f $scaffoldsFile $sgaScaffFile > $newScaffLog");
	getLostScafNameFromFasta($scaffoldsFile,$sgaScaffFile,$sgaLostFile);	
	execCmd("$getseq -list $sgaLostFile -f $scaffoldsFile -o $scaffLostFile");
	execCmd("cat $newScaffFile $scaffLostFile > $allScaffFile");
}

sub checkMap2qc ($$$) {
	my ($help,$scaffoldsFile,$genomeSize) = @_;
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus map2qc -f assembly.fa -e estimate_size -s tags_coordinates.txt [options] [-h]\n";
	$synopsis .= "\t-f <string>\tassembly.fa: assembly file (FASTA)\n";
	$synopsis .= "\t-e <int>\testimate_size: genome estimate size (bp)\n";
	$synopsis .= "\t-s <string>\ttags_coordinates.txt: sorted file according to mapping position of tags \n\n";
	$synopsis .= "OPTIONAL PARAMETERS:\n";		
	$synopsis .= "\t-p <string>\tprefix for output files (default: magus)\n\n";
	$synopsis .= "EXAMPLE:\n\tmagus map2qc -f Arabidopsis.fa -e 120000000 -s tags_coordinates.txt -p Arabido\n\n";
#	$synopsis .= "This module gives some quality metrics about new assembly.\n";
#	$synopsis .= "It takes in entry magus file prefix_ordered_tags.txt and the old assembly fasta file.\n\n";
	
	if ($help) {
		die $synopsis;
	}
	if ($scaffoldsFile eq "") {
		die "$synopsis";
	}
	if ($genomeSize == 0 ) {
		die "$synopsis";
	} 
}

sub addNmetrics ($$$$$) {
	my ($values,$name,$ref,$log,$file) = @_;
	my $name50 = 0;
	my $name75 = 0;
	my $name90 = 0;
	my $sum = 0;
	my @tab = ();
		
		
		
	foreach ( sort { $b <=> $a } @$values ) {
		$sum += $_;
		foreach my $i (1..100) {	
			if ((!defined $tab[$i-1]) && ($sum > $ref*$i/100)) {
				$tab[$i-1] = $_;
			}
		}	
	}
	
	open (GRAPH, ">$file")
		or die "Cannot open statistics file : $file : $!\n";
	
	foreach my $i (1..100) {
		if (!defined $tab[$i-1]) {
			$tab[$i-1] = 0;
		}		
		print GRAPH "$i $tab[$i-1]\n";
	}
	
	close GRAPH;
	
	print $log "${name}50 size: $tab[49] pb\n";
	print $log "${name}75 size: $tab[74] pb\n";
	print $log "${name}90 size: $tab[89] pb\n";
}

sub makeGraphs ($$$$$) {
	my ($graph,$An,$AnA,$AnG,$RPath) = @_;	
	
	my $R =  $RPath."R";
	
	my $cmd = "echo 'An <- read.csv(file=\"$An\",head=FALSE,sep=\" \")
AnA <- read.csv(file=\"$AnA\",head=FALSE,sep=\" \")
AnG <- read.csv(file=\"$AnG\",head=FALSE,sep=\" \")
png(\"$graph\",width=600,height=450)
plot(AnG,type=\"l\",col=\"green\",main=\"Anchored values\",xlab=\"X\",ylab=\"bp\")
lines(AnA,type=\"l\",col=\"blue\")
lines(An,type=\"l\",col=\"red\")
legend(\"topright\",c(\"An\",\"AnA\",\"AnG\"),col=c(\"red\",\"blue\",\"green\"),lty=1, cex=0.8,)
dev.off()
' | $R --slave > /dev/null"; # redirect because of null device message with dev.off()
	execCmd($cmd); 
}

sub map2qc ($$$$$$$) {
	my ($orderedTagsFile,$scaffoldsFile,$statsFile,$genomeSize,$prefix,$RPath,$fastalengthPath) = @_;
	my ($lastScaff,$lastPos,$lastTag,$lastRank,$lastMol) = ("",-1,"",-1,"");
	my %tagsCounter = ();
	my @anchoredSegLength = ();
	my $isLastGood = 0;
	
	# scaffoldId position tag rank molecule
	open (FILE, "<$orderedTagsFile")
		or die "Cannot open orderedTagsFile : $orderedTagsFile : $!\n";
	
	while (<FILE>) {
		chomp $_;
		
		if ($_ =~ /BAD MAPPPED TAG/) {
			next;
		}
		
		my ($scaff,$position,$tag,$rank,$molecule) = split (/\s+/, $_);
		$tagsCounter{$scaff}{'alignedTags'} += 1;
		
		# initialize
		if (!defined $tagsCounter{$scaff}{'goodPairs'}) {
			$tagsCounter{$scaff}{'goodPairs'} = 0;
		}
		if (!defined $tagsCounter{$scaff}{'validSize'}) {
			$tagsCounter{$scaff}{'validSize'} = 0;
		}
		if (!defined $tagsCounter{$scaff}{'anchorSegLength'}) {
			$tagsCounter{$scaff}{'anchorSegLength'} = 0;
		}
		
		if ($scaff eq $lastScaff) {
			my $tagPair = $lastTag."_".$tag;	
			# tags not mapped have been excluded, so ranks are continuous
			if ( ($rank == $lastRank) || ($rank+1 == $lastRank) || ($rank == $lastRank+1 ) ) {
				my $distancePair = $position - $lastPos;
				$tagsCounter{$scaff}{'goodPairs'} += 1; 
				$tagsCounter{$scaff}{'validSize'} += $distancePair;	
				$tagsCounter{$scaff}{'anchorSegLength'} += $distancePair;	
				$isLastGood = 1;
			}	
			else {
				$tagsCounter{$scaff}{'badPairs'}{$tagPair} = "$scaff $lastPos $position $lastTag $tag $lastRank $rank $molecule";
				$isLastGood = 0;
				if ($tagsCounter{$scaff}{'anchorSegLength'} != 0) {
					push @anchoredSegLength, $tagsCounter{$scaff}{'anchorSegLength'};
				}
				$tagsCounter{$scaff}{'anchorSegLength'} = 0;
			}
		}
		else {
			$isLastGood = 0;
			if ($tagsCounter{$scaff}{'anchorSegLength'} != 0) {
				push @anchoredSegLength, $tagsCounter{$lastScaff}{'anchorSegLength'}; # 0 if only one tag
			}
		}
		# change values for next line
		($lastScaff,$lastPos,$lastTag,$lastRank,$lastMol) = ($scaff,$position,$tag,$rank,$molecule);
	}
	close FILE;
	
	# last value
	if ($isLastGood && ($tagsCounter{$lastScaff}{'anchorSegLength'} != 0)) {
		push @anchoredSegLength, $tagsCounter{$lastScaff}{'anchorSegLength'};
	}
	
	# assembly size
	my $assemblySize = `fastalength $scaffoldsFile |  awk '{ sum+=\$1; } END { print sum; }'`;
	chomp $assemblySize;
	
	# get all scaff and length
	my $scaffoldsLengths = getScaffoldslengths($scaffoldsFile,$fastalengthPath);
	my $anchoredAssemblySize = 0;
	my $anchoredScaffolds = 0;
	my $anchoredScaffoldsLength = 0;
	my $nonAnchoredScaffolds = 0;
	my $nonAnchoredScaffoldsLength = 0;
	my $tagsMapped = 0;
	my $tagsCouples = 0;
	
	foreach my $scaff (keys %$scaffoldsLengths) {
		if (!defined $tagsCounter{$scaff}) {
			$nonAnchoredScaffolds++;
			$nonAnchoredScaffoldsLength += $$scaffoldsLengths{$scaff};
		}
		else {
			$anchoredAssemblySize += $tagsCounter{$scaff}{'validSize'};
			$anchoredScaffolds++;
			$anchoredScaffoldsLength += $$scaffoldsLengths{$scaff};
			$tagsMapped += $tagsCounter{$scaff}{'alignedTags'};
			$tagsCouples += $tagsCounter{$scaff}{'goodPairs'};
		}
	}
	
	my $allScaffolds = $anchoredScaffolds+$nonAnchoredScaffolds;
	
	open (my $handleMetrics, ">$statsFile")
		or die "Cannot open statistics file : $statsFile : $!\n";
	
	print $handleMetrics "\nassembly size: $assemblySize\n";
	print $handleMetrics "TA: amount of tags aligned on assembly : $tagsMapped\n";
	print $handleMetrics "TC: consistent couple of tags : $tagsCouples\n";

	my $plotFile = "${prefix}_quality_metrics.png";
	my $AnFile = "${prefix}_An.csv";
	my $AnAFile = "${prefix}_AnA.csv";
	my $AnGFile = "${prefix}_AnG.csv";

	addNmetrics (\@anchoredSegLength,"An",$anchoredAssemblySize,$handleMetrics,$AnFile); # N Anchored vs all anchored
	addNmetrics (\@anchoredSegLength,"AnA",$assemblySize,$handleMetrics,$AnAFile); # N Anchored vs all assembly
	addNmetrics (\@anchoredSegLength,"AnG",$genomeSize,$handleMetrics,$AnGFile); # N Anchored vs genome size
	
	print $handleMetrics "anchored scaffolds: $anchoredScaffolds (".int($anchoredScaffolds*100/$allScaffolds)."%)\n";
	print $handleMetrics "cumulative anchored scaffolds size: $anchoredScaffoldsLength (".int($anchoredScaffoldsLength*100/$assemblySize)."%)\n";
	print $handleMetrics "number of scaffolds without tags: $nonAnchoredScaffolds (".int($nonAnchoredScaffolds*100/$allScaffolds)."%)\n";
	print $handleMetrics "cumulative scaffolds size without tags: $nonAnchoredScaffoldsLength (".int($nonAnchoredScaffoldsLength*100/$assemblySize)."%)\n";
	print $handleMetrics "mean of tags by scaffold: ".int($tagsMapped/$allScaffolds)."\n\n";
	
	close $handleMetrics;
	
	makeGraphs($plotFile,$AnFile,$AnAFile,$AnGFile,$RPath);
}

sub checkAll ($$$$$$$) {
	my ($help,$mapFile,$tagsFile,$checkedTagsFile,$scaffoldsFile,$readsData,$genomeSize) = @_;
	
	my $synopsis = "\nVersion 1.0.1\n\nUsage: magus all -w wgpFile -t tags.bam -f assembly.fa -e estimate_size -b file.bam,m,sd,s [options] [-h]\n";
	$synopsis .= "\t -w <string>\twgpFile: WGP data\n";
	$synopsis .= "\t -t <string>\ttags.bam: tags alignment on the assembly (BAM)\n";
	$synopsis .= "\t -f <string>\tassembly.fa: assembly file (FASTA)\n";
	$synopsis .= "\t -e <int>\testimate_size: genome estimate size (bp)\n";
	$synopsis .= "\t -b <string>\tfile.bam,m,sd,s: paired reads alignment (BAM), library median size (bp), library standart deviation (bp), reads size (bp)\n\n";;
	$synopsis .= "OPTIONAL PARAMETERS:\n";
	$synopsis .= "\t -s <string>\tsorted file according to mapping position of tags";
	#$synopsis .= "\t(Format : #scaffoldId     position        tagId   rank    contigBac)\n";
	$synopsis .= "\t(default : prefix_tags_coordinates.txt)\n";
	$synopsis .= "\t -a <string>\tanchored tags on assembly";
	#$synopsis .= "\t(Format : #groupID      scaffoldId      minimum rank    maximum rank     number of tags)\n";
	$synopsis .= "\t(default : prefix_anchored_assembly.txt)\n";
	$synopsis .= "\t -l <string>\toutput file containing links between contigs/scaffolds"; 
	#$synopsis .= "\t\t(Format : id1_id2 where id is a scaffold name)\n";
	$synopsis .= "\t(default : prefix_map_links.txt)\n";
	$synopsis .= "\t -c <string>\tfile containing links in DE format";
	$synopsis .= "\t(default : prefix_validated_map_links.de)\n";
	$synopsis .= "\t -m <string>\tBin path\t(default: \$PATH)\n";
	$synopsis .= "\t -v <string>\tpath to samtools\t(default: \$PATH)\n";
	$synopsis .= "\t -r <string>\tpath to R\t(default: \$PATH)\n";
	$synopsis .= "\t -q <string>\tpath to fastalength\t(default: \$PATH)\n";
	$synopsis .= "\t -z <string>\tpath to sga\t(default: \$PATH)\n";
	$synopsis .= "\t -g <string>\tpath to getseq\t(default: \$PATH)\n";
	$synopsis .= "\t -p <string>\tprefix for output files\t(default: magus)\n";
	$synopsis .= "\t -h\t\tthis help\n\n";
	$synopsis .= "EXAMPLE:\n\tmagus all -w tagsWgp.out -t mapping.bam -f Arabidopsis.fa -b mapping_library1.bam,3500,600,101";
	$synopsis .=  " -b mapping_library2.bam,6000,1000,151 -b mapping_library3.bam,9000,1400,251 -e 120000000 \n";
	$synopsis .= "\tmagus all -w tagsWgp.out -t mapping.bam -f Arabidopsis.fa -b mapping_library1.bam,3500,600,101 -e 120000000 -s tags_coordinates.txt -a anchoring_file.txt -l links_file.txt -c links.de";
	$synopsis .=  " -v /env/bin/PathToSamtools -r /env/bin/PathToR -q /env/bin/PathToFstlg -z /env/bin/PathToSga -g /env/bin/PathToGetseq\n";
	$synopsis .= "\tmagus all -w tagsWgp.out -t mapping.bam -f Arabidopsis.fa -b mapping_library1.bam,3500,600,101 -e 120000000";
	$synopsis .=  " -m /env/bin/ -p Arabido\n\n";
	
	if ($help) {
		die $synopsis;
	}
	if ($mapFile eq "") {
		die "$synopsis";
	}
	if ($tagsFile eq "") {
		die "$synopsis";
	}
	if ($scaffoldsFile eq "") {
		die "$synopsis";
	}
	if ($genomeSize == 0 ) {
		die "$synopsis";
	} 
	if ($checkedTagsFile eq "") {
		die "$synopsis";
	}
	if (scalar(@$readsData) == 0 ) {
		die "$synopsis";
	} 
	#checkWgp2map($help,$mapFile,$tagsFile);
	#checkMap2qc($help,$scaffoldsFile,$genomeSize);
	#checkMap2links($help,$checkedTagsFile);
	#checkPairs2links($help,$scaffoldsFile,$readsRef);
	#checkLinks2scaf($help,$scaffoldsFile);
}

return 1;
