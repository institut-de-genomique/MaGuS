package Magus::AlignTags;

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
use Moose;
use Magus::MatchBam;

extends 'Magus::AlignmentBam';

has '+matches' => ( 
	       isa => 'HashRef[Magus::MatchBam]', # key is tag id 
	       default => sub { {} } # initialize it   
	       );	       

has 'mapLinks' => (
		   is => 'rw',
	       isa => 'ArrayRef[Str]',  
	       default => sub { [] } # initialize it   
	       ); 

sub fill ($) { 
	my ($self,$bamFile) = @_;
	my %hash =();
	my $samtoolsPath = $self->binPath()."samtools";
	# do not keep segment unmapped
	open (SAM, "$samtoolsPath view -F 4 $bamFile |")
		or die "[ERROR] :: fill :: Cannot open $bamFile : $!\n";
		
	while (<SAM>) {
		my $sam = $_;
		chomp;
		# keep only single hit
		if ($_ =~ /\tX0:i:1\t/ ) {
			
			my @data = split(/\t/,$_);
	
			# create Match set
			my $scaffoldId = $data[2];
			my $tag = $data[0];
			my $position = $data[3];
			
			if ($scaffoldId =~ /_/) {
				die "[DIE] Fasta header can not contains \'_\' character !!!\n";
			}
						
			my $obj = new Magus::MatchBam();
			$obj->scaffoldId($scaffoldId);
			$obj->tag($tag);
			$obj->position($position);
			#$obj->sam($sam); not need at this time
			
			$hash{$tag} = $obj;
						
		}		
	}	
	close SAM;
	$self->matches(\%hash);
}

sub getBadTags ($$) {
	my ($self,$order,$orderedTagsFile) = @_;
	my @badTags = ();
	
	open (ORDER, ">$orderedTagsFile")
		or die "Cannot open $orderedTagsFile : $!\n";
		
	print ORDER "#scaffoldId\tposition\ttagId\trank\tcontigBac\n";
	
	foreach my $scaff (keys %$order) {
		my $hash = $$order{$scaff};
		my $last = "";
		my $lastLine = "";
		my $single = "";
		my $singleLine = "";
		
		foreach my $position (sort {$a <=> $b} keys %$hash) { 
			my $molecule = $$hash{$position}{'molecule'};
			my $tagId = $$hash{$position}{'tag'};
			# scaffoldId position tag rank molecule
			$lastLine = "$scaff\t$position\t$tagId\t$$hash{$position}{'rank'}\t$molecule\n";
			
			if ($last eq "") {
				print ORDER $lastLine;
				$last = $molecule;
				next;
			}
			
		  	if ($last eq $molecule) { # same
		  		print ORDER $lastLine;
		  		if ($single ne "") {
		  			print ORDER "[BAD MAPPPED TAG] : $singleLine";
		  			my @tab = split(/\t/,$singleLine);
		  			push @badTags,  $tab[2];
		  			$single = "";
		  		}	
		  	}
		  	else {
		  		if ($single eq "") { # do not know
		  			$single = $molecule;
		  			$singleLine = $lastLine;
		  		}
		  		elsif ($single eq $molecule) { # change
		  			print ORDER $lastLine;
		  			print ORDER $singleLine;
		  			$single = "";
		  			$last = $molecule;
		  		}	
		  		else { # exclude
		  			print ORDER $lastLine;
		  			print ORDER "[BAD MAPPPED TAG] : $singleLine";
		  			my @tab = split(/\t/,$singleLine);
		  			push @badTags, $tab[2];
		  			$single = $molecule;
		  		}
		  	}
			
		}
	}
	close ORDER;
	return \@badTags;
}

sub addTagsData ($$) {
	my ($self,$tagsMap,$orderedTagsFile) = @_;
	my %order = ();
	
	foreach my $match (values %{$self->matches()}) {
		my $tagId = $match->tag();
		my $scaffoldId = $match->scaffoldId();
		
		# add rank and molecule
		if (defined $$tagsMap{$tagId}) {			
			my $rank = $$tagsMap{$tagId}->rank();
			my $moleculeId = $$tagsMap{$tagId}->moleculeId();
			$match->rank($rank);
			$match->moleculeId($moleculeId);
			
			# for quality statistics
			$$tagsMap{$tagId}->isMapped(1);
			
			# get order on scaffold
			my $position = $match->position();
			
			# keep only single tag for one position
			if (defined $order{$scaffoldId}{$position}) {
				delete ${$self->matches()}{$tagId};
				my $ref = $order{$scaffoldId}{$position};
				if (defined ${$self->matches()}{$ref}) {
					delete ${$self->matches()}{$ref};
				}
			}
			else {
				$order{$scaffoldId}{$position}{'tag'} = $tagId;
				$order{$scaffoldId}{$position}{'rank'} = $rank;
				$order{$scaffoldId}{$position}{'molecule'} = $moleculeId;
			}
		}
		else {
			delete ${$self->matches()}{$tagId};
		}		
	}
	
	# clean bad mapped tags
	my $badTagsRef = $self->getBadTags(\%order,$orderedTagsFile);		
	foreach my $tagId (@$badTagsRef) {
		delete ${$self->matches()}{$tagId};
	}
}

sub getAnchors () {
	my ($self) = @_;
	my %anchor = ();	
	
	foreach my $match (values %{$self->matches()}) {
		my $scaffoldId = $match->scaffoldId();
		my $moleculeId = $match->moleculeId();
		my $rank = $match->rank();
		
		# count number of links
		if (defined $anchor{$moleculeId}{$scaffoldId}{'tags'}) {
			$anchor{$moleculeId}{$scaffoldId}{'tags'} += 1;
		}
		else {
			$anchor{$moleculeId}{$scaffoldId}{'tags'} = 1;
		}
		# get minimum rank
		if (defined $anchor{$moleculeId}{$scaffoldId}{'min'}) {
			if ($rank < $anchor{$moleculeId}{$scaffoldId}{'min'}) {
				$anchor{$moleculeId}{$scaffoldId}{'min'} = $rank;
			}
		}
		else {
			$anchor{$moleculeId}{$scaffoldId}{'min'} = $rank;
		} 
		# get maximum rank
		if (defined $anchor{$moleculeId}{$scaffoldId}{'max'}) {
			if ($rank > $anchor{$moleculeId}{$scaffoldId}{'max'}) {
				$anchor{$moleculeId}{$scaffoldId}{'max'} = $rank;
			}
		}
		else {
			$anchor{$moleculeId}{$scaffoldId}{'max'} = $rank;
		} 	
	}
	return \%anchor;
}



sub printAnchors ($) {
	my ($self,$anchorsFile) = @_;
	
	open (ANCHORS, ">$anchorsFile")
		or die "Cannot open $anchorsFile : $!\n";
		
	print ANCHORS "#contigBac\tscaffoldId\tminimum rank\tmaximum rank\t number of tags\n";	
			
	my $anchor = $self->getAnchors();
		
	foreach my $moleculeId (sort {$a <=> $b} keys %$anchor) { 
		my $hash = $$anchor{$moleculeId};
		foreach my $scaff (sort { $$hash{$a}{'min'} <=> $$hash{$b}{'min'}}  keys %$hash) {
			print ANCHORS "$moleculeId\t$scaff\t$$hash{$scaff}{'min'}\t$$hash{$scaff}{'max'}\t$$hash{$scaff}{'tags'}\n";
		}
	}
	close ANCHORS;
}




sub saveLinks ($) {
	my ($self,$links) = @_;		
	my $size = $#$links;
	if ($size > 0) { # more than one scaffold
		for (my $i = 0; $i <= $size-1; $i++) {
			if ($$links[$i] ne $$links[$i+1]) {			
				# push 
				push @{$self->mapLinks()}, "$$links[$i]_$$links[$i+1]";
			}
		}
	}
}

sub createLinks ($) {
	my ($self,$anchorsFile) = @_;
	
	open (ANCHORS, "<$anchorsFile")
		or die "Cannot open $anchorsFile : $!\n";
	
	my $lastMoleculeId = "";
	my $id = "";
	my $min = 0;
	my $max = 0;
	my @links = ();
	
	##  list putative links 
	while (<ANCHORS>) {
		chomp $_;
		
		if ($_ =~ /^#/) {
			next; # exclude description line
		}
		
		my ($moleculeId,$scaffId,$scaffMin,$scaffMax,$nbLinks) = split (/\s+/,$_);	
				
		## first
		if ($lastMoleculeId eq "") {
			#initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
			next;
		}
		
		if ($moleculeId != $lastMoleculeId) {
			push @links, $id;
			# save data
			$self->saveLinks(\@links);
			# delete old data
			@links = ();
			# re-initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
		}
		else {	
			push @links, $id; # simple approach
		}
	}
	close ANCHORS;
}

sub createLinksNotIncluded ($) {
	my ($self,$anchorsFile) = @_;
	
	open (ANCHORS, "<$anchorsFile")
		or die "Cannot open $anchorsFile : $!\n";
	
	my $lastMoleculeId = "";
	my $id = "";
	my $min = 0;
	my $max = 0;
	my @links = ();
	
	##  list putative links 
	while (<ANCHORS>) {
		chomp $_;
		
		if ($_ =~ /^#/) {
			next; # exclude description line
		}
		
		my ($moleculeId,$scaffId,$scaffMin,$scaffMax,$nbLinks) = split (/\s+/,$_);	
			
		## first
		if ($lastMoleculeId eq "") {
			#initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
			next;
		}
		
		if ($moleculeId != $lastMoleculeId) {
			push @links, $id;
			# save data
			$self->saveLinks(\@links);
			# delete old data
			@links = ();
			# re-initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
		}
		else {			
			if ( $scaffMin > $min ) {
				if ($scaffMax > $max) { # not included
					push @links, $id;
					$id = $scaffId;
					$min = $scaffMin;
					$max = $scaffMax;
				}
			}
			else {	# same min
				if ($scaffMax > $max) { # keep max end
					$max = $scaffMax;
					$id = $scaffId;
				}
			}
		}
	}
	close ANCHORS;
}

sub createLinksNoNoise ($) {
	my ($self,$anchorsFile) = @_;
	
	open (ANCHORS, "<$anchorsFile")
		or die "Cannot open $anchorsFile : $!\n";

	my $lastMoleculeId = "";
	my $id = "";
	my $min = 0;
	my $max = 0;
	my @links = ();
	
	##  list putative links 
	while (<ANCHORS>) {
		chomp $_;
		
		if ($_ =~ /^#/) {
			next; # exclude description line
		}
		
		my ($moleculeId,$scaffId,$scaffMin,$scaffMax,$nbLinks) = split (/\s+/,$_);	
		
		# remove noise
		if ($scaffMax-$scaffMin >= $nbLinks) { 
			next;
		}
				
		## first
		if ($lastMoleculeId eq "") {
			#initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
			next;
		}
		
		if ($moleculeId != $lastMoleculeId) {
			push @links, $id;
			# save data
			$self->saveLinks(\@links);
			# delete old data
			@links = ();
			# re-initialize
			$lastMoleculeId = $moleculeId;
			$id = $scaffId;
			$min = $scaffMin;
			$max = $scaffMax;
		}
		else {			
			if ( $scaffMin > $min ) { 
				push @links, $id;
				$id = $scaffId;
				$min = $scaffMin;
				$max = $scaffMax;
			}
			else {	# same min
				if ($scaffMax > $max) { # keep max end
					$max = $scaffMax;
					$id = $scaffId;
				}
			}
		}
	}
	close ANCHORS;
}

sub printLinks ($) {
	my ($self,$handle) = @_;
	my %hash = ();
	foreach (@{$self->mapLinks()}) {
		if (!defined $hash{$_}) { # only once
			$hash{$_} = 1;
			print $handle "$_\n";
		}	
	}
}

sub printData () {
	print "NOT YET IMPLENTED \n";
}

return 1;