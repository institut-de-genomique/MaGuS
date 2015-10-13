package Magus::AlignReadsMP;

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
use Magus::CloneLink;

extends 'Magus::AlignmentBam';

has '+matches' => ( # surcharge
	       isa => 'HashRef[ArrayRef[Magus::MatchBam]]', # key is read id 
	       default => sub { {} } # initialize it   
	       );    
       
has 'matchesRead1' => ( 
		   is => 'rw',
	       isa => 'HashRef[ArrayRef[Magus::MatchBam]]', # key is read id 
	       default => sub { {}} # initialize it   
	       );  

has 'matchesRead2' => ( 
		   is => 'rw',
	       isa => 'HashRef[ArrayRef[Magus::MatchBam]]', # key is read id 
	       default => sub { {} } # initialize it   
	       );              
       
has 'readsByScaffold' => ( 
	       is => 'rw',
	       isa => 'HashRef[HashRef[Int]]', # key is scaffoldId value is template
	       default => sub { {} } # initialize it  
	       );
	       
has 'scaffoldsByRead' => ( 
	       is => 'rw',
	       isa => 'HashRef[HashRef[Int]]', # key is template value is scaffoldId
	       default => sub { {} } # initialize it  
	       );	       

has 'samByLink' => ( 
	       is => 'rw',
	       isa => 'HashRef[HashRef[Str]]', # key is scaffoldIds then bank data
	       default => sub { {} } # initialize it  
	       );	       

sub cleanMatches() {
	my ($self) = @_;	
	foreach my $key (keys %{$self->matchesRead1()}) {
		delete ${$self->matchesRead1()}{$key};
	}
	foreach my $key (keys %{$self->matchesRead2()}) {
		delete ${$self->matchesRead2()}{$key};
	}
}

sub fill ($$$$) {
	my ($self,$readsFile,$bankMean,$bankDev,$readsSize) = @_;
	my $samtoolsPath = $self->binPath()."samtools";
	my $bank = 	$bankMean."_".$bankDev."_".$readsSize;
	
	open (SAM, "$samtoolsPath view -F 4 $readsFile |")
		or die "[ERROR] :: fill :: Cannot open $readsFile : $!\n";
	
	# DEBUG open (SAM, $readsFile)
	# DEBUG		or die "[ERROR] :: fill :: Cannot open $readsFile : $!\n";
			
	while (<SAM>) {
		my $sam = $_;
		chomp;
		my @data = split(/\t/,$_);
		if ( ($data[6] eq "=") || ($data[6] eq "*") ) { # same scaffold or undefined
			next; 
		}
		my $scaffold1 = $data[2];	
		my $scaffold2 = $data[6];
		
		if ( ($scaffold1 =~ /_/) || ($scaffold2 =~ /_/) ) {
				die "[DIE] Fasta header can not contains \'_\' character !!!\n";
		}
		
		my $key = $scaffold1."_".$scaffold2;
		push @{${$self->samByLink()}{$key}{$bank}}, $sam; # memory !!!		
	}
}

sub createOneLinkMatches ($$$) {
	my ($self,$scaff1,$scaff2) = @_;
	my $hash = $self->samByLink();
	my $key = $scaff1."_".$scaff2;	
	
	foreach (keys %{$$hash{$key}}) {		
		my ($bankMean,$bankDev,$readsSize) = split(/_/,$_);
	
		foreach my $sam (@{$$hash{$key}{$_}}) {			
			chomp $sam;
			my @data = split(/\t/,$sam);
			my $template = $data[0];
			my $scaffold1 = $data[2];
			my $position1 = $data[3];
			my $segment = ( ($data[1] & 64) != 64 ) ? 1 : 2; # read1 or read2
			my $orientation1 = ( ($data[1] & 16) != 16 ) ? "+" : "-"; # reverse complemented

			my $obj = new Magus::MatchBam();
			$obj->scaffoldId($scaffold1);
			$obj->position($position1);
			$obj->orientation($orientation1);
			$obj->template($template);
			$obj->segment($segment);
			$obj->bankMean($bankMean);
			$obj->bankDev($bankDev);
			$obj->readsSize($readsSize);
			#$obj->sam($sam);
		
			# store data in object			
			if ($segment == 1) {
				push @{${$self->matchesRead1()}{$template}}, $obj;
			}
			elsif ($segment == 2) {
				push @{${$self->matchesRead2()}{$template}}, $obj;
			}
			else {
				die "Cannot recognize read number 1 or 2 ? $segment \n";
			}		
		}
	}
}

sub createCloneLinks ($$$$) {
	my ($self,$scaff1,$scaff2,$length1,$length2) = @_;
	my @links = ();
	$self->createOneLinkMatches($scaff1,$scaff2);
	$self->createOneLinkMatches($scaff2,$scaff1); # other direction 
	my $matchesByRead1 = $self->matchesRead1();
	my $matchesByRead2 = $self->matchesRead2();
	
	foreach my $template (keys %$matchesByRead1) {	
		
		my $link = new Magus::CloneLink();
		
		foreach my $m1 (@{$$matchesByRead1{$template}}) {	
			foreach my $m2 (@{$$matchesByRead2{$template}}) {		
				my $link = new Magus::CloneLink();
			
				if ( $m1->scaffoldId() eq $m2->scaffoldId() ) {
					next;
				}					
				elsif ( ($m1->scaffoldId() eq $scaff1) && ($m2->scaffoldId() eq $scaff2) ) {
					$m1->length($length1);
					$m2->length($length2);
					$link->match1($m1);
					$link->match2($m2);
					$link->sizeGap();
					push @links, $link;		
				} 
				elsif ( ($m1->scaffoldId() eq $scaff2) && ($m2->scaffoldId() eq $scaff1) ) {
					$m1->length($length2);
					$m2->length($length1);				
					$link->match1($m2);
					$link->match2($m1);
					$link->sizeGap();
					push @links, $link;
				}
				else {
					die "[ERROR] ".$m1->scaffoldId()." and ".$m2->scaffoldId()." not equal to $scaff1 or $scaff2 !\n";
 				}	
			}
		}
	}	
	
	$self->cleanMatches();
	return \@links;	
}

sub printData () {
	print "NOT YET IMPLENTED \n";
}

return 1;