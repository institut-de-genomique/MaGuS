package Magus::Connection;

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
use Magus::CloneLink;

has 'orientation'  => ( # -+ -- ++ +-
            is => 'rw',
            isa => 'Str',
            );

has 'scaffoldId1' => ( 
            is => 'rw',
            isa => 'Str',
            );
            
has 'scaffoldId2' => (
            is => 'rw',
            isa => 'Str',
            );
            
has 'links' => (
            is => 'rw',
            isa => 'ArrayRef[Magus::CloneLink]',
            default => sub { [] } # initialize it   
            );

has 'validLinksNb' => (
	       is => 'rw',
	       isa => 'Int'
	       );

has 'gapMean' => (
	       is => 'rw',
	       isa => 'Int'
	       );
	            
has 'gapStdDev' => (
	       is => 'rw',
	       isa => 'Int'
	       );	
	
	
sub getMean ($) {
	my ($self,$ref) = @_;
  	my $sum = 0;
  	my $sqrdif=0;
  	foreach (@$ref) {
  		$sum += $_;
  	}
  	my $mean = $sum/scalar(@$ref);
	return (int($mean));
}


sub getStats ($) {
	my ($self,$ref) = @_;
   	if (scalar(@$ref) == 0 ) {
   		warn "[WARN] : you need one value to have a mean !\n";
   		return "0","0";
   	}
   	my $sqrdif=0;
	my $mean = $self->getMean($ref);
	foreach (@$ref){
    	$sqrdif+= ($mean-$_)**2;
   	}
 	my $sd = sqrt($sqrdif/scalar(@$ref));
 	return $mean,int($sd);
} 
	
sub sizeGap () {
	my ($self) = @_;
	# mean if gap and standard dev
	my @gapSizes = ();
	my $validLinksNb = 0;
	my $filename = '';
	#open (LINKS_DETAILS, ">>links_details.txt");
	open (LINKS_DETAILS, ">>output_coord_links/$filename");
	foreach my $link ( @{$self->links()} ) {		
		if ($link->isValid()) {
			$validLinksNb++;
			$filename = $self->scaffoldId1."_".$self->scaffoldId2."_".$self->orientation();
			open (LINKS_DETAILS, ">>output_coord_links/$filename");
			print LINKS_DETAILS $self->scaffoldId1,"\t",$self->scaffoldId2,"\t",$link->match1()->length(),"\t",$link->match2()->length(),"\t";
			if ($self->orientation() eq "-+") {
				print LINKS_DETAILS "+\t+\t";
			}	
			if ($self->orientation() eq "--") {
				print LINKS_DETAILS "+\t-\t";
			}
			if ($self->orientation() eq "+-") {
				print LINKS_DETAILS "-\t-\t";
			}
			if ($self->orientation() eq "++") {
				print LINKS_DETAILS "-\t+\t";
			}
			print LINKS_DETAILS $link->match1()->position(),"\t",$link->match2()->position(),"\t",$link->gap(),"\n";

			push @gapSizes, $link->gap();
		}
		close (LINKS_DETAILS);
	}
	$self->validLinksNb($validLinksNb);
	if (scalar(@gapSizes) != 0) {
		my ($mean,$sd) = $self->getStats(\@gapSizes);
		$self->gapMean($mean);
		$self->gapStdDev($sd);	
	}	
	else {
		$self->gapMean(0);
		$self->gapStdDev(0);	
	}

}	
	    
sub printDeFormat ($) {
	my ($self,$file) = @_;
	open (CONNECTIONS, ">>$file")
		or die "[ERROR] :: printDeFormat :: Cannot open $file : $!\n";
	if ($self->orientation() eq "-+") {
		print CONNECTIONS $self->scaffoldId1." ".$self->scaffoldId2()."+,".$self->gapMean().",".$self->validLinksNb().",".$self->gapStdDev()."\n";
	}	
	if ($self->orientation() eq "--") {
		print CONNECTIONS $self->scaffoldId1." ".$self->scaffoldId2()."-,".$self->gapMean().",".$self->validLinksNb().",".$self->gapStdDev()."\n";
	}
	if ($self->orientation() eq "+-") {
		print CONNECTIONS $self->scaffoldId1." ; ".$self->scaffoldId2()."+,".$self->gapMean().",".$self->validLinksNb().",".$self->gapStdDev()."\n";
	}
	if ($self->orientation() eq "++") {
		print CONNECTIONS $self->scaffoldId1." ; ".$self->scaffoldId2()."-,".$self->gapMean().",".$self->validLinksNb().",".$self->gapStdDev()."\n";
	}
}	    
	       
return 1;
