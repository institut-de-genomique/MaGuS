package Magus::Map;

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
use Data::Dumper;
use Magus::Tag;

has 'type' => (
            is => 'rw',
            isa => 'Str',            							
			default => 'map'
		);
	          
has 'tags'	=> (
			is => 'rw',
	       	isa => 'HashRef[Magus::Tag]', # key is tagId 
	       	default => sub { {} } # initialize it
			);	  
	 
sub createMapData ($$$$) {
	my ($self,$moleculeId,$tagsRef,$ranksRef,$tagsMapped) = @_;
	my $i = 0;
	my $rankCount = 0;
	my $lastRank = 0;
	my $hash = $self->tags();
	
	if ( scalar(@$tagsRef) != scalar(@$ranksRef) ) {
		die "[ERROR] :: createMapData :: parsing map data : number of tags and ranks are different ! \n";
	}
	
	foreach my $tagId (@$tagsRef) {		
		# if not mapped next
		if (!defined $$tagsMapped{$tagId}) {
			$i++;
			next;
		}
		
		my $rank = $$ranksRef[$i];
		#convert rank into integer
		if ($rank != $lastRank) { 
			$rankCount++;
		}
		# create Tag object
		my $obj = new Magus::Tag();
		$obj->moleculeId($moleculeId);
		$obj->tag($tagId);
		$obj->rank($rankCount);
		
		$$hash{$tagId} = $obj;
		
		$lastRank = $rank;		
		$i++;
	} 
	$self->tags($hash);
}

sub printMapData() {
	my ($self) = @_;
	my %h = keys(%{$self->molecules()});
	foreach my $key (sort keys(%h) ) {
			my $obj = $h{$key}; 
			$obj->printTags();
	}
}

return 1;	       