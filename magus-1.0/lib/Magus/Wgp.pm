package Magus::Wgp;

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

extends 'Magus::Map';

has '+type' => ( 
				default => 'wgp' 
			);

sub fill ($$) {
	my ($self,$mapFile,$tagsMapped) = @_;
	my $moleculeId = "";
	my @ranks = ();
	my @tags = ();
	my $exclude = 0;
	my $tagIn = 0;
	my $rankIn = 0;
	
	# keep only tags mapped
	my %mapped = ();
	foreach (keys %$tagsMapped) {
		$mapped{$_} = 1;
	}
	
	open (MAP, $mapFile)
		or die "[ERROR] :: fill :: Cannot open $mapFile : $!\n";
		
	while (<MAP>) {
		chomp;
		if ($_ =~ /(\d+)/) {
			$moleculeId = $1;
			$tagIn = 1;
		}
		elsif ($tagIn) {
			@tags = split(/\s+/, $_);
			$tagIn = 0;
			$rankIn = 1;
		}
		elsif ($rankIn) {
			@ranks = split(/\s+/, $_);				
			$self->createMapData($moleculeId,\@tags,\@ranks,\%mapped);
			$rankIn = 0;
		}
		# do not take non placed tags on contig 0
		elsif ($exclude) {
			next;
		}
		else {
			die "[ERROR] :: fill :: parsing map file \n\n $_\n\n";
		}
	}	
}



return 1;