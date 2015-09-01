package Magus::CloneLink;

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
use Magus::Match;

has 'isValid' => (
            is => 'rw',
            isa => 'Bool' 
            );

has 'orientation' => ( # -+ -- ++ +-
            is => 'rw',
            isa => 'Str',
            );

has 'match1' => (
            is => 'rw',
            isa => 'Magus::Match'
            );

has 'match2' => (
	       is => 'rw',
	       isa => 'Magus::Match'
	       );

has 'gap' => (
	       is => 'rw',
	       isa => 'Int'
	       ); 
	       
has 'bankMean' => (
	       is => 'rw',
	       isa => 'Int'
	       ); 
	            
has 'bankDev' => (
	       is => 'rw',
	       isa => 'Int'
	       );
	
has 'readsSize' => (
	       is => 'rw',
	       isa => 'Int'
	       );   	
	
sub addSizes() {
	my ($self) = @_;
	
	if ( $self->match1()->bankMean() == $self->match2()->bankMean() ) {
		$self->bankMean($self->match1()->bankMean());
	}
	else {
		die "[ERROR] read1 and read2 must come from the same bank !\n";
	}
	if ( $self->match1()->bankDev() == $self->match2()->bankDev() ) {
		$self->bankDev($self->match1()->bankDev());
	}
	else {
		die "[ERROR] read1 and read2 must come from the same bank !\n";
	}
	if ( $self->match1()->readsSize() == $self->match2()->readsSize() ) {
		$self->readsSize($self->match1()->readsSize());
	}
	else {
		die "[ERROR] read1 and read2 must come from the same bank !\n";
	}
}	
	
sub checkStatus () {
	my ($self) = @_;
	if ( ($self->gap() > -200) && ( $self->gap() < $self->bankMean()+3*$self->bankDev() ) ) {
		$self->isValid(1);	
	}	
	else {
		$self->isValid(0);
	}
}	
	
sub sizeGap () {
	my ($self) = @_;
	
	if ( (!defined $self->match1()) || !defined $self->match2()) {
		die "[ERROR] sizeGap match1 and match2 need to be defined !\n";
	}
	
	$self->addSizes();
	
	my $pos1 = 0;
	my $pos2 = 0;
	
	my $or = $self->match1()->orientation().$self->match2()->orientation();
	$self->orientation($or);
	
	if ($self->match1()->orientation() eq "+") {
		$pos1 = $self->match1()->position() + $self->readsSize();
	}
	else {
		$pos1 = $self->match1()->length() - $self->match1()->position();
	}
	
	if ($self->match2()->orientation() eq "+") {
		$pos2 = $self->match2()->position() + $self->readsSize();
	}
	else {
		$pos2 = $self->match2()->length() - $self->match2()->position();
	}	

	my $gap = $self->bankMean() - $pos1 - $pos2;
	$self->gap($gap);
	$self->checkStatus();
}		    
	       
return 1;