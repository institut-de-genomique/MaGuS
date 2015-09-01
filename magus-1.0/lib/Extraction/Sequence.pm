package Extraction::Sequence;

use strict;
use warnings;

sub new {
    my ( $proto, $name, $seq, $phase, $hash_ref ) = @_;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    if ( !defined $phase ) { $phase = 1; }
    if ( $phase != 1 && $phase != -1 ) {
        die "[Sequence] Error : - new - Unknown phase : $phase.\n";
    }
    $self->{NAME}     = $name;
    $self->{SEQ}      = $$seq;
    $self->{HASH_REF} = $hash_ref;
    $self->{PHASE}    = $phase;
    $self->{_ARRAY}   = undef;
    return $self;
}

# Virtual function to implement in sub-classes
sub array {
    my $self = shift;
    warn "[Sequence] Error - array - : no array method defined in sub class\n";
    die;
}

sub string {
    my $self = shift;
    warn
      "[Sequence] Error - string - : no string method defined in sub class\n";
    die;
}

sub elemEqual {
    my $self = shift;
    warn
      "[Sequence] Error - string - : no string method defined in sub class\n";
    die;
}

sub getComplementary {
    my $self = shift;
    warn
      "[Sequence] Error - string - : no string method defined in sub class\n";
    die;
}

sub _cacheArray {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) { $self->{_ARRAY} = $arg; }
    else                { return $self->{_ARRAY}; }
}

sub name {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) { $self->{NAME} = $arg; }
    else                { return $self->{NAME}; }
}

sub frame {
    my ($self, $arg) = @_;
    if(defined $arg) { $self->{PHASE} = $arg; }
    else { return $self->{PHASE}; }
}

sub seq {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) { $self->{SEQ} = $arg; }
    else                { return $self->{SEQ}; }
}

sub getLength {
    my $self     = shift;
    my $refarray = $self->array();
    return scalar(@$refarray);
}

sub _verifDebFin {
    my ( $self, $deb, $fin, $length ) = @_;
    if ( !defined $length ) { $length = $self->getLength(); }
    if ( $deb < 0 ) {
        warn "[Sequence] Warning : - verifDebFin - ", $self->name(),
          " : start position adjusted to 0 : $deb.\n";
        $deb = 0;
    }
    if ( $fin < 0 ) {
        warn "[Sequence] Warning : - verifDebFin - ", $self->name(),
          " : end position adjusted to $length : $fin.\n";
        $fin = $length;
    }
    if ( $deb > $length ) {
        warn "[Sequence] Warning : - verifDebFin - ", $self->name(),
          " : start position adjusted to ", $length, " : $deb.\n";
        $deb = $length;
    }
    if ( $fin > $length ) {
        warn "[Sequence] Warning : - verifDebFin - ", $self->name(),
          " : end position adujsted to ", $length, " : $fin.\n";
        $fin = $length;
    }
    return ( $deb, $fin );
}

sub getPortion {
    my ( $self, $deb, $fin ) = @_;
    $deb--;
    my ( $debT, $finT ) = $self->_verifDebFin( $deb, $fin );
    my $refarray        = undef;
    if ( defined $self->_cacheArray() ) { $refarray = $self->_cacheArray(); }
    else { $refarray = $self->array(); $self->_cacheArray($refarray); }
    my @slice = splice( @$refarray, $debT, $finT - $debT );

    # Correction BUG - jmaury AT genoscope.cns.fr
    # 20/01/2009 - splice modifie le tableau, il faut donc le restaurer
    # par la suite. On fait un 2eme splice avec comme taille 0.
    splice( @$refarray, $debT, 0, @slice );
    return ( $debT, $finT, $self->string( \@slice ) );
}

sub getFrag {
    my ( $self, $deb, $fin, $softname ) = @_;
    if ( !defined $softname ) { $softname = 1; }
    my ( $newDeb, $newFin, $seq ) = $self->getPortion( $deb, $fin );
    my $newName;
    if ($softname) { $newName = $self->name(); }
    else {
        $newName = $self->name()
          . " extract from $deb($newDeb) to $fin($newFin) (length ";
        $newName .= ( $fin - $deb + 1 ) . ").";
    }
    return $self->new( $newName, \$seq, $self->{PHASE} );
}

sub _getRev {
    my ($self) = @_;
    return CORE::reverse $self->seq();
}

sub formatSeq {
    my ( $self, $colonne, $prefix ) = @_;
    if ( !defined $prefix ) { $prefix = ""; }
    return $self->_fct_fold( $prefix, $colonne );
}

sub formatSeq2fh {
    my ( $self, $fh, $colonne, $prefix ) = @_;
    if ( !defined $prefix ) { $prefix = ""; }
    $self->_fct_fold_fh( $fh, $prefix, $colonne );
}

sub getOcc {
    my ( $self, $elem ) = @_;
    my @pos      = ();
    my $refarray = $self->array();
    for ( my $cpt = 0 ; $cpt < scalar(@$refarray) ; $cpt++ ) {
        if ( $self->elemEqual( $refarray->[$cpt], $elem ) ) {
            push( @pos, $cpt + 1 );
        }
    }
    return \@pos;
}

sub replace {
    my ( $self, $newval, $pos ) = @_;
    my $refarray                = $self->array();
    foreach my $i (@$pos) {
        if ( $i - 1 < scalar(@$refarray) ) { $refarray->[ $i - 1 ] = $newval; }
    }
    $self->seq( $self->string($refarray) );
}

# Compatibility with previous version of this library
sub FormatSeq {
    my ( $self, $colonne ) = @_;
    warn
"[Sequence] Warning - FormatSeq - : Deprecated method, use formatSeq instead of FormatSeq\n";
    return $self->formatSeq($colonne);
}

sub IUPACtoN {
    my $self    = shift;
    my $newSeq  = $self->{SEQ};
    $newSeq     =~ tr/RYMKWSBDHVrymkwsbdhv/N/;
    return $newSeq;
}


1;
