package Extraction::QSequence;
use warnings;
use strict;
use Text::Wrap;

use Extraction::Sequence;

if( $] >= 5.010000 ){
    use base( "Extraction::Sequence" );
}
else {
    our @ISA = ( "Extraction::Sequence" );
}


sub new {
    my ( $proto, $name, $seq ) = @_;
    my $class = ref($proto) || $proto;
    if ( !ref($seq) ) {
        return new( $class, $name, \$seq );
    }
    $$seq =~ s/\n/ /g;
    $$seq =~ s/^\s*//g;
    my $self = new Extraction::Sequence( $name, $seq );
    bless $self, $class;
    return $self;
}

sub array {
    my $self = shift;
    my @array = split( /\s+/, $self->seq() );
    return \@array;
}

sub string {
    my ( $self, $refarray ) = @_;
    return join( " ", @$refarray );
}

sub elemEqual {
    my ( $self, $elem1, $elem2 ) = @_;
    return ( $elem1 == $elem2 );
}

sub replace {
    my ( $self, $newval, $pos ) = @_;
    my @qualtab = split( /\s+/, $self->seq() );
    foreach (@$pos) {
        if ( $_ - 1 < scalar(@qualtab) ) { $qualtab[ $_ - 1 ] = $newval; }
    }
    $self->seq( join( " ", @qualtab ) );
}

# Formate une chaine sur long_lign colonnes
#-------------------------------------------
sub _fct_fold {
    my ( $self, $prefix, $long_lign ) = @_;

    if ( !defined $long_lign ) { $long_lign = 20; }

    my $i   = 0;
    my $res = "";

    my @qualtab = split( /\s+/, $self->seq() );

    my ( $cpt, $first ) = ( -1, 1 );
    for ( $i = 0 ; $i < scalar(@qualtab) ; $i++ ) {
        if ( $cpt == $long_lign ) {
            $res .= "\n" . $prefix . sprintf( "%2s", $qualtab[$i] );
            $cpt = -1;
        }
        else {
            if ( !$first ) { $res .= " " . sprintf( "%2s", $qualtab[$i] ); }
            else { $res .= sprintf( "%2s", $qualtab[$i] ); $first = 0; }
        }
        $cpt++;
    }
    return $res . "\n";
}

1;
