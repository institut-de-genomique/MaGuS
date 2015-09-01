################################################################################
# *
# * ExtraSeq::FSequence.pm
# *
# * This software is a computer program whose purpose is to do stuff
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
# * requirements in conditions enabling the security of their systems
# * and/or data to be ensured and,  more generally, to use and operate it
# * in the same conditions as regards security.
# *
# * The fact that you are presently reading this means that you have had
# * knowledge of the CeCILL license and that you accept its terms.
# *
################################################################################

=head1 NAME

ExtraSeq::FSequence - Class for the a dna sequence

=head1 SYNOPSIS

use ExtraSeq::FSequence;

=over

=item # Create a new dna sequence object

 $seq = new ExtraSeq::FSequence("my_sequence", "acgatcgatcgt");

=item # Load annotations wich are reference data

 $ens = new Compare::EnsSegElem($file_annot, $field_nom_annot, $field_deb_annot,
                                $field_fin_annot, $field_strand_annot, $filtre_annot);
 $ens->loadFile();
 $kpp->annotations($ens);

=back

=head1 DESCRIPTION

The FSequence class supply methods to manipulate dna sequences
This class is included in ExtraSeq package, and inherit ExtraSeq::Sequence class.

=cut

package Extraction::FSequence;

use strict;
use warnings;
use Text::Wrap;
use Exporter;

if( $] >= 5.010000 ){
    use base( "Extraction::Sequence" );
}
else {
    our @ISA = ( "Extraction::Sequence" );
}

=head2 new

    Title   : new
    Usage   : my $seq = Extraction::FSequence("my_sequence", "atctatctatgatcttcgagtggcta");
    Function: Instantiate a new FSequence object
    Returns : Extraction::FSequence object
    Args    : * string : sequence name
              * string : dna sequence
=cut

sub new {
    my ( $proto, $name, $seq ) = @_;
    my $class = ref($proto) || $proto;
    if ( !ref($seq) ) {
        return new( $class, $name, \$seq );
    }
    $$seq       =~ s/\n//g;
    my $self    = new Extraction::Sequence( $name, $seq );
    bless $self, $class;
    return $self;
}

=head2 array

    Title   : array
    Usage   : my $array = $seq->array();
    Function: return a listref of the dna sequence
    Returns : a listref which which contain the dna sequence, one base per cell
    Args    : none
=cut

sub array {
    my $self = shift;
    my @array = split( //, $self->seq() );
    return \@array;
}

=head2 array

    Title   : string
    Usage   : my $str = $seq->string(\@my_array);
    Function: transform the given list in a string
    Returns : string
    Args    : a listref
=cut

sub string {
    my ( $self, $refarray ) = @_;
    return join( "", @$refarray );
}

sub elemEqual {
    my ( $self, $elem1, $elem2 ) = @_;
    return ( $elem1 eq $elem2 );
}

sub getLength {
    my $self = shift;
    return length( $self->seq() );
}

sub getPortion {
    my ( $self, $deb, $fin ) = @_;
    my $len = $self->getLength();
    $deb--;
    my ( $debT, $finT ) = $self->_verifDebFin( $deb, $fin, $len );
    return ( $debT, $finT, substr( $self->seq(), $debT, $finT - $debT ) );
}

sub getGC {
    my $self = shift;
    my $seq  = $self->seq();
    my ( $pgc, $nbw ) = ( 0, 0 );

    my $len = $self->getLength();
    my @occ = ( $seq =~ /[GC]/ig );

    return ( ( scalar(@occ) / $len ) * 100 );
}

sub getWinGC {
    my ( $self, $win, $pas ) = @_;
    if ( !defined $win ) { $win = 100; }
    if ( !defined $pas ) { $pas = 1; }
    my $seq = $self->seq();
    my @pgc;

    my $len = $self->getLength();
    if ( $win > $len ) {
        warn
"[Sequence] Warning : - getGC - The window size is greater than the sequence length\n";
        return ();
    }

    my $pos = $len - $win;
    for ( my $i = 0 ; $i < $pos ; $i += $pas ) {
        my $frag = substr( $seq, $i, $win );
        my @occ_gc = ( $frag =~ /[GC]/ig );
        my @occ_at = ( $frag =~ /[AT]/ig );
        if ( ( scalar(@occ_gc) + scalar(@occ_at) ) / $win < 0.75 ) { next; }
        push @pgc, $i + 1 . ":" . ( scalar(@occ_gc) / $win ) * 100;
    }

    return @pgc;
}

sub mask {
    my ( $self, $deb, $fin, $char ) = @_;
    if ( !defined $char ) { $char = "N"; }
    my ( $start, $stop ) = $self->_verifDebFin( $deb, $fin );
    my $long = $stop - $start + 1;
    my $mask = $char x $long;
    my $s    = $self->seq();
    substr( $s, $start - 1, $long ) = $mask;
    $self->seq($s);
    return 1;
}

sub getComplementary {
    my $self = shift;
    my $seq  = $self->seq();
    $seq =~ tr/acgtnACGTN/tgcanTGCAN/;
    my $r = CORE::reverse $seq;
    return $r;
}

sub revcomp {
    my $self = shift;
    $self->seq( $self->getComplementary() );
    $self->{PHASE} = ( $self->{PHASE} == 1 ) ? -1 : 1;
    return $self->seq();
}

sub getTranslation {
    my $self = shift;
    my $code = shift || 1;
    return $self->to_aa( 1, $code );
}

sub translation {
    my $self = shift;
    my $code = shift || 1;
    $self->seq( $self->getTranslation($code) );
    return $self->seq();
}

sub mutation {
    my ( $self, $pos, $newbase ) = @_;
    substr( $self->{SEQ}, $pos - 1, 1, $newbase );
}

sub insertion {
    my ( $self, $pos, $insertion ) = @_;
    substr( $self->{SEQ}, $pos - 1, 0, $insertion );
}

sub mut_wyw {
    my ( $self, $pos, $insertion, $overlap ) = @_;
    substr( $self->{SEQ}, $pos - 1, $overlap, $insertion );
}

sub load_new_seq {
    my ( $self, $seq ) = @_;
    $self->{SEQ} = $seq;
}

sub deletion {
    my ( $self, $pos, $length ) = @_;
    if ( !defined $length ) { $length = 1; }
    substr( $self->{SEQ}, $pos - 1, $length, "" );
}

# Formate une chaine sur long_lign colonnes
#-------------------------------------------
sub _fct_fold {
    my ( $self, $prefix, $long_lign ) = @_;

    if ( !defined $long_lign ) { $long_lign = 60; }

    my $longstring;
    $longstring = $self->getLength();
    my $i   = 0;
    my $seq = $self->seq();
    my $res = "";

    for ( $i = 0 ; $i < $longstring ; $i += $long_lign ) {
        $res .= $prefix . substr( $seq, $i, $long_lign );
        if ( $i + $long_lign < $longstring ) { $res .= "\n"; }
    }
    return $res . "\n";
}

sub _fct_fold_fh {
    my ( $self, $fh, $prefix, $long_lign ) = @_;

    if ( !defined $long_lign ) { $long_lign = 60; }

    my $i          = 0;
    my $longstring = $self->getLength();
    my $seq        = $self->seq();

    $fh->autoflush(1);
    for ( $i = 0 ; $i < $longstring ; $i += $long_lign ) {
        print $fh $prefix, substr( $seq, $i, $long_lign ), "\n";
    }
}

sub to_aa {
    my ( $self, $phase, $code ) = @_;
    if ( !defined $phase ) { $phase = 1; }
    if ( !defined $code )  { $code  = 1; }
    my $nt     = $self->array();
    my $taille = $self->getLength();
    my ( $result, $deb, $fin, $offset ) = ( "", $phase - 1, $taille, 2 );
    my %code = (
        "GCA" => "A",
        "GCC" => "A",
        "GCG" => "A",
        "GCT" => "A",
        "TGC" => "C",
        "TGT" => "C",
        "GAC" => "D",
        "GAT" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "TTC" => "F",
        "TTT" => "F",
        "GGA" => "G",
        "GGC" => "G",
        "GGG" => "G",
        "GGT" => "G",
        "CAC" => "H",
        "CAT" => "H",
        "ATA" => "I",
        "ATC" => "I",
        "ATT" => "I",
        "AAA" => "K",
        "AAG" => "K",
        "TTA" => "L",
        "TTG" => "L",
        "CTA" => "L",
        "CTC" => "L",
        "CTG" => "L",
        "CTT" => "L",
        "ATG" => "M",
        "AAC" => "N",
        "AAT" => "N",
        "CCA" => "P",
        "CCC" => "P",
        "CCG" => "P",
        "CCT" => "P",
        "CAA" => "Q",
        "CAG" => "Q",
        "AGA" => "R",
        "AGG" => "R",
        "CGA" => "R",
        "CGC" => "R",
        "CGT" => "R",
        "CGG" => "R",
        "AGC" => "S",
        "AGT" => "S",
        "TCA" => "S",
        "TCC" => "S",
        "TCT" => "S",
        "TCG" => "S",
        "ACA" => "T",
        "ACC" => "T",
        "ACT" => "T",
        "ACG" => "T",
        "GTA" => "V",
        "GTC" => "V",
        "GTT" => "V",
        "GTG" => "V",
        "TGG" => "W",
        "TAC" => "Y",
        "TAT" => "Y",
        "TAG" => "*",
        "TAA" => "*",
        "TGA" => "*",
    );    # standard genetic code

    if ( defined $code && $code != 1 ) {
        if ( $code == 6 ) {    # ciliates genetic code
            $code{"TAG"} = "Q";
            $code{"TAA"} = "Q";
        }
        elsif ( $code == 5 ) {    # mitochondrial insect genetic code
            $code{"AGA"} = "S";
            $code{"AGG"} = "S";
            $code{"ATA"} = "M";
            $code{"TGA"} = "W";
        }
	elsif($code == 12) # The Alternative Yeast Nuclear Code
        {
            $code{"CTG"} = "S";
        }
    }

    for ( my $i = $deb ; $i + $offset < $fin ; $i += 3 ) {
        my $aa = "";
        $aa = $nt->[$i] . $nt->[ $i + 1 ] . $nt->[ $i + 2 ];
        $result .= ( defined $code{uc($aa)} ) ? $code{uc($aa)} : "X";
    }
    return $result;
}

# Compatibility with previous version of this library
sub toComplementary {
    my $self = shift;
    warn
"[Sequence] Warning - toComplementary - : Deprecated method, use revcomp instead of toComplementary\n";
    return $self->revcomp();
}

sub getCompFrag {
    my $self = shift;
    warn
"[Sequence] Error - getCompFrag - : This method was removed, use revcomp and getFrag\n";
    die;
}

1;
