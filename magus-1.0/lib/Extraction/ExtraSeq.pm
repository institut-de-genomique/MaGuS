package Extraction::ExtraSeq;

use strict;
use warnings;
use Extraction::FSequence;
use Extraction::QSequence;

my $DEBUG = 0;

use constant READONLY  => 0;
use constant READWRITE => 1;

use constant FASTA_TYPE   => 0;
use constant QUALITY_TYPE => 1;
my %TYPE_NAME = ( 0 => "fasta", 1 => "quality" );
use constant FIRST_INDEX => 0;

sub new {
    my ( $proto, $bank, $cache, $type ) = @_;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    $self->{BANKFILE}           = $bank;
    $self->{HANDLE}             = undef;
    $self->{SEQ2HASH}           = ();
    $self->{SEQ2INDEX}          = ();
    $self->{FIRSTSEQ}           = "";
    $self->{NBSEQ}              = 0;
    $self->{CACHE}              = ();
    $self->{USE_CACHE}          = $cache || 1;
    $self->{MAX_SIZE_CACHE}     = 100000000;
    $self->{CURRENT_CACHE_SIZE} = 0;
    $self->{TYPE}               = $type || FASTA_TYPE;
    $self->{HOMOGENEOUS_WIDTH}  = 0;
    $self->{TOTAL_LENGTH}       = 0;
    $self->{VERBOSE}            = 0;
    return $self;
}

sub bankfile {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{BANKFILE} = $arg; }
    else                { return $self->{BANKFILE}; }
}

sub max_size_cache {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{MAX_SIZE_CACHE} = $arg; }
    else                { return $self->{MAX_SIZE_CACHE}; }
}

sub firstseq {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{FIRSTSEQ} = $arg; }
    else                { return $self->{FIRSTSEQ}; }
}

sub nbseq {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{NBSEQ} = $arg; }
    else                { return $self->{NBSEQ}; }
}

sub total_length {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{TOTAL_LENGTH} = $arg; }
    else                { return $self->{TOTAL_LENGTH}; }
}

sub usecache {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{USE_CACHE} = $arg; }
    else                { return $self->{USE_CACHE}; }
}

sub verbose {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{VERBOSE} = $arg; }
    else                { return $self->{VERBOSE}; }
}

sub type {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) {
        if ( !defined $TYPE_NAME{$arg} ) {
            warn "[ExtraSeq] Error : Unknown database type $arg !\n";
            die;
        }
        $self->{TYPE} = $arg;
    }
    else { return $self->{TYPE}; }
}

sub homogeneous_width {
    my $self = shift;
    my $arg  = shift;
    if ( defined $arg ) { $self->{HOMOGENEOUS_WIDTH} = $arg; }
    else                { return $self->{HOMOGENEOUS_WIDTH}; }
}

sub isHomogeneousWidth {
    my $self = shift;
    if ( defined $self->{HOMOGENEOUS_WIDTH} > 0 ) { return 1; }
    return 0;
}

sub spaceInCache {
    my $self = shift;
    return $self->{MAX_SIZE_CACHE} - $self->{CURRENT_CACHE_SIZE};
}

sub _putSeqInCache {
    my ( $self, $seq, $name, $flag_RO_RW ) = @_;
    my $size = length($$seq);
    warn "[ExtraSeq] Debug : Put in cache seq : $name : $size\n" if $DEBUG;
    ( $self->{CACHE} )->{$name} =
      { SEQ => $self->_sequenceFactory( $name, $seq ), SIZE => $size };
    $self->{CURRENT_CACHE_SIZE} = $self->{CURRENT_CACHE_SIZE} + $size;

    # Correction BUG - jmaury AT genoscope.cns.fr
    # 05/06/2008 - ne pas renvoyer un pointeur sur une sequence du cache,
    # mais plutot une copie de l'objet sequence du cache
    if ( $flag_RO_RW == READONLY ) {
        return ( $self->{CACHE} )->{$name}->{SEQ};
    }
    else { return $self->_sequenceFactory( $name, $seq ); }
}

sub putSeqInCache {
    my ( $self, $seq, $name, $flag_RO_RW ) = @_;
    my $size = length($$seq);
    if ( $size + $self->{CURRENT_CACHE_SIZE} <= $self->{MAX_SIZE_CACHE} ) {
        return $self->_putSeqInCache( $seq, $name, $flag_RO_RW );
    }
    return undef;
}

sub forceSeqInCache {
    my ( $self, $seq, $name, $flag_RO_RW ) = @_;
    my $size = length($seq);
    warn "[ExtraSeq] Debug : Force to cache seq : $name : $size\n" if $DEBUG;
    if ( $size > $self->{MAX_SIZE_CACHE} ) {
        $self->clearCache();
        return $self->_putSeqInCache( $seq, $name, $flag_RO_RW );
    }
    if ( $size + $self->{CURRENT_CACHE_SIZE} > $self->{MAX_SIZE_CACHE} ) {
        $self->libereCache($size);
        return $self->_putSeqInCache( $seq, $name, $flag_RO_RW );
    }
    return $self->_putSeqInCache( $seq, $name, $flag_RO_RW );
}

sub seqInCache {
    my ( $self, $seqName, $flag_RO_RW ) = @_;
    if ( $self->usecache() ) {
        if ( exists( ( $self->{CACHE} )->{$seqName} ) ) {
            warn "[ExtraSeq] Debug : sequence $seqName found in cache\n"
              if $DEBUG;

           # Correction BUG - jmaury AT genoscope.cns.fr
           # 05/06/2008 - ne pas renvoyer un pointeur sur une sequence du cache,
           # mais plutot une copie de l'objet sequence du cache
            if ( $flag_RO_RW == READONLY ) {
                return ( 1, ( $self->{CACHE} )->{$seqName}->{SEQ} );
            }
            else {
                my $seq = ( $self->{CACHE} )->{$seqName}->{SEQ};
                return ( 1, $self->_sequenceFactory( $seqName, $seq->seq() ) );
            }
        }
    }
    warn "[ExtraSeq] Debug : sequence $seqName not found in cache\n" if $DEBUG;
    return ( 0, undef );
}

sub libereCache {
    my ( $self, $quantity ) = @_;
    my $tot = 0;
    warn "[ExtraSeq] Debug : must free : $quantity\n" if $DEBUG;
    while ( my ( $clef, $valeur ) = each %{ $self->{CACHE} } ) {
        my $size = $valeur->{SIZE};
        $self->{CURRENT_CACHE_SIZE} = $self->{CURRENT_CACHE_SIZE} - $size;
        delete( ( $self->{CACHE} )->{$clef} );
        warn "[ExtraSeq] Debug :   delete sequence : $clef : $size\n" if $DEBUG;
        $tot += $size;
        if ( $tot >= $quantity ) { return 1; }
    }
}

sub clearCache {
    my $self = shift;
    $self->{CACHE}              = ();
    $self->{CURRENT_CACHE_SIZE} = 0;
    warn "[ExtraSeq] Debug : clear cache...\n" if $DEBUG;
}

sub loadBank {
    my ( $self, $nomaj, $min_size ) = @_;
    if ( !defined $nomaj )    { $nomaj    = 1; }
    if ( !defined $min_size ) { $min_size = -1; }
    open( BANK, "<" . $self->bankfile() ) || die $self->bankfile(),
      " : Fichier inexistant\n";
    my (
        $precSeq, $endLine, $first,     $cpt, $verifType,
        $name,    $width,   $sameWidth, $len, $totlen
    ) = ( undef, 0, 1, FIRST_INDEX, 0, undef, 0, 1, 0, 0 );
    while (<BANK>) {
        if ( $_ =~ /^>(\S+)/o ) {
            my %seq = ();
            warn "[ExtraSeq] Debug : load sequence $1...\n" if $DEBUG;
            $seq{DEBUT} = tell(*BANK);
            $seq{INDEX} = $cpt++;
            $self->nbseq( $self->nbseq() + 1 );
            if ( defined $precSeq ) {
                $precSeq->{FIN}       = tell(*BANK) - length($_) - 2;
                $precSeq->{LONG2READ} = $precSeq->{FIN} - $precSeq->{DEBUT} + 1;
                $precSeq->{LENGTH}    = $len;
                if ( $len < $min_size ) {
                    delete $self->{SEQ2HASH}{ $precSeq->{NAME} };
                    delete $self->{SEQ2INDEX}[ $precSeq->{INDEX} ];
                    $seq{INDEX} = $seq{INDEX} - 1;
                    $cpt--;
                }
                else {
                    $totlen += $len;
                }
                if ( $sameWidth == -1 ) { $sameWidth = 1; }
            }
            if ($nomaj) { $name = $1; }
            else        { ( $name = $1 ) =~ tr/a-z/A-Z/; }
            if ($first) { $first = 0; $self->{firstSeq} = $name; }
            $seq{NAME} = $name;
            $self->{SEQ2HASH}{$name} = $precSeq = \%seq;
            $self->{SEQ2INDEX}[ $seq{INDEX} ] = \%seq;
            $len = 0;
            warn "[ExtraSeq] Debug : sequence $1 loaded.\n" if $DEBUG;
        }
        else {
            if ( !$verifType ) {
                if ( $_ =~ /^[a-zA-Z]+/ ) {
                    if ( $self->{TYPE} != FASTA_TYPE ) {
                        warn "[ExtraSeq] Warning : Change database type from ",
                          $TYPE_NAME{ $self->{TYPE} },
                          " to ", $TYPE_NAME{ (FASTA_TYPE) }, " [",
                          $self->bankfile(), "]\n"
                          if $self->verbose();
                        $self->{TYPE} = FASTA_TYPE;
                    }
                }
                elsif ( $_ =~ /^\s*[0-9]+/ ) {
                    if ( $self->{TYPE} != QUALITY_TYPE ) {
                        warn "[ExtraSeq] Warning : change database type from ",
                          $TYPE_NAME{ $self->{TYPE} },
                          " to ", $TYPE_NAME{ (QUALITY_TYPE) }, " [",
                          $self->bankfile(), "]\n"
                          if $self->verbose();
                        $self->{TYPE} = QUALITY_TYPE;
                    }
                }
                else {
                    warn "[ExtraSeq] Warning : can't detect the database type\n"
                      if $self->verbose();
                }
                $verifType = 1;
            }
            $len += length($_) - 1;
            if ( $sameWidth == -1 ) { $sameWidth = 0; }
            else {
                if ( $width == 0 )          { $width     = length($_); }
                if ( $width != length($_) ) { $sameWidth = -1; }
            }
        }
        $endLine = tell(*BANK) - 2;
    }
    if ( $sameWidth == -1 || $sameWidth == 1 ) {
        warn
"[ExtraSeq] Debug : Homogeneous width detected for that bank : $width\n"
          if $DEBUG;
        $self->homogeneous_width($width);
    }
    if ( defined $precSeq ) {
        $precSeq->{FIN}       = $endLine;
        $precSeq->{LONG2READ} = $precSeq->{FIN} - $precSeq->{DEBUT} + 1;
        $precSeq->{LENGTH}    = $len;
        if ( $len < $min_size ) {
            delete $self->{SEQ2HASH}{ $precSeq->{NAME} };
            delete $self->{SEQ2INDEX}[ $precSeq->{INDEX} ];
        }
        else {
            $totlen += $len;
        }
    }
    $self->total_length($totlen);
    close(BANK);
}

sub isSeq {
    my ( $self, $seqname, $nomaj ) = @_;
    if ( !defined $seqname ) { return 0; }
    if ( !defined $nomaj ) { $nomaj = 0; }
    unless ( defined $nomaj && $nomaj ) { $seqname =~ tr/a-z/A-Z/; }
    my $seq = $self->{SEQ2HASH}{$seqname};
    if   ( !defined $seq ) { return 0; }
    else                   { return 1; }
}

sub _sequenceFactory {
    my ( $self, $name, $seq ) = @_;
    my $type = $self->type();
    if ( $type == FASTA_TYPE ) {
        return new Extraction::FSequence( $name, $seq );
    }
    if ( $type == QUALITY_TYPE ) {
        return new Extraction::QSequence( $name, $seq );
    }
}

sub _fileRead {
    my ( $self, $start, $len, $tampon ) = @_;
    my $handle = $self->{HANDLE};
    if ( !defined $handle ) {
        open( $handle, "<" . $self->bankfile() )
          || die "[ExtraSeq] Error - fileRead - : Can't open file ",
          $self->bankfile(), " in [r] mode\n";
        $self->{HANDLE} = $handle;
    }
    seek( $handle, $start, 0 );
    read( $handle, $$tampon, $len );

    #close($handle);
}

sub _getSeq {
    my ( $self, $seqName, $exact, $nomaj ) = @_;
    my $seq;
    if ( !defined $exact ) { $exact = 0; }
    if ( !defined $nomaj ) { $nomaj = 1; }
    if ( defined $nomaj && !$nomaj ) { $seqName =~ tr/a-z/A-Z/; }
    if ( $seqName =~ m/\*/ ) {
        while ( my ( $clef, $valeur ) = each %{ $self->{SEQ2HASH} } ) {
            if ( $clef =~ m/$seqName/ ) { $seq = $valeur; }
        }
        if ( !defined $seq ) {
            die "[ExtraSeq] Error - getSeq - : $seqName : Unknown sequence\n";
        }
    }
    else {
        $seq = $self->{SEQ2HASH}{$seqName};
    }
    if ( !defined $seq ) {
        if ($exact) {
            die "[ExtraSeq] Error - getSeq - : $seqName : Unknown sequence\n";
        }
        else { $seq = $self->_getSeq( $seqName . ".*", 0, $nomaj ); }
    }
    return $seq;
}

sub getSeqIndex {
    my ( $self, $seqName ) = @_;
    my $seq = $self->_getSeq($seqName);
    return $seq->{INDEX};
}

sub wrapper_getSeq {
    my ( $self, $flag_RO_RW, $seqName, $exact, $nomaj ) = @_;
    if ( $self->usecache() ) {
        my ( $ret, $s ) = $self->seqInCache( $seqName, $flag_RO_RW );
        if ($ret) { return $s; }
    }

    #if (!defined $exact) { $exact=0; }
    #if (!defined $nomaj) { $nomaj=1; }
    my $seq = undef;

    #if ($exact) { $seq = $self->_getSeq($seqName, 1, $nomaj); }
    #else { $seq = $self->_getSeq($seqName, 0, $nomaj); }
    $seq        = $self->_getSeq( $seqName, $exact, $nomaj );
    my $tampon  = "";
    $self->_fileRead( $seq->{DEBUT}, $seq->{LONG2READ}, \$tampon );
    $seqName    =~ s/\*//g;
    if ( $self->usecache() ) {
        return $self->forceSeqInCache( \$tampon, $seqName, $flag_RO_RW );
    }
    return $self->_sequenceFactory( $seqName, \$tampon );
}

sub getSeq_readonly {
    my ( $self, $seqName, $exact, $nomaj ) = @_;
    return $self->wrapper_getSeq( READONLY, $seqName, $exact, $nomaj );
}

sub getSeq {
    my ( $self, $seqName, $exact, $nomaj ) = @_;
    if ( !defined $seqName ) {
        warn "[ExtraSeq] Error - getSeq - Empty sequence name.\n";
        exit 1;
    }
    return $self->wrapper_getSeq( READWRITE, $seqName, $exact, $nomaj );
}

sub getSeqByIndex {
    my ( $self, $index, $nomaj ) = @_;
    if ( $index > $self->{NBSEQ} + 1 ) {
        die "[ExtraSeq] Error - getSeqByIndex - : index($index) > nbseq(",
          $self->{NBSEQ} + 1, ").";
    }
    my $seq = $self->{SEQ2INDEX}[$index];
    if ( $self->usecache() ) {
        my ( $ret, $s ) = $self->seqInCache( $seq->{NAME}, READONLY );
        if ($ret) { return $s; }
    }
    my $tampon = "";
    $self->_fileRead( $seq->{DEBUT}, $seq->{LONG2READ}, \$tampon );
    if ( $self->usecache() ) {
        return $self->forceSeqInCache( \$tampon, $seq->{NAME}, READONLY );
    }
    return $self->_sequenceFactory( $seq->{NAME}, \$tampon );
}

sub getFirstSeq {
    my $self = shift;
    my $seq;
    $seq = $self->{SEQ2HASH}{ $self->{firstSeq} };
    if ( !defined $seq ) {
        die "[ExtraSeq] Error - getSeq - : ", $self->{firstSeq},
          " : Unknown sequence\n";
    }
    my $tampon = "";
    $self->_fileRead( $seq->{DEBUT}, $seq->{LONG2READ}, \$tampon );
    return $self->_sequenceFactory( $self->{firstSeq}, \$tampon );
}

sub getFrag {
    my ( $self, $nom, $deb, $fin, $softname ) = @_;

 #    if($self->homogeneous_width()) {
 #	my $seq = $self->_getSeq($nom);
 #	my $tampon = "";
 #	$self->_fileRead($seq->{DEBUT}+($deb+int($deb/$self->homogeneous_width()))+2,
 #			 ($fin-$deb+1) + int(($fin-$deb+1)/$self->homogeneous_width())+2,
 #			 \$tampon);
 #	$nom =~ s/\*//g;
 #	return $self->_sequenceFactory($nom, \$tampon);
 #    }
    my $seq = $self->getSeq_readonly($nom);
    return $seq->getFrag( $deb, $fin, $softname );
}

sub getCompFrag {
    my ( $self, $nom, $deb, $fin ) = @_;
    my $seq = $self->getSeq($nom);
    $seq->revcomp();
    return $seq->getFrag( $deb, $fin );
}

sub getLength {
    my ( $self, $nom ) = @_;
    if ( $self->{TYPE} != FASTA_TYPE ) {
        my $seq = $self->getSeq_readonly($nom);
        return $seq->getLength();
    }
    else {
        my $seq = $self->{SEQ2HASH}{$nom};
        if ( !defined $seq ) {
            die "[ExtraSeq] Error - getSeq - : $nom : Unknown sequence\n";
        }
        return $seq->{LENGTH};
    }
}

sub getAllSeqName {
    my $self  = shift;
    my $sort  = shift || 1;
    my @names = ();
    foreach my $nom ( keys %{ $self->{SEQ2HASH} } ) {
        push @names, $nom;
    }
    if ( !$sort ) { return \@names; }
    my @sort = sort { $a cmp $b } @names;
    return \@sort;
}

sub getAllSeqNameOriginalOrder {
    my $self  = shift;
    my @names = ();
    for ( my $i = FIRST_INDEX ; $i < $self->nbseq() ; $i++ ) {
        my $seq = $self->{SEQ2INDEX}[$i];
        push( @names, $seq->{NAME} );
    }
    return \@names;
}

# Compatibility with previous version of this library
sub getSeqFrag {
    my $self = shift;
    warn
"[ExtraSeq] Error - getSeqFrag - : This method was removed, use Sequence::getFrag or ExtraSeq::getFrag\n";
    die;
}

1;
