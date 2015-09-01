
=head1 NAME

ToolBox::FileTools::FileTools -

=head1 AUTHORS

2008, Genoscope - CNS, Odile Rogier, orogier@genoscope.cns.fr

=head1 SYNOPSIS

use ToolBox::FileTools::FileTools;

=head1 DESCRIPTION

FileTools provides useful functions to manage files

=cut

package ToolBox::FileTools::FileTools;

use strict;
use vars qw(@ISA @EXPORT);
use Exporter;
use File::Temp qw/ tempfile tempdir /;
use File::Path;
use FileHandle;

@ISA = qw( Exporter );
@EXPORT =
  qw(openFile create_dir create_tmp_dir create_tmp_file create_LSF_tmp_file create_PFM_tmp_file create_LSF_tmp_dir);

# Open a file in different modes -----------------------------------------------
sub openFile {
    my ( $file, $mode, $default ) = @_;
    if ( !defined $default ) {
        $default = ( $mode eq "r" ) ? "-" : ">-";
    }
    my $fh;
    if ( defined $default && $file eq $default ) {
        $fh = new FileHandle($file);
    }
    else {
        $fh = new FileHandle( $file, $mode );
    }
    if ( !defined $fh ) {
        warn("Impossible d'ouvrir le fichier $file en mode [$mode].");
    }
    return $fh;
}

# Creates a specified directory if it does not exist and returns 0 if the mkdir command failed
sub create_dir {
    my $dir = shift;
    unless ( -d $dir ) {
        unless ( mkpath($dir) ) {
            if ($@) {
                print "Couldn't create $dir: $@";
            }
            return (0);
        }
    }
    return (1);
}

# Create a new temporary directory ---------------------------------------------
sub create_tmp_dir {
    my $prefix     = shift || "tmp_";
    my $tmpDirName = "${prefix}XXXXXXXXXX";
    my $dirname    = tempdir($tmpDirName);
    if ( !defined $dirname || $dirname eq "" ) {
        warn "[Error] : Can't create tempdir $dirname $!.\n";
        exit(1);
    }
    system("chmod -R 0775 $dirname");
    return $dirname;
}

# Create a new temporary directory to be used for bsubs for Etnas --------------
sub create_LSF_tmp_dir {
    my $prefix     = shift || "tmp_";
    my $tmpDirName = "/env/cns/tmp/${prefix}XXXXXXXXXX";
    my $dirname    = tempdir($tmpDirName);
    if ( !defined $dirname || $dirname eq "" ) {
        warn "[Error] : Can't create tempdir $dirname $!.\n";
        exit(1);
    }
    system("chmod -R 0775 $dirname");
    return $dirname;
}

# Create a new temporary file --------------------------------------------------
sub create_tmp_file {
    my $prefix = shift || "jobs_";
    my $suffix = shift || ".dat";
    my $tmpFileName = "/tmp/${prefix}XXXXXXXXXX";
    my ( $fh_output, $filename ) = tempfile( $tmpFileName, SUFFIX => $suffix );
    if ( !defined $fh_output ) {
        warn "[Error] : Can't create tempfile $filename $!.\n";
        exit(1);
    }
    return ( $fh_output, $filename );
}

# Create temporary file to be used for bsubs on Etnas --------------------------
sub create_LSF_tmp_file {
    my $prefix = shift || "jobs_";
    my $suffix = shift || ".dat";
    my $tmpFileName = "/env/cns/tmp/${prefix}XXXXXXXXXX";
    my ( $fh_output, $filename ) = tempfile( $tmpFileName, SUFFIX => $suffix );
    if ( !defined $fh_output ) {
        warn "[Error] : Can't create LSF tempfile $filename $!.\n";
        exit(1);
    }
    return ( $fh_output, $filename );
}

# Create a temporary file for the mutation detection platform ------------------
sub create_PFM_tmp_file {
    my $prefix = shift || "jobs_";
    my $suffix = shift || ".dat";
    my $tmpFileName =
      "/env/cns/proj/Plateforme_mutation/tmp_pfm/${prefix}XXXXXXXXXX";
    my ( $fh_output, $filename ) = tempfile( $tmpFileName, SUFFIX => $suffix );
    if ( !defined $fh_output ) {
        warn "[Error] : Can't create LSF tempfile $filename $!.\n";
        exit(1);
    }
    return ( $fh_output, $filename );
}

=head1 LIST OF METHODS

=cut

=head2 create_LSF_tmp_file()

=over 4

=item Function :
Creates a tempory file to '/env/cns/tmp' (read from any node of the
 cluster) and returns a filehandle and file name.

=item Usage    : 
my ($tmp_fh, $tmp_name) = &create_LSF_tmp_file();

=item Returns  : filehandle, string

=item Args     :

=over 4

=item -prefix (String)  : optional (default is 'jobs_')

=item -suffix (String)  : optional (default is '.dat')

=back

=back

=cut

1;
