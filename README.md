# MaGuS
MaGuS (Map-GUided Scaffolding) is a scaffolder and a reference-free evaluator of assembly quality. It uses a draft genome assembly, a genome map, and high-throughput sequencing paired-end data. MaGuS provides quality metrics and performed a new scaffolding. It has been succesfully tested on the Arabidopsis genome with Illumina reads and a Whole-Genome Profiling (WGP) map.

MaGuS run the five following steps :

- wgp2map: create a map genome file based on WGP data.
- map2qc: analyse the assembly quality based on colinearity between the assembly and the map.
- map2links: create links (map-links) between scaffolds .
- pairs2links: use NGS data to validate the map-links, orient the scaffolds and estimate the gaps size and create a '.de'.
- links2scaf: output the new final assembly in fasta format.

The use of MaGuS is not restricted to WGP data, other map types can be used. However they have to be formatted in the MaGuS map format (see below).

MaGuS is distributed open-source under CeCILL FREE SOFTWARE LICENSE. Check out http://www.cecill.info/ for more information about the contents of this license.

MaGuS website http://www.genoscope.cns.fr/magus

Contact : magus [a] genoscope [.] cns [.] fr

PRE-REQUISITES
--------------

  - A Linux based operating system.
  - Perl 5.10.1 or higher installed.
  - Perl GetOpt::Long module (http://search.cpan.org/dist/Getopt-Long/)
  - Perl FileHandle module (http://search.cpan.org/dist/FileHandle/)
  - Perl Moose module (http://search.cpan.org/dist/Moose/)
  - Perl Switch module (http://search.cpan.org/dist/Switch/)
  - samtools 0.1.19 
  - SGA 0.10.13 
  - R 3.0.2 (cran.r-project.org)

DEPENDENCIES
------------
fastalength:  program from exonerate availaible under the LGPL license (http://www.ebi.ac.uk/~guy/exonerate/). 


INSTALLATION
------------

  1. Clone this GitHub repository
  2. Modify if needed the Perl interpreter that has been set to : /usr/bin/perl
  3. Add MaGuS libraries in $PATH (i.e. `PERL5LIB=./lib/:$PERL5LIB`)
  4. To test the program download this example dataset
  `wget http://www.genoscope.cns.fr/externe/magus/datasets/MaGuS_example_arabido.tar.gz`
  5. Untar/unzip the archive :
  `tar -zxvf MaGuS_example_arabido.tar.gz`
  6. Run MaGuS on the example data set :
```
$ magus all -wgp tagsWgp.out -tags tags.bam -reads mp1.bam,5350,1000,76 -reads mp2.bam,5350,1000,76 -scaff assembly.fa -prefix arabido -genome 119667750
```

RUNNING EXAMPLE
--------------

The directory "MaGuS_example_arabido" contains the example files to test MaGuS.

### Inputs

-wgp WGP sorted tags

-tags tags alignment on the assembly of Arabidopsis (BAM).

-reads paired-reads alignments on draft assembly (BAM), mean size of paired-reads fragment (bp), standard deviation in bp of the paired-reads fragment size (bp), reads length.

-scaff genome assembly (FASTA)

-prefix output files prefix.

-genome genome size  (bp).

N.B: Several paired-end libraries of different fragment size can be used simultanously, each library has to be added with the -reads options 


### Outputs

MaGuS runs on arabidopsis data in less than 20 minutes et needs 1,3 GB of memory.

### Options

Each step can be executed separately as follow:
```
MaGuS step1 -option1 -option2  
```
### More informations

```
magus -h (more help for each module whith -h, i.e. : magus wgp2map -h)
```

Download the documentation http://www.genoscope.cns.fr/externe/magus/magus-1.0.pdf.

ACKNOWLEDGMENTS
---------------
Carole Dossat, Jean-Marc Aury and Amin Madoui - MaGuS's authors

This work was financially supported by the Genoscope,
Institut de Genomique, CEA and Agence Nationale de la
Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).


