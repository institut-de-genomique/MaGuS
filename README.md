# MaGuS
MaGuS (Map-Guided Scaffolding) is an efficient map-guided scaffolder to improve assembly and a powerful reference-free evaluator of quality assembly.

This modular tool uses a draft genome assembly, a genome map and high-throughput sequencing paired-end data.

MaGuS provides new quality metrics and perfomed a new scaffolding of the assembly increasing the continuity by creating new links in low-covered regions and highly repeated regions where other usual scaffolding methods lack consistency.

MaGuS is distributed open-source under CeCILL FREE SOFTWARE LICENSE. Check out http://www.cecill.info/ for more information about the contents of this license.

MaGuS home on the web is http://www.genoscope.cns.fr/magus

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
The following binary fastalength comes from the exonerate software availaible under the LGPL license. 
See http://www.ebi.ac.uk/~guy/exonerate/ for details.
The following binary getseq is a Genoscope tool.

INSTALLATION
------------

  1. Clone this GitHub repository
  3. Modify if needed the Perl interpreter that has been set to : /usr/bin/perl
  4. To test the program download this example dataset
  `wget http://www.genoscope.cns.fr/externe/magus/datasets/MaGuS_example_arabido.tar.gz`
  5. Untar/unzip the archive :
  `tar -zxvf MaGuS_example_arabido.tar.gz`
  6. Run MaGuS :
```
$ magus all -wgp tagsWgp.out -tags tags.bam -reads mp1.bam,5350,1000,76 -reads mp2.bam,5350,1000,76 -scaff assembly.fa -prefix arabido -genome 119667750
```

RUNNING MaGuS
--------------
This part describes the different steps required to run MaGuS

The directory "MaGuS_example_arabido" provides all files needed tu run MaGuS.

### Inputs

- tagsWgp.out : WGP tags provided in Keygene format.

- tags.bam : tags alignments on dr'aft assembly previously build in bam format.

- mp1.bam and mp2.bam : paired reads alignments on dr'aft assembly previously build in bam format.

- 5350 : mean size of paired reds bank in bp.

- 1000 : standard deviation of the paired reads bank.

- 76 : reads length.

- assembly.fa : fasta file of draft assembly previously build.

- arabido : prefix for output files.

- 119667750 : genome size in bp.


Warning : MaGuS doesn't support character _ on fasta header files.


### Results

MaGuS runs on arabidopsis data in less than 20 minutes et needs 1,3 GB of memory.

### Options

You can run each step separatly with options :

- wgp2map
- map2qc
- map2links
- pairs2links
- links2scaf

### More informations

```
magus -h
```

Download the documentation http://www.genoscope.cns.fr/externe/magus/magus-1.0.pdf.

ACKNOWLEDGMENTS
---------------
Carole Dossat, Jean-Marc Aury and Amin Madoui - MaGuS's authors

This work was financially supported by the Genoscope,
Institut de Genomique, CEA and Agence Nationale de la
Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).


