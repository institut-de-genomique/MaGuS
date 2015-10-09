# MaGuS
MaGuS (Map-GUided Scaffolding) is a scaffolder and a reference-free evaluator of assembly quality. It uses a draft genome assembly, a genome map, and high-throughput sequencing paired-end data. It has been succesfully tested on the Arabidopsis genome with Illumina reads and a Whole-Genome Profiling (WGP) map.

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

There are two ways to run MaGuS. The most common way to run it is:
```
 magus all -w wgpFile -t tags.bam -f assembly.fa -e estimate_size -b file.bam,m,sd,s
```

User can also choose to run MaGuS pipeline step by step as:

step1:
```
magus wg2map -w wgpFile -t tags.bam
```
step2:
```
magus map2qc -f assembly.fa -e estimate_size -s tags_coordinates.txt
```
step3:
```
magus map2links -a anchoring_file.txt
```
step4:
```
magus pairs2links -f assembly.fa -l links_file.txt -b file.bam,m,sd,s
```
step5:
```
magus links2scaf -f assembly.fa -c links.de
```

N.B: Several mapped paired-end libraries (BAM file) of different fragment size can be used simultanously with the -b option for each one
example: -reads pairs1.bam,mean1,sd1,length1 -reads pairs2.bam,mean2,sd2,length2 ...

### Options

##### Options for all (wgp2map-map2links-pairs2links-links2scaf-map2qc)
```
         -w <string>    wgpFile: WGP data
         -t <string>    tags.bam: tags alignment on the assembly (BAM)
         -f <string>    assembly.fa: assembly file (FASTA)
         -e <int>       estimate_size: genome estimate size (bp)
         -b <string>    file.bam,m,sd,s: paired reads alignment (BAM), library median size (bp), library standart deviation (bp), reads size (bp)

OPTIONAL PARAMETERS:
         -s <string>    sorted file according to mapping position of tags       (default : prefix_tags_coordinates.txt)
         -a <string>    anchored tags on assembly       (default : prefix_anchored_assembly.txt)
         -l <string>    output file containing links between contigs/scaffolds  (default : prefix_map_links.txt)
         -c <string>    file containing links in DE format      (default : prefix_validated_map_links.de)
         -m <string>    Bin path        (default: $PATH)
         -v <string>    path to samtools        (default: $PATH)
         -r <string>    path to R       (default: $PATH)
         -q <string>    path to fastalength     (default: $PATH)
         -z <string>    path to sga     (default: $PATH)
         -g <string>    path to getseq  (default: $PATH)
         -p <string>    prefix for output files (default: magus)
         -h             this help

```
##### Options for wgp2map
        -w <string>     wgpFile: WGP data
        -t <string>     tags.bam: tags alignment on the assembly (BAM)

OPTIONAL PARAMETERS:
        -p <string>     prefix for output files (default : magus)
        -h              this help




### Output

- map2links output
-prefix_map_links.txt: the list of the scaffold links inferred from the genome map, (scaf1_scaf2)

-prefix_ordered_tags.txt: the tags sorted by position on the scaffolds

  col 1: scaffold Id     
  col 2: position
  col 3: tagId
  col 4: rank
  col 5: group ID
  
-prefix_anchorage.txt: the position of the scaffolds on the genome map, contains 5 columns
  col 1: group ID                     
  col 2: scaffold Id
  col 3: minimum tag rank
  col 4: maximum tag rank
  col 5: number of tags
  
- pairs2links output

-prefix_validated_map_links.de: map-links validated by paired reads in .de format (SGA specific format)

-prefix_unvalidated_map_links.txt: list of map-links not validated by the paired reads


- links2scaf output

-prefix_all_scaffolds.fa: the final assembly (FASTA)

-prefix_scaffolds.fa: SGA scaffolding output (FASTA)

-prefix_sga.scaf: scaffolding information (SGA specific format) 

-prefix_sga_scaf_lost.txt: list of scaffolds not considered by SGA

-prefix_lost_scaffolds.fa: scaffold not taken by SGA during the scaffolding step (FASTA)


- map2qc output

-prefix_An.csv: the Anx vlaues for x = 1 to 100%

-prefix_AnA.csv: the AnAx values for x = 1 to 100%

-prefix_AnG.csv: the AnGx vlaues for x = 1 to 100%

-prefix_quality_metrics.png: Anx, AnAx and AnGx plots

-prefix_quality_metrics.txt: summary quality metrics of Anx, AnAx and AnGx values for x=0.5, x=0.75 and x=0.9


### More informations

```
magus -h (more help for each module whith -h, i.e. : magus wgp2map -h)
```

Download the documentation http://www.genoscope.cns.fr/externe/magus/magus-1.0.pdf.

ACKNOWLEDGMENTS
---------------
Carole Dossat, Jean-Marc Aury and Mohammed-Amin Madoui - MaGuS's authors

This work was financially supported by the Genoscope,
Institut de Genomique, CEA and Agence Nationale de la
Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).


