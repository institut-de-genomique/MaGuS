# MaGuS
MaGuS (Map-GUided Scaffolding) is a scaffolder and a reference-free evaluator of assembly quality. It uses a draft genome assembly, a genome map, and high-throughput sequencing paired-end data. It has been succesfully tested on the Arabidopsis genome with Illumina reads and a Whole-Genome Profiling (WGP) map.

MaGuS run the five following steps :
```
- wgp2map: create a map genome file based on WGP data.
- map2qc: analyse the assembly quality based on colinearity between the assembly and the map.
- map2links: create links (map-links) between scaffolds .
- pairs2links: use NGS data to validate the map-links, orient the scaffolds and estimate the gaps size and create a '.de'.
- links2scaf: output the new final assembly in fasta format.
```
The use of MaGuS is not restricted to WGP data, other map types can be used. However they have to be formatted in the MaGuS map format (see below).

MaGuS is distributed open-source under CeCILL FREE SOFTWARE LICENSE. Check out http://www.cecill.info/ for more information about the contents of this license.

MaGuS website http://www.genoscope.cns.fr/externe/magus

Contact : magus [a] genoscope [.] cns [.] fr

PRE-REQUIREMENTS
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

  1. Download the .zip file MaGuS-master: 
  `wget https://github.com/institut-de-genomique/MaGuS/archive/master.zip`
  2. Unzip it:  
  `unzip Arabido_data.tar.gz`
  3. Modify if needed the Perl interpreter that has been set to : /usr/bin/perl
  4. Download the example dataset available on the website http://www.genoscope.cns.fr/externe/magus or
  `wget http://www.genoscope.cns.fr/externe/magus/datasets/MaGuS/Arabido/Arabido_data.tar.gz`
  5. Untar/unzip it:  
  `tar -zxvf MaGuS_example_arabido.tar.gz`
  6. Add MaGuS libraries in $PATH (i.e. `PERL5LIB=$(pwd)/MaGuS-master/magus-1.0/lib/:$PERL5LIB`)
  7. Run MaGuS on the example data set :
```
$ magus all -w Arabido_data/tagsWgp.out -t Arabido_data/mapped_tag.bam -b Arabido_data/mp_map1_2.bam,5414,1000,76 -b Arabido_data/mp_map2_2.bam,5414,1000,76 -f Arabido_data/Arabido.fa -p Arabido -e 119667750
```

RUNNING MaGuS
--------------

There are two ways to run MaGuS. The most common way is:
```
 magus all -w wgpFile -t tags.bam -f assembly.fa -e estimate_size -b file.bam,m,sd,s
```


User can also choose to run MaGuS pipeline step by step as:

- step1:
```
magus wg2map -w wgpFile -t tags.bam
```
- step2:
```
magus map2qc -f assembly.fa -e estimate_size -s tags_coordinates.txt
```
- step3:
```
magus map2links -a anchoring_file.txt
```
- step4:
```
magus pairs2links -f assembly.fa -l links_file.txt -b file.bam,m,sd,s
```
- step5:
```
magus links2scaf -f assembly.fa -c links.de
```

OPTIONS
--------------

- ##### Options for all (wgp2map-map2qc-map2links-pairs2links-links2scaf)

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

- ##### Options for wgp2map

        -w <string>     wgpFile: WGP data
        -t <string>     tags.bam: tags alignment on the assembly (BAM)

        OPTIONAL PARAMETERS:
        -p <string>     prefix for output files (default : magus)
        -h              this help

- ##### Options for map2qc

        -f <string>     assembly.fa: assembly file (FASTA)
        -e <int>        estimate_size: genome estimate size (bp)
        -s <string>     tags_coordinates.txt: sorted file according to mapping position of tags

        OPTIONAL PARAMETERS:
        -p <string>     prefix for output files (default: magus)
        -h              this help

- ##### Options for map2links

        -a <string>     anchoring_file.txt: anchored assembly (wgp2map output)

        OPTIONAL PARAMETERS:
        -l <string>     output file containing links between contigs/scaffolds  (default: prefix_map_links.txt)
        -m <string>     Bin path        (default: $PATH)
        -p <string>     prefix for output files (default: magus)
        -h              this help

- ##### Options for pairs2links

        -f <string>     assembly.fa: assembly file (FASTA)
        -l <string>     links_file.txt: file containing links between contigs/scaffolds
        -b <string>     file.bam,m,sd,s: paired reads alignment (BAM), library median size (bp), library standart deviation (bp), reads size (bp)

        OPTIONAL PARAMETERS:
        -v <string>     path to samtools        (default: $PATH)
        -p <string>     prefix for output files (default: magus)
        -h              this help


Several mapped paired-end libraries (BAM file) can be used simultanously with the -b option for each one of them. Example:

magus all -w tagsWgp.out -t mapping.bam -f Arabidopsis.fa -b mapping_library1.bam,3500,600,101 -b mapping_library2.bam,6000,1000,151 -b mapping_library3.bam,9000,1400,251 -e 120000000

- ##### Options for links2scaf	
        -f <string>     assembly.fa: assembly file (FASTA)
        -c <string>     links.de: file containing links in DE format

        OPTIONAL PARAMETERS:
        -z <string>     path to sga     (default: $PATH)
        -g <string>     path to getseq  (default: $PATH)
        -p <string>     prefix for output files (default: magus)
        -h              this help

OUTPUT
--------------

- wgp2map output

-${prefix}_tags_coordinate.txt: File containing anchored tags on assembly sorted by mapping position

  col 1: scaffold ID  
  col 2: position    
  col 3: tagId  
  col 4: rank	
  col 5: group ID 
  
-${prefix}_anchored_assembly.txt: MaGuS format File of anchored scaffolds on the genome map

  col 1: group ID                
  col 2: scaffold ID	
  col 3: minimum tag rank	
  col 4: maximum tag rank	
  col 5: number of tags	

- map2qc output	
-prefix_An.csv: Anx values for x = 1 to 100%  
-prefix_AnA.csv: AnAx values for x = 1 to 100%  
-prefix_AnG.csv: AnGx vlaues for x = 1 to 100%  
-prefix_quality_metrics.png: Anx, AnAx and AnGx plots 
-prefix_quality_metrics.txt: summary quality metrics of Anx, AnAx and AnGx values for x=0.5, x=0.75 and x=0.9 

- map2links output	
-${prefix}_map_links.txt: list of map-links between scaffolds (scaf1_scaf2)

- pairs2links output	
-${prefix}_validated_map_links.de: map-links validated by paired reads in .de format (SGA specific format)  
-${prefix}_unvalidated_map_links.txt: list of map-links not validated by the paired reads	
-${prefix}_besst_validated_map_links.log: links statistics  

- links2scaf output	
-${prefix}_final_scaffolds.fa: the final assembly (FASTA)	
-${prefix}_scaffolds.fa: SGA scaffolding output (FASTA) 
-${prefix}_unused_scaffolds.fa: scaffold not used for scaffolding step (FASTA)  
-${prefix}_sga_scaf_unused.txt: list of scaffolds not considered by SGA 
-${prefix}_sga.scaf: scaffolding information (SGA specific format)  
-${prefix}_sga_scaffold2fasta.log: scaffolding statistics (SGA specific format) 
-${prefix}_sga_scaffold.log: scaffolding statistics/informations (SGA specific format)  


### More informations
--------------
If you have questions about MaGuS, you may ask them to amadoui [at] genoscope [.] cns [.] fr, cdossat [at] genoscope [.] cns [.] fr, ldagata [at] genoscope [.] cns [.] fr and jmaury [at] genoscope [.] cns [.] fr . You may also create an issue to ask questions on github website: https://github.com/institut-de-genomique/MaGuS/issues. 


ACKNOWLEDGMENTS
---------------
Carole Dossat, Jean-Marc Aury and Mohammed-Amin Madoui - MaGuS's authors

This work was financially supported by the Genoscope,
Institut de Genomique, CEA and Agence Nationale de la
Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).


