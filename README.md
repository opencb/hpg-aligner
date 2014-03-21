HPG Aligner README

Welcome to HPG Aligner !

HPG Aligner is an ultrafast and highly sensitive Next-Generation Sequencing (NGS) read mapping.

COMPONENTS
----------

 hpg-aligner - The executable file to map DNA/RNA sequences


CONTACT
------- 
  You can contact any of the following developers:

    * Joaquin Tarraga (jtarraga@cipf.es)
    * Hector Martinez (martineh@cipf.es)
    * Ignacio Medina (imedina@cipf.es)


BUILDING
--------

  Scons is used as building system, to build the library download 'develop' branch (for the most recent library) and execute:

    $ scons


RUNNING
-------

  For command line options invoke:

    $ ./bin/hpg-aligner -h



  In order to map DNA sequences:

    1) Create the HPG Aligner index based on suffix arrays by invoking:

      $ ./bin/hpg-aligner build-sa-index -g <ref-genome-fasta-file> -i <index-directory>

    2) Map by invoking:

      $ ./bin/hpg-aligner dna -i <index-directory> -f <fastq-file> -o <output-directory>



  In order to map RNA sequences:

    1) Create the HPG Aligner index based on Burrows-Wheeler Transform by invoking:

      $ ./bin/hpg-aligner build-bwt-index -g <ref-genome-fasta-file> -i <index-directory> -r <compression-ratio>

    2) Map by invoking:

      $ ./bin/hpg-aligner rna -i <index-directory> -f <fastq-file> -o <output-directory>




DOCUMENTATION
-------------


  You can find more info about HPG project at:

  http://wiki.opencb.org/projects/hpg/doku.php
