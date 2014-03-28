HPG Aligner README

Welcome to HPG Aligner !

HPG Aligner is an ultrafast and highly sensitive Next-Generation Sequencing (NGS) read mapping.

COMPONENTS
----------

 hpg-aligner - The executable file to map DNA/RNA sequences


CONTACT
------- 
  You can contact any of the following developers:

    * Joaquín Tárraga (jtarraga@cipf.es)
    * Héctor Martínez (martineh@cipf.es)
    * Ignacio Medina (imedina@cipf.es)


DOWNLOAD and BUILDING
---------------------

  HPG Aligner has been opened to the community and released in GitHub, so you can download by invoking the following commands:

    $ git clone https://github.com/opencb/hpg-aligner.git
    Cloning into 'hpg-aligner'...
    remote: Reusing existing pack: 1441, done.
    remote: Total 1441 (delta 0), reused 0 (delta 0)
    Receiving objects: 100% (1441/1441), 1.85 MiB | 122 KiB/s, done.
    Resolving deltas: 100% (882/882), done.
    $ cd hpg-aligner
    $ git submodule update --init
    Submodule 'lib/hpg-libs' (https://github.com/opencb/hpg-libs.git) registered for path 'lib/hpg-libs'
    Cloning into 'lib/hpg-libs'...
    remote: Reusing existing pack: 7735, done.
    remote: Total 7735 (delta 0), reused 0 (delta 0)
    Receiving objects: 100% (7735/7735), 26.82 MiB | 79 KiB/s, done.
    Resolving deltas: 100% (4430/4430), done.
    Submodule path 'lib/hpg-libs': checked out '962f531ef0ffa2a6a665ae6fba8bc2337c4351a9'

  For the most recent HPG Aligner version, choose the 'develop' Git branch:

    $ git checkout develop

  Finally, use Scons to build the HPG Aligner applilcation:

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
