HPG Aligner README

Welcome to HPG Aligner !

HPG Aligner is an ultrafast and highly sensitive Next-Generation Sequencing (NGS) read mapping.

COMPONENTS
----------

 hpg-aligner - The executable file to map DNA/RNA sequences


CONTACT
------- 
  You can contact any of the following developers:

    * Héctor Matínez (martneh@uji.es)


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

  For the most recent HPG Aligner version, choose the 'develop-multi' Git branch (both for the hpg-libs and the HPG Aligner):

    $ cd lib/hpg-libs
    $ git checkout develop-multi
    $ cd ../..
    $ git checkout develop-multi

  For build HPG Aligner yo need Git, GCC, Scons and Tcmalloc library. If some of they are not install in your system, you can try for Ubuntu and Debian:

    $ apt-get install git 
    $ apt-get install gcc
    $ apt-get install scons
    $ Tcmalloc library, http://code.google.com/p/gperftools/
	
  And you can try for Centos and Fedora:
    
    $ yum install git 
    $ yum install gcc
    $ yum install scons
    $ Tcmalloc library, http://code.google.com/p/gperftools/
    
  Finally, use Scons to build the HPG Aligner application:

    $ scons
  
RUNING
-------
	
  For run HPG Multi Aligner Framework , you must select first 'dna/rna' mode, next with '-c' command, select the command line of the mapper with '%I' for input file and '%O' for output file or output path. After that, select the input file with '-f' and output with '-o', if the mapper needs a temporal file to write the partial output, you can indicate it with the option '--tmp-file', and if you want to write the partial output in a temporal path (hard disk of the nodes) you can indicate it with the option '--tmp-path'. For RNA mode, you must indicate with '-i' option the path for BWT_INDEX or SA_INDEX of HPG Aligner.
  
  Example, run with DNA, Bowtie and two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner dna -c "bowtie2 -p 16 -S %O -x /work/user/genomes/bowtie/hs.73 %I" -f /work/user/datasets/40M_300nt_r0.001.bwa.read1.fastq -o /work/user/final-ouput --tmp-path /tmp/patial-output --tmp-file partial-alignments.sam

  Example, run with RNA, Tophat and two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work/genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output

  If you want process with the second phase to improve the alignments results, you can activate this with '--second-phase' option. The reads no mapped or incompleted mapped will be remapped with HPG Aligner, if you want other mapper you can indicate this with '--second-command' option.

  Example, run with RNA, Tophat, second phase and two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work//genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output --second-phase

  Example, run with RNA, Tophat, second phase with other mapper and two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work//genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output --second-phase --second-command "python mapsplice.py -c /work/user/gneomes/map_splice/ -x /work/user/genomes/map_splice/hs.73 -1 %I -p 16 -o %O" 


  For generate HPG SA Index:

    $ ./bin/hpg-aligner build-sa-index -g /hpg/genome/human/GRCH_37_ens73.fa -i /tmp/sa-index-human73/ 

  For generate HPG BWT Index:

    $ ./bin/hpg-aligner build-bwt-index -g /home/user/Homo_sapiens.fa -i /home/user/INDEX/  -r 8

ATENTION: 
     Mvapich2 has an incompatibility with tcmalloc. You can install mvapich2 with "--disable-registration-cache" option or before you run HPG Aligner mpi do:

    $ export MV2_USE_LAZY_MEM_UNREGISTER=0


  If you compile MPI HPG Multi Aligner Framework without tcmalloc the perferomance will be degraded

