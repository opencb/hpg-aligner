MPI HPG Aligner generic framework README
========================================

Welcome to MPI HPG Aligner generic framework!

HPG Aligner is an ultra fast and highly sensitive Next-Generation
Sequencing (NGS) DNA-Seq and RNA-Seq aligner. MPI HPG Aligner generic
framework is an alignment framework that can be leveraged to
coordinately run any *single-node* DNA-Seq or RNA-Seq aligner, taking
advantage of the resources of a cluster, without having to modify any
portion of the original software.


Components
----------

* ``hpg-aligner`` - An ultra fast and highly sensitive NGS DNA-Seq and
  RNA-Seq aligner.
* ``hpg-multialigner`` - An MPI generic-aligner framework.


Contact
-------

If required, please contact Héctor Martínez <martineh@uji.es> for
assistance.


Installation
------------

HPG Aligner is open source and can be freely downloaded from
GitHub. To download the MPI HPG Aligner generic framework with
overlapped input and output version, please issue the following
commands:

    $ git clone https://github.com/opencb/hpg-aligner.git
    $ cd hpg-aligner
    $ git submodule update --init
    $ cd lib/hpg-libs
    $ git checkout develop-multi-fifo-output
    $ cd ../..
    $ git checkout develop-multi-fifo-output

To build MPI HPG Aligner generic framework, the next software should
be installed first: GCC, Scons, and the TCMALLOC library. To install
GCC and Scons on Ubuntu and Debian, please use:

    $ sudo apt-get install gcc scons
	
On CentOs and Fedora:
    
    $ sudo yum install gcc scons

To install the TCMALLOC library, please follow the instructions at
<https://github.com/gperftools/gperftools/>. Once installed, please
update the TCMALLOC include and library absolute paths defined in the
``hpg-aligner/SConstruct`` file.

Once GCC, Scons, and TCMALLOC are installed, MPI HPG Aligner generic
framework can be built using the following command on the
``hpg-aligner`` directory:

    $ scons


Usage
-----

For running the MPI HPG Aligner generic framework, you must use
``mpirun`` to execute ``/bin/hpg-multialigner`` with the following
options:

1. Which kind of sequencing should be performed: ``dna`` or ``rna``.

2. The command line of the mapper to be used on the first phase
   (``-c`` option).The mapper command line must use ``%I`` for the
   input file and ``%O`` for the output file or path.

3. A temporary file path, if the mapper requires
   one (``--tmp-file`` option).

4. Select the input file (``-f`` option).

5. Select the output (``-o`` option).

6. Optionally, a path to the temporary nodes output (``--tmp-path``
   option). For increased performance, this path should be local to the
   nodes.

7. Optionally, a second phase that will try to map those reads no
   mapped by the aligner called on the first phase can be activated
   (``--second-phase`` option). By default, ``hpg-aligner`` will be
   used on this second phase.

8. The BWT INDEX or SA INDEX path must be indicated (`-i` option).

9. In case a different aligner is to be used on the second phase, its
   command line should be indicated (``--second-command`` option).


### Usage examples

DNA-Seq, Bowtie, two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner dna -c "bowtie2 -p 16 -S %O -x /work/user/genomes/bowtie/hs.73 %I" -f /work/user/datasets/40M_300nt_r0.001.bwa.read1.fastq -o /work/user/final-ouput --tmp-path /tmp/patial-output --tmp-file partial-alignments.sam


RNA-Seq, Tophat, two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work/genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output


RNA-Seq, Tophat first, HPG Aligner next, two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work//genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output --second-phase


RNA-Seq, Tophat first, MapSplice next, two nodes:

    $ mpirun -np 3 -hosts compute-0,compute-1,compute-0 ./bin/hpg-multialigner rna -c "tophat2 --no-convert-bam -p 16 --no-sort-bam -o %O /work//genomes/bowtie/hs.73 %I" -i /work/user/bwt-index/ -f /work/user/datasets/10M_100nt_r0.001.rna.fastq -o /work/user/tophat.out --tmp-path /work/user/partial-output --second-phase --second-command "python mapsplice.py -c /work/user/gneomes/map_splice/ -x /work/user/genomes/map_splice/hs.73 -1 %I -p 16 -o %O" 


Generate HPG SA Index:

    $ ./bin/hpg-aligner build-sa-index -g /hpg/genome/human/GRCH_37_ens73.fa -i /tmp/sa-index-human73/ 


Generate HPG BWT Index:

    $ ./bin/hpg-aligner build-bwt-index -g /home/user/Homo_sapiens.fa -i /home/user/INDEX/  -r 8


MVAPICH2 considerations
-----------------------

Please, be aware that MVAPICH2 and TCMALLOC does not play nice
together by default. In order to use MVAPICH2 with MPI HPG Aligner
generic framework, that relies on TCMALLOC, either install MVAPICH2
with the ``--disable-registration-cache`` option, or execute the
following command before calling ``mpirun``:

	$ export MV2_USE_LAZY_MEM_UNREGISTER=0

Not doing any of the previous could lead to inconsistent results.
