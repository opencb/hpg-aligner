aligner
=======

Welcome to NGS data aligner

'bioinfo-libs' is a C library with bioinformatic data format parsers and BWT and SW aligners. It's part of HPG project.


COMPONENTS
This library consist of:
 * several data format parsers: gff, bam, vcf, bed, ...
 * HPC implementation of two widely used sequence aligners: Burrows-Whheler Transform (BWT) and Smith-Waterman (SW)


CONTACT
You can contact any of the following developers:
 * Joaquin Tarraga (jtarraga@cipf.es)
 * Hector Martinez (martineh@cipf.es)
 * Ignacio Medina (imedina@cipf.es)


BUILDING
Scons is used as building system, to build the library download 'next' branch (for the most recent library) and execute:

 $ scons

GCC4.4+ is the only requirement.


DOCUMENTATION
You can find more info about HPG project and bioinfo-libs at:

 http://wiki.opencb.org/projects/hpg/doku.php
