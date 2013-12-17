/*
 * alig_launcher.c
 *
 *  Created on: Dec 13, 2013
 *      Author: rmoreno
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <argtable2.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <sys/stat.h>

int align_launch(char *reference, char *bam, char *output);

int alig_bam(int argc, char **argv)
{

    struct arg_file *refile = arg_file0("r",NULL,"<reference>","reference genome compressed file (dna_compression.bin)");
    struct arg_file *infile = arg_file0("b",NULL,"<input>","input BAM file");
    struct arg_file *outfile = arg_file0("o",NULL,"<output>","output recalibrated BAM file, default:\"output.bam\"");
    struct arg_lit  *help    = arg_lit0("h","help","print this help and exit");
    struct arg_lit  *version = arg_lit0(NULL,"version","print version information and exit");
    struct arg_end  *end     = arg_end(20);

    void* argtable[] = {refile,infile,outfile,help,version,end};
    const char* progname = "hpg-bam a";
    int nerrors;
    int exitcode=0;

    /* verify the argtable[] entries were allocated sucessfully */
    if (arg_nullcheck(argtable) != 0)
    {
        /* NULL entries were detected, some allocations must have failed */
        printf("%s: insufficient memory\n",progname);
        exitcode=1;
		goto exit;
    }

    /* set any command line default values prior to parsing */
    outfile->filename[0]="output.bam";

    /* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
    {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        printf("This program process BAM files to generate local realigned BAM files.\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        exitcode=0;
        goto exit;
    }

    /* special case: '--version' takes precedence error reporting */
    if (version->count > 0)
    {
        printf("'%s' genomics.\n",progname);
        printf("2013, Raúl Moreno - OpenCB\n");
        exitcode=0;
        goto exit;;
    }

	/* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto exit;
    }

    /* special case: uname with no command line options induces brief help */
    if (argc==1)
    {
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=0;
        goto exit;
    }

    /* Incorrect case: Null or more reference files */
    if (refile->count != 1)
    {
		printf("Please, specify one reference file.\n");
		exitcode = 1;
		goto exit;
	}

    //Incorrect case: No input files
	if (infile->count == 0)
	{
		printf("Please, specify one input files.\n");
		return 1;
	}

	//Incorrect case: Input file dont exists
	if (infile->count > 0)
	{
		FILE *f = NULL;
		f = fopen(infile->filename[0], "r");
		if(f)
		{
			fclose(f);
		}
		else
		{
			printf("Inexistent input file \"%s\".\n", infile->filename[0]);
			return 1;
		}
	}

	//Incorrect case: More than one output files
	if (outfile->count > 1)
	{
		printf("Please, specify only one output file.\n");
		return 1;
	}

	//Incorrect case: No reference in phase 1
	if (refile->count == 0)
	{
		printf("Please, specify one reference file with -r.\n");
		return 1;
	}

	//Incorrect case: Reference file dont exists
	if (refile->count > 0)
	{
		FILE *f = NULL;
		f = fopen(refile->filename[0], "r");
		if(f)
		{
			fclose(f);
		}
		else
		{
			printf("Inexistent reference file \"%s\".\n", refile->filename[0]);
			return 1;
		}
	}

    /* normal case: realignment */
	{
		exitcode = align_launch(refile->filename[0], infile->filename[0], outfile->filename[0]);
	}

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int
align_launch(char *reference, char *bam, char *output)
{
	char *dir, *base, *bamc, *outputc, *infofilec, *datafilec;

	assert(reference);
	assert(bam);
	assert(output);

	init_log();

	//Obtain reference filename and dirpath from full path
	dir = strdup(reference);
	dir = dirname(dir);
	base = strrchr(reference, '/');

	//Obtain data from bam
	bamc = strdup(bam);
	printf("Reference dir: %s\n", dir);
	printf("Reference name: %s\n", base);
	alig_bam_file2(bamc, base, dir);
	free(bamc);
	free(dir);

	stop_log();

	return 0;
}



