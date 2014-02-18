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
#include <errno.h>

#include "aux/timestats.h"
#include "aligner/alig.h"

int align_launch(char *reference, char *bam, char *output, int threads);

int alig_bam(int argc, char **argv)
{

    struct arg_file *refile = arg_file0("r",NULL,"<reference>","reference genome compressed file (dna_compression.bin)");
    struct arg_file *infile = arg_file0("b",NULL,"<input>","input BAM file");
    struct arg_file *outfile = arg_file0("o",NULL,"<output>","output recalibrated BAM file, default:\"output.bam\"");
    struct arg_lit  *help    = arg_lit0("h","help","print this help and exit");
    struct arg_lit  *version = arg_lit0(NULL,"version","print version information and exit");
    struct arg_end  *end     = arg_end(20);

    void* argtable[] = {refile,infile,outfile,help,version,end};
    const char* progname = "hpg-bam realign v"ALIG_VER;
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
		exitcode = align_launch((char *)refile->filename[0], (char *)infile->filename[0], (char *)outfile->filename[0], 1);
	}

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int
align_launch(char *reference, char *bam, char *output, int threads)
{
	char *dir, *base, *bamc, *outputc, *infofilec, *datafilec, *sched;
	int err;
	double init_time, end_time;

	assert(reference);
	assert(bam);
	assert(output);

	init_log();
	LOG_FILE("realign.log","w");
	LOG_LEVEL(LOG_ERROR_LEVEL);

	//Set schedule if not defined
	setenv("OMP_SCHEDULE", "static", 0);
	sched = getenv("OMP_SCHEDULE");

	//Time measures
#ifdef D_TIME_DEBUG

	char filename[100];
	char intaux[20];
	char cwd[1024];

	//Initialize stats
	if(time_new_stats(20, &TIME_GLOBAL_STATS))
	{
		printf("ERROR: FAILED TO INITIALIZE TIME STATS\n");
	}


	if (getcwd(cwd, sizeof(cwd)) != NULL)
	{
		printf("Current working dir: %s\n", cwd);
	}
	else
	{
		perror("WARNING: getcwd() dont work");
	}

	strcpy(filename, cwd);
	strcat(filename,"/stats/");
	/*if(sched)
		strcat(filename,sched);
	else
	{
		printf("ERROR: Obtaining OMP_SCHEDULE environment value\n");
	}*/

	//Create stats directory
	printf("Creating stats directory: %s\n", filename);
	err = mkdir(filename, S_IRWXU);
	if(err != 0 && errno != EEXIST)
	{
		perror("WARNING: failed to create stats directory");
	}
	else
	{
		//strcat(filename,"_");
		//sprintf(intaux, "%d", MAX_BATCH_SIZE);
		//strcat(filename, intaux);
		strcat(filename, "_");
		sprintf(intaux, "%d", threads);
		strcat(filename, intaux);
		strcat(filename, "_.stats");

		//Initialize stats file output
		if(time_set_output_file(filename, TIME_GLOBAL_STATS))
		{
			printf("ERROR: FAILED TO INITIALIZE TIME STATS FILE OUTPUT\n");
		}

		printf("STATISTICS ACTIVATED, output file: %s\n\n", filename);
	}

#endif

	//System info
	printf("------------\n");
	printf("System info:\n");
#ifdef __SSE2__	//SSE Block
	printf("SSE2 Activated\n");
#endif
#ifdef __SSSE3__
	printf("Extended SSE3 Activated\n");
#endif
	printf("------------\n");

	//Obtain reference filename and dirpath from full path
	dir = strdup(reference);
	dir = dirname(dir);
	base = strrchr(reference, '/');

	//Obtain data from bam
	bamc = strdup(bam);
	printf("Reference dir: %s\n", dir);
	printf("Reference name: %s\n", base);
#ifdef D_TIME_DEBUG
	init_time = omp_get_wtime();
#endif
	alig_bam_file(bamc, base, dir);
#ifdef D_TIME_DEBUG
	end_time = omp_get_wtime();
	time_add_time_slot(D_SLOT_TOTAL, TIME_GLOBAL_STATS, end_time - init_time);
#endif
	free(bamc);
	free(dir);

	//Print times
#ifdef D_TIME_DEBUG
	double min, max, mean;

	//Print time stats
	printf("----------------------------\nTIME STATS: \n");

	printf("\n====== General times ======\n");
	time_get_mean_slot(D_SLOT_TOTAL, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_TOTAL, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_TOTAL, TIME_GLOBAL_STATS, &max);
	printf("Total time to realign -> %.2f s - min/max = %.2f/%.2f\n",
			mean, min, max);

	time_get_mean_slot(D_SLOT_INIT, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_INIT, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_INIT, TIME_GLOBAL_STATS, &max);
	printf("Time used to initialize aligner -> %.2f s - min/max = %.2f/%.2f\n",
			mean, min, max);

	/*time_get_mean_slot(D_SLOT_PROCCESS, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_PROCCESS, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_PROCCESS, TIME_GLOBAL_STATS, &max);
	printf("Time used to proccess every read -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);*/

	printf("\n====== Haplotype ======\n");
	time_get_mean_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &max);
	printf("Time used for read region extraction -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &max);
	printf("Time used for read haplotype extraction -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &max);
	printf("Time used for read realign per haplotype -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);


	printf("\n====== I/O ======\n");
	time_get_mean_slot(D_SLOT_READ, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_READ, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_READ, TIME_GLOBAL_STATS, &max);
	printf("Time used for read alignment (mean) -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_WRITE, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_WRITE, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_WRITE, TIME_GLOBAL_STATS, &max);
	printf("Time used for write alignment (mean) -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

#endif

	stop_log();

	return 0;
}



