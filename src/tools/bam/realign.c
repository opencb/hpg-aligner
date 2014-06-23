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
#include "aux/aux_common.h"
#include "bfwork/bam_file_ops.h"
//#include "bfwork/bfwork.h"
//#include "aligner/alig.h"

int align_launch(const char *reference, const char *bam, const char *output, const char *stats_path, int recalibration);
//ERROR_CODE wander_bam_file(char *bam_path, char *ref_path, char *outbam);

int alig_bam(int argc, char **argv)
{
	char *refc, *bamc, *outputc, *statsc;

	struct arg_file *recal = arg_lit0(NULL,"recalibrate","performs a base quality recalibration");
    struct arg_file *refile = arg_file0("r",NULL,"<reference>","reference genome compressed file (dna_compression.bin)");
    struct arg_file *infile = arg_file0("b",NULL,"<input>","input BAM file");
    struct arg_file *outfile = arg_file0("o",NULL,"<output>","output processed BAM file");
    struct arg_file *stats = arg_file0("s",NULL,"<stats>","folder to store timing stats");
    struct arg_lit  *help    = arg_lit0("h","help","print this help and exit");
    struct arg_lit  *version = arg_lit0(NULL,"version","print version information and exit");
    struct arg_end  *end     = arg_end(20);

    void* argtable[] = {recal,refile,infile,outfile,stats,help,version,end};
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
    outfile->filename[0]=NULL;

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

	bamc = NULL;
	if(infile->count > 0)
		bamc = strdup(infile->filename[0]);

	outputc = NULL;
	if(outfile->count > 0)
		outputc = strdup(outfile->filename[0]);

	refc = NULL;
	if(refile->count > 0)
		refc = strdup(refile->filename[0]);

	statsc = NULL;
	if(stats->count > 0)
		statsc = strdup(stats->filename[0]);

    /* normal case: realignment */
	{
		exitcode = align_launch(	refc, bamc, outputc, statsc,recal->count);
	}

	if(bamc)
		free(bamc);
	if(outputc)
		free(outputc);
	if(refc)
		free(refc);
	if(statsc)
		free(statsc);

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int
align_launch(const char *reference, const char *bam, const char *output, const char *stats_path, int recalibration)
{
	int err, threads;
	char *sched;
	double init_time, end_time;

	assert(reference);
	assert(bam);

	init_log();
	LOG_FILE("realign.log","w");
	LOG_VERBOSE(1);
	LOG_LEVEL(LOG_WARN_LEVEL);

	//Set schedule if not defined
	setenv("OMP_SCHEDULE", "static", 0);
	sched = getenv("OMP_SCHEDULE");
	threads = omp_get_max_threads();

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

	if(recalibration)
	{
		printf("Executing realign and recalibration\n");

		//Realign and recalibration
		alig_recal_bam_file(bam, reference, NULL, NULL, output, 200, NULL);
	}
	else
	{
		printf("Executing realign\n");

		//Only realign
		alig_bam_file(bam, reference, output, stats_path);
	}

	stop_log();

	return 0;
}

/*int
dummy_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int i;

	assert(fwork);
	assert(region);

	LOG_INFO_F("Processing over region %d:%d-%d with %d reads\n", region->chrom + 1, region->init_pos + 1, region->end_pos +1, region->size);

	for(i = 0; i < region->size; i++)
	{
		assert(region->reads[i]);
	}

	//Update user data
	size_t *reads;
	bfwork_lock_user_data(fwork, (void **)&reads);
	*reads += region->size;
	bfwork_unlock_user_data(fwork);

	//Get local user data
	size_t *regions;
	bfwork_local_user_data(fwork, (void **)&regions);
	if(regions == NULL)
	{
		//Local data is not initialized
		regions = (size_t *)malloc(sizeof(size_t));
		*regions = 0;
		//printf("%d - POP\n", omp_get_thread_num());

		//Set local data
		bfwork_local_user_data_set(fwork, regions);
	}

	//Increase local user data
	*regions += 1;

	//usleep(10 * region->size); //1 us per alignment

	return NO_ERROR;
}
*/
/*void
reduce_reads_dummy(void *data, void *dest)
{
	size_t *s_data, *s_dest;
	assert(data);
	assert(dest);

	//Cast pointers
	s_data = (size_t *)data;
	s_dest = (size_t *)dest;

	//Add
	*s_dest += *s_data;
}*/

/*static inline ERROR_CODE
wander_bam_file_dummy(bam_fwork_t *fwork, const char *in, const char *out, const char *ref)
{
	bfwork_context_t context;
	size_t reads_proc = 0;
	size_t regions_proc = 0;

	assert(fwork);
	assert(in);
	assert(out);
	assert(ref);

	//Init wandering
	bfwork_init(fwork);

	//Create realigner context
	bfwork_context_init(&context,
			(int (*)(void *, bam_region_t *, bam1_t *))realign_wanderer,
			(int (*)(void *, bam_region_t *))dummy_processor);

#ifdef D_TIME_DEBUG
	//Init timing
	bfwork_context_init_timing(&context, "dummy");
#endif

	//Set context user data
	bfwork_context_set_user_data(&context, &reads_proc);

	//Configure wanderer
	bfwork_configure(fwork, in, out, ref, &context);

	//Run wander
	bfwork_run(fwork);

	//Print user data
	printf("\nREDUCED DATA, Reads processed: %d\n", reads_proc);

	//Print reduced local user data
	bfwork_context_local_user_data_reduce(&context, &regions_proc, reduce_reads_dummy);
	printf("\nREDUCED LOCAL DATA, Regions processed: %d\n", regions_proc);

	//Free local user data
	bfwork_context_local_user_data_free(&context, NULL);

#ifdef D_TIME_DEBUG
	//Print times
	bfwork_context_print_times(&context);

	//Destroy wanderer time
	bfwork_context_destroy_timing(&context);
#endif

	//Destroy wanderer
	bfwork_destroy(fwork);

	//Destroy context
	bfwork_context_destroy(&context);

	return NO_ERROR;
}*/

/*ERROR_CODE
wander_bam_file(char *bam_path, char *ref_path, char *outbam)
{
	//Times
	double times;

	//Wanderer
	bam_fwork_t fwork;

	assert(bam_path);
	assert(ref_path);

	//Dummy wandering
	//wander_bam_file_dummy(&fwork, bam_path, outbam, ref_path);

	return NO_ERROR;
}*/

