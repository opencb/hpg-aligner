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
#include "wanderer/wanderer.h"

int align_launch(char *reference, char *bam, char *output, int threads);
ERROR_CODE wander_bam_file(char *bam_path, char *ref_name, char *ref_path, char *outbam);

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
	LOG_LEVEL(LOG_INFO_LEVEL);
	//LOG_LEVEL(LOG_ERROR_LEVEL);

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
	//alig_bam_file(bamc, base, dir, output);
	wander_bam_file(bamc, base, dir, output);
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

	printf("\n====== Times per alignment ======\n");
	time_get_mean_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, &max);
	printf("Time used for region extraction -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, &max);
	printf("Time used for haplotype extraction -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, &max);
	printf("Time used for realign -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	printf("\n====== Iterations ======\n");

	time_get_mean_slot(D_SLOT_IT_PROCESS, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_IT_PROCESS, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_IT_PROCESS, TIME_GLOBAL_STATS, &max);
	printf("Time used for process -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_IT_READ, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_IT_READ, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_IT_READ, TIME_GLOBAL_STATS, &max);
	printf("Time used for read -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_IT_WRITE, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_IT_WRITE, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_IT_WRITE, TIME_GLOBAL_STATS, &max);
	printf("Time used for write -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);


	printf("\n====== I/O ======\n");
	time_get_mean_slot(D_SLOT_ALIG_READ, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_ALIG_READ, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_ALIG_READ, TIME_GLOBAL_STATS, &max);
	printf("Time used for read alignment (mean) -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

	time_get_mean_slot(D_SLOT_ALIG_WRITE, TIME_GLOBAL_STATS, &mean);
	time_get_min_slot(D_SLOT_ALIG_WRITE, TIME_GLOBAL_STATS, &min);
	time_get_max_slot(D_SLOT_ALIG_WRITE, TIME_GLOBAL_STATS, &max);
	printf("Time used for write alignment (mean) -> %.2f us - min/max = %.2f/%.2f\n",
			mean*1000000.0, min*1000000.0, max*1000000.0);

#endif

	stop_log();

	return 0;
}

int
realign_wanderer(bam_wanderer_t *wanderer, bam_region_t *current_region, bam1_t *read)
{
	int i, err;
	int processed = 0;
	size_t length;

	//Current region
	int reg_valid = 0;
	size_t init_pos = SIZE_MAX;
	size_t end_pos = SIZE_MAX;
	size_t aux_init_pos;
	size_t aux_end_pos;
	size_t read_pos;

	assert(wanderer);
	assert(current_region);
	assert(read);

	//Filter read
	if(filter_read(read, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK))
	{
		//Read is not valid for process
		return WANDER_READ_FILTERED;
	}

	//Get read position
	read_pos = read->core.pos;

	//Has region an init position?
	/*if(current_region->init_pos == SIZE_MAX)
	{
		current_region->chrom = read->core.tid;
		current_region->init_pos = read_pos;
	}*/

	//Inside this region?
	if(current_region->end_pos != SIZE_MAX)
	{
		if(	current_region->chrom != read->core.tid
				|| current_region->end_pos < read->core.pos)
		{
			LOG_INFO_F("READ %d:%d\n",
					read->core.tid + 1, read->core.pos + 1);
			//Not in window region
			return WANDER_REGION_CHANGED;
		}
	}

	//Get interval for this alignment
	err = region_get_from_bam1(read, &aux_init_pos, &aux_end_pos);
	if(err)
	{
		LOG_ERROR_F("Trying to get region from invalid read: %s\n", bam1_qname(read));
		return INVALID_INPUT_BAM;
	}

	//This alignment have an interval?
	if(aux_end_pos != SIZE_MAX)
	{
		//Interval found

		//Update region chrom
		current_region->chrom = read->core.tid;

		//Update region start position
		if(current_region->init_pos == SIZE_MAX || current_region->init_pos > aux_init_pos)
		{
			current_region->init_pos = aux_init_pos;
		}

		//Update region end position
		if(current_region->end_pos == SIZE_MAX || current_region->end_pos < aux_end_pos)
		{
			current_region->end_pos = aux_end_pos;
		}
	}

	return NO_ERROR;
}

int
realign_processor(bam_wanderer_t *wanderer, bam_region_t *region)
{
	int i;

	assert(wanderer);
	assert(region);

	LOG_INFO_F("Processing over window %d:%d-%d with %d reads\n", region->chrom + 1, region->init_pos + 1, region->end_pos +1, region->size);

	for(i = 0; i < region->size; i++)
	{
		assert(region->reads[i]);
	}

	return NO_ERROR;
}

ERROR_CODE
wander_bam_file(char *bam_path, char *ref_name, char *ref_path, char *outbam)
{
	//Files
	bam_file_t *bam_f;
	bam_file_t *out_bam_f;
	genome_t* ref;
	int bytes;

	//Wanderer
	bam_wanderer_t wanderer;

	assert(bam_path);
	assert(ref_name);
	assert(ref_path);
	assert(outbam);

	//Open bam
	{
		printf("Opening BAM from \"%s\" ...\n", bam_path);
		bam_f = bam_fopen(bam_path);
		assert(bam_f);
		printf("BAM opened!...\n");
	}

	//Open reference genome
	{
		printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
		ref = genome_new(ref_name, ref_path);
		assert(ref);
		printf("Reference opened!...\n");
	}

	//Create new bam
	{
		printf("Creating new bam file in \"%s\"...\n", outbam);
		//init_empty_bam_header(orig_bam_f->bam_header_p->n_targets, recal_bam_header);
		out_bam_f = bam_fopen_mode(outbam, bam_f->bam_header_p, "w");
		assert(out_bam_f);
		bam_fwrite_header(out_bam_f->bam_header_p, out_bam_f);
		out_bam_f->bam_header_p = NULL;
		printf("New BAM initialized!...\n");
	}

	LOG_VERBOSE(1);
	//LOG_LEVEL(LOG_WARN_LEVEL);

	//Init wandering
	bwander_init(&wanderer);

	//Configure wanderer
	bwander_configure(&wanderer, bam_f, out_bam_f,
			(int (*)(void *, bam_region_t *, bam1_t *))realign_wanderer,
			(int (*)(void *, bam_region_t *))realign_processor);

	//Run wander
	bwander_run(&wanderer);

	//Destroy wanderer
	bwander_destroy(&wanderer);

	//Free context
	//printf("\nDestroying context...\n");
	//fflush(stdout);
	//alig_destroy(&context);

	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	printf("BAM closed.\n");

	printf("\nClosing reference file...\n");
	genome_free(ref);
	printf("Reference closed.\n");

	printf("Closing \"%s\" BAM file...\n", outbam);
	bam_fclose(out_bam_f);
	printf("BAM closed.\n");

	printf("Realignment around indels DONE.\n");

	return NO_ERROR;
}

