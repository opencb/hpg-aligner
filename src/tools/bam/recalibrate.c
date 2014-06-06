/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <argtable2.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <sys/stat.h>
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "wanderer/wanderer.h"
#include "recalibrate/recal_config.h"
#include "recalibrate/recal_structs.h"
#include "recalibrate/bam_recal_library.h"

#define RECALIBRATE_COLLECT 			0x01
#define RECALIBRATE_RECALIBRATE 	0x02

ERROR_CODE wander_bam_file_(uint8_t flags, char *bam_path, char *ref_name, char *ref_path, char *data_file, char *info_file, char *outbam, int cycles);

int mymain(	int full,
			int p1,
			int p2,
			int cycles,
			int threads,
			const char *reference, 
			int refcount, 
			const char *input, 
			int incount,
			const char *output, 
			int outcount,
			const char *datafile, 
			int datacount,
			const char *infofile,
			int infocount,
			const char *compfile,
			int compcount)
{	
	char *dir, *base, *inputc, *outputc, *infofilec, *datafilec;
	char *sched;
	char cwd[1024];
	
	//Incorrect case: No input files
    if (incount == 0)
    {
		printf("Please, specify one input files.\n");
		return 1;
	}

    //Incorrect case: Input file dont exists
    if (incount == 1)
    {
    	FILE *f = NULL;
    	f = fopen(input, "r");
    	if(f)
    	{
    		fclose(f);
    	}
    	else
    	{
    		printf("Inexistent input file \"%s\".\n", input);
    		return 1;
    	}
    }
	
	//Incorrect case: More than one output files
    if (outcount > 1)
    {
		printf("Please, specify only one output file.\n");
		return 1;
	}
	
	//Incorrect case: If only phase 2, then datafile must be specified
    if (!p1 && p2 && datacount == 0)
    {
		printf("Please, when using \"-2\" argument, specify valid datafile.\n");
		return 1;
	}
	
	//Default case: No phase selected
    if (!full && !p1 && !p2 && !compcount)
    {
		printf("No recalibration phase specified or comparation, executing full recalibration.\n");
		full = 1;
	}

	//If full recalibration, execute two phases
	if (full)
	{
		p1 = 1;
		p2 = 1;
	}
	
	//Incorrect case: No reference in phase 1
	if (p1 && refcount == 0)
	{
		printf("Please, specify reference file with phase 1.\n");
		return 1;
	}

	//Incorrect case: Reference file dont exists
	if (p1 && refcount > 0)
	{
		FILE *f = NULL;
		f = fopen(reference, "r");
		if(f)
		{
			fclose(f);
		}
		else
		{
			printf("Inexistent reference file \"%s\".\n", reference);
			return 1;
		}
	}

	//Set schedule if not defined
	setenv("OMP_SCHEDULE", "static", 0);
	sched = getenv("OMP_SCHEDULE");

	//Time measures
	#ifdef D_TIME_DEBUG

	char filename[100];
	char intaux[20];

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
			perror("WARNING: getcwd() dont work\n");
		}

		strcpy(filename, cwd);
		strcat(filename,"/stats/");
		if(sched)
			strcat(filename,sched);
		else
		{
			printf("ERROR: Obtaining OMP_SCHEDULE environment value\n");
		}

		//Create stats directory
		mkdir(filename, S_IRWXU);

		strcat(filename,"_");
		sprintf(intaux, "%d", MAX_BATCH_SIZE);
		strcat(filename, intaux);
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

	#endif
	
	//Printf proc caps
	printf_proc_features();
	#ifdef __MMX__
	printf("Using MMX features\n");
	#endif
	#ifdef __SSE__
	printf("Using SSE features\n");
	#endif
	#ifdef __SSE2__
	printf("Using SSE2 features\n");
	#endif
	#ifdef __SSE3__
	printf("Using SSE3 features\n");
	#endif

	init_log();

	//Set num of threads
	//omp_set_num_threads(threads);
	printf("Threading with %d threads and %s schedule\n", omp_get_max_threads(), sched);

	//Execute phase 1
	if (p1)
	{
		printf("\n===================\nExecuting phase 1.\n===================\n\n");
		
		//Obtain reference filename and dirpath from full path
		dir = strdup(reference);
		dir = dirname(dir);
		base = strrchr(reference, '/');
		
		//Print data
		infofilec = NULL;
		if(infocount > 0)
		{
			infofilec = strdup(infofile);
		}
		
		//Save data file
		datafilec = NULL;
		if(datacount > 0)
		{
			datafilec = strdup(datafile);
		}

		//Obtain data from bam
		inputc = strdup(input);
		printf("Reference dir: %s\n", dir);
		printf("Reference name: %s\n", base);
		//recal_get_data_from_file(inputc, base, dir, data);
		wander_bam_file_(RECALIBRATE_COLLECT, inputc,base, dir, datafilec, infofilec, NULL, cycles);
		free(inputc);
		if(datafilec)
			free(datafilec);
		if(infofilec)
			free(infofilec);
	}
	
	//Execute phase 2
	if (p2)
	{
		printf("\n===================\nExecuting phase 2.\n===================\n\n");
		
		//Get datafile name
		datafilec = NULL;
		if(datacount > 0)
		{
			datafilec = strdup(datafile);
		}

		inputc = strdup(input);
		outputc = NULL;
		if(output)
			outputc = strdup(output);

		//Recalibrate bam
		//recal_recalibrate_bam_file(inputc, data, outputc);
		wander_bam_file_(RECALIBRATE_RECALIBRATE, inputc, NULL, NULL, datafilec, NULL, outputc, cycles);

		free(inputc);
		if(outputc)
			free(outputc);
		if(datafilec)
			free(datafilec);
	}
	
	
	//Print times
	#ifdef D_TIME_DEBUG
		double min, max, mean;

		//Print time stats
		printf("----------------------------\nTIME STATS: \n");
			
		//Times from phase 1
		if (p1)
		{	
			printf("\n-------\n");
			printf("PHASE 1\n");
			printf("-------\n");

			printf("\n====== Total time to collect ======\n");
			time_get_min_slot(D_SLOT_PH1_COLLECT_BAM, TIME_GLOBAL_STATS, &min);
			printf("TIME USED FOR BAM COLLECTION -> %.2f s\n", min);

			printf("\n====== Deltas proccess ======\n");
			time_get_min_slot(D_SLOT_PROCCESS_DELTAS, TIME_GLOBAL_STATS, &min);
			printf("Time used for deltas proccess -> %.2f ms\n", min*1000.0);

			printf("\n====== Batch iteration ======\n");
			time_get_mean_slot(D_SLOT_PH1_ITERATION, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH1_ITERATION, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH1_ITERATION, TIME_GLOBAL_STATS, &max);
			printf("Time used for iteration (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
					mean*1000.0, min*1000.0, max*1000.0);

			printf("\n====== Inside iteration ======\n");
			time_get_mean_slot(D_SLOT_PH1_READ_BATCH, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH1_READ_BATCH, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH1_READ_BATCH, TIME_GLOBAL_STATS, &max);
			printf("(P) Time used for batch read (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
					mean*1000.0, min*1000.0, max*1000.0);

			time_get_mean_slot(D_SLOT_PH1_COLLECT_BATCH, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH1_COLLECT_BATCH, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH1_COLLECT_BATCH, TIME_GLOBAL_STATS, &max);
			printf("(P) Time used for batch collect (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
					mean*1000.0, min*1000.0, max*1000.0);

			time_get_mean_slot(D_SLOT_PH1_COLLECT_REDUCE_DATA, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH1_COLLECT_REDUCE_DATA, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH1_COLLECT_REDUCE_DATA, TIME_GLOBAL_STATS, &max);
			printf("Time used for reduce data (mean) -> %.2f micros - min/max = %.2f/%.2f\n",
					mean*1000000.0, min*1000000.0, max*1000000.0);

			printf("\n====== Alignment ======\n");
			time_get_mean_slot(D_SLOT_PH1_COLLECT_ALIG, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH1_COLLECT_ALIG, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH1_COLLECT_ALIG, TIME_GLOBAL_STATS, &max);
			printf("Time used for alig collect (mean) -> %.2f micros - min/max = %.2f/%.2f\n",
					mean*1000000.0, min*1000000.0, max*1000000.0);
		}	
		
		if (p2)
		{
			printf("\n-------\n");
			printf("PHASE 2\n");
			printf("-------\n");

			printf("\n====== Total time to recalibrate ======\n");
			time_get_min_slot(D_SLOT_PH2_RECALIBRATE, TIME_GLOBAL_STATS, &min);
			printf("TIME USED FOR BAM RECALIBRATION -> %.2f s\n", min);

			printf("\n====== Batch iteration ======\n");
			time_get_mean_slot(D_SLOT_PH2_ITERATION, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH2_ITERATION, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH2_ITERATION, TIME_GLOBAL_STATS, &max);
			printf("Time used for iteration (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
					mean*1000.0, min*1000.0, max*1000.0);

			printf("\n====== Inside iteration ======\n");
			time_get_mean_slot(D_SLOT_PH2_READ_BATCH, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH2_READ_BATCH, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH2_READ_BATCH, TIME_GLOBAL_STATS, &max);
			printf("(P) Time used for batch read (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
								mean*1000.0, min*1000.0, max*1000.0);

			time_get_mean_slot(D_SLOT_PH2_PROCCESS_BATCH, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH2_PROCCESS_BATCH, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH2_PROCCESS_BATCH, TIME_GLOBAL_STATS, &max);
			printf("(P) Time used for batch proccess (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
								mean*1000.0, min*1000.0, max*1000.0);

			time_get_mean_slot(D_SLOT_PH2_WRITE_BATCH, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH2_WRITE_BATCH, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH2_WRITE_BATCH, TIME_GLOBAL_STATS, &max);
			printf("(P) Time used for batch write (mean) -> %.2f ms - min/max = %.2f/%.2f\n",
								mean*1000.0, min*1000.0, max*1000.0);

			printf("\n====== Alignment proccess ======\n");
			//Print recalibrate stats
			time_get_mean_slot(D_SLOT_PH2_RECAL_ALIG, TIME_GLOBAL_STATS, &mean);
			time_get_min_slot(D_SLOT_PH2_RECAL_ALIG, TIME_GLOBAL_STATS, &min);
			time_get_max_slot(D_SLOT_PH2_RECAL_ALIG, TIME_GLOBAL_STATS, &max);
			printf("Time used for alig recalibration (mean) -> %.2f micros - min/max = %.2f/%.2f\n",
					mean*1000000.0, min*1000000.0, max*1000000.0);
		}
			
		//Free memory from stats
		//time_destroy_stats(&TIME_GLOBAL_STATS);
	#endif
	
	//Compare
	if(compcount)
	{
		if(p1 || p2)
		{
			compare_bams_qual(output, compfile, 76);
		}
		else
		{
			compare_bams_qual(input, compfile, 76);
		}
	}
	
	stop_log();

	return 0;
}


int recalibrate_bam(int argc, char **argv)
{
	int thr;

    //struct arg_lit  *recurse = arg_lit0("R",NULL,                       "recurse through subdirectories");
    //struct arg_int  *repeat  = arg_int0("k","scalar",NULL,              "define scalar value k (default is 3)");
    //struct arg_str  *defines = arg_strn("D","define","MACRO",0,argc+2,  "macro definitions");
    //struct arg_file *outfile = arg_file0("o",NULL,"<output>",           "output file (default is \"-\")");
    //struct arg_file *infiles = arg_filen(NULL,NULL,NULL,1,argc+2,       "input file(s)");
    //struct arg_end  *end     = arg_end(20);
    
    struct arg_lit  *full = arg_lit0("f",NULL,"execute full recalibration");
    struct arg_lit  *phase1 = arg_lit0("1",NULL,"execute only first part of recalibration");
    struct arg_lit  *phase2 = arg_lit0("2",NULL,"execute only second part of recalibration, need data file");
    struct arg_int  *cycles  = arg_int0(NULL,"cycles",NULL,"define number of cycles of bams");
    struct arg_int  *threads  = arg_int0(NULL,"num-threads",NULL,"define number of threads to use");
    struct arg_file *refile = arg_file0("r",NULL,"<reference>","reference genome compressed file (dna_compression.bin)");
    struct arg_file *infile = arg_file0("b",NULL,"<input>","input BAM file");
    struct arg_file *outfile = arg_file0("o",NULL,"<output>","output recalibrated BAM file, default:\"output.bam\"");
    struct arg_file *datafile = arg_file0("d",NULL,"<data>","data file containing recalibration information");
    struct arg_file *infofile = arg_file0("i",NULL,"<info>","info file (human readable) containing recalibration information");
    struct arg_file *compfile = arg_file0("c",NULL,"<compfile>","bam file to compare with recalibrated bam");
    struct arg_lit  *help    = arg_lit0("h","help","print this help and exit");
    struct arg_lit  *version = arg_lit0(NULL,"version","print version information and exit");
    struct arg_end  *end     = arg_end(20);
    
    void* argtable[] = {full,phase1,phase2,cycles,threads,refile,infile,outfile,datafile,infofile,compfile,help,version,end};
    const char* progname = "hpg-bam recalibrate";
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
        printf("This program process BAM files to generate recalibrated BAM files.\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        exitcode=0;
        goto exit;
    }
    
    /* special case: '--version' takes precedence error reporting */
    if (version->count > 0)
    {
        printf("'%s' genomics.\n",progname);
        printf("April 2013, Raúl Moreno\n");
        exitcode=0;
        goto exit;
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
    
    /* Incorrect case: Null or more reference files with -F -1 -2*/
    if (refile->count != 1 && (full->count > 0 || phase1->count > 0 || phase2->count > 0))
    {
		printf("Please, specify one reference file.\n");
		exitcode = 1;
		goto exit;
	}
	
	/* Incorrect case: No cycles*/
    if (cycles->count == 0)
    {
		printf("Please, specify number of cycles with -c.\n");
		exitcode = 1;
		goto exit;
	}
    
    //Default case: 1 thread
	if (threads->count == 0 || threads->ival[0] <= 0)
	{
		thr = 1;
	}
	else
	{
		thr = threads->ival[0];
	}


    /* normal case: take the command line options at face value */
    exitcode = mymain(	full->count,
						phase1->count,
						phase2->count,
						cycles->ival[0],
						thr,
						refile->filename[0], 
						refile->count, 
						infile->filename[0], 
						infile->count, 
						outfile->filename[0], 
						outfile->count,
						datafile->filename[0],
						datafile->count,
						infofile->filename[0],
						infofile->count,
						compfile->filename[0],
						compfile->count);

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int
recalibrate_wanderer(bam_wanderer_t *wanderer, bam_region_t *region, bam1_t *read)
{
	assert(wanderer);
	assert(region);
	assert(read);

	//Filter read
	if(filter_read(read, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK))
	{
		//Read is not valid for process
		return WANDER_READ_FILTERED;
	}

	//Update region bounds
	if(region->init_pos > read->core.pos)
	{
		region->init_pos = read->core.pos;
		region->chrom = read->core.tid;
	}
	if(region->end_pos < read->core.pos)
	{
		region->end_pos = read->core.pos;
	}

	return NO_ERROR;
}

int
recalibrate_collect_processor(bam_wanderer_t *wanderer, bam_region_t *region)
{
	int err, i;
	recal_info_t *data;
	recal_data_collect_env_t *collect_env;
	bam1_t *read;
	size_t *cycles;

	//Get data
	bwander_local_user_data(wanderer, (void **)&data);
	if(data == NULL)
	{
		//Local data is not initialized
		data = (recal_info_t *)malloc(sizeof(recal_info_t));

		//Lock cycles
		bwander_lock_user_data(wanderer, (void **)&cycles);
		recal_init_info(*cycles, data);
		bwander_unlock_user_data(wanderer);

		//Set local data
		bwander_local_user_data_set(wanderer, data);
	}

	//Initialize get data environment
	collect_env = (recal_data_collect_env_t *) malloc(sizeof(recal_data_collect_env_t));
	recal_get_data_init_env(data->num_cycles, collect_env);

	//Obtain data from all reads in region
	for(i = 0; i < region->size; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Get data
		omp_set_lock(&wanderer->reference_lock);
		recal_get_data_from_bam_alignment(read, wanderer->reference, data, collect_env);
		omp_unset_lock(&wanderer->reference_lock);
	}

	//Destroy environment
	recal_get_data_destroy_env(collect_env);

	return NO_ERROR;
}

void
reduce_data(void *data, void *dest)
{
	recal_info_t *data_ptr, *dest_ptr;
	assert(data);
	assert(dest);

	//Cast pointers
	data_ptr = (recal_info_t *)data;
	dest_ptr = (recal_info_t *)dest;

	//Combine
	recal_reduce_info(dest_ptr, data_ptr);
}

void
destroy_data(void *data)
{
	recal_info_t *aux;

	assert(data);

	aux = (recal_info_t *)data;
	recal_destroy_info(aux);
}

static inline ERROR_CODE
wander_bam_file_recalibrate(bam_wanderer_t *wanderer, bam_file_t *in_bam, bam_file_t *out_bam, genome_t *ref, const char *in_data_f, const char *out_data_f, const char *out_info_f, int cycles)
{
	//Data
	recal_info_t *data;
	U_CYCLES aux_cycles;

	assert(wanderer);
	assert(in_bam);
	assert(ref);
	assert(cycles != 0);

	//Init wandering
	bwander_init(wanderer);

	//Configure wanderer
	bwander_configure(wanderer, in_bam, out_bam, ref,
			(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
			(int (*)(void *, bam_region_t *))recalibrate_collect_processor);


	printf("Cycles: %d\n",cycles);

	//Load data file
	if(in_data_f)
	{
		recal_load_recal_info(in_data_f, data);
	}
	else
	{
		//Create new data
		aux_cycles = cycles;
		data = (recal_info_t *)malloc(sizeof(recal_info_t));
		recal_init_info(aux_cycles, data);
	}

	//Set user data
	int cycles_param = cycles;
	bwander_set_user_data(wanderer, &cycles_param);

	//Run wander
	bwander_run(wanderer);

	//Reduce data
	bwander_local_user_data_reduce(wanderer, data, reduce_data);
	bwander_local_user_data_free(wanderer, destroy_data);

	//Delta processing
	recal_calc_deltas(data);

	//Save data file
	if(out_data_f)
	{
		recal_save_recal_info(data, out_data_f);
	}

	//Save info file
	if(out_info_f)
	{
		recal_fprint_info(data, out_info_f);
	}

	//Free data memory
	recal_destroy_info(data);
	free(data);

	//Destroy wanderer
	bwander_destroy(wanderer);

	return NO_ERROR;
}

ERROR_CODE
wander_bam_file_(uint8_t flags, char *bam_path, char *ref_name, char *ref_path, char *data_file, char *info_file, char *outbam, int cycles)
{
	//Files
	char *in_data_f = NULL;
	char *out_data_f = NULL;
	bam_file_t *bam_f = NULL;
	bam_file_t *out_bam_f = NULL;
	genome_t* ref = NULL;
	int bytes;

	//Data
	recal_info_t data;

	//Times
	double times;

	//Wanderer
	bam_wanderer_t wanderer;

	assert(bam_path);
	assert(ref_name);
	assert(ref_path);

#ifdef D_TIME_DEBUG
	times = omp_get_wtime();
#endif

	//Get phases
	if((flags & RECALIBRATE_RECALIBRATE) && !(flags & RECALIBRATE_COLLECT))
	{
		//Second phase only, needs data file
		assert(data_file);
	}

	//Open bam
	{
		printf("Opening BAM from \"%s\" ...\n", bam_path);
		bam_f = bam_fopen(bam_path);
		assert(bam_f);
		printf("BAM opened!...\n");
	}

	//Open reference genome, if collect phase is present
	if(flags & RECALIBRATE_COLLECT)
	{
		printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
		ref = genome_new(ref_name, ref_path);
		assert(ref);
		printf("Reference opened!...\n");
	}

	//Create new bam
	if(outbam != NULL && (flags & RECALIBRATE_RECALIBRATE))
	{
		printf("Creating new bam file in \"%s\"...\n", outbam);
		//init_empty_bam_header(orig_bam_f->bam_header_p->n_targets, recal_bam_header);
		out_bam_f = bam_fopen_mode(outbam, bam_f->bam_header_p, "w");
		assert(out_bam_f);
		bam_fwrite_header(out_bam_f->bam_header_p, out_bam_f);
		out_bam_f->bam_header_p = NULL;
		printf("New BAM initialized!...\n");
	}

#ifdef D_TIME_DEBUG
	times = omp_get_wtime() - times;
	time_add_time_slot(D_FWORK_INIT, TIME_GLOBAL_STATS, times);
#endif

	//Setup data file pointers
	if(data_file != NULL)
	{
		//If recalibrate data file
		if((flags & RECALIBRATE_RECALIBRATE) && (flags & RECALIBRATE_COLLECT))
		{
			//Full recalibration
			in_data_f = data_file;
			out_data_f = data_file;
		}
		else
		{
			if(flags & RECALIBRATE_COLLECT)
			{
				//Only collect
				out_data_f = data_file;
			}
			else
			{
				//Only recalibrate
				in_data_f = data_file;
			}
		}
	}

	//Recalibrate wandering
	wander_bam_file_recalibrate(&wanderer, bam_f, out_bam_f, ref, in_data_f, out_data_f, info_file, cycles);

	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	printf("BAM closed.\n");

	if(ref != NULL)
	{
		printf("\nClosing reference file...\n");
		genome_free(ref);
		printf("Reference closed.\n");
	}

	if(outbam != NULL)
	{
		printf("Closing \"%s\" BAM file...\n", outbam);
		bam_fclose(out_bam_f);
		printf("BAM closed.\n");
	}

	return NO_ERROR;
}
