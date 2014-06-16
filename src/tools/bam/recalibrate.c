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
#include "bfwork/bfwork.h"
#include "recalibrate/recal_config.h"
#include "recalibrate/recal_structs.h"
#include "recalibrate/bam_recal_library.h"

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
	char *refc, *inputc, *outputc, *infofilec, *datafilec;
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
	LOG_FILE("recalibration.log","w");
	LOG_VERBOSE(1);
	LOG_LEVEL(LOG_WARN_LEVEL);

	//Set num of threads
	//omp_set_num_threads(threads);
	printf("Threading with %d threads and %s schedule\n", omp_get_max_threads(), sched);

	//Obtain reference filename and dirpath from full path
	refc = NULL;
	if(refcount > 0)
		refc = strdup(reference);

	//Print data
	infofilec = NULL;
	if(infocount > 0)
		infofilec = strdup(infofile);
	
	//Save data file
	datafilec = NULL;
	if(datacount > 0)
		datafilec = strdup(datafile);

	//Input BAM
	inputc = NULL;
	if(incount > 0)
		inputc = strdup(input);

	//Output BAM
	outputc = NULL;
	if(output)
		outputc = strdup(output);

	//Recalibrate
	if(p1 && p2)
	{
		printf("Full recalibration\n");
		recal_bam_file(RECALIBRATE_COLLECT | RECALIBRATE_RECALIBRATE, inputc, refc, datafilec, infofilec, outputc, cycles);
	}
	else
	{
		if(p1)
		{
			printf("Phase 1 recalibration\n");
			recal_bam_file(RECALIBRATE_COLLECT, inputc, refc, datafilec, infofilec, outputc, cycles);
		}
		else
		{
			printf("Phase 2 recalibration\n");
			recal_bam_file(RECALIBRATE_RECALIBRATE, inputc, refc, datafilec, infofilec, outputc, cycles);
		}
	}
	if(inputc)
		free(inputc);
	if(outputc)
		free(outputc);
	if(refc)
		free(refc);
	if(datafilec)
		free(datafilec);
	if(infofilec)
		free(infofilec);
	
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
    outfile->filename[0]=NULL;
    
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
