import os

# Initialize the environment with path variables, CFLAGS, and so on
hpglib_path = '#lib/hpg-libs/c/'
third_party_path = '#/lib/hpg-libs/third_party'
third_party_hts_path = '#/lib/hpg-libs/third_party/htslib'
third_party_samtools_path = '#/lib/hpg-libs/third_party/samtools'

#system_include = '/usr/include'
#system_libs = '/usr/lib' 

#extrae_include = '/home/hmartinez/opt/extrae-2.5.1/include'
#extrae_libs    = '/home/hmartinez/opt/extrae-2.5.1/lib'

#other_libs = '/home/hmartinez/opt/lib/'
#other_include = '/home/hmartinez/opt/include/'

compiler = ARGUMENTS.get('compiler', 'gcc')
build_tools = [ 'default', 'packaging' ]

env = Environment(tools = build_tools,
      		  CC = compiler,
                  CFLAGS = '-Wall -std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -D_REENTRANT',
		  CPPPATH = ['#', '#src', '#include', "#src/tools/bam/", hpglib_path + 'src', third_party_path, third_party_samtools_path, third_party_hts_path, '/usr/include', '/usr/local/include', '/usr/include/libxml2', '/usr/lib/openmpi/include'],
		  LIBPATH = [hpglib_path + 'build', third_party_hts_path, third_party_samtools_path, '/usr/lib', '/usr/local/lib'],
                  LIBS = ['curl', 'dl', 'gsl', 'gslcblas', 'm', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])


if os.environ.has_key('C_INCLUDE_PATH'):
   for dir in os.getenv('C_INCLUDE_PATH').split(':'):
       env.Append(CPPPATH=[dir])
elif os.environ.has_key('C_PATH'):
   for dir in os.getenv('C_PATH').split(':'):
       env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
   for dir in os.getenv('LIBRARY_PATH').split(':'):
       env.Append(LIBPATH=[dir])
elif os.environ.has_key('LD_LIBRARY_PATH'):
   for dir in os.getenv('LD_LIBRARY_PATH').split(':'):
       env.Append(LIBPATH=[dir])

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

if int(ARGUMENTS.get('verbose', '0')) == 1:
    env['CFLAGS'] += ' -D_VERBOSE'

if int(ARGUMENTS.get('timing', '0')) == 1:
    env['CFLAGS'] += ' -D_TIMING'

env['objects'] = []

# Compile dependencies
SConscript(['lib/hpg-libs/SConstruct'])

#SConscript(['%s/SConscript' % bioinfo_path,
#            '%s/SConscript' % commons_path
#            ], exports = ['env', 'debug', 'compiler'])

envprogram = env.Clone()
envprogram['CFLAGS'] += ' -DNODEBUG -mssse3 -DD_TIME_DEBUG'

aligner = envprogram.Program('#bin/hpg-aligner',
             source = [Glob('src/*.c'),
		       Glob('src/tools/bam/aux/*.c'),
	     	       Glob('src/tools/bam/bfwork/*.c'),
		       Glob('src/tools/bam/recalibrate/*.c'),
	     	       Glob('src/tools/bam/aligner/*.c'),
		       Glob('src/build-index/*.c'),
		       Glob('src/dna/clasp_v1_1/*.c'),
		       Glob('src/dna/*.c'),
	               Glob('src/rna/*.c'),
	               Glob('src/bs/*.c'),
	               Glob('src/sa/*.c'),
		       "%s/sam_opts.o" % third_party_samtools_path,
		       "%s/bam_sort.o" % third_party_samtools_path,
		       "%s/bam_index.o" % third_party_samtools_path,
                      "%s/build/libhpg.a" % hpglib_path,
                      "%s/libbam.a" % third_party_samtools_path,
                      "%s/libhts.a" % third_party_hts_path
                      ]
           )

bam  = envprogram.Program('#bin/hpg-bam',
             source = [Glob('src/tools/bam/*.c'), 
	     	       Glob('src/tools/bam/aux/*.c'),
	     	       Glob('src/tools/bam/bfwork/*.c'),
	     	       Glob('src/tools/bam/recalibrate/*.c'),
	     	       Glob('src/tools/bam/aligner/*.c'),
		       "%s/sam_opts.o" % third_party_samtools_path,
		       "%s/bam_sort.o" % third_party_samtools_path,
		       "%s/bam_index.o" % third_party_samtools_path,
                       "%s/build/libhpg.a" % hpglib_path,
                       "%s/libbam.a" % third_party_samtools_path,
                       "%s/libhts.a" % third_party_hts_path
                      ]
           )

fastq = envprogram.Program('#bin/hpg-fastq',
             source = [Glob('src/tools/fastq/*.c'), 
                       "%s/build/libhpg.a" % hpglib_path
                      ]
           )

#Depends(aligner, bam, fastq)

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
