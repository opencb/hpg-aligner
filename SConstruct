
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/hpg-libs/bioinfo-libs'
commons_path = '#lib/hpg-libs/common-libs'
#math_path = '#libs/math'

system_include = '/usr/include'
system_libs = '/usr/lib' 

extrae_include = '/home/hmartinez/opt/extrae-2.5.1/include'
extrae_libs    = '/home/hmartinez/opt/extrae-2.5.1/lib'

other_libs = '/home/hmartinez/opt/lib/'
other_include = '/home/hmartinez/opt/include/'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-Wall -std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -D_REENTRANT',
                  CPPPATH = ['#', '#src', '#src/tools/bam', bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path, system_include, '%s/libxml2' % system_include ],
                  LIBPATH = [commons_path, bioinfo_path, system_libs],
                  LIBS = ['xml2', 'm', 'z', 'curl', 'dl', 'bioinfo', 'common'],
                  LINKFLAGS = ['-fopenmp'])

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3'

if int(ARGUMENTS.get('verbose', '0')) == 1:
    env['CFLAGS'] += ' -D_VERBOSE'

if int(ARGUMENTS.get('timing', '0')) == 1:
    env['CFLAGS'] += ' -D_TIMING'

env['objects'] = []

# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

envprogram = env.Clone()
envprogram['CFLAGS'] += ' -DNODEBUG -mssse3 -DD_TIME_DEBUG'

bams = envprogram.Program('#bin/hpg-bam',
             source = [Glob('src/tools/bam/*.c'), 
	     	       Glob('src/tools/bam/aux/*.c'),
	     	       Glob('src/tools/bam/bfwork/*.c'),
	     	       Glob('src/tools/bam/recalibrate/*.c'),
	     	       Glob('src/tools/bam/aligner/*.c'),
                       "%s/libbioinfo.a" % bioinfo_path,
                       "%s/libcommon.a" % commons_path
                      ]
           )

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
		       "%s/libcommon.a" % commons_path,
		       "%s/libbioinfo.a" % bioinfo_path
                      ]
           )
Depends(aligner, bams)

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
