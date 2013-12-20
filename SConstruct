
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/hpg-libs/bioinfo-libs'
commons_path = '#lib/hpg-libs/common-libs'
#math_path = '#libs/math'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp',
                  CPPPATH = ['#', '#src', '#src/tools/bam/recalibrator/include', bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path ],
                  LIBPATH = [commons_path, bioinfo_path],
                  LIBS = ['m', 'z', 'curl'],
                  LINKFLAGS = ['-fopenmp'])

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

# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

env.Program('#bin/hpg-aligner',
             source = [Glob('src/*.c'),
		       Glob('src/build-index/*.c'),
		       Glob('src/dna/*.c'),
	               Glob('src/rna/*.c'),
	               Glob('src/bs/*.c'),
	               Glob('src/sa/*.c'),
		       "%s/libcommon.a" % commons_path,
		       "%s/libbioinfo.a" % bioinfo_path
                      ]
           )

env.Program('#bin/hpg-bam',
             source = [Glob('src/tools/bam/*.c'),
	     	       Glob('src/tools/bam/recalibrator/*.c'),
	     	       Glob('src/tools/bam/recalibrator/aux_library/*.c'),
	     	       Glob('src/tools/bam/recalibrator/bam_recal/*.c'),
	     	       Glob('src/tools/bam/recalibrator/timestats/*.c'),
                       "%s/libbioinfo.a" % bioinfo_path,
                       "%s/libcommon.a" % commons_path
                      ]
           )

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
