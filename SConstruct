
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/hpg-libs/bioinfo-libs'
commons_path = '#lib/hpg-libs/common-libs'
#math_path = '#libs/math'
bam_path = "#src/tools/bam"
aux_path = "#src/tools/bam/aux"
recal_path = "#src/tools/bam/recalibrate"
alig_path = "#src/tools/bam/aligner"

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp',
                  CPPPATH = ['#', '#src', bam_path, bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path ],
                  LIBPATH = [commons_path, bioinfo_path],
                  LIBS = ['m', 'z', 'curl'],
                  LINKFLAGS = ['-fopenmp'])

if int(ARGUMENTS.get('debug', '1')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3'

env['objects'] = []

# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

env.Program('#bin/hpg-bam',
             source = [Glob('src/tools/bam/*.c'), 
	     	       Glob('%s/*.c' % aux_path),
	     	       Glob('%s/*.c' % recal_path),
	     	       Glob('%s/*.c' % alig_path),
                       "%s/libbioinfo.a" % bioinfo_path,
                       "%s/libcommon.a" % commons_path
                      ]
           )

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
