import argparse
import os
import filecmp
 
#------------------------------------------------
# parameters handling
#------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--branch", help="branch name")
parser.add_argument("--outdir", help="output directory to store results", required=True)
parser.add_argument("--fqdir", help="input directory containing the FastQ files", required=True)
parser.add_argument("--bindir", help="input directory containing executable files or scripts", required=True)
parser.add_argument("--bwadir", help="input directory where the BWA MEM executable is located", required=True)
parser.add_argument("--bwaindex", help="BWA MEM index", required=True)
parser.add_argument("--hpgdir", help="input directory where the HPG Aligner executable is located", required=True)
parser.add_argument("--hpgindex", help="HPG Aligner index directory", required=True)
parser.add_argument("--cores", help="number of cores", required=True)
args = parser.parse_args()

if args.verbosity:
    print "verbosity turned on"

#------------------------------------------------
# plot results function using R script
#------------------------------------------------

def plot_results(bwamem_sam, hpg_sam, image, title):
    print "\tplotting results"
    cmd = args.bindir + "/wgsim_eval.pl alneval -ag20 " + bwamem_sam + " > " + bwamem_sam + ".eval"
    os.system(cmd)

    cmd = args.bindir + "/wgsim_eval.pl alneval -ag20 " + hpg_sam + " > " + hpg_sam + ".eval"
    os.system(cmd)

    cmd = "Rscript " + args.bindir + "/r-plot.mapq.R " + bwamem_sam + ".eval " + hpg_sam + ".eval " + image + " \"" + title + "\""
    os.system(cmd)
    return;

#------------------------------------------------
# run BWA MEM
#------------------------------------------------

def run_bwa_mem(outdir, fq1, fq2, sam):
    print "\trunning BWA MEM"
    cmd = args.bwadir + "/bwa mem"
    cmd += " -t " + args.cores
    cmd += " " + args.bwaindex + " " + fq1 + " " + fq2 + " > " + outdir + "/" + sam
    os.system(cmd)
    #time ~/soft/bwa/bwa mem -t 8 ~/data/GRCh38/bwa-index/Homo_sapiens.GRCh38.dna.fa ~/data/GRCh38/sim-datasets/pair1.1M.125nt.fa  ~/data/GRCh38/sim-datasets/pair2.1M.125nt.fa  > out/out.1M.125nt.sam
    return;

#------------------------------------------------
# run HPG Aligner
#------------------------------------------------

def run_hpg_aligner():
    print "\trunning HPG Aligner"
    return;

#------------------------------------------------
# run programs for this commit id
#------------------------------------------------

def run_tests(id):
    outdir = args.outdir + "/" + id
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # run simultated dataset #1
    print "running tests for simulated dataset #1:"
    run_bwa_mem(outdir, args.fqdir + "/sim1.pair1.fq", args.fqdir + "/sim1.pair2.fq", "bwamem.sim1.sam")
    run_hpg_aligner()
    plot_results(outdir + "/bwamem.sim1.sam", outdir + "/bwamem.sim1.sam", outdir + "/sim1.pdf", "Simulated 1")

    # run real dataset #1
#    print "running tests for real dataset #1:"
#    run_bwa_mem()
#    run_hpg_aligner()
#    plot_results()

    # run real dataset #2
#    print "running tests for real dataset #2:"
#    run_bwa_mem()
#    run_hpg_aligner()
#    plot_results()
    
    return;

#------------------------------------------------
# main function
#------------------------------------------------

last_id_filename = "/tmp/last.id"
curr_id_filename = "/tmp/curr.id"

# get the last commit ID
os.system('cd ~/appl/hpg-aligner; git log | head -1 | cut -d " " -f2 > ' + last_id_filename)

# compare the last ID to the current ID
cmp = filecmp.cmp(curr_id_filename, last_id_filename) 
if cmp:
    print "same commit IDs, nothing to do"
else:
    # save the last ID and run the benchmarking
    os.system('cp ' + last_id_filename + ' ' + curr_id_filename)
    id = (open(last_id_filename)).readline().split('\n', 1)[0]
    print "different commit IDs, save the last ID and run benchmarking in folder " + id
    run_tests(id)
    
print "done!"
