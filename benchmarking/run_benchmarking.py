#!/usr/bin/env python

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
    print "Running BWA MEM: " + cmd
    os.system(cmd)
    #time ~/soft/bwa/bwa mem -t 8 ~/data/GRCh38/bwa-index/Homo_sapiens.GRCh38.dna.fa ~/data/GRCh38/sim-datasets/pair1.1M.125nt.fa  ~/data/GRCh38/sim-datasets/pair2.1M.125nt.fa  > out/out.1M.125nt.sam
    return;

#------------------------------------------------
# run HPG Aligner
#------------------------------------------------

def run_hpg_aligner(outdir, fq1, fq2, sam):
    print "\trunning HPG Aligner"
    out = outdir + "/hpg-aligner"
    os.system("mkdir " + out)
    cmd = args.hpgdir + "/bin/hpg-aligner dna --report-best"
    cmd += " -t " + args.cores
    cmd += " -i " + args.hpgindex + " -f " + fq1 + " -j " + fq2 + " -o " + out
    print "Running HPG Aligner: " + cmd
    os.system(cmd)
    os.system("cp " + out + "/alignments.sam " + outdir + "/" + sam)
    return;

#------------------------------------------------
# run programs for this commit id
#------------------------------------------------

def run_tests(id):
    if id is None:
        outdir = args.outdir
    else:
        outdir = args.outdir + "/" + args.branch + "/" + id

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print "output results will be stored at folder: " + outdir

    # run simultated dataset #1
    print "running tests for simulated dataset #1:"
    #run_bwa_mem(outdir, args.fqdir + "/sim1.pair1.fq", args.fqdir + "/sim1.pair2.fq", "bwamem.sim1.sam")
    run_hpg_aligner(outdir, args.fqdir + "/sim1.pair1.fq", args.fqdir + "/sim1.pair2.fq", "hpgaligner.sim1.sam")
    plot_results(outdir + "/bwamem.sim1.sam", outdir + "/hpgaligner.sim1.sam", outdir + "/sim1.pdf", "Simulated 1")

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
# get_commit_id
#------------------------------------------------

def get_commit_id(name):
    res = None
    last_id_filename = args.outdir + "/" + name + "/last.id"
    curr_id_filename = args.outdir + "/" + name + "/curr.id"

    # get the last commit ID
    os.system('cd ' + args.hpgdir + '; git log | head -1 | cut -d " " -f2 > ' + last_id_filename)

    # compare the last ID to the current ID
    cmp = filecmp.cmp(curr_id_filename, last_id_filename) 
    if cmp:
        res = None
    else:
        # save the last ID and run the benchmarking
        os.system('cp ' + last_id_filename + ' ' + curr_id_filename)
        res = (open(last_id_filename)).readline().split('\n', 1)[0]
    
    return res;

#------------------------------------------------
# main function
#------------------------------------------------

id = None;

if args.branch:
    id = get_commit_i(args.branch)
    if (id is None):
        print "no changes for the branch '" + args.branch + "' (same commit IDs): nothing to do"
    else:
        print "new commit for branch '" + args.branch + "' (different commit IDs): save the last commit ID and run benchmarking"
        run_tests(id)
else:
    print "run benchmarking directly (without branch specified)"
    run_tests(id)

print "done!"
