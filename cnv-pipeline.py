#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 18:16:57 2021
@author: MShahrokh
"""
import sys, getopt
import optparse
from cnv_sidefuncs import *
from cnv_plots import *
import datetime
# A wrapper to generate copy number index values

def main(argv):
    parser = optparse.OptionParser()
    parser.add_option('-b', '--targeted_bam_dir', action="store", dest="targetedBAM_dir",help="Directory of Targeted-only BAM files", default=None)
    parser.add_option('-B', '--all_bam_dir', action="store", dest="allBAM_dir",help="Directory of BAM files contraining off-target", default=None)
    parser.add_option('-s', '--selector_file', action="store", dest="selectorFile",help="Selector file", default=None)
    parser.add_option('-g', '--genome_fasta', action="store", dest="genomeFASTA",help="Genome FASTA file", default=None)
    parser.add_option('-G', '--genome_bed', action="store", dest="genomeBED",help="Genome BED file", default=None)
    parser.add_option('-m','--norm_cohort1_t', action="store", dest="normalPath1_t",help="Directory of Targeted-only BAM files [for normal cohort I]", default=None)
    parser.add_option('-M','--norm_cohort1_w', action="store", dest="normalPath1_w",help="Directory of BAM files contraining off-target [for normal cohort I]", default=None)
    parser.add_option('-n','--norm_cohort2_t', action="store", dest="normalPath2_t",help="Directory of Targeted-only BAM files [for normal cohort II]", default=None)
    parser.add_option('-N','--norm_cohort2_w', action="store", dest="normalPath2_w",help="Directory of BAM files contraining off-target [for normal cohort II]", default=None)
    parser.add_option('-p','--norm_cohort1_t_b', action="store", dest="normalPath1_t_b",help="Path to the targeted background (cohort I)", default=None)
    parser.add_option('-P','--norm_cohort1_w_b', action="store", dest="normalPath1_w_b",help="Path to the off-target background (cohort I)", default=None)
    parser.add_option('-q','--norm_cohort2_t_b', action="store", dest="normalPath2_t_b",help="Path to the targeted background (cohort II)", default=None)
    parser.add_option('-Q','--norm_cohort2_w_b', action="store", dest="normalPath2_w_b",help="Path to the off-target background (cohort II)", default=None)
    parser.add_option('-z', '--binsize', type="int", action="store", dest="binSize",help="Bin Size", default=1E5)
    parser.add_option('-L', '--maxLength', type="int", action="store", dest="maxLength",help="Maximum cfDNA fragment size", default=1000)
    parser.add_option('-i', '--sampleInfo', action="store", dest="sampleInfo",help="Tab delimited file- three columns (sampleName barcode_deduped samtools_deduped)", default=None)
    parser.add_option('-a', '--annotateBED', action="store", dest="annotBED",help="Annotation BED file", default=None)
    parser.add_option('-r', '--multiprocess',type="int", action="store", dest="multiprocess",help="Number of processes for multi-processing of coverage calculations", default=16)
    parser.add_option('-o', '--outdir_path', action="store", dest="outdir",help="Output directory for the final files", default=None)
    parser.add_option('-x', '--filter_bed', action="store", dest="filterBED",help="Regions to exclude (e.g., centromeres)", default=None)
    parser.add_option('-c', '--backgroundDir', action="store", dest="bgdir",help="Directory to four background files", default=None)
    parser.add_option('-v', '--calculateCovLogical_t', type="int", action="store", dest="calccov_t",help="Set to 0 if coverage files of target samples already EXIST", default=1)
    parser.add_option('-u', '--calculateCovLogical_b', type="int", action="store", dest="calccov_b",help="Set to 0 if coverage files of background samples already EXIST", default=1)
    parser.add_option('-f', '--lowresolution', type="int", action="store", dest="lowres",help="Set to 1 if selector is bigger than 1M (e.g., exome data)", default=0)
    old_stdout = sys.stdout
    options, args = parser.parse_args()
    print(options.calccov_t)
    if options.outdir==None:
        options.outdir = options.targetedBAM_dir   
    else:
        os.makedirs(options.outdir,exist_ok=True)
    if options.lowres == 1:
        options.lowres = True
    else:
        options.lowres = False
    intermediateDir = options.outdir+"/intermediatefiles/"
    os.makedirs(intermediateDir,exist_ok=True)
    log_file = open(intermediateDir+"message.log","w")
    sys.stdout = log_file
    ##
    if options.bgdir!=None:
        options.normalPath1_t_b = options.bgdir+"/background-ontarget-cohort1.txt"
        options.normalPath1_w_b = options.bgdir+"/background-offtarget-cohort1.txt"
        options.normalPath2_t_b = options.bgdir+"/background-ontarget-cohort2.txt"
        options.normalPath2_w_b = options.bgdir+"/background-offtarget-cohort2.txt"
    if options.normalPath1_t_b==None and options.normalPath1_t==None:
        raise ValueError("There must be either a path for the initial normalization background or directory of bam files to build one from")
    if options.normalPath1_w_b==None and options.normalPath1_w==None:
        raise ValueError('"There must be either a path for the initial normalizatio background or directory of bam files to build one from"')

    if options.normalPath2_t_b==None and options.normalPath2_t==None:
        raise ValueError("There must be either a path for the final Z-score background or directory of bam files to build one from")
    if options.normalPath2_w_b==None and options.normalPath2_w==None:
        raise ValueError('"There must be either a path for the final Z-score background or directory of bam files to build one from"')
    
    ## Off-target coverage calculations
    binbedname = intermediateDir + os.path.basename(options.genomeBED) + ".binned."+str(options.binSize)+".bed"
    find_antitarget(genome_size=options.genomeBED,savename=binbedname, binsize = options.binSize)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tGenome is binned and stored in this file: "+binbedname)
    binBED_GC = add_gc_col(binbedname,genome = options.genomeFASTA, extend = None)
    paddedSel_file = padd_selector(bedname=options.selectorFile, paddsize=500, outdir=intermediateDir)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: coverage of off-target fragments for target samples")
    if options.calccov_t==1:
        calc_coverage_antitarget_parallel(bamdir=options.allBAM_dir,binbed=binBED_GC,bedname = paddedSel_file, maxLength=options.maxLength, multiprocess=options.multiprocess)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of off-target fragments for target samples")

    # On-target coverage calculations [returning the file names of GC-corrected depth]
    bedname_new = prepare_bed_base(bedname = options.selectorFile, genome = options.genomeFASTA, outdir=intermediateDir, lowres = options.lowres, binsz = 60)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: coverage of on-target fragments for target samples")
    if options.calccov_t==1:
        calc_coverage_ontarget_parallel(bamdir=options.targetedBAM_dir,bedname = bedname_new,bedfile = options.selectorFile,genome_size = options.genomeBED, lowres = options.lowres, multiprocess=options.multiprocess) 
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of on-target fragments for target samples")
    ## CREATE BACKGROUND *IF NOT PROVIDED*
    if options.normalPath1_t_b==None:
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: Coverage of on-target fragments for control cohort I samples")
        if options.calccov_b==1:
            calc_coverage_ontarget_parallel(bamdir=options.normalPath1_t,bedname = bedname_new,bedfile = options.selectorFile, genome_size = options.genomeBED, lowres = options.lowres, multiprocess=options.multiprocess) 
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of on-target fragments for control cohort I samples")
        calc_bg_stats(options.normalPath1_t,savename = options.normalPath1_t+"/background-ontarget-cohort1.txt",countzero=True)
        options.normalPath1_t_b = options.normalPath1_t+"/background-ontarget-cohort1.txt"
        
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: CNV index and annotations for on-target (control cohort I)")
    if options.normalPath1_w_b==None:
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: Coverage of off-target fragments for control cohort I samples")
        if options.calccov_b==1:
            calc_coverage_antitarget_parallel(bamdir=options.normalPath1_w,binbed=binBED_GC,bedname = paddedSel_file, maxLength=options.maxLength, multiprocess=options.multiprocess) 
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of off-target fragments for control cohort I samples")
        calc_bg_stats(options.normalPath1_w,savename = options.normalPath1_w+"/background-offtarget-cohort1.txt",countzero=True)
        options.normalPath1_w_b = options.normalPath1_w+"/background-offtarget-cohort1.txt"
        
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: CNV index and annotations for off-target (control cohort I)")
    ## First generate the off-target bin background
    cov_preparation(options.targetedBAM_dir,options.normalPath1_t_b,numrow=5000)
    cov_preparation(options.allBAM_dir,options.normalPath1_w_b,numrow=100)
    ## In the following, we annotate the CNV info using the user provided genomic-coordinates of a set of pre-defined regions
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: CNV index and annotations for on-target and off-target (target samples)")
    annotate_index(options.targetedBAM_dir,options.annotBED,options.filterBED)
    annotate_index(options.allBAM_dir,options.annotBED,options.filterBED)
     
    ##### This is where we incorporate the information from the 2nd control cohort ****
    if options.normalPath2_t_b==None:
        if options.normalPath2_t!=options.normalPath1_t:
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: Coverage of on-target fragments for control cohort II samples")
            if options.calccov_b==1:
                calc_coverage_ontarget_parallel(bamdir=options.normalPath2_t,bedname = bedname_new,bedfile = options.selectorFile, genome_size = options.genomeBED,lowres = options.lowres, multiprocess=options.multiprocess) 
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of on-target fragments for control cohort II samples")
        cov_preparation(options.normalPath2_t,options.normalPath1_t_b,numrow=5000)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: CNV index and annotations for on-target (control cohort II)")
        annotate_index(options.normalPath2_t,options.annotBED,options.filterBED)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: CNV index and annotations for on-target (control cohort II)")
        calc_bg_stats(covdir=options.normalPath2_t,savename = options.normalPath2_t+"/background-ontarget-cohort2.txt", str2use = "gc.corrected.norm.log.std.index",str2srch = ".GCcorrected.index.annotated")
        options.normalPath2_t_b = options.normalPath2_t+"/background-ontarget-cohort2.txt"

    if options.normalPath2_w_b==None:
        if options.normalPath2_w!=options.normalPath1_w:
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: Coverage of off-target fragments for control cohort II samples")
            if options.calccov_b==1:
                calc_coverage_antitarget_parallel(bamdir=options.normalPath2_w,binbed=binBED_GC,bedname = paddedSel_file, maxLength=options.maxLength, multiprocess=options.multiprocess) 
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: Coverage of off-target fragments for control cohort II samples")
        cov_preparation(options.normalPath2_w,options.normalPath1_w_b,numrow=100)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: CNV index and annotations for off-target (control cohort II)")
        annotate_index(options.normalPath2_w,options.annotBED,options.filterBED)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: CNV index and annotations for off-target (control cohort II)")
        calc_bg_stats(covdir=options.normalPath2_w,savename = options.normalPath2_w+"/background-offtarget-cohort2.txt", str2use = "gc.corrected.norm.log.std.index",str2srch = ".GCcorrected.index.annotated")
        options.normalPath2_w_b = options.normalPath2_w+"/background-offtarget-cohort2.txt"
    ## In the following, we take the last step by computing the stats for final index transformation using the 2nd cohort of controls
    ## Final Z-score generation!
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tCalculating: summarizing by pre-defined regions for CNV index and Z-scores...")
    index2z(options.targetedBAM_dir,options.allBAM_dir,options.normalPath2_t_b,options.normalPath2_w_b,options.annotBED,sampleinfo=options.sampleInfo,outdir=options.outdir)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tDone: summarizing by pre-defined regions for CNV index and Z-scores...")
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\tGenerating ideograms...")
    ideogram(options.outdir,options.genomeBED)
    sys.stdout = old_stdout
    log_file.close()
if __name__ == "__main__":
   main(sys.argv[1:])    
   
   
   
   