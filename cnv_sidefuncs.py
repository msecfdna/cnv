#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 13:31:42 2021

@author: MShahrokh
"""
import sys
import numpy as np
import pandas as pd
import re # For string manipulations
import os # this is to run shell commands via the python environment
import time
import statsmodels.api as sm # This is for smoothing
from scipy.interpolate import interp1d
from os.path import exists
import random
import subprocess
from multiprocessing import Pool
from functools import partial
from pandas import DataFrame
from scipy.stats import uniform
from scipy.stats import randint
import matplotlib.pyplot as plt
from cnv_plots import *
import scipy.signal as signal
# we first open a bed file for the selector (on-target)
# a set of fixed regions for the off-target-

def sample_cov_base(bamname,bedname,bedfile, genome_size,mapq=20):
    outname = re.sub('\.bam$', '.base.coverage', bamname, count=1)
    command_str = "samtools view " + bamname + " -q "+str(mapq) 
    command_str = command_str + " | awk '{if ($2==99 || $2==163) print $3\"\t\"$4-1\"\t\"$4+$9}' | bedtools intersect -a stdin -b  " + bedfile + "   "
    command_str = command_str + "| bedtools genomecov -i stdin -g " 
    command_str = command_str+genome_size + " -dz | awk '{print $1\"\t\"$2\"\t\"$2+1\"\t\"$3}' > " + outname
    subprocess.run(command_str,shell=True)
    bedin_baselevel = pd.read_csv(bedname,sep="\t",names=['#chr','start','end','gc'],skiprows=1)
    covfile = pd.read_csv(outname,sep="\t",names=['#chr','start','end','f_depth'])
    covfile_merged = pd.merge(bedin_baselevel,covfile,on=['#chr','start','end'],how = 'left')
    covfile_merged.to_csv(outname,index=False,sep="\t")
    return(outname)
def prepare_bed_base(bedname, genome, outdir,extend = 2500, lowres = False, binsz=60): # This is for base level info.
    bedbasename = os.path.basename(bedname)
    if lowres:
        bedname_base = outdir+"/"+re.sub("\.bed$",".lowresbase.bed",bedbasename)
        if not os.path.exists(bedname_base):
            convert_bed_bin(bedname, binsz = binsz, save = bedname_base)
        outname = add_gc_col(bedname_base,genome=genome, extend = extend)
        
    else:
        bedname_base = outdir+"/"+re.sub("\.bed$",".base.bed",bedbasename)
        if not os.path.exists(bedname_base):
            convert_bed_base(bedname, save = bedname_base)
        outname = add_gc_col(bedname_base,genome=genome, extend = extend)
    return(outname)
def prepare_bed(bedname,genome,extend=2500): # This is for the bins alone! 
    add_gc_col(bedname,genome=genome,extend=extend)

def frag_cov(samrow,start,end):
    v = list(np.arange(float(samrow[3])-1,float(samrow[3])+float(samrow[8])+1))
    return(v)

#start = time.time();frag_cov_bam(bamname,bedname);end = time.time();print(end-start)
def padd_selector(bedname,paddsize,outdir):
    bedbasename = os.path.basename(bedname)
    paddedname = outdir+"/"+re.sub('\.bed$','.padded.'+str(paddsize)+".bed",bedbasename)
    commandstr = "cat "+bedname + " | awk '{print $1\"\t\"$2-500\"\t\"$3+500}' | bedtools merge -i stdin | bedtools sort -i stdin > "+paddedname
    subprocess.run(commandstr,shell=True)
    return(paddedname)
def sample_cov_targetbin(bamname,bedname,maxLength=1000,mapq = 20): # Coverage across bins [non-baselevel] 
    # *** This function expects 7 columns with a 4th column with name for the region to aggregate on! 
    dat = pd.read_csv(bedname,sep="\t",header="infer") 
    dat = dat.set_axis(['#chr', 'start', 'end','gc'], axis=1, inplace=False)
    chrs = ["chr"+str(i) for i in range(1,23)]
    chrs.extend(["chrX","chrY"])
    #beddat = pd.read_csv(bedname,sep="\t",header="infer",skiprows=1) open(bedname)
    outname = re.sub('\.bam$', '.base.coverage', bamname, count=1)
    outname_temp = re.sub('\.bam$', '.base.coverage.temp', bamname, count=1)
    cc = 0
    #lines =  f.readlines()[1:]
    for ch in chrs:
        cc+=1
        #fields = line.split()
        command_str = "samtools view -F3084 -q "+ str(mapq)+ " " + bamname + " "+ ch
        command_str = command_str + " | awk '{if (($2==99 || $2==163)"+ "&& ($9<=" + str(maxLength)+" && $9>=0)) print $3\"\t\"$4-1\"\t\"$4+$9}' " + "  > " + outname_temp
        # print("******************")
        # print(command_str)
        # print("******************")
        subprocess.run(command_str,shell=True,stdout=subprocess.DEVNULL)
        if os.stat(outname_temp).st_size != 0:
            if cc==1:
                command_str2 = "cat " + bedname + "| awk '{if ($1==\"" +ch+ "\") print $0}' " + "| bedtools coverage -b " + outname_temp + " -a stdin " + " > " + outname 
                subprocess.run(command_str2,shell=True,stdout=subprocess.DEVNULL)  
            else:
                command_str2 = "cat " + bedname + "| awk '{if ($1==\"" +ch+ "\") print $0}' "+ "| bedtools coverage -b " + outname_temp + " -a stdin " + " >> " + outname 
                subprocess.run(command_str2,shell=True,stdout=subprocess.DEVNULL)  
    os.remove(outname_temp)
    fnam = pd.read_csv(outname,sep="\t",header=None) 
    dat_ = dat
    dat_["cov"] = fnam.iloc[:,4]
    dat_ = dat_.set_axis(['#chr', 'start', 'end','gc','f_depth'], axis=1, inplace=False)
    dat_.to_csv(outname,index=False,sep="\t")
    return(outname)

def sample_cov_antitarget(bamname,binbed,targetExclude,maxLength=1000,mapq = 20): # Coverage across bins [non-baselevel] 
    # *** This function expects 7 columns with a 4th column with name for the region to aggregate on! 
    dat = pd.read_csv(binbed,sep="\t",header="infer") 
    dat = dat.set_axis(['#chr', 'start', 'end','gc'], axis=1, inplace=False)
    chrs = ["chr"+str(i) for i in range(1,23)]
    chrs.extend(["chrX","chrY"])
    #beddat = pd.read_csv(bedname,sep="\t",header="infer",skiprows=1) open(bedname)
    outname = re.sub('\.bam$', '.bin.coverage', bamname, count=1)
    outname_temp = re.sub('\.bam$', '.bin.coverage.temp', bamname, count=1)
    cc = 0
    #lines =  f.readlines()[1:]
    for ch in chrs:
        cc+=1
        #fields = line.split()
        command_str = "samtools view -F3084 -q "+ str(mapq)+ " " + bamname + " "+ ch
        command_str = command_str + " | awk '{if (($2==99 || $2==163)"+ "&& ($9<=" + str(maxLength)+" && $9>=0)) print $3\"\t\"$4-1\"\t\"$4+$9}' | bedtools subtract -a stdin -b " + targetExclude + " -A > " + outname_temp
        # print("******************")
        # print(command_str)
        # print("******************")
        subprocess.run(command_str,shell=True,stdout=subprocess.DEVNULL)
        if os.stat(outname_temp).st_size != 0:
            if cc==1:
                command_str2 = "cat " + binbed + "| awk '{if ($1==\"" +ch+ "\") print $0}' " + "| bedtools coverage -b " + outname_temp + " -a stdin " + " > " + outname 
                subprocess.run(command_str2,shell=True,stdout=subprocess.DEVNULL)  
            else:
                command_str2 = "cat " + binbed + "| awk '{if ($1==\"" +ch+ "\") print $0}' "+ "| bedtools coverage -b " + outname_temp + " -a stdin " + " >> " + outname 
                subprocess.run(command_str2,shell=True,stdout=subprocess.DEVNULL)  
    os.remove(outname_temp)
    fnam = pd.read_csv(outname,sep="\t",header=None) 
    dat_ = dat
    dat_["cov"] = fnam.iloc[:,4]
    dat_ = dat_.set_axis(['#chr', 'start', 'end','gc','f_depth'], axis=1, inplace=False)
    dat_.to_csv(outname,index=False,sep="\t")
    return(outname)

def add_gc_col(bedname,genome, extend = None, removeZero = True):
    outname = re.sub(".bed",".gc.bed",bedname)
    outnameTemp = re.sub(".bed",".gc.temp.bed",bedname)
    #print(os.path.exists(outname))
   # if not os.path.exists(outname):
    if extend!=None:
        dat = pd.read_csv(bedname,sep="\t",header=0) 
        dat2 = dat
        dat2.iloc[:,1] = dat.iloc[:,1]-extend
        dat2.iloc[:,2] = dat.iloc[:,2]+extend
        dat2["start_orig"] = dat2.iloc[:,1]+extend
        dat2["end_orig"] = dat2.iloc[:,2]-extend
        dat2.iloc[:,1] = [x if x>0 else 0 for x in dat2.iloc[:,1]]
        dat2.to_csv(outnameTemp,index=False,sep="\t",header=False)
        infile = outnameTemp
    else:
        infile = bedname
    command_str = "bedtools nuc -fi " + genome +" -bed " + infile + ">" + outname
    subprocess.run(command_str, shell = True)
    dat = pd.read_csv(outname,sep="\t") 
    
    #print(dat)
    if extend!=None:
        datfinal = dat[["#1_usercol","4_usercol","5_usercol","7_pct_gc"]]
    else:
        datfinal = dat[["#1_usercol","2_usercol","3_usercol","5_pct_gc"]]
    if removeZero:
        datfinal = datfinal[datfinal.iloc[:,3]!=0] # removing zero GC contents
    datfinal.to_csv(outname,index=False,sep="\t")
    return(outname)
def convert_bed_bin(bedin, binsz = 60, save = None): # Convert a tile-level selector to a base-level bed file
    dat = pd.read_csv(bedin,sep="\t",header=None) 
    nrow = dat.shape[0]
    mydf = pd.DataFrame()
    c = 0
    st_s = []
    ed_s = []
    ch_s = []
    for j in range(nrow):
        cur_row= dat.iloc[j,:]
        ch = cur_row[0]
        st = cur_row[1]
        ed = cur_row[2]
        numtiles = np.ceil((ed-st)/binsz)
        if (numtiles==1):
            st_s.extend([float(st)])
            ed_s.extend([float(st)])
            ch_s.extend([ch]*(int(1)))
        else:
            rmd = np.remainder(ed-st,binsz)
            ovlp = np.ceil((binsz -rmd )/(numtiles-1))
            if rmd==0:
                ovlp = 0
            stp = binsz - ovlp 
            newstarts = np.arange(st,st+numtiles*stp,int(stp))
            newends = newstarts+binsz
            newends[int(numtiles-1)] = ed
            st_s.extend((newstarts))
            ed_s.extend((newends))
            ch_s.extend([ch]*(int(numtiles)))
        
    mydf["#chr"] = ch_s
    mydf["start"] = [int(x) for x in st_s]
    mydf["end"] = [int(x) for x in ed_s]
    if save!=None:
        mydf.to_csv(save,index=False,sep="\t")
    else:
        return(mydf)
        
def convert_bed_base(bedin, save = None): # Convert a tile-level selector to a base-level bed file
    dat = pd.read_csv(bedin,sep="\t",header=None) 
    nrow = dat.shape[0]
    mydf = pd.DataFrame()
    c = 0
    st_s = []
    ed_s = []
    ch_s = []
    for j in range(nrow):
        cur_row= dat.iloc[j,:]
        ch = cur_row[0]
        st = cur_row[1]
        ed = cur_row[2]
        st_s.extend((range(st,ed)))
        ed_s.extend((range(st+1,ed+1)))
        ch_s.extend([ch]*(ed-st))
    mydf["#chr"] = ch_s
    mydf["start"] = st_s
    mydf["end"] = ed_s
    if save!=None:
        mydf.to_csv(save,index=False,sep="\t")
    else:
        return(mydf)
def files_in_dir(srchdir, str2srch = ".bam"):
    files = []
    for file in os.listdir(srchdir):
        if file.endswith(str2srch):
            files.append(os.path.join(srchdir, file))
    return(files)

def gc_correct(covfile,lowess_frac = 0.1, autosome = True, minFCov = 5):
    lowess = sm.nonparametric.lowess
    dat = pd.read_csv(covfile,sep="\t") 
    chrs = dat.iloc[:,0]
    y = dat.iloc[:,4]+0.1
    zerocov = np.where(y<=(minFCov+0.1)) # where are the regions with less than N fragments! 
    if autosome:
        y = y / np.nanmedian(y[np.invert(chrs.isin(["chrX","chrY"]))])
    else:
        y = y / np.nanmedian(y)
    x = dat.iloc[:,3]
    x_autosome = x[np.invert(chrs.isin(["chrX","chrY"]))]
    y_autosome =  y[np.invert(chrs.isin(["chrX","chrY"]))]
    zerocov_autosome = np.array(np.intersect1d(zerocov,np.where(np.invert(chrs.isin(["chrX","chrY"])))))
    ylog = np.log(y)
    for index in sorted(list(zerocov_autosome),reverse=True):
        del y_autosome[index]
    for index in sorted(list(zerocov_autosome),reverse=True):
        del x_autosome[index]
    ylog_autosome = np.log(y_autosome)
    order=np.argsort(x_autosome)
    ylog_ordered = ylog_autosome.iloc[order]
    x_ordered = x_autosome.iloc[order]
    ylog_rolled = ylog_ordered.rolling(15,min_periods=3).median()
    x_rolled = x_ordered.rolling(15,min_periods=3).median()
    z = lowess(ylog_rolled, x_rolled, frac= lowess_frac)
    f = np.interp(x, z[:,0], z[:,1])
    ynew = np.exp(ylog - f + np.nanmedian(ylog_autosome))
    ynew = [0 if x<0 else x for x in ynew]
    dat["gc.corrected"] = ynew
    dat.loc[zerocov]["gc.corrected"] = np.nan
    outname = covfile + ".GCcorrected"
    dat.to_csv(outname,index=False,sep="\t")
    return(outname)
    # bamfiles = files_in_dir(bamdir, str2srch = ".bam")
    
    # cov_files_gccorrected = []
    # for bf in bamfiles:
    #     baselevel_fname = sample_cov_base(bamname=bf,bedname=bedname,bedfile=bedfile,genome_size=genome_size)
    #     baselevel_gc_fname = gc_correct(baselevel_fname)
    #     cov_files_gccorrected.append(baselevel_gc_fname)
    # return(cov_files_gccorrected)

def calc_coverage_ontarget_single(bf,bedname,bedfile, genome_size, lowres = False):
    if lowres:
        baselevel_fname = sample_cov_targetbin(bamname=bf,bedname = bedname)
    else:
        baselevel_fname = sample_cov_base(bamname=bf,bedname=bedname,bedfile = bedfile, genome_size=genome_size)
    baselevel_gc_fname = gc_correct(baselevel_fname)
        
def calc_coverage_ontarget_parallel(bamdir,bedname,bedfile,genome_size, lowres = False, multiprocess = 16):
    bamfiles = files_in_dir(bamdir, str2srch = ".bam")
    pool = Pool(processes = multiprocess) 
    pool.map(partial(calc_coverage_ontarget_single,bedname=bedname,bedfile=bedfile,genome_size=genome_size, lowres = lowres), bamfiles) 
    pool.close()

def calc_coverage_antitarget(bamdir,binbed,bedname,maxLength=1000):
    # calculation of anti-target coverage
    bamfiles = files_in_dir(bamdir, str2srch = ".bam")
   # binbe = prepare_bed(bedname=binbed,genome=genome)
    cov_files_gccorrected = []
    for bf in bamfiles:
        rgnlevel_fname = sample_cov_antitarget(bf,binbed=binbed,targetExclude=bedname,maxLength=maxLength)
        #rgnlevel_fname = sample_cov_rgns(bamname=bf,bedname=bedname_new)
        rgnlevel_gc_fname = gc_correct(rgnlevel_fname)
        cov_files_gccorrected.append(rgnlevel_gc_fname)
    return(cov_files_gccorrected)
def calc_coverage_antitarget_single(bf,binbed,bedname,maxLength=1000):
        rgnlevel_fname = sample_cov_antitarget(bamname=bf,binbed=binbed,targetExclude=bedname,maxLength=maxLength)
        rgnlevel_gc_fname = gc_correct(rgnlevel_fname)

def calc_coverage_antitarget_parallel(bamdir,binbed,bedname,maxLength=1000, multiprocess = 16):
    # calculation of anti-target coverage
    bamfiles = files_in_dir(bamdir, str2srch = ".bam")
    pool = Pool(processes=multiprocess) 
    pool.map(partial(calc_coverage_antitarget_single,binbed=binbed,bedname=bedname,maxLength=maxLength), bamfiles)
    pool.close()
def sumzero(x):
    return(sum(x==0)/len(x))
def calc_bg_stats(covdir,savename = None, str2use = "gc.corrected",str2srch = ".GCcorrected",countzero = False):
    # background stats calculations
    covfiles = files_in_dir(covdir, str2srch = str2srch)
    allrows = pd.DataFrame()
    allrows_count = pd.DataFrame()
    bg_info = pd.DataFrame()
    fcounter = 0
    for covf in covfiles:
        fcounter+=1
        dat = pd.read_csv(covf,sep="\t",header=0)
        #print(dat)
        allrows[covf] = dat[str2use]
        if countzero:
            allrows_count[covf] = dat["f_depth"]
    dat=dat.rename(columns = {'chr':'#chr'})
    bg_info["#chr"] = dat["#chr"]
    bg_info["start"] = dat["start"]
    bg_info["end"] = dat["end"]
    bg_info["gc"] = dat["gc"]
    bg_info["f_depth"] = np.nan
    
    if countzero:
        numzeros = (allrows_count == 0).sum(axis=1)/allrows_count.shape[1]
        bg_info[str2use+".sumzero"] = numzeros
    bg_info[str2use+".mu"] = allrows.mean(axis=1,skipna=True)
    bg_info[str2use+".median"] = allrows.median(axis=1,skipna=True)
    bg_info[str2use+".std"] = allrows.std(axis=1,skipna=True)
    if savename!=None:
        bg_info.to_csv(savename,index=False,sep="\t")
    return(bg_info)
def find_antitarget(genome_size, savename, binsize = 1E5):
    commandstr0 = "bedtools makewindows -g " + str(genome_size) + " -w " + str(binsize)
    commandstr1 = commandstr0 + "| awk '{print $1\"\t\"$2\"\t\"$3}' "+" > " + savename
    os.system(commandstr1)
    
def sample_base_std(dat,numrow=5000,str2calc = "gc.corrected.norm.log"):
    # sample-specific standard-deviation calculations
    # numrow must be set to 100 for off-target regions
    # numrow must be set to 5000 for on-target regions
    iters = 1000 # Number of iterations to select random bases (or bins) to calcualte STD from!
    chrs = ["chr"+str(i) for i in range(1,23)]
    #dat = pd.read_csv(covfile,sep="\t",header="infer")
    stds_ = []
    for i in range(iters):
        curr_chr = random.choice(chrs)
        dat_select = dat[dat.iloc[:,0]==curr_chr]
        if dat_select.shape[0]<numrow:
            std_i = np.nanstd(dat_select[str2calc])
        else:
            curr_pos = random.choice(list(dat_select.index))
            dat_select2 = dat_select[abs(dat_select.index-curr_pos)<=numrow/2]
            std_i = np.nanstd(dat_select2[str2calc])
        stds_.append(std_i)
    return(np.nanmedian(stds_))

def nz2controls(covfile,bgdat,str2calc = "gc.corrected",str2norm = "gc.corrected.median"):
    # Normalization to median across controls (a column is added)
    #
    dat = pd.read_csv(covfile,sep="\t")
    dat[str2calc+".norm"] = dat[str2calc]/bgdat[str2norm]
    return(dat)
def log_vals(datin,str2calc = "gc.corrected.norm"): # Log-transforamtion (a column is added) 
    # This function expects a column of str2calc[default: gc.corrected.norm] in the coverage file]
    #dat = pd.read_csv(covfile,sep="\t")
    datin[str2calc+".log"] = np.log2(datin[str2calc]+1E-6)
    return(datin)

def cov_preparation(covdir,bgfile,numrow=5000):
    rmvzerodepth = 0.33
    highCoV = 1
    covfiles = files_in_dir(covdir, str2srch = ".GCcorrected")
    bgdat = pd.read_csv(bgfile,sep="\t",header="infer")
    for covfile in covfiles:
        # Normalization to the normal cohort I:
        dat = nz2controls(covfile,bgdat,str2calc = "gc.corrected",str2norm = "gc.corrected.median") 
        # Taking the log
        dat_logged = log_vals(dat,str2calc = "gc.corrected.norm")
        # calculating the sample level STD
        dat_logged.loc[(bgdat['gc.corrected.sumzero']>rmvzerodepth),"gc.corrected.norm.log"] = np.nan
        dat_logged.loc[(((bgdat['gc.corrected.std'])/abs(bgdat['gc.corrected.mu']))>highCoV),"gc.corrected.norm.log"] = np.nan
        sample_std = sample_base_std(dat_logged,numrow=numrow,str2calc = "gc.corrected.norm.log")
        # Dividing the normalized log-tranformed coverage by the STD
        dat_logged["gc.corrected.norm.log.std"] = dat_logged["gc.corrected.norm.log"]/sample_std
        dat_logged=dat_logged.rename(columns = {'chr':'#chr'})
        #dat_logged.loc[(bgdat['gc.corrected.sumzero']>rmvzerodepth),"gc.corrected.norm.log.std"] = 'nan'
        dat_logged.to_csv(covfile+".index",index=False,sep="\t")
def IQRV(df):
    return((df.quantile(0.75,axis=1)-df.quantile(0.25,axis=1))/df.quantile(0.5,axis=1))
def IQRV_row(x):
    return((np.quantile(x,0.75)-np.quantile(x,0.25))/abs(np.quantile(x,0.5)))
def annotate_index(covdir,annotBED,filterRgns=None):
    annot_dat = pd.read_csv(annotBED,sep="\t",names=['#chr','start','end','name'],header=0)
    annot_dat['name2'] = annot_dat['#chr']+"_" + annot_dat['name']
    covfiles = files_in_dir(covdir, str2srch = ".GCcorrected.index")
    for covfile in covfiles:
        covfileTemp = covfile+".annotTemp"
        covfileannotate = covfile+".annotated"
        command_str = "bedtools intersect -b " + covfile  + " -a " + annotBED +" -wao " + " | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$8\"\t\"$13}' " #+ covfileTemp
        if filterRgns!=None:
            command_str = command_str + "| bedtools intersect -a stdin -b "+filterRgns + " -c > " + covfileTemp
        else:
            command_str = command_str + " > " + covfileTemp
        subprocess.run(command_str,shell=True)
        dat = pd.read_csv(covfileTemp,sep="\t",header=None)
        # print("*")
        # print(dat)
        dat = dat.reset_index(drop=True)
        #dat=dat.rename(columns = {0:'#chr',1:"start",2:"end",3:"name"})
        dat['name2'] = dat.iloc[:,0] + "_" + dat.iloc[:,3]
        
        dat.loc[(dat[4]=="."),4] = np.nan
        dat[4] = dat[4].astype(float)
        dat.loc[(dat[5]=="."),5] = np.nan
        if filterRgns!=None:
            dat.loc[(dat[6]==1),5] = np.nan
        dat[5] = dat[5].astype(float)
        
        # print("***")
        # print(dat)
        dat_grouped = dat.groupby(by="name2").agg(np.nanmedian)
        dat_grouped2 = dat.groupby(by="name2").agg(IQRV_row)
        
        dat_grouped['#chr'] = [x.split('_')[0] for x in dat_grouped.index]
        dat_grouped['name'] = [x.split('_')[1] for x in dat_grouped.index]
        #print("*******")
        #print(dat_grouped)
        dat_grouped_ = pd.DataFrame({'#chr': list(dat_grouped['#chr']),
                                     'start' : list(dat_grouped.iloc[:,0]),
                                     'end' : list(dat_grouped.iloc[:,1]),
                                     'name': list(dat_grouped['name']), 
                                     'gc': list(dat_grouped[4]),
                                     'gc.corrected.norm.log.std.index':list(dat_grouped[5]),
                                     'IQR_median': list(dat_grouped2[4])},
                                    columns = ['#chr','start','end','name','gc','gc.corrected.norm.log.std.index','IQR_median'])
        dat_grouped_.reset_index()
        #os.remove(covfileTemp)
        dat_grouped_merged = pd.merge(annot_dat,dat_grouped_,on=['#chr','start','end','name'],how = 'left')
        dat_grouped_merged.to_csv(covfileannotate,sep="\t",index=False)

def index2z(covdirT,covdirW,covdirT_bg,covdirW_bg,annotBED,sampleinfo=None,outdir=None):
    rmvzerodepth = 0.33
    targeted_bg = pd.read_csv(covdirT_bg,sep="\t",names=['#chr','start','end','name','gc',
                                                       'gc.corrected.norm.log.std.index.mu',
                                                       'gc.corrected.norm.log.std.index.median',
                                                       'gc.corrected.norm.log.std.index.std'],header=0)
    targeted_bg = targeted_bg.reset_index()
    antitarget_bg = pd.read_csv(covdirW_bg,sep="\t",names=['#chr','start','end','name','gc',
                                                       'gc.corrected.norm.log.std.index.mu',
                                                       'gc.corrected.norm.log.std.index.median',
                                                       'gc.corrected.norm.log.std.index.std'],header=0)
    antitarget_bg = antitarget_bg.reset_index()
    annot_dat = pd.read_csv(annotBED,sep="\t",names=['#chr','start','end','name'],header=0)
    #annot_dat['name2'] = annot_dat['#chr']+"_"+annot_dat['name']
    covfilesT = files_in_dir(covdirT, str2srch = ".GCcorrected.index.annotated")
    covfilesW = files_in_dir(covdirW, str2srch = ".GCcorrected.index.annotated")
    for covfile in covfilesW:
        dat = pd.read_csv(covfile,sep="\t")
        #dat = dat.reset_index()
        dat["gc.corrected.norm.log.std.index"+".z"] = (dat["gc.corrected.norm.log.std.index"]-antitarget_bg['gc.corrected.norm.log.std.index.mu'])/antitarget_bg['gc.corrected.norm.log.std.index.std']
        #dat["gc.corrected.norm.log.std.index"+".z"][antitarget_bg['gc.corrected.norm.log.std.index.sumzero']>rmvzerodepth] = 'nan'
        dat["gc.corrected.norm.log.std.index"+".stdnorm"] = (dat["gc.corrected.norm.log.std.index"])/antitarget_bg['gc.corrected.norm.log.std.index.std']

        dat.to_csv(covfile+".zScore",sep="\t",index=False)
    for covfile in covfilesT:
        dat = pd.read_csv(covfile,sep="\t")
       # dat = dat.reset_index()
        dat["gc.corrected.norm.log.std.index"+".z"] = (dat["gc.corrected.norm.log.std.index"]-targeted_bg['gc.corrected.norm.log.std.index.mu'])/targeted_bg['gc.corrected.norm.log.std.index.std']
        dat["gc.corrected.norm.log.std.index"+".stdnorm"] = (dat["gc.corrected.norm.log.std.index"])/targeted_bg['gc.corrected.norm.log.std.index.std']

        dat.to_csv(covfile+".zScore",sep="\t",index=False)
    if sampleinfo!=None:
        sampleinfo_dat = pd.read_csv(sampleinfo,sep="\t",names=['sample','targeted','wg'],header=0)
        nsamples = sampleinfo_dat.shape[0]
        for index, sinf in sampleinfo_dat.iterrows():
            samplename = sinf[0]
            sample_t = re.sub('\.bam$', ".base.coverage.GCcorrected.index.annotated.zScore", sinf[1], count=1)
            sample_w = re.sub('\.bam$', ".bin.coverage.GCcorrected.index.annotated.zScore", sinf[2], count=1)
            dat_t = pd.read_csv(sample_t,sep="\t")
            dat_w = pd.read_csv(sample_w,sep="\t")
            w_t_raw = dat_t['IQR_median']
            w_w_raw = dat_w['IQR_median']
            w_t_final = 1 - w_t_raw/(w_t_raw+w_w_raw)
            w_w_final = 1-w_t_final
            dat_t['gc.corrected.norm.log.std.index.z.off'] = dat_w['gc.corrected.norm.log.std.index.z']
            
            dat_t['gc.corrected.norm.log.std.index.stdnorm.off'] = dat_w['gc.corrected.norm.log.std.index.stdnorm']
            dat_t['gc.off'] = dat_w['gc']
            dat_t['gc.corrected.norm.log.std.index.z.Final'] = (dat_t['gc.corrected.norm.log.std.index.z.off']+dat_t['gc.corrected.norm.log.std.index.z'])/np.sqrt(2)
            dat_t['gc.corrected.norm.log.std.index.zWeighted.Final'] = (w_w_final*dat_t['gc.corrected.norm.log.std.index.z.off']+w_t_final*dat_t['gc.corrected.norm.log.std.index.z'])/np.sqrt(w_t_final**2+w_w_final**2)
            dat_t['gc.corrected.norm.log.std.index.stdnorm.Final'] = (dat_t['gc.corrected.norm.log.std.index.stdnorm.off']+dat_t['gc.corrected.norm.log.std.index.stdnorm'])/np.sqrt(2)
            dat_t['gc.corrected.norm.log.std.index.stdnormWeighted.Final'] = (w_w_final*dat_t['gc.corrected.norm.log.std.index.stdnorm.off']+w_t_final*dat_t['gc.corrected.norm.log.std.index.stdnorm'])/np.sqrt(w_t_final**2+w_w_final**2)
            dat_t['Weight_on'] = w_t_final
            dat_t['Weight_off'] = w_w_final
            dat_t.loc[(dat_t['gc.corrected.norm.log.std.index.z.Final'].isnull()),'gc.corrected.norm.log.std.index.z.Final'] = dat_t['gc.corrected.norm.log.std.index.z.off']
            dat_t.loc[(dat_t['gc.corrected.norm.log.std.index.stdnorm.Final'].isnull()),'gc.corrected.norm.log.std.index.stdnorm.Final'] = dat_t['gc.corrected.norm.log.std.index.stdnorm.off']
            dat_t.loc[(dat_t['gc.corrected.norm.log.std.index.stdnormWeighted.Final'].isnull()),'gc.corrected.norm.log.std.index.stdnormWeighted.Final'] = dat_t['gc.corrected.norm.log.std.index.stdnorm.off']
            dat_t.loc[(dat_t['gc.corrected.norm.log.std.index.zWeighted.Final'].isnull()),'gc.corrected.norm.log.std.index.zWeighted.Final'] = dat_t['gc.corrected.norm.log.std.index.z.off']
            
            dat_t = dat_t.rename(columns = {'gc.corrected.norm.log.std.index.stdnorm':'gc.corrected.norm.log.std.index.stdnorm.on',
                                            'gc.corrected.norm.log.std.index.z':'gc.corrected.norm.log.std.index.z.on',
                                            'gc':'gc.on'})
            
            sub2save = ["#chr","start","end","name","gc.on","gc.off","gc.corrected.norm.log.std.index.stdnorm.on","gc.corrected.norm.log.std.index.z.on",
                        "gc.corrected.norm.log.std.index.stdnorm.off","gc.corrected.norm.log.std.index.z.off",
                        "gc.corrected.norm.log.std.index.stdnorm.Final","gc.corrected.norm.log.std.index.z.Final","Weight_on","Weight_off",
                        "gc.corrected.norm.log.std.index.stdnormWeighted.Final","gc.corrected.norm.log.std.index.zWeighted.Final"]
            dat2save = dat_t[sub2save]
            if outdir!=None:
                dat2save.to_csv(outdir+"/"+samplename+".cnvZscores",sep="\t",index=False)
            else:
                dat2save.to_csv(covdirT+"/"+samplename+".cnvZscores",sep="\t",index=False)
            
def ideogram(covdir,genomeBED):
    covfiles = files_in_dir(covdir, str2srch = ".cnvZscores")
    
    for covfile in covfiles:
        savename = re.sub('.cnvZscores', '', covfile, count=1)
        plot_ideogram_individual(covfile, genomeBED, "gc.corrected.norm.log.std.index.z.off", "gc.corrected.norm.log.std.index.z.on", 'CNV index (Z-score space)', savename)
        savename_combined = re.sub('.cnvZscores', '.combined', covfile, count=1)
        plot_ideogram_individual_singlescore(covfile=covfile, genomeBED=genomeBED,str2plot= "gc.corrected.norm.log.std.index.z.Final", ylabel = "Stouffer's Z-score", savename=savename_combined)
        savename2 = re.sub('.cnvWeightedZscores', '', covfile, count=1)
        plot_ideogram_individual(covfile=covfile, genomeBED=genomeBED, stroff = "gc.corrected.norm.log.std.index.stdnorm.off", stron = "gc.corrected.norm.log.std.index.stdnorm.on", ylabel = "Combined Index", savename=savename2)
        savename2_combined = re.sub('.cnvZscores','.weighted.combined', covfile, count=1)
        plot_ideogram_individual_singlescore(covfile=covfile, genomeBED=genomeBED, str2plot="gc.corrected.norm.log.std.index.zWeighted.Final", ylabel = "Weighted Stouffer's Z-score", savename=savename2_combined)
# def fdr(p_vals):
#     from scipy.stats import rankdata
#     ranked_p_values = rankdata(p_vals)
#     fdr = p_vals * len(p_vals) / ranked_p_values
#     fdr[fdr > 1] = 1
#     return fdr
def threshold(covdir,zCol = "gc.corrected.norm.log.std.index.zWeighted.Final", qThresh = 0.99):
    covfiles = files_in_dir(covdir, str2srch = ".cnvZscores")
    mydf = pd.DataFrame()
    for covfile in covfiles:
        dat = pd.read_csv(covfile,sep="\t")
        mydf[covfile] = dat[zCol]
    threshes_amp = mydf.quantile(qThresh,axis=1)
    threshes_del = mydf.quantile(1-qThresh,axis=1)
    threshdf = pd.DataFrame()
    threshdf["#chr"] = dat["#chr"]
    threshdf["start"] = dat["start"]
    threshdf["end"] = dat["end"]
    threshdf["name"] = dat["name"]
    threshdf["threshold_amp"] = list(threshes_amp)
    threshdf["threshold_del"] = list(threshes_del)
    threshdf.to_csv(covdir+"/event.thresholds.txt",sep="\t",index=False)
    return(covdir+"/event.thresholds.txt")
def callCNV_fdr(covdir,zCol = "gc.corrected.norm.log.std.index.zWeighted.Final"):
    import scipy
    import statsmodels.stats.multitest
    covfiles = files_in_dir(covdir, str2srch = ".cnvZscores")
    for covfile in covfiles:
        dat = pd.read_csv(covfile,sep="\t")
        zscores = dat[[zCol]]
        p_values = scipy.stats.norm.sf(abs(zscores))*2
        fdr_vals = [statsmodels.stats.multitest.fdrcorrection(x)[1][0] for x in p_values]
        zscores2 = [x for x in zscores.iloc[:,0]]
        for FDR in [0.01,0.05,0.1]:
            binarycalls = [np.sign(y) if x<=FDR else 0 for x,y in zip(fdr_vals,zscores2)]
            dat["CNVcalls"+"FDR"+str(FDR*100)] = binarycalls
            dat.loc[(dat[zCol].isnull()),"CNVcalls"+"FDR"+str(FDR*100)] = np.nan
        dat.to_csv(covfile+".CNVcallsFDR",sep="\t",index=False)
def callCNV_normal(covdir,thresholdFile, zCol = "gc.corrected.norm.log.std.index.zWeighted.Final"):
    import scipy
    import statsmodels.stats.multitest
    covfiles = files_in_dir(covdir, str2srch = ".cnvZscores")
    thresholds = pd.read_csv(thresholdFile,sep="\t")
    for covfile in covfiles:
        dat = pd.read_csv(covfile,sep="\t")
        dat_merged = pd.merge(thresholds,dat,on=['#chr','start','end','name'],how = 'left')
        ternaryCalls = [1 if x>y else -1 if x<=z else 0 for x,y,z in zip(dat_merged[zCol],dat_merged["threshold_amp"],dat_merged["threshold_del"])]
        dat_merged["CNVcalls_threshold"] = ternaryCalls
        dat_merged.loc[(dat_merged[zCol].isnull()),"CNVcalls_threshold"] = np.nan
        dat_merged.to_csv(covfile+".CNVcallsThreshold",sep="\t",index=False)            
def mergeFiles_fdr(covdir,FDR = 0.01):
    covfiles = files_in_dir(covdir, str2srch = ".CNVcallsFDR")
    allcnvcalls = pd.DataFrame()
    for covfile in covfiles:
        dat = pd.read_csv(covfile,sep="\t")
        samname = re.sub(".cnvZscores.CNVcallsFDR","",os.path.basename(covfile),count=1)
        dat_amp = dat.loc[dat["CNVcallsFDR"+str(FDR*100)]==1]
        dat_amp_sel = dat_amp[["#chr","name"]]
        dat_amp_sel["sample"] = samname
        dat_amp_sel["Z"] = dat_amp["gc.corrected.norm.log.std.index.zWeighted.Final"]
        dat_amp_sel["event"] = "amplification"
        dat_del = dat.loc[dat["CNVcallsFDR"+str(FDR*100)]==-1]
        dat_del_sel = dat_del[["#chr","name"]]
        dat_del_sel["sample"] = samname
        dat_del_sel["Z"] = dat_del["gc.corrected.norm.log.std.index.zWeighted.Final"]
        dat_del_sel["event"] = "deletion"
        allcnvcalls = allcnvcalls.append(dat_amp_sel)
        allcnvcalls = allcnvcalls.append(dat_del_sel)
    allcnvcalls.to_csv(covdir+"/allsamples.calledevents.byFDR"+str(100*FDR)+".txt",sep="\t",index=False)   
def mergeFiles_threshold(covdir):
    covfiles = files_in_dir(covdir, str2srch = ".CNVcallsThreshold")
    allcnvcalls = pd.DataFrame()
    for covfile in covfiles:
        dat = pd.read_csv(covfile,sep="\t")
        samname = re.sub(".cnvZscores.CNVcallsThreshold","",os.path.basename(covfile),count=1)
        dat_amp = dat.loc[dat["CNVcalls_threshold"]==1]
        dat_amp_sel = dat_amp[["#chr","name"]]
        dat_amp_sel["sample"] = samname
        dat_amp_sel["Z"] = dat_amp["gc.corrected.norm.log.std.index.zWeighted.Final"]
        dat_amp_sel["event"] = "amplification"
        dat_del = dat.loc[dat["CNVcalls_threshold"]==-1]
        dat_del_sel = dat_del[["#chr","name"]]
        dat_del_sel["sample"] = samname
        dat_del_sel["Z"] = dat_del["gc.corrected.norm.log.std.index.zWeighted.Final"]
        dat_del_sel["event"] = "deletion"
        allcnvcalls = allcnvcalls.append(dat_amp_sel)
        allcnvcalls = allcnvcalls.append(dat_del_sel)
    allcnvcalls.to_csv(covdir+"/allsamples.calledevents.byThreshold.txt",sep="\t",index=False)   
    