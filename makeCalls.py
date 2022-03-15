#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:37:28 2021

@author: MShahrokh
"""

import numpy as np
import pandas as pd
import sys, getopt
import optparse
from cnv_sidefuncs import *
def main(argv):
    parser = optparse.OptionParser()
    parser.add_option('-t', '--thresholds', action="store", dest="thresholdFile",help="Path to the thresholds if already calculated", default=None)
    parser.add_option('-r', '--results', action="store", dest="results_dir",help="Path to the results", default=None)
    parser.add_option('-c', '--normalPath', action="store", dest="normalPath",help="Path to the .cnvZscore files of a set of normals to calculate the thresholds from", default=None)
    options, args = parser.parse_args()
    if options.thresholdFile!=None:
        callCNV_normal(options.results_dir,options.thresholdFile)
        mergeFiles_threshold(options.results_dir)
    if options.normalPath!=None and options.thresholdFile==None:
        options.thresholdFile = threshold(options.normalPath,zCol = "gc.corrected.norm.log.std.index.zWeighted.Final", qThresh = 0.99)
        callCNV_normal(options.results_dir,options.thresholdFile)
        mergeFiles_threshold(options.results_dir)
    callCNV_fdr(options.results_dir,zCol = "gc.corrected.norm.log.std.index.zWeighted.Final")
    mergeFiles_fdr(options.results_dir,FDR = 0.01)
    mergeFiles_fdr(options.results_dir,FDR = 0.05)
    mergeFiles_fdr(options.results_dir,FDR = 0.1)
if __name__ == "__main__":
  main(sys.argv[1:])    
   