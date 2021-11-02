#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 23:19:41 2021

@author: MShahrokh
"""

# This script is for cell-of-origin classification of DLBCL
# As input, it requires copy number information, somatic point mutations and translocations.
# Features are based on Esfahani, et al., Blood 2019. 
# The script is also capable of producing labels for the NCI (EZB, etc) & DFCI labels (C1-C5)


def coo_classifier(featuremat, featureweights, mode = "COO"):
    # The classification mode can be COO, NCI (EZB, BN2, etc.), or DFCI (C1-C5). 
    