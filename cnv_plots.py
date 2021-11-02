#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 12:48:23 2021

@author: MShahrokh
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_ideogram_individual(covfile, genomeBED, stroff, stron ,ylabel, savename):
    savefile = savename + ".ideogram.png"
    gnm = pd.read_csv(genomeBED,header=None,sep="\t")
    chrs = ["chr"+str(i) for i in range(1,23)]
    gnm = gnm[gnm.iloc[:,0].isin(chrs)]
    gsize = sum(gnm.iloc[:,1])
    gnm['fraction'] = np.cumsum(gnm.iloc[:,1]/gsize)
    gnm['cumulative'] = gnm['fraction'] - gnm.iloc[:,1]/gsize/2
    gnm['size2'] = np.cumsum(gnm.iloc[:,1])-gnm.iloc[:,1]
    
    gnm =  gnm.rename(columns = {0:'#chr', 1:"size"})
    df = pd.read_csv(covfile,header=0,sep="\t")
    df_merged = pd.merge(df,gnm,on=['#chr'],how = 'left')
    center = list(df_merged["start"])+list(df_merged["end"])
    position = [float(x)+y for x,y,z in zip(center,df_merged['size2'],df_merged['size'])]
    position2 = [float(x)+y+50000 for x,y,z in zip(center,df_merged['size2'],df_merged['size'])]
    df_1 = pd.DataFrame({'#chr': list(df_merged['#chr']), 
                         'position' : [a /gsize for a, b, c in zip(position, list(df_merged['fraction']),list(df_merged['size2']) )],
                         'name': list(df_merged['name']),
                         'zScore': list(df_merged[stroff]),
                         'fraction' : list(df_merged['fraction']),
                         'cumulative':list(df_merged['cumulative']),
                         'type': 4})
    df_2 = pd.DataFrame({'#chr': list(df_merged['#chr']),
                         'position' : [a  /gsize for a, b, c in zip(position2, list(df_merged['fraction']),list(df_merged['size2']))],
                         'name': list(df_merged['name']),'zScore': list(df_merged[stron]),
                         'fraction' : list(df_merged['fraction']),
                         'cumulative':list(df_merged['cumulative']),
                         'type': 10})
    dfs = [df_1,df_2]
    df2plot = pd.concat(dfs)
    df2plot["#chr"] = df2plot["#chr"].astype('category')
    df2plot["#chr"] = df2plot["#chr"].cat.set_categories(['chr%i' % i for i in range(1,23)], ordered=True)
    df2plot = df2plot.sort_values('#chr')
    df_grouped = df2plot.groupby(('#chr'))
    
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    fig.suptitle(savename, fontsize=16)
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        #group.plot(kind='scatter', x='position', y='zScore',color=colors[num % len(colors)], ax=ax, s=group['type'])
        plt.scatter(x=group['position'], y=group['zScore'],color=colors[num % len(colors)], s=list(group['type']))
        plt.axvline(x=group['fraction'].iloc[0],linestyle='--')
        plt.axhline(y=1.959964,linestyle='dotted')
        plt.axhline(y=-1.959964,linestyle='dotted')
        x_labels.append(name)
        x_labels_pos.append((group['cumulative'].iloc[0]))
        #x_labels_pos.append((group['position']))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels,rotation=45)
    ax.set_xlim([0, 1])
    ax.set_ylim([-20, 20])
    
    # x axis label
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylabel)
    plt.savefig(savefile)
        # show the graph
    #plt.show()


def plot_ideogram_individual_singlescore(covfile, genomeBED, str2plot,ylabel, savename):
    savefile = savename + ".ideogram.png"
    gnm = pd.read_csv(genomeBED,header=None,sep="\t")
    chrs = ["chr"+str(i) for i in range(1,23)]
    gnm = gnm[gnm.iloc[:,0].isin(chrs)]
    gsize = sum(gnm.iloc[:,1])
    gnm['fraction'] = np.cumsum(gnm.iloc[:,1]/gsize)
    gnm['cumulative'] = gnm['fraction'] - gnm.iloc[:,1]/gsize/2
    gnm['size2'] = np.cumsum(gnm.iloc[:,1])-gnm.iloc[:,1]
    
    gnm =  gnm.rename(columns = {0:'#chr', 1:"size"})
    df = pd.read_csv(covfile,header=0,sep="\t")
    df_merged = pd.merge(df,gnm,on=['#chr'],how = 'left')
    center = list(df_merged["start"])+list(df_merged["end"])
    position = [float(x)+y for x,y,z in zip(center,df_merged['size2'],df_merged['size'])]
    position2 = [float(x)+y+50000 for x,y,z in zip(center,df_merged['size2'],df_merged['size'])]
    df_1 = pd.DataFrame({'#chr': list(df_merged['#chr']), 
                         'position' : [a /gsize for a, b, c in zip(position, list(df_merged['fraction']),list(df_merged['size2']) )],
                         'name': list(df_merged['name']),
                         'zScore': list(df_merged[str2plot]),
                         'fraction' : list(df_merged['fraction']),
                         'cumulative':list(df_merged['cumulative']),
                         'type': 4})
      
    dfs = df_1
    df2plot = dfs
    df2plot["#chr"] = df2plot["#chr"].astype('category')
    df2plot["#chr"] = df2plot["#chr"].cat.set_categories(['chr%i' % i for i in range(1,23)], ordered=True)
    df2plot = df2plot.sort_values('#chr')
    df_grouped = df2plot.groupby(('#chr'))
    
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    fig.suptitle(savename, fontsize=16)
    ax = fig.add_subplot(111)
    colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        #group.plot(kind='scatter', x='position', y='zScore',color=colors[num % len(colors)], ax=ax, s=group['type'])
        plt.scatter(x=group['position'], y=group['zScore'],color=colors[num % len(colors)], s=list(group['type']))
        plt.axvline(x=group['fraction'].iloc[0],linestyle='--')
        plt.axhline(y=1.959964,linestyle='dotted')
        plt.axhline(y=-1.959964,linestyle='dotted')
        x_labels.append(name)
        x_labels_pos.append((group['cumulative'].iloc[0]))
        #x_labels_pos.append((group['position']))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels,rotation=45)
    ax.set_xlim([0, 1])
    ax.set_ylim([-20, 20])
    
    # x axis label
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylabel)
    plt.savefig(savefile)
        # show the graph
    #plt.show()