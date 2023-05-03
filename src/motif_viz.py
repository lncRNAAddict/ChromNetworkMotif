#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:47:28 2023
@author: soibamb
"""
import numpy as np
import pandas as pd
import networkx as nx
import os
import time 
import re
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def motif_code_to_colors_3(code_string,broad_states):
        colormap=[0,0,0]
        num_broad_states = len(broad_states)
    
        # B-A-C motif
        if re.search(r"[^0]",code_string[0:num_broad_states]):
            code_substring_1 = code_string[0:num_broad_states]
            code_substring_2 = code_string[num_broad_states:]
            broad_index = 0
            start = 0
            for code in code_substring_1:
                code_ = int(code)
                if code_  > 0:
                    colormap[start:start+code_] = [broad_index for i in range(start,start+code_)]
                    start = start + code_
                broad_index = broad_index + 1
            
            broad_index = code_substring_2.index('1')
            tmp = colormap[1]
            #print(colormap)
            colormap[1] = broad_index
            colormap[2] = tmp
         
        # B-A-C-B motif
        else: 
            code_substring = code_string[num_broad_states:]
            broad_index = 0
            start = 0
            for code in code_substring:
                code_ = int(code)
                if code_  > 0:
                    colormap[start:start+code_] = [broad_index for i in range(start,start+code_)]
                    start = start + code_
                broad_index = broad_index + 1
        colormap_state = [broad_states[i] for i in colormap]
        return colormap,colormap_state

def motif_code_to_colors_4(code_string,broad_states):
    
        num_broad_states = len(broad_states)
        flag_1 = re.findall(r"[1234]",code_string[0:num_broad_states])
        flag_2 = re.findall(r"[1234]",code_string[num_broad_states:2*num_broad_states])
        flag_3 = re.findall(r"[1234]",code_string[2*num_broad_states:])
        #print(code_string)
        if len(flag_1) == 0 and len(flag_2) == 0 and len(flag_3) > 0:
            colormap = motif_4_type1(code_string,num_broad_states)
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
        
        elif len(flag_1) == 0 and len(flag_2) > 0 and len(flag_3) > 0:
            colormap= motif_4_type2(code_string,num_broad_states)
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
        
        elif len(flag_1) == 0 and len(flag_2) > 0  and len(flag_3) == 0:
            colormap= motif_4_type3(code_string,num_broad_states)  
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
        
        elif len(flag_1) > 0 and len(flag_2) > 0 and len(flag_3) > 0:
            colormap= motif_4_type4(code_string,num_broad_states)   
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
        
        elif len(flag_1) > 0 and len(flag_2) > 0 and len(flag_3) == 0:
            colormap= motif_4_type5(code_string,num_broad_states)
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
        
        elif len(flag_1) > 0 and len(flag_2) == 0 and len(flag_3) > 0:
            colormap= motif_4_type6(code_string,num_broad_states)
            colormap_state = [broad_states[i] for i in colormap ]
            return colormap,colormap_state
  
        '''   
        elif len(flag_1) == 0 and len(flag_2) > 0 and len(flag_3) > 0:
            colormap= motif_4_type7(code_string,num_broad_states)
            colormap_state = [broad_states[i] for i in colormap ]
            #print(colormap,colormap_state)
            return colormap,colormap_state
     '''
  
def motif_4_type1(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring = code_string[2*num_broad_states:]
    broad_index = 0
    start = 0
    for code in code_substring:
         code_ = int(code)
         if code_  > 0:
              colormap[start:start+code_] = [broad_index for i in range(start,start+code_)]
              start = start + code_
         broad_index = broad_index + 1
    
    return colormap
        
def motif_4_type2(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring_1 = code_string[num_broad_states:2*num_broad_states]
    code_substring_2 = code_string[2*num_broad_states:]
    
    x_2 = [i for i, letter in enumerate(code_substring_1) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_1) if letter == '1']
    if len(x_2) > 0:
        colormap[1] = x_2[0]
        colormap[3] = x_2[0]
    if len(x_1) > 0:
        colormap[1] = x_1[0]
        colormap[3] = x_1[1]
        
    x_2 = [i for i, letter in enumerate(code_substring_2) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_2) if letter == '1']
    if len(x_2) > 0:
        colormap[0] = x_2[0]
        colormap[2] = x_2[0]
    if len(x_1) > 0:
        colormap[0] = x_1[0]
        colormap[2] = x_1[1]  
    
    return colormap
        
def motif_4_type3(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring= code_string[num_broad_states:2*num_broad_states]
    #print(code_substring)
    broad_index = 0
    start = 0
    for code in code_substring:
        code_ = int(code)
        #print(code_)
        if code_  > 0:
            colormap[start:start+code_] = [broad_index for i in range(start,start+code_)]
            start = start + code_
        broad_index = broad_index + 1
    #print(colormap)
    return colormap

def motif_4_type4(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring_1 = code_string[0:num_broad_states]
    code_substring_2 = code_string[num_broad_states:2*num_broad_states]
    code_substring_3 = code_string[2*num_broad_states:]
   
    x_1 = [i for i, letter in enumerate(code_substring_1) if letter == '1']
    colormap[1] = x_1[0]

    x_2 = [i for i, letter in enumerate(code_substring_2) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_2) if letter == '1']
    if len(x_2) > 0:
        colormap[0] = x_2[0]
        colormap[3] = x_2[0]
    if len(x_1) > 0:
        colormap[0] = x_1[0]
        colormap[3] = x_1[1]  
        
    x_3 = [i for i, letter in enumerate(code_substring_3) if letter == '1']
    colormap[2] = x_3[0]    
    
    return colormap

def motif_4_type5(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring_1 = code_string[0:num_broad_states]
    code_substring_2 = code_string[num_broad_states:2*num_broad_states]
    
    x_2 = [i for i, letter in enumerate(code_substring_1) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_1) if letter == '1']
    if len(x_2) > 0:
        colormap[0] = x_2[0]
        colormap[3] = x_2[0]
    if len(x_1) > 0:
        colormap[0] = x_1[0]
        colormap[3] = x_1[1]
        
    x_2 = [i for i, letter in enumerate(code_substring_2) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_2) if letter == '1']
    if len(x_2) > 0:
        colormap[1] = x_2[0]
        colormap[2] = x_2[0]
    if len(x_1) > 0:
        colormap[1] = x_1[0]
        colormap[2] = x_1[1]
    
    return colormap
'''
def motif_4_type6(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring_2 = code_string[num_broad_states:2*num_broad_states]
    code_substring_3 = code_string[2*num_broad_states:]
    
    x_2 = [i for i, letter in enumerate(code_substring_2) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_2) if letter == '1']
    if len(x_2) > 0:
        colormap[0] = x_2[0]
        colormap[1] = x_2[0]
    if len(x_1) > 0:
        colormap[0] = x_1[0]
        colormap[1] = x_1[1]
        
    x_2 = [i for i, letter in enumerate(code_substring_3) if letter == '2']
    x_1 = [i for i, letter in enumerate(code_substring_3) if letter == '1']
    if len(x_2) > 0:
        colormap[2] = x_2[0]
        colormap[3] = x_2[0]
    if len(x_1) > 0:
        colormap[2] = x_1[0]
        colormap[3] = x_1[1]  
    
    return colormap
'''
def motif_4_type6(code_string,num_broad_states):
    colormap=[0,0,0,0]
    code_substring_1= code_string[0:num_broad_states]
    code_substring_3= code_string[2*num_broad_states:]
    #print(code_substring)
    broad_index = 0
    start = 0
    for code in code_substring_1:
        code_ = int(code)
        #print(code_)
        if code_  > 0:
            colormap[start:start+code_] = [broad_index for i in range(start,start+code_)]
            start = start + code_
        broad_index = broad_index + 1 
    
    x_1 = [i for i, letter in enumerate(code_substring_3) if letter == '1']
    colormap[3] = x_1[0]
    return colormap


def generate_heatmap_motif_4(motif_code_string_list,broad_states,output):
    colors = [(1, 0, 0), (1, 1, 0), (0, 0, 1)] # R -> yellow -> B
    cmap_name = 'state_list'

    num_broad_states = len(broad_states)
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=num_broad_states)

    x = [motif_code_to_colors_4(code_string,broad_states)[0] for code_string in motif_code_string_list]
    x_state = [motif_code_to_colors_4(code_string,broad_states)[1] for code_string in motif_code_string_list]
    
    fig, ax = plt.subplots()
    im=ax.imshow(np.array(x), vmin=0,vmax=(num_broad_states-1),cmap=cmap,aspect='auto')
    ax.set_xticks(range(0,4))
    ax.set_xticklabels(['A', 'B','C','D'])
    ax.set_yticks(range(0,len(motif_code_string_list)))
    #ax.set_yticklabels(motif_code_string_list)
    ax.set_yticklabels([])
    cb = fig.colorbar(im,ticks=range(0,num_broad_states))
    cb.ax.set_yticklabels(broad_states) 
    plt.savefig(output+".png")
    #plt.imshow()
    #print(x) 
    #print(x_state)
    
    return pd.DataFrame(x_state,columns=['A', 'B','C','D'])
    
def generate_heatmap_motif_3(motif_code_string_list,broad_states,output):

    colors = [(1, 0, 0), (1, 1, 0), (0, 0, 1)] # R -> y -> B
    cmap_name = 'state_list'

    num_broad_states = len(broad_states)
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=num_broad_states)
    
    x = [motif_code_to_colors_3(code_string,broad_states)[0] for code_string in motif_code_string_list]
    x_state = [motif_code_to_colors_3(code_string,broad_states)[1] for code_string in motif_code_string_list]
    #x_state = [motif_code_to_colors_4(code_string,broad_states) for code_string in motif_code_string_list]
    fig, ax = plt.subplots()
    im=ax.imshow(np.array(x), vmin=0,vmax=(num_broad_states-1),cmap=cmap,aspect='auto')
    ax.set_xticks(range(0,3))
    ax.set_xticklabels(['A', 'B','C'])
    print(broad_states)
    ax.set_yticks(range(0,len(motif_code_string_list)))
    #ax.set_yticklabels(motif_code_string_list)
    ax.set_yticklabels([])
    
    cb = fig.colorbar(im,ticks=range(0,num_broad_states))
    cb.ax.set_yticklabels(broad_states) 
    plt.savefig(output+".png")
    #plt.imshow()
    #print(x)
    #print(x_state)
    #print(pd.DataFrame(x_state,columns=['A', 'B','C']))
    return pd.DataFrame(x_state,columns=['A', 'B','C'])

def motif_results_to_plots_3(motif_results_df,broad_states,p_value =0.05,min_freq = 1, output="motif_heatmap"):
    motif_results_df_ = motif_results_df.loc[motif_results_df['p_value'] < p_value].copy(deep=True)
    motif_results_df_ = motif_results_df_.loc[motif_results_df_['real_count'] > motif_results_df_['random_mean_count']]
    motif_results_df_ = motif_results_df_.loc[motif_results_df_['real_count'] > min_freq]
    print(motif_results_df_)
    motif_results_df_['motif_type']=motif_results_df_['motifs'].str.split("_",expand=True).iloc[:,0]
    motif_results_df_['motif_code']=motif_results_df_['motifs'].str.split("_",expand=True).iloc[:,1]
    
    uniq_motifs = list(set(motif_results_df_['motif_type'].to_list()))
    dframe=list()
    
    for motif in uniq_motifs:
        #print(motif)
        motif_results_df_tmp = motif_results_df_.loc[motif_results_df_['motif_type'] == motif].copy(deep=True)
        #print(motif_results_df_tmp)
        d = generate_heatmap_motif_3(motif_results_df_tmp['motif_code'].to_list(),broad_states,output+"."+motif)
        d['p_value'] = motif_results_df_tmp['p_value'].to_list()
        #d['z_score'] = motif_results_df_tmp['z_score'].to_list()
        d['real_count'] = motif_results_df_tmp['real_count'].to_list()
        d['random_mean_count'] = motif_results_df_tmp['random_mean_count'].to_list()
        d.insert(0,'motif_type',motif_results_df_tmp['motif_type'].to_list())
        d.insert(0,'chromatin_code',motif_results_df_tmp['motif_code'].to_list())
        dframe.append(d)
    result_df = pd.concat([dframe[i] for i in range(len(dframe))], axis=0)
    result_df['chromatin_code']= result_df['chromatin_code'].astype(str)
    return(result_df)
  

def motif_results_to_plots_4(motif_results_df,broad_states,p_value =0.05,min_freq=1,output="motif_heatmap"):
    motif_results_df_ = motif_results_df.loc[motif_results_df['p_value'] < p_value].copy(deep=True)
    motif_results_df_ = motif_results_df_.loc[motif_results_df_['real_count'] > motif_results_df_['random_mean_count']]
    motif_results_df_ = motif_results_df_.loc[motif_results_df_['real_count'] > min_freq]
    motif_results_df_['motif_type']=motif_results_df_['motifs'].str.split("_",expand=True).iloc[:,0]
    motif_results_df_['motif_code']=motif_results_df_['motifs'].str.split("_",expand=True).iloc[:,1]
    #print(motif_results_df_)
    uniq_motifs = list(set(motif_results_df_['motif_type'].to_list()))

    dframe=list()
    for motif in uniq_motifs:
        print(motif)
        motif_results_df_tmp = motif_results_df_.loc[motif_results_df_['motif_type'] == motif].copy(deep=True)
        #print(motif_results_df_tmp)
        d = generate_heatmap_motif_4(motif_results_df_tmp['motif_code'].to_list(),broad_states,output+"."+motif)
        d['p_value'] = motif_results_df_tmp['p_value'].to_list()
        #d['z_score'] = motif_results_df_tmp['z_score'].to_list()
        d['real_count'] = motif_results_df_tmp['real_count'].to_list()
        d['random_mean_count'] = motif_results_df_tmp['random_mean_count'].to_list()
       
        #print(d.iloc[:,0:6])
        d.insert(0,'motif_type',motif_results_df_tmp['motif_type'].to_list())
        d.insert(0,'chromatin_code',motif_results_df_tmp['motif_code'].to_list())
        dframe.append(d)
    result_df = pd.concat([dframe[i] for i in range(len(dframe))], axis=0)
    result_df['chromatin_code']= result_df['chromatin_code'].astype(str)
    return result_df
#motif_code_string_list =['000000400000','000000000400','000000001111','000000000121']
#d = pd.read_csv('GM12878.motif.results.csv')
#motif_results_to_plots_3(d,['Active','Weak','Poised','Repressed'],p_value =0.05)
#generate_heatmap_motif_4(motif_code_string_list,['Active','Weak','Poised','Repressed'])



'''
def motif_visualization(motifs_df,broad_states,motif_size=3,pvalue=0.05):
    motifs_df_1 = motifs_df.loc[motifs_df['p_value']<pvalue,:]
    
    if motif_size == 3:
        
        
        
    else:
        
        print("Motif size has to be 3 or 4")
    
  '''  
'''  
colormap = motif_code_to_colors('00000300',['a','b','c','d'])
print(colormap)
colormap = motif_code_to_colors('01100010',['a','b','c','d'])
print(colormap)
colormap = motif_code_to_colors('00000111',['a','b','c','d'])
print(colormap)
'''
