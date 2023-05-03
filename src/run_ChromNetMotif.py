#import utilities_parallel
import numpy as np
import pandas as pd
from optparse import OptionParser
import time
import random
import sys
from utilities_parallel import *
from motif_viz import *

'''
Parsing input section
'''

parser = OptionParser()
parser.add_option("-g", "--network", action="store", type="string", dest="chromatin_state_network_file")
parser.add_option("-o", "--output", action="store", type="string", dest="output_prefix",default="test")
parser.add_option("-n",type="int", dest="num_random",default = 100)
parser.add_option("-m",type="int", dest="motif_size",default = 3)
parser.add_option("-p",type="int", dest="num_processes",default = 1)
parser.add_option("-t",type="float", dest="p_value",default = 0.05)
parser.add_option("-f",type="int", dest="min_freq",default = 1)




(options, args) = parser.parse_args()
print("==============================================================")
print("chromatin state network file:",options.chromatin_state_network_file)
#print("chromatin state file:",options.chromatin_file)
print("output prefix:",options.output_prefix)
print("no. of random networks:",options.num_random)
print("motif size:",options.motif_size)
print("no. of processes requested:",options.num_processes)
print("p-value:",options.p_value)
print("minimum frequency of motif in real network:",options.num_processes)
print("==============================================================")
chromatin_state_network_file = options.chromatin_state_network_file
#network_file = options.network_file
#chromatin_file = options.chromatin_file
output_prefix = options.output_prefix
motif_size = options.motif_size
num_random = options.num_random
num_processes = options.num_processes
p_value = options.p_value
min_freq = options.min_freq


#(gtrie_network_df,chromatin_net_df) = prepare_network_files(chromatin_state_network_file)
#network_file = output_prefix +'.gtrie.txt'
#chromatin_file = output_prefix+'.chromatin.state.net.csv'
#gtrie_network_df.to_csv(network_file,sep=' ',header=False,index=False)
#chromatin_net_df.to_csv(chromatin_file,sep=',',header=True,index=False)

print("==============================================================")
chromatin_state_df=pd.read_csv(chromatin_state_network_file,sep=",")
broad_states = list(set(chromatin_state_df.iloc[:,2].to_list() + chromatin_state_df.iloc[:,3].to_list()))
print("Chromatin States detected ")
print(broad_states)
print("==============================================================")


if __name__ == '__main__':

 print("Detecting subgraphs in real network ")
 real_results = get_subgraphs(chromatin_state_network_file,motif_size,output_prefix,broad_states)
 #print(real_results)
 
 print("Running permutations on random graphs ")
 start_time = time.time()
 network_file = output_prefix +'.gtrie.txt'
 chromatin_file = output_prefix + '.chromatin.state.net.csv'
 
 results = get_random_subgraphs_parallel(network_file,motif_size,broad_states,chromatin_file,num_random,num_processes)
 results_df = pd.concat(results,ignore_index=True)
 results_df.to_csv(output_prefix+".random.motifs.csv",index=False)
 command = 'rm *random.subgraphs.txt' 
 os.system(command)
 command = 'rm *random.results.txt' 
 os.system(command) 
 command = 'rm *random.txt' 
 os.system(command)
 print("Total time taken --- %s seconds ---" % (time.time() - start_time))
 




 # compare real and random networks
 d = compare_random_real_motifs(real_results,results_df,num_random)
 #print(d)
 d.to_csv(output_prefix + ".motif.results.csv",index=False)
 
 d = pd.read_csv(output_prefix + ".motif.results.csv")
 print(broad_states)
 if motif_size  == 3:
     results_df = motif_results_to_plots_3(d,broad_states,p_value,min_freq,output_prefix)
     results_df.to_csv(output_prefix + ".motif.results.csv",index=False)
 if motif_size == 4:
     results_df = motif_results_to_plots_4(d,broad_states,p_value,min_freq,output_prefix)
     results_df.to_csv(output_prefix + ".motif.results.csv",index=False)
     


 #visualize motifs
 #generate_heatmap_motif_3(motif_code_string_list,broad_states)
 
#import time
#start_time = time.time()
#get_random_subgraphs('GM12878.gtrie.network.csv',3,1,10,'GM12878.network.new.csv')
#print("--- %s seconds ---" % (time.time() - start_time))

#get_random_subgraphs(network_file,motif_size,1,num_random,chromatin_file)
#get_subgraphs('GM12878.gtrie.network.csv',4,'GM12878.4nodes','GM12878.network.new.csv')

