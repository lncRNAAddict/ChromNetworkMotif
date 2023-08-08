import pandas as pd
import numpy as np
import networkx as nx
import os
import time 
from joblib import Parallel, delayed, parallel_backend
import multiprocessing as mp
import random

def prepare_network_files(chromatin_state_network_file):
    #print("preparing gtrie network file for extracting subgraphs\n")
    chromatin_state_net_df=pd.read_csv(chromatin_state_network_file,sep=",")
    gtrie_network = chromatin_state_net_df.iloc[:,[0,1]].copy(deep=True)
    chromatin_state_network = chromatin_state_net_df.copy(deep=True)
    chromatin_locs = list(set(chromatin_state_net_df.iloc[:,0].to_list() + chromatin_state_net_df.iloc[:,1].to_list()))
    chromatin_locs_dict = {}
    i = 1
    for locs in chromatin_locs:
        chromatin_locs_dict[str(i)] = locs
        #print(i)
        gtrie_network[gtrie_network == locs] = i
        chromatin_state_network[chromatin_state_network == locs] = i
        i=i+1
    gtrie_network['weight'] = 1
    return(gtrie_network,chromatin_state_network,chromatin_locs_dict)
        
    

'''
Extract all subgraphs of size = k

G_file = name of the network file in gtrie format
size = size of the subgraph (motif) to be extracted
output = prefix to the name of the output files
chromatin_state_file = File containing chromatin states of nodes in G_file
'''
#('GM12878.gtrie.network.csv',4,'GM12878','GM12878.network.csv')
def get_subgraphs(chromatin_state_network_file,size,output,broad_states):
    
    (gtrie_network_df,chromatin_net_df,chromatin_locs_dict) = prepare_network_files(chromatin_state_network_file)
    G_file = output +'.gtrie.txt'
    chromatin_state_file = output + '.chromatin.state.net.csv'
    gtrie_network_df.to_csv(G_file,sep=' ',header=False,index=False)
    chromatin_net_df.to_csv(chromatin_state_file,sep=',',header=True,index=False)
    
    command = 'gtrieScanner -s ' +str(size) + ' esu -g ' + G_file + ' -oc ' + output + ".subgraphs.txt -o " + G_file+".results.txt"
    os.system(command + " > " + G_file +".tmp.console.txt")
    
    #df = get_chromatin_states(chromatin_matrix_file)
    subgraph_df = pd.read_csv(output + ".subgraphs.txt",sep=" ",header=None)
    #print(subgraph_df)
    #states = ['TssA','TssFlnk','TssFlnkU','TssFlnkD','Tx','TxWk','EnhG1','EnhG2','EnhA1','EnhA2','EnhWk','ZNF/Rpts','Het','TssBiv','EnhBiv','ReprPC','ReprPCWk','Quies']
    
    d = pd.read_csv(G_file,sep=' ',header=None)
    d.columns=['node1','node2','edge']
    G = nx.from_pandas_edgelist(d,'node1','node2')
    chromatin_state_df=pd.read_csv(chromatin_state_file,sep=",")
    #(chromatin_state_motif_df,chromatin_state_motif_df_occurences)=get_chromatin_state_motifs(subgraph_df,G,chromatin_state_df,size,broad_states).to_frame(name='count').reset_index()
    (chromatin_state_motif_df,chromatin_state_motif_df_occurences)=get_chromatin_state_motifs(subgraph_df,G,chromatin_state_df,size,broad_states) #.to_frame(name='count').reset_index()
    chromatin_state_motif_df=chromatin_state_motif_df.to_frame(name='count').reset_index()
    #print(chromatin_state_motif_df_occurences)
    chromatin_state_motif_df.to_csv(output+'.motifs.txt',header=False,index=False)
    if size == 3:
      chromatin_state_motif_df_occurences[['node1','node2','node3']] = chromatin_state_motif_df_occurences['nodes'].str.split("#",expand=True)
      chromatin_state_motif_df_occurences['node1']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node1'].to_list() ]
      chromatin_state_motif_df_occurences['node2']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node2'].to_list() ]
      chromatin_state_motif_df_occurences['node3']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node3'].to_list() ]
    if size == 4:
      chromatin_state_motif_df_occurences[['node1','node2','node3','node4']] = chromatin_state_motif_df_occurences['nodes'].str.split("#",expand=True)
      chromatin_state_motif_df_occurences['node1']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node1'].to_list() ]
      chromatin_state_motif_df_occurences['node2']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node2'].to_list() ]
      chromatin_state_motif_df_occurences['node3']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node3'].to_list() ]        
      chromatin_state_motif_df_occurences['node4']=[chromatin_locs_dict[node] for node in  chromatin_state_motif_df_occurences['node4'].to_list() ]        

    
    chromatin_state_motif_df_occurences = chromatin_state_motif_df_occurences.drop(columns=['nodes'])
    chromatin_state_motif_df_occurences.to_csv(output+'.motifs.locations.csv',header=True,index=False)
    return chromatin_state_motif_df
    #algorithm
    '''
    repeat for 1000 times
    randomize the network
    Extract all subgraphs of size = k
    For each type of chromatin state marked motif, get count

    For each type of chromatin state marked motif, get p-value, and normalized z-score
    '''
    #G_random = nx.double_edge_swap(G,5*G.number_of_edges(),10*G.number_of_edges())

import random
def get_one_random_subgraphs(G,size,chromatin_state_df,broad_states,method='I'):
        run_no = round(random.uniform(1,1000), 5)
        output_file=str(run_no)+".random.motifs.txt"
        #output_file=str(row.iloc[0,0])+".random.motifs.txt"
        #if os.path.isfile(output_file):
        #  os.system('rm ' + output_file) 
        #run_no = row['col1']
        output=str(run_no)

        G_random = nx.double_edge_swap(G,5*G.number_of_edges(),10*G.number_of_edges())  
            
        random_G_df = nx.to_pandas_edgelist(G_random)
        random_G_df['edge']=1
        
        random_G_df.to_csv(output+'random.txt',sep=" ",header=None,index=None)
        command = 'gtrieScanner -s ' + str(size) + ' esu -g  ' + output +'random.txt -oc '+output+'.random.subgraphs.txt -o ' +output+'.random.results.txt'
        os.system(command + " > random.tmp.console.txt")
        random_subgraph_df = pd.read_csv(output+".random.subgraphs.txt",sep=" ",header=None)
        
        chromatin_state_df_ = chromatin_state_df.copy(deep=True)
        
        if method == 'I':
            chromatin_state_df_ = chromatin_state_df_
            #print("randomization method  = I")
            
        elif method == 'II':
            #print("randomization method  = II")
            d1=dict(zip(chromatin_state_df_.iloc[:,0],chromatin_state_df_.iloc[:,2]))
            d2=dict(zip(chromatin_state_df_.iloc[:,1], chromatin_state_df_.iloc[:,3]))
            d2.update(d1)
            states = list(d2.values())
            #for k in range(0,10):
            random.shuffle(states)
            d2 = dict(zip(d2,states))
            chromatin_state_df_.iloc[:,2]=[d2[node] for node in chromatin_state_df_.iloc[:,0].to_list()]
            chromatin_state_df_.iloc[:,3]=[d2[node] for node in chromatin_state_df_.iloc[:,1].to_list()]
        
        else:
            print('randomization method options should be I or II')
            exit()
            
        random_motif_df,tmp_ = get_chromatin_state_motifs(random_subgraph_df,G_random,chromatin_state_df_,size,broad_states)#.to_frame(name='count').reset_index()
        random_motif_df = random_motif_df.to_frame(name='count').reset_index()
        random_motif_df['iteration']=run_no
        #if not os.path.isfile(output_file):
        #   random_motif_df.to_csv(output_file,index=False,header=False)
        #else:
        #   random_motif_df.to_csv(output_file,mode='a',index=False,header=False)
        #os.system('rm ' + output +'random.txt')
        #os.system('rm ' + output +'random.subgraphs.txt')
        return(random_motif_df)


def get_random_subgraphs_parallel(G_file,motif_size,broad_states,chromatin_state_file,num_random,n_processors=1,rand_method = 'I'):
# Step 1: Init multiprocessing.Pool()
  
  
  print("Number of cores detected = ",mp.cpu_count())
  print("Number of processes requested = ",n_processors)
  if(mp.cpu_count() < n_processors):
    n_processors = mp.cpu_count()
    print("number of processes used = ", n_processors)

  pool = mp.Pool(n_processors)

  d = pd.read_csv(G_file,sep=' ',header=None)
  d.columns=['node1','node2','edge']
  G = nx.from_pandas_edgelist(d,'node1','node2')

  chromatin_state_df=pd.read_csv(chromatin_state_file,sep=",")
  #print(chromatin_state_df)
  results = Parallel(n_jobs=n_processors)(delayed(get_one_random_subgraphs)(G,motif_size,chromatin_state_df,broad_states,rand_method) for row in range(num_random))
  pool.close()
  #print(results)

  
  return results  
    
'''
nswap = number of double edges swap to perform
max_tries = integer (optional) Maximum number of attempts to swap edges
'''
def get_random_subgraphs(G_file,size,start,end,chromatin_state_file):
    #big_random_motif=[]
    output_file=str(start)+"-"+str(end)+".random.motifs.txt"
    if os.path.isfile(output_file):
       os.system('rm ' + output_file) 
    for i in range(start,end):
        #print(i)
        output=str(i)
        d = pd.read_csv(G_file,sep=' ',header=None)
        d.columns=['node1','node2','edge']
        G = nx.from_pandas_edgelist(d,'node1','node2')
	    #eprint(i)
        #random_G_df = randomize_df(d,10)
        #G_random = nx.from_pandas_edgelist(random_G_df,'node1','node2')
        G_random = nx.double_edge_swap(G,G.number_of_edges(),10*G.number_of_edges())
        #for t in range(0,10):
        #G_random = nx.double_edge_swap(G_random,int(G_random.number_of_edges()/2))
        random_G_df = nx.to_pandas_edgelist(G_random)
        random_G_df['edge']=1
        random_G_df.to_csv(output+'random.txt',sep=" ",header=None,index=None)
        command = 'gtrieScanner -s ' +str(size) + ' esu -g  ' + output +'random.txt -oc '+output+'random.subgraphs.txt -o random.results.txt'
        os.system(command)
        random_subgraph_df = pd.read_csv(output+"random.subgraphs.txt",sep=" ",header=None)
        random_motif_df,tmp_=get_chromatin_state_motifs(random_subgraph_df,G_random,chromatin_state_file,size).to_frame(name='count').reset_index()
        random_motif_df['iteration']=i
        #print(random_motif_df)
        if not os.path.isfile(output_file):
           random_motif_df.to_csv(output_file,index=False,header=False)
        else:
           random_motif_df.to_csv(output_file,mode='a',index=False,header=False)
        os.system('rm ' + output +'random.txt')
        os.system('rm ' + output +'random.subgraphs.txt')


def get_chromatin_state_motifs_one_row(row,G,d2,broad_states,motif_size,num_broad_states):
       chromatin_motif = np.zeros(shape=(motif_size-1,num_broad_states),dtype='int')
       #print(chromatin_motif)
       #print("===========================\n")
       #print(motif_size,num_broad_states)
       motif=[]
       #print(row)
       H = G.subgraph(row.iloc[1:].to_list())
       for node in H.nodes():
              motif.append(H.degree[node])
              #print("node",node)
              state_index = broad_states.index(d2[node])
              chromatin_motif[H.degree[node]-1,state_index] = chromatin_motif[H.degree[node]-1,state_index] + 1 
       
       #numpy 2d array to string
       chromatin_motif=np.reshape(chromatin_motif,(chromatin_motif.shape[0]*chromatin_motif.shape[1]))
       
       
       chromatin_motif=''.join([str(i) for i in chromatin_motif])
       
       motif=''.join([str(i) for i in sorted(motif,reverse=True)])
       #return [motif,chromatin_motif]
       return motif+":"+chromatin_motif+":" + "#".join([str(node) for node in H.nodes()])
    
def get_chromatin_state_motifs(subgraph_df,G,chromatin_state_df,motif_size,broad_states):
   
    
    d1=dict(zip(chromatin_state_df.iloc[:,0],chromatin_state_df.iloc[:,2]))
    d2=dict(zip(chromatin_state_df.iloc[:,1], chromatin_state_df.iloc[:,3]))
    d2.update(d1)
    #print(d2)
    #broad_states = list(set(chromatin_state_df.iloc[:,2].to_list() + chromatin_state_df.iloc[:,3].to_list()))
    #print("itendified unqiues states",broad_states)
    num_broad_states = len(broad_states)
    
    #motif_=[]
    #chromatin_motif_=[] 
    #print(broad_states)
    
    result = subgraph_df.apply(get_chromatin_state_motifs_one_row,args=(G,d2,broad_states,motif_size,num_broad_states),axis=1)
    #print(result)
    motif_ = [motif.split(":")[0] for motif in result]
    chromatin_motif_ = [motif.split(":")[1] for motif in result]
    chromatin_motif_nodes = [motif.split(":")[2] for motif in result]
    
       #chromatin_motif_df=pd.DataFrame(columns=['motif','chromatin_state'])
    d = {'motif':motif_ , 'chromatin_motif':chromatin_motif_}
    chromatin_motif_df=pd.DataFrame(data=d)
    chromatin_motif_df_nodes = chromatin_motif_df.copy(deep=True)
    chromatin_motif_df_nodes['nodes'] = chromatin_motif_nodes
    chromatin_motif_df=chromatin_motif_df.groupby(['motif', 'chromatin_motif']).size()
       #chromatin_motif_df['motif']=motif_
       #chromatin_motif_df['chromatin_motif']=chromatin_motif_
    return (chromatin_motif_df,chromatin_motif_df_nodes)

'''   
def get_chromatin_state_motifs(subgraph_df,G,chromatin_state_file,motif_size):
   
    chromatin_state_df=pd.read_csv(chromatin_state_file,sep=",",header=None)
    d1=dict(zip(chromatin_state_df.iloc[:,0],chromatin_state_df.iloc[:,2]))
    d2=dict(zip(chromatin_state_df.iloc[:,1], chromatin_state_df.iloc[:,3]))
    d2.update(d1)
    
    broad_states = list(set(chromatin_state_df.iloc[:,2].to_list() + chromatin_state_df.iloc[:,3].to_list()))
    print("itendified unqiues states",broad_states)
    num_broad_states = len(broad_states)
    
    motif_=[]
    chromatin_motif_=[] 
    for i in range(0,subgraph_df.shape[0]):
       #states = ['TssA','TssFlnk','TssFlnkU','TssFlnkD','Tx','TxWk','EnhG1','EnhG2','EnhA1','EnhA2','EnhWk','ZNF/Rpts','Het','TssBiv','EnhBiv','ReprPC','ReprPCWk','Quies']
       #broad_states =['active','weak','poised','repressed']

       chromatin_motif = np.zeros(shape=(motif_size-1,num_broad_states),dtype='int')
       
       motif=[]
       H = G.subgraph(subgraph_df.iloc[i,1:].to_list())
       for node in H.nodes():
              motif.append(H.degree[node])
              state_index = broad_states.index(d2[node])
              chromatin_motif[H.degree[node]-1,state_index] = chromatin_motif[H.degree[node]-1,state_index] + 1 
       
       #numpy 2d array to string
       chromatin_motif=np.reshape(chromatin_motif,(chromatin_motif.shape[0]*chromatin_motif.shape[1]))
       chromatin_motif=''.join([str(i) for i in chromatin_motif])
       chromatin_motif_.append(chromatin_motif)
       #print(chromatin_motif)
       motif=''.join([str(i) for i in sorted(motif,reverse=True)])
       motif_.append(motif)
       #f=open(output+".motifs.txt","w")
       #f.write(motif+","+chromatin_motif+"\n")
       #f.close()
       
       #chromatin_motif_df=pd.DataFrame(columns=['motif','chromatin_state'])
    d = {'motif':motif_ , 'chromatin_motif':chromatin_motif_}
    chromatin_motif_df=pd.DataFrame(data=d)
    chromatin_motif_df=chromatin_motif_df.groupby(['motif', 'chromatin_motif']).size()
       #chromatin_motif_df['motif']=motif_
       #chromatin_motif_df['chromatin_motif']=chromatin_motif_
    return(chromatin_motif_df)
'''


'''
compare real and random motifs occurences
compute p-value, z-score
'''

import statistics
def compare_random_real_motifs(real_motifs_df,random_motifs_df,num_random):
    #real_motifs_df = pd.read_csv(real_motifs_file,dtype={'motif': 'str','chromatin_motif':'str','count':'int64'})
    #real_motifs_df.columns=['motif','chromatin_motif','count']
    real_motifs_df['group']=real_motifs_df['motif'].astype('str') + "_" + real_motifs_df['chromatin_motif'].astype('str')
    #random_motifs_df = pd.read_csv(random_motifs_file,dtype={'motif': 'str','chromatin_motif':'str','count':'int64','iteration':'int64'})
    #random_motifs_df.columns=['motif','chromatin_motif','count','iteration']
    random_motifs_df['group']=random_motifs_df['motif'].astype('str') +"_" + random_motifs_df['chromatin_motif'].astype('str')

    p_value = []
    z_score = [] 
    real_freq = []
    motifs = []
    random_mean_freq =[]
    for group in real_motifs_df['group']:
        motifs.append(group)
        real_count = real_motifs_df[real_motifs_df['group']==group]['count'].to_list()[0]
        #print(real_count)
        tmp = random_motifs_df[random_motifs_df['group'] == group].copy()
        if(tmp.shape[0] > 0):
           m=statistics.mean(tmp['count'].to_list())
           random_mean_freq.append(m)
           if(tmp.shape[0] > 1):
              sd=statistics.stdev(tmp['count'].to_list())
           else:
              sd=0
           tmp_=tmp.loc[tmp['count']>=real_count]
           #print(tmp.shape)
           p_value.append(tmp_.shape[0]/num_random)
           if sd > 0:
             z_score.append((real_count-m)/sd)
           else:
             z_score.append('NA')
        else:
           p_value.append(0)
           z_score.append('NA')
           random_mean_freq.append(0)
        real_freq.append(real_count)   
    d=pd.DataFrame(columns=['motifs','p_value','real_count','random_mean_count'])
    d['motifs']=motifs
    d['p_value']=p_value
    #d['z_score']=z_score
    d['real_count']=real_freq
    d['random_mean_count']=random_mean_freq
    return(d)



#(A,B) = prepare_network_files('test.csv')


