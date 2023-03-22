# ChromNetworkMotif
ChromNetMotif is a python tool to extract enriched chromatin-state marked network motifs (of size 3 and 4).

## External dependencies
- `gtrieScanner`


## Installation and Dependencies

PySmooth has been tested with Python 3.8.12 version. It should work with Python >= 3.0 version. We recommend installing the anaconda python distributon. Download anaconda python distribution from https://www.anaconda.com/products/distribution and install following the instructions provided.

PySmooth depends on the following python libraries. These libraries are already included in the anaconda distribution. Therefore, you do not need to install them.

- `numpy`
- `Pandas`
- `Sklearn`
- `matplotlib`

You can simply download the following scripts from `PySmooth` GitHub page and put them in a single folder. 

- `utilities.py`
- `smooth.py`
- `ImputeMissingGenotype.py`
- `run_smooth.py`

`PySmooth` can be executed by running the script `run_smooth.py`

## Running `run_chromnetmotif.py`

### Input Network File format

There is no header in the file. Each row in the file represents an undirected edge in the network. 

The network file MUST have the following 3 columns. Two consecutive columns must be separated by single space.

- Column 1: Node 1 in the edge. Must be integer.
- Column 2: Node 2 in the edge. Must be integer.
- Column 3: edge weight. Since `ChromNetworkMotif` only handles undirected and unweighted networks, the third column is always set to `1`. 

A screeshot of a portion of an example input file is shown below

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/network_example_File.jpg)


### Input Chromatin state File format

There is a header in the file. The header is marked as `from_node,to_node,from_broad_state,to_broad_state`
Each row in the file represents an undirected edge in the network with additional information about the chromatin state of the nodes in the edge.

The network file MUST have the following 4 columns. Two consecutive columns must be separated by comma.

- Column 1: Node 1 in the edge. Must be integer.
- Column 2: Node 2 in the edge. Must be integer.
- Column 3: Chromatin state of Node 1. Must be a string.
- Column 4: Chromatin state of Node 2. Must be a string.

A screeshot of a portion of an example input file is shown below

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/chromatin_state_file.jpg)


ChromNetworkMotif takes the following arguments

- `-g` or `--network`: Name of the input network file. This MUST be provided.
- `-c` or `--chromatin`: Name of the chromatin state file. This MUST be provided.
- '-o' or '--output': Prefix to name of output files to be generated. If not provided, default is `test`.
- `-n` : number of random networks to be generated fpr motif enrichment testing. Default is 100.
- `-m` : number of nodes in the motif or motif size to be extracted. Default is 3.
- `-p` : number of processors to be used. Default is 1.

### Outputs

- <output>.motifs.txt: This file contains the 
