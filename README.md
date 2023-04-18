# ChromNetMotif
`ChromNetMotif` is a python tool to extract significantly enriched network motifs (of size 3 and 4) marked by chromatin states. The nodes in the network are marked by chromatin states.

## External dependencies

`ChromNetMotif` has been tested with Python 3.8.12 version. It should work with Python >= 3.0 version. We recommend installing the anaconda python distributon. Download anaconda python distribution from https://www.anaconda.com/products/distribution and install following the instructions provided.

PySmooth depends on the following python libraries. These libraries are already included in the anaconda distribution. Therefore, you do not need to install them.

- `numpy`
- `Pandas`
- `Sklearn`
- `matplotlib`


`ChromNetMotif` requires that the following external tool  is already installed.

- `gtrieScanner`

`gtrieScanner` is a tool that extract subgraphs using the `g-trie` data structure. This tool can be installed from the github site https://github.com/ComplexNetworks-DCC-FCUP/gtrieScanner.

After installation, add the following lines in the `~/.bash_profile` or `~/.bashrc` file.

`PATH=$PATH:<path to installation folder of gtrieScanner>`

`export PATH`


## Installation and Dependencies


You can simply download the following scripts from `PySmooth` GitHub page and put them in a single folder. 

- `utilities.py`
- `smooth.py`
- `ImputeMissingGenotype.py`
- `run_smooth.py`


## Usage

`ChromNetMotif` is executed using the python script `run_ChromNetMotif.py` in the python command line.

`run_ChromNetkMotif.py` takes the following arguments

- `-g` or `--network`: Name of the input chromatin state network file. This MUST be provided.
- '-o' or `--output`: Prefix to name of output files to be generated. If not provided, default is `test`.
- `-n` : number of random networks to be generated fpr motif enrichment testing. Default is 500.
- `-m` : number of nodes in the motif or motif size to be extracted. Default is 3. It can only handle motif of size 3 and 4. 
- `-p` : number of processors to be used. Default is 1.
- `-t` : p-value threshold to detect statistically enriched motifs. Default value is 0.05.

First, change working directory to the folder where the `ChromNetworkMotif` scripts are stored. You can do that by simply typing the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

`cd <path to where ChromNetworkMotif scripts are stored>`

Once the working directory is set, shown below is an example of running `ChromNetworkMotif`.

`python run_chromnetmotif.py -g <path to the network file>/my_network_file.csv  -o my_output -n 500 -m 3 -p 4 -t 0.05`

The code above will extract chromatin state marked motifs of size 4 and generate all output files with prefix `my_output`. Computation will be done using 4 processors in parallel. P-value of 0.05 will be used to identify significant motifs and for visualization.
  
### Input Chromatin state network file format

The chromatin state network file contains a header, which should read as `from_node,to_node,from_broad_state,to_broad_state`.
Each row in the file represents an undirected edge in the network with additional information about the chromatin state of the nodes in the edge.

The network file MUST have the following 4 columns. Two consecutive columns must be separated by comma.

- Column 1: Node 1 in the edge. Must be in chromosome:start-end format.
- Column 2: Node 2 in the edge. Must be in chromosome:start-end format.
- Column 3: Chromatin state of Node 1. Must be a string.
- Column 4: Chromatin state of Node 2. Must be a string.

A screeshot of a portion of an example input file is shown below

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/chromatin_state_file.jpg)

### Outputs

- `<output>.motif.results.txt`: This file contains the chromatin-state marked motifs present in the network. This first column represent the chromatin code and motif code for the chromatin-state marked motif. Combination of these two columns represent one unique chromatin-state marked motif. The next 3 or four columns represent the chromatin state of the nodes in the motif. The subsequenct columns indicate p-value, z-score, motif frequency in the network, mean motif frequency in random networks. 
An example file is shown below.

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_results.JPG)
- `<output>.<motif_code>.png`: For each `motif code`, one `png` file to vizualize the chromatin-state marked motifs will be generated. For example, for motif of size size 3 with the motif code `222`, an example `png` file may look like as shown below. Motif code `222` represents the motif of size 3, where every node is connected to the other nodes. Each row in the heatmap represents a chromatin-state marked motif. In the example heatmap shown below, there are four possible chromatin states: `weak`, `repressed`, `poised`, and `active` which are color coded. 

![chromatin-state marked motifs for motif of size 3 with motif code `222`](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/hela_motif_3.222.png)
- `<output>.motifs.locations.txt`: This file contains all instances of each chromatin-state marked motifs in the network. motif code, chromatin motif code, the interacting chromatin nodes. 

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_location.JPG)
