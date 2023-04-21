# ChromNetMotif

`ChromNetMotif` is a tool that offers an easy-to-use Python command line interface to extract chromatin-state marked motifs from a chromatin interaction network data. The tool can extract occurrences and frequencies of different chromatin-state marked motifs but, also computes p-values and z-scores to quantify the statistical enrichment of the motifs in the network. The statistical metrics are computed by comparing to several random networks the tools generates internally. Visualization files are also generated that allow the user to interpret the motifs. Because the statistical enrichment calculation can be computationally time intensive, ChromNetMotif contains feature that allows analysis on random networks in parallel in a multicore processor environment. 

## External dependencies

`ChromNetMotif` has been tested with Python 3.8.12 version and runs on `macOS` and `Linux` environments. It should work with Python >= 3.0 version. We recommend installing the anaconda Python distributon. Download anaconda python distribution from https://www.anaconda.com/products/distribution and install following the instructions provided.

`ChromNetMotif` depends on the following python libraries. 

- `Numpy`
- `Pandas`
- `Matplotlib`
- `Joblib`: For multicore processing. 
- `networkx`: For using network-related data structures and algorithms. 

The libraries `Numpy`, `Pandas`, and `Matplotlob` are already included in the anaconda distribution. Therefore, you do not need to install them. Installation instructions for `Joblib` and `networkx` can be found in https://joblib.readthedocs.io/en/latest/ and https://networkx.org/, respectively.

`ChromNetMotif` requires that the following external tool  is already installed.

- `gtrieScanner`

`gtrieScanner` is a tool that extract subgraphs using the `g-trie` data structure. This tool can be installed from the github site https://github.com/ComplexNetworks-DCC-FCUP/gtrieScanner.

After installation, add the following lines in the `~/.bash_profile` or `~/.bashrc` file.

`PATH=$PATH:<path to installation folder of gtrieScanner>`

`export PATH`

## Installation


You can simply download the following scripts from `ChromNetMotif` GitHub page and put them in a single folder. 

- `utilities_parallel.py`
- `motif_viz.py`
- `run_ChromNetMotif.py`


## Usage

`ChromNetMotif` is executed using the python script `run_ChromNetMotif.py` in the python command line.

`run_ChromNetkMotif.py` takes the following arguments

- `-g` or `--network`: Name of the input chromatin state network file. This MUST be provided.
- '-o' or `--output`: Prefix to name of output files to be generated. If not provided, default is `test`.
- `-n` : number of random networks to be generated fpr motif enrichment testing. Default is 500.
- `-m` : number of nodes in the motif or motif size to be extracted. Default is 3. It can only handle motif of size 3 and 4. 
- `-p` : number of processors to be used. Default is 1.
- `-t` : p-value threshold to detect statistically enriched motifs. Default value is 0.05.
- `-f` : minimum frequency in the real network for a chromatin-state marked motif to be counted as statistically enriched. Default value is 1. There might be some chromatin-state marked motifs which occur at very instances in the network. The user has the choice to disregard such motifs using this threshold.

First, change working directory to the folder where the `ChromNetworkMotif` scripts are stored. You can do that by simply typing the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

`cd <path to where ChromNetworkMotif scripts are stored>`

Once the working directory is set, shown below is an example of running `ChromNetworkMotif`.

`python run_ChromNetMotif.py -g <path to the network file>/my_network_file.csv  -o my_output -n 500 -m 3 -p 4 -t 0.05`

The code above will extract chromatin-state marked motifs of size 4 and generate all output files with prefix `my_output`. Computation will be done using 4 processors in parallel. P-value of 0.05 will be used to identify significant motifs selected for visualization.
  
### Input Chromatin state network file format

The chromatin state network file shoule contains a header that reads as `from_node,to_node,from_broad_state,to_broad_state`.
Each row in the file represents an undirected edge in the network with additional information about the chromatin state of the nodes in the edge.

The network file MUST have the following 4 columns. Two consecutive columns must be separated by comma.

- Column 1: Node 1 in the edge. Must be in `chromosome:start-end format`.
- Column 2: Node 2 in the edge. Must be in `chromosome:start-end format`.
- Column 3: Chromatin state of Node 1. Must be a string.
- Column 4: Chromatin state of Node 2. Must be a string.

A screeshot of a portion of an example input file is shown below


`Example chromatin state network file`

![](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/chromatin_state_file.JPG)


### Outputs

`ChromNetMotif` generates three main output files: `<my_output>.motif.results.txt`, `<output>.<motif_code>.png`, and  `<output>.motifs.locations.csv`.  They are described below.


- `<my_output>.motif.results.txt`: This file contains all the extracted chromatin-state marked motifs present in the network. An example file is shown below.

`Example chromatin-state marked motif (size = 3) results file `

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_results.JPG)

The second column in the `<my_output>.motif.results.txt` file represents the motif code for the chromatin-state marked motif. The first column in the `<my_output>.motif.results.txt` file represents the chromatin code that describes the chromatin states of the nodes involved in the motif. Combination of the first two columns in the file represent one unique chromatin-state marked motif. The next 3 or 4 columns represent the chromatin states of the nodes in the motif. The subsequenct columns indicate p-value, z-score, motif frequency in the actual network, mean motif frequency in random networks. Note that all chromatin-state marked motifs (whether they are enriched enriched or not are present in this file).

`ChromNetMotif` uses two unique motif codes for motifs of size 3 and 6 unique codes for motif of size 4 as shown below. 

`Motif codes for unique motifs of size 3 and 4`

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_code.JPG)

- `<output>.motifs.locations.txt`: This file contains all instances of each chromatin-state marked motifs in the network. motif code, chromatin motif code, the interacting chromatin nodes. 

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_location.JPG)


### Chromatin-state marked motifs visualizations

- `<output>.<motif_code>.png`: For each `motif code`, one `png` file to vizualize the chromatin-state marked motifs will be generated. Two example output images corresponding to two motif codes of motif size = 3 are shown below. Each row in the heatmap represents a chromatin-state marked motif. In the example heatmap shown below, there are four possible chromatin states: `weak`, `repressed`, `poised`, and `active` which are color coded. In this file only the motifs which are statistically enriched in the 

- For motif of size size 3 with the motif code `211`, there are three chromatin-state marked motifs.



![](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_size_3_heatmaps.JPG)


## Example networks
The folder `Data` contains a network file.

