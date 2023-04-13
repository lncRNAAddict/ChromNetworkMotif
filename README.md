# ChromNetMotif
`ChromNetMotif` is a python tool to extract significantly enriched network motifs (of size 3 and 4) marked by chromatin states. The nodes in the network are marked by chromatin states.

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

## Running `run_ChromNetMotif.py`

### Input Chromatin state network file format

There is a header in the file. The header is marked as `from_node,to_node,from_broad_state,to_broad_state`
Each row in the file represents an undirected edge in the network with additional information about the chromatin state of the nodes in the edge.

The network file MUST have the following 4 columns. Two consecutive columns must be separated by comma.

- Column 1: Node 1 in the edge. Must be in chromosome:start-end format.
- Column 2: Node 2 in the edge. Must be in chromosome:start-end format.
- Column 3: Chromatin state of Node 1. Must be a string.
- Column 4: Chromatin state of Node 2. Must be a string.

A screeshot of a portion of an example input file is shown below

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/chromatin_state_file.jpg)


ChromNetkMotif takes the following arguments

- `-g` or `--network`: Name of the input network file. This MUST be provided.
- '-o' or `--output`: Prefix to name of output files to be generated. If not provided, default is `test`.
- `-n` : number of random networks to be generated fpr motif enrichment testing. Default is 100.
- `-m` : number of nodes in the motif or motif size to be extracted. Default is 3. It can only handle motif of size 3 and 4. 
- `-p` : number of processors to be used. Default is 1.

First, change working directory to the folder where the `ChromNetworkMotif` scripts are stored. You can do that by simply typing the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

`cd <path to where ChromNetworkMotif scripts are stored>`

Once the working directory is set, shown below is an example of running `ChromNetworkMotif`.

`python run_chromnetmotif.py -g <path to the network file>/my_network_file.csv -c <path to the network file>/my_chromatin_file.csv -o <path to output folder>/my_output -n 500 -m 4 -p 4`

The code above will extract chromatin state marked motifs of size 4 and generate all output files with prefix `my_output` in the folder `<path to output folder>`. P-values and Z-scores are assigned based on 500 randomized networks. Computation will be done using 4 processors in parallel.
  


### Outputs

- `<output>.motif.results.txt`: This file contains the chromatin-state marked motifs present in the network. This first column represent the chromatin code and motif code for the chromatin-state marked motif. Combination of these two columns represent one unique chromatin-state marked motif. The next 3 or four columns represent the chromatin state of the nodes in the motif. The subsequenct columns indicate p-value, z-score, motif frequency in the network, mean motif frequency in random networks. 
An example file is shown below.

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/motif_results.JPG)
- `<output>.<motif_code>.png`: For each `motif code`, one `png` file to vizualize the chromatin-state marked motifs will be generated. For example, for motif of size size 3 with the motif code `222`, an example `png` file may look like.

![chromatin-state marked motifs for motif of size 3 with motif code `222`](https://github.com/lncRNAAddict/ChromNetworkMotif/blob/main/Figures/hela_motif_3.222.png)
- `<output>.motif.results.txt`
