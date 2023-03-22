# ChromNetworkMotif
ChromNetMotif is a python tool to extract enriched chromatin-state marked network motifs (of size 3 and 4).

## External dependencies
- `gtrieScanner`
## Input



4438 4439 1

4436 4439 1

4434 4439 1

4435 4439 1

4436 4437 1

4444 4445 1

4440 4445 1

4434 4445 1

4435 4445 1


# PySmooth

Here, we present a python implementation of SMOOTH (van Os,H. et al. (2005) ) called `PySmooth` which offers an easy-to-use command line interface and solves the drawbacks of SMOOTH. PySmooth reads the input genotype file and identifies singletons based on the algorithm described in SMOOTH with some modifications to allow four genotype codes, and flexible parameters. Unlike SMOOTH which doesn’t correct the singletons and missing data, PySmooth corrects genotype errors using a k-nearest algorithm. At each step, PySmooth generates summary files and visualizations that can be inspected by the user for further interpretation.


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

