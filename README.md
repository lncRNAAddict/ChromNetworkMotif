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

## Running `run_smooth.py`

### Input Network File format

There is no header in the file. Each row in the file represents an undirected edge in the network.

The network file MUST have the following 3 columns. Two consecutive columns must be separated by single space.

- Column 1: Node 1 in the edge. Must be integer.
- Column 2: Node 2 in the edge. Must be integer.
- Column 3: edge weight. Since `ChromNetworkMotif` only handles undirected and unweighted networks, the third column is always set to `1`. 

A screeshot of a portion of an example input file is shown below

![Example Input Genotype File](https://github.com/lncRNAAddict/PySmooth/blob/main/example/GenotypeInput.PNG)

### Running PySmooth

PySmooth takes the following arguments

- `-i` or `--input`: Name of the input genotype file. This MUST be provided.
- '-o' or '--output': Prefix to name of output files to be generated. If not provided, default is `test`.
- `-c` or `--chr`: list of chromosome names to perform analysis on. Names should be separated by comma (e.g `chr1,chr2,chr3`). Default is to run through all the chromosomes in the genotype file.
- `-l` or `--lower`: Lowest threshold for identifying singletons. Default is 0.70.
- `-u`or `--upper`: Highest threshold for identifying singletons. Default is 0.98.
- `-g` or `--gap`: PySmooth iteratively identifies singletons starting with the highest threshold till the lowest threshold. This parameter is used to decreased the threshold at each iteration. Default is 0.02.
- `-k` : number of nearest neighbors to be used to assign correct genotype to singleton or missing item. Default value is 30.

First, change working directory to the folder where the `PySmooth` scripts are stored. You can do that by simply typing the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

`cd <path to where PySmooth scripts are stored>`

Once the working directory is set, shown below are two examples of running `PySmooth`.

`python run_smooth.py -i <path to the genotype file>/my_genotype_file.csv`

The code above will analyze  each chromosome detected and generate all output files with prefix `test` in the folder `<path to the genotype file>`
  
`python run_smooth.py -i <path to the genotype file>/my_genotype_file.csv -o <path to output folder>/my_output -c chr1 -l 0.80 -u 0.98 -g 0.02`

The code above will analyze for chromosome `chr1`and generate all output files with prefix `my_output` in the folder `<path to output folder>`.

### Outputs

For each chromosome, PySmooth Generates the following outputs.

- Three summary csv files: `<output>_<chr>.stats.csv`, `<output>_<chr>_singletons_stats.csv`, and `<output>_<chr>_imputed_stats.csv` that contain `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file, after singleton detection, and after error correction. Examples are shown below.

- Three bar plot png files: `<output>_<chr>.stats.png`, `<output>_<chr>_singletons_stats.png`, and `<output>_<chr>_imputed_stats.png` bar plot files that contains `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file, after singleton detection, and after error correction, respectively. Example images are shown below.

![alt text](https://github.com/lncRNAAddict/ChromNetworkMotif/Figures/Slide3.PNG)

- Three heatmap files: `<output>_<chr>.heatmap.png`, `<output>_<chr>_singletons_heatmap.png`, and `<output>_<chr>_imputed_heatmap.png` that visualize a color-coded image of different genotype codes in the original file, after singleton detection, and after error correction, respectively. Example images are shown below.

![alt text](https://github.com/lncRNAAddict/PySmooth/blob/main/example/Slide2.PNG)

- `<output>_<chr>_singletons.csv`: genotype file with singleton detected. Singletons are marked as `S`. 
- `<output>_<chr>_imputed.csv`: genotype file after error correction.




### References
van Os,H. et al. (2005) SMOOTH: a statistical method for successful removal of genotyping errors from high-density genetic linkage data. Theor. Appl. Genet., 112, 187–94.

