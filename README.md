# ChromNetworkMotif
ChromNetMotif is a python tool to extract enriched chromatin-state marked network motifs (of size 3 and 4).

## External dependencies
- `gtrieScanner`
## Input

- `network file` with three columns, columns are separated by a single space. Each row in the file represents an undirected edge in the network, the first and second columns represent the two nodes in the edge, the third column represents the edge weight. Since `ChromNetworkMotif` only handles undirected and unweighted networks, the third column is always set to `1`. The nodes are indicated by integer values in the first and second column.

4438 4439 1

4436 4439 1

4434 4439 1

4435 4439 1

4436 4437 1

4444 4445 1

4440 4445 1

4434 4445 1

4435 4445 1
