# Python code for running evolutionary dynamics

The following Python codes are included in this folder:

- `bd-on-hypergraph.py` simulates the birth-death process under model 1 on an arbitrary connected hypergraph. It returns the fixation probability for a single mutant for each fitness value, $r$.
- `bd-on-one-mode-projection.py` simulates the birth-death process under model 1 on the weighted one-mode projection of an arbitrary connected hypergraph. It returns the fixation probability for a single mutant for each fitness value, $r$.

Here we provide an artificial hypergraph data as an example:

- `artificial_edgelist.txt` the first column represents the node ID. The second column represents the ID of a hyperedge which the node belongs to.
