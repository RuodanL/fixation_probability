# Running the birth-death process under model 1

The following Python codes are included in this folder:

- `bd-on-hypergraph.py` simulates the birth-death process on an arbitrary connected hypergraph under model 1. It returns the fixation probability for a single mutant for each fitness value, $r$.
- `bd-on-one-mode-projection.py` simulates the birth-death process on the weighted one-mode projection of an arbitrary connected hypergraph under model 1. It returns the fixation probability for a single mutant for each fitness value, $r$.

Here we provide an artificial hypergraph data, `artificial_edgelist.txt`, as an example. Each row represents the membership of a node to a hyperedge. The first column represents the node ID. It runs from 0. The second column represents the ID of a hyperedge which the node belongs to. It also runs from 0.

All the other files contain the results presented in the article.
