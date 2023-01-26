# Fixation probability on hypergraphs

Here we provide the code that accompanies the following paper:

R. Liu, N. Masuda.
Fixation dynamics on hypergraphs.
Preprint on arXiv: [arXiv:2301.05343](https://arxiv.org/abs/2301.05343)

## Jupyter notebook to create figures

```
figure_X.ipynb
```
contains the code to produce fig X in the manuscript, where $X = 2, 4, 5$, and E.6.
The folder [notebooks/data/](https://github.com/RuodanL/fixation_probability/tree/main/notebooks/data) contains the numerical results based on which Figures 2, 4, 5, and E.6 are produced.

## Python code for running evolutionary dynamics

The following Python codes are also included:

- bd-on-hypergraph.py -- to simulate the birth-death process on an arbitrary connected hypergraph under model 1 and get the fixation probability of a single mutant for each fitness value of $r$.
- bd-on-one-mode-projection.py -- to simulate the birth-death process on the weighted one-mode projection of an arbitrary connected hypergraph and get the fixation probability of a single mutant for each fitness value of $r$.
