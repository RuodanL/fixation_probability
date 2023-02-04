# Fixation probability on hypergraphs

Here we provide the code that accompanies the following paper:

R. Liu, N. Masuda.
Fixation dynamics on hypergraphs.
Preprint on arXiv: [arXiv:2301.05343](https://arxiv.org/abs/2301.05343)

When you use our code for your publications, please cite this paper.

## Jupyter notebook to create figures

```
figure_X.ipynb
```
contains the code to produce Figure X in the manuscript, where $X = 2, 4, 5$, and E.6.
The folder [notebooks/data/](https://github.com/RuodanL/fixation_probability/tree/main/notebooks/data) contains the numerical results based on which Figures 2, 4, 5, and E.6 are produced.

## Python code for running evolutionary dynamics

The following Python codes are also included:

- `bd-on-hypergraph.py` to simulate the birth-death process on an arbitrary connected hypergraph under model 1 and get the fixation probability for a single mutant for each fitness value, $r$.
- `bd-on-one-mode-projection.py` to simulate the birth-death process on the weighted one-mode projection of an arbitrary connected hypergraph and get the fixation probability for a single mutant for each fitness value, $r$.

## Installation and Dependencies

The `pentapy` package can be installed via [pip](https://pypi.org/project/pentapy/). On Windows, you can install [WinPython](https://winpython.github.io) to get Python and pip. Then, run

```
pip install pentapy
```
There are pre-built wheels for Linux, MacOS and Windows for most Python versions (2.7, 3.4-3.7).

The `SciPy` package can be installed via pip from [PyPI](https://pypi.org/project/scipy/) as follows.

```
pip install scipy
```
Alternatively, SciPy is part of the [Anaconda](https://docs.continuum.io/anaconda/) distribution and can be installed with Anaconda or Miniconda as follows:

```
conda install scipy
```

