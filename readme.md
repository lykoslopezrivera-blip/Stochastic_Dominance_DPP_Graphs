# Stochastic Dominance in Discrete Parking Processes on Graphs

Code accompanying the paper **"Stochastic Dominance in Discrete Parking Processes on Graphs"** by Robinson Higuita-Diaz, Lykos López and Alejandro Roldán-Correa.

The repository contains the exact computations used in the paper to compare the *sequential* and *static* models of the discrete parking process (unfriendly seating problem) on arbitrary finite simple graphs: enumeration of all graphs up to a given size, exact distributions of the number of occupied vertices, stochastic dominance checks, counterexample search, and the numerical evidence behind the conjecture on graphs of maximum degree at most 3.

## Abstract

In a recent study, it was established that the number of occupied sites in the sequential discrete parking process first-order stochastically dominates its static counterpart on a path graph (also referred to as a line graph). In this paper, we explore whether this property holds for general graph topologies. We show that the stochastic dominance extends naturally to cycle graphs. Furthermore, an exhaustive computational analysis confirms that the dominance holds for all 31 connected simple graphs with fewer than 6 vertices. For graphs with 6 or more vertices, we demonstrate that the property does not hold universally by constructing explicit counterexamples. We also show that the property is inherited by disjoint unions when all connected components share it, and we exhibit mixed unions in which stochastic comparability can be preserved or entirely lost. Finally, supported by algorithmic computations, we formulate a conjecture proposing that stochastic dominance is preserved for all graphs with maximum degree at most 3.

## Authors

- Robinson Higuita-Diaz
- Lykos López
- Alejandro Roldán-Correa


## Repository structure

Stochastic_Dominance_DPP_Graphs/
├── GeneralArrangement/          # Exhaustive analysis of arbitrary graphs
│   ├── Notebooks/
│   │   └── StructureAndExamples.ipynb
│   └── Source/
│       ├── Calc_PerGraph.py
│       ├── GeneralSearch.py
│       └── Graph.py
├── Conjecture/                  # Evidence for the Δ(G) ≤ 3 conjecture
│   ├── Notebooks/
│   │   └── verification.ipynb
│   └── Source/
│       ├── RandomGraph.py
│       └── Verify_General_D_1.py
├── requirements.txt
├── .gitignore
└── README.md


## Getting started

### Requirements

- Python 3.12
- Jupyter Notebook or JupyterLab
- Python packages listed in `requirements.txt` (NumPy, Matplotlib, NetworkX)
- **nauty** (only for `GeneralSearch.py`), which provides the `geng` executable

