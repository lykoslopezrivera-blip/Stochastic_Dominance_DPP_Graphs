# Stochastic Dominance in Discrete Parking Processes on Graphs

Code accompanying the paper **"Stochastic Dominance in Discrete Parking Processes on Graphs"** by Robinson Higuita-Diaz, Lykos López and Alejandro Roldán-Correa.

The repository contains the exact computations used in the paper to compare the *sequential* and *static* models of the discrete parking process (unfriendly seating problem) on arbitrary finite simple graphs: enumeration of all graphs up to a given size, exact distributions of the number of occupied vertices, stochastic dominance checks, counterexample search, and the numerical evidence behind the proof of the main theorem in the paper and the conjecture on graphs of maximum degree at most 3.
## Authors

- Robinson Higuita-Diaz
- Lykos López
- Alejandro Roldán-Correa


## Repository structure
```text
Stochastic_Dominance_DPP_Graphs/
├── GeneralArrangement/
│   ├── Notebooks/
│   │   └── StructureAndExamples.ipynb
│   └── Source/
│       ├── Calc_PerGraph.py
│       ├── GeneralSearch.py
│       └── Graph.py
├── Conjecture/
│   ├── Notebooks/
│   │   └── verification.ipynb
│   └── Source/
│       ├── RandomGraph.py
│       └── Verify_General_D_1.py
├── LinearArrangement/
│   ├── Estimates_EN.ipynb
│   └── r_n_monoticity.ipynb
├── requirements.txt
├── .gitignore
└── README.md
```


### Requirements

- Python 3.12
- Jupyter Notebook or JupyterLab
- Python packages listed in `requirements.txt` (NumPy, Matplotlib, NetworkX)
- **nauty** (only for `GeneralSearch.py`), which provides the `geng` executable

