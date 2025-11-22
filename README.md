# Engineering_holography

The code is used to find the optimal graph for two instance of the hyperbolic pentagon code [Pastawski et al.,
JHEP 2015:149 (2015)]:

1) Small instance of 16 qubits (Main_16.py):
   
     * Option 1 (Value='all'): Explore all non-ismorphic graphs using the library gsc* and print the ones with minimal edges.
   
     * Option 2 (Value='ch'): Select the particular set of Hadamard gates to obtain the graph state from the manuscript.  
    
     
2) Bigger instance of 36 qubits (Main_36.py)
   
     * Option 1 (Value='opt'): Found a better graph by applying Hadamard gates in different four pentagons plaquettes a couple of times.
   
     * Option 2 (Value='ch'): Select the particular set of Hadamard gates to obtain the graph state from the manuscript.  

In both cases the program also gives you the files GraphBM.pdf and GraphAM.Pdf. The first shows the layout of the pentagons before the meaurement and the second the resulting graph after the measurement (without optimizing).


*The library gsc is a modified version of https://github.com/sammorley-short/gsc which works for python3.

The preprint of the work is in https://arxiv.org/abs/2209.08954.


## DAGI validation experiment

`dagi_validation.py` automates the feasibility study outlined in
“DAGI Validation with Stabilizer Graph Codes.pdf”. It reconstructs the 16-qubit
graph state using the manuscript’s Hadamard choices, then runs three families of
information-theoretic checks against bulk qubit `15`:

- **Recovery wedge**: the five-neighbour set `[0, 1, 2, 12, 14]` that retrieves
  two bits of bulk information with a pure fifth-order (`f_5`) contribution.
- **Outside-wedge control**: a disjoint five-qubit set `[3, 4, 5, 6, 7]` that
  yields zero mutual information and zero `f_k`, confirming no leakage outside
  the entanglement wedge.
- **Fragment-growth controls**: smaller subsets (e.g., `[0,1,2]`, `[0,2,14]`,
  `[0,1,2,12]`) that all produce zero mutual information and `f_k` until the
  full recovery set is reached, supporting the absence of partial access before
  threshold.

To reproduce the calculation and store `dagi_results.json`:

```bash
python -m pip install -r requirements.txt
python dagi_validation.py
```

The script prints the total mutual information, the individual `f_k`
contributions, and writes per-scenario tables to `dagi_results.json` for further
analysis or plotting. The stored mutual-information and `f_k` tables are emitted
as sorted lists of subsets and values to remain JSON-friendly.

