# DAGI Validation Run Report

This report captures the latest execution of `dagi_validation.py`, which rebuilds the 16-qubit stabilizer graph state and runs three information-theoretic probes around bulk qubit 15: the recovery wedge, a disjoint control set outside the wedge, and fragment-growth controls that remain below the recovery threshold.

## Execution

Command (re-run in the container on 2025-11-22 UTC):
```bash
python -m pip install -r requirements.txt
python dagi_validation.py
```

Environment notes:
- Python 3.12
- Dependencies pinned in `requirements.txt` to make the calculation reproducible (`networkx`, `numpy`, `matplotlib`, `pynauty`, `sympy`).

## Results by scenario

- **Recovery wedge** (`[0, 1, 2, 12, 14]`): mutual information **2 bits**, with sole contribution from **`f_5 = 2.0`** and synergy ratio **R_{>=3} = 1.0** at order 5.
- **Outside-wedge control** (`[3, 4, 5, 6, 7]`): mutual information **0**, and all `f_k` values **0**, demonstrating no leakage of bulk information beyond the wedge.
- **Fragment-growth controls**:
  - 3-qubit subsets `[0,1,2]` and `[0,2,14]`: mutual information **0**, all `f_k` **0**.
  - 4-qubit subset `[0,1,2,12]`: mutual information **0**, all `f_k` **0**.
  - These confirm bulk information appears only once the full five-qubit wedge is included.

Stored mutual-information and Möbius tables for each scenario live in `dagi_results.json`, with subsets serialized as sorted `{subset, value}` entries for easy parsing.

## Plot regeneration

Running `python dagi_plots.py` converts the saved JSON tables into three figures stored in the repository root:

- `mutual_information_curve.png` (mutual information versus fragment size, including the outside-wedge control)
- `moebius_contributions_stack.png` (`f_k` contributions for the recovery wedge)
- `higher_order_fraction.png` (synergy ratio trajectory)

The current branch already includes the freshly generated PNGs alongside `dagi_results.json`, so the visualizations are available without re-running the simulation or plotting steps.
