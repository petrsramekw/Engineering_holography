# DAGI Validation Run Report

This report captures the latest execution of `dagi_validation.py`, which rebuilds the 16-qubit stabilizer graph state, computes mutual information between bulk qubit 15 and its five-neighbour recovery set, and evaluates the higher-order synergy ratio.

## Execution

Command (re-run in the container on 2025-11-22 UTC):
```bash
python -m pip install -r requirements.txt
python dagi_validation.py
```

Environment notes:
- Python 3.12
- Dependencies pinned in `requirements.txt` to make the calculation reproducible (`networkx`, `numpy`, `matplotlib`, `pynauty`, `sympy`).

## Key outputs

- Bulk qubit: **15**
- Fragment set: **[0, 1, 2, 12, 14]**
- Total mutual information I(bulk : fragments): **2 bits**
- Möbius aggregated contributions `f_k`: **{1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 2.0}**
- Running synergy ratio R_{>=3}: **[(1, 0.0), (2, 0.0), (3, 0.0), (4, 0.0), (5, 1.0)]**
- Stored tables are serialized as JSON lists of `{subset, value}` entries sorted by subset size to simplify downstream parsing.

All intermediate tables and values are stored in `dagi_results.json` for reuse.
