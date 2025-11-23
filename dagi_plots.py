"""Generate DAGI validation plots from saved JSON results.

This script reads `dagi_results.json` in the repository root, extracts the
recovery-wedge and outside-wedge control experiments, and produces three PNG
figures:

- mutual_information_curve.png
- moebius_contributions_stack.png
- higher_order_fraction.png

It makes no external calls and only depends on the pinned Python requirements.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt

RESULTS_PATH = Path("dagi_results.json")
PLOT_FILENAMES = {
    "mutual_information": "mutual_information_curve.png",
    "fk": "moebius_contributions_stack.png",
    "synergy_ratio": "higher_order_fraction.png",
}


class ExperimentNotFound(ValueError):
    """Raised when the expected experiment label is missing."""


def load_results(path: Path) -> Dict:
    with path.open() as f:
        return json.load(f)


def normalize_experiments(data: Dict) -> Dict[str, Dict]:
    """Return a mapping from label to experiment payload.

    Supports both the list-of-experiments representation written by
    ``dagi_validation.py`` and the mapping-style example in the plotting
    instructions (e.g., ``{"wedge": {...}, "outside_control": {...}}``).
    """

    experiments = data.get("experiments", {})
    if isinstance(experiments, dict):
        return {str(k): v for k, v in experiments.items()}

    normalized: Dict[str, Dict] = {}
    for exp in experiments:
        label = exp.get("label")
        if label:
            normalized[str(label)] = exp
    return normalized


def find_experiment(data: Dict, *labels: str) -> Dict:
    experiments = normalize_experiments(data)
    for label in labels:
        match = experiments.get(label)
        if match is not None:
            return match
    raise ExperimentNotFound(
        f"Experiment not found; looked for any of: {', '.join(labels)}"
    )


def aggregate_by_size(subsets: Iterable[Dict[str, object]]) -> Dict[int, float]:
    """Reduce subset-valued measurements to a size->max mapping.

    The stored mutual information and f-values list one entry per subset. For
    plotting by fragment size we take the maximum value for each size, which is
    sufficient for the symmetric HaPPY state (nonzero only once per size).
    """

    by_size: Dict[int, List[float]] = defaultdict(list)
    for entry in subsets:
        subset = entry.get("subset")
        value = float(entry.get("value", 0))
        if subset is None:
            continue
        by_size[len(subset)].append(value)

    return {size: max(values) if values else 0.0 for size, values in by_size.items()}


def mutual_information_by_size(exp: Dict) -> Dict[int, float]:
    mi = exp.get("mutual_information", {})
    if isinstance(mi, dict):
        return {int(k): float(v) for k, v in mi.items()}
    return aggregate_by_size(mi)


def extract_synergy_ratio(exp: Dict) -> List[Tuple[int, float]]:
    ratios = exp.get("synergy_ratio", [])
    # Ratios may be stored as list of pairs or a mapping.
    if isinstance(ratios, dict):
        return sorted(((int(k), float(v)) for k, v in ratios.items()), key=lambda t: t[0])
    return sorted(((int(size), float(val)) for size, val in ratios), key=lambda t: t[0])


def plot_mutual_information(wedge: Dict, control: Dict, output: Path) -> None:
    wedge_mi = mutual_information_by_size(wedge)
    control_mi = mutual_information_by_size(control)

    sizes = sorted(set(wedge_mi.keys()) | set(control_mi.keys()))
    wedge_vals = [wedge_mi.get(s, 0.0) for s in sizes]
    control_vals = [control_mi.get(s, 0.0) for s in sizes]

    plt.figure(figsize=(5, 4))
    plt.plot(sizes, wedge_vals, marker="o", label="Wedge fragment")
    plt.plot(sizes, control_vals, marker="s", label="Outside control")
    plt.xlabel("Fragment size")
    plt.ylabel("Mutual information I(bulk : fragment) [bits]")
    plt.title("Mutual Information vs Fragment Size")
    plt.grid(True, linewidth=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def plot_fk(wedge: Dict, output: Path) -> None:
    fk = wedge.get("fk", {})
    k_vals = sorted(int(k) for k in fk.keys())
    f_vals = [float(fk[str(k)]) for k in k_vals]

    plt.figure(figsize=(5, 4))
    plt.bar(k_vals, f_vals)
    plt.xlabel("Order k")
    plt.ylabel("Total contribution f_k [bits]")
    plt.title("Irreducible information contributions f_k (wedge)")
    plt.xticks(k_vals)
    plt.grid(axis="y", linewidth=0.3)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def plot_synergy_ratio(wedge: Dict, control: Dict, output: Path) -> None:
    wedge_ratios = extract_synergy_ratio(wedge)
    control_ratios = extract_synergy_ratio(control)

    sizes = sorted({size for size, _ in wedge_ratios} | {size for size, _ in control_ratios})
    wedge_map = {size: val for size, val in wedge_ratios}
    control_map = {size: val for size, val in control_ratios}

    wedge_vals = [wedge_map.get(s, 0.0) for s in sizes]
    control_vals = [control_map.get(s, 0.0) for s in sizes]

    plt.figure(figsize=(5, 4))
    plt.plot(sizes, wedge_vals, marker="o", label="Wedge fragment")
    plt.plot(sizes, control_vals, marker="x", label="Outside control")
    plt.xlabel("Fragment size")
    plt.ylabel(r"Synergy ratio $R_{\geq 3}$")
    plt.title("Synergy Ratio vs Fragment Size")
    plt.ylim(-0.05, 1.05)
    plt.grid(True, linewidth=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main() -> None:
    results = load_results(RESULTS_PATH)
    wedge = find_experiment(results, "recovery_wedge", "wedge")
    control = find_experiment(results, "outside_wedge_control", "outside_control")

    plot_mutual_information(wedge, control, Path(PLOT_FILENAMES["mutual_information"]))
    plot_fk(wedge, Path(PLOT_FILENAMES["fk"]))
    plot_synergy_ratio(wedge, control, Path(PLOT_FILENAMES["synergy_ratio"]))

    print("Generated plots:")
    for name in PLOT_FILENAMES.values():
        print(f" - {name}")


if __name__ == "__main__":
    main()
