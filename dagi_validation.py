"""DAGI validation experiment using the stabilizer graph code output.

This script reconstructs the 16-qubit graph state from ``Main_16`` using the
Hadamard choices from the manuscript (gates 2, 6, 10, and 14). It then computes
subregion entropies, mutual informations with a selected bulk qubit, performs a
Möbius inversion to extract ``f_k`` synergy terms, and reports the running
synergy ratio R_{>=3}.

The calculation relies on the graph-state entropy identity
S(A) = rank_2(Γ_{A,\bar{A}}), where Γ is the adjacency matrix and rank is taken
modulo 2. All outputs are written to ``dagi_results.json`` for reuse.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Mapping, Sequence, Tuple

import networkx as nx
import numpy as np

import Functions as f
import LC_explore


@dataclass
class GraphStateData:
    adjacency: np.ndarray
    bulk_nodes: List[int]
    boundary_nodes: List[int]
    bulk_target: int
    fragment_set: List[int]


def build_graph_state() -> GraphStateData:
    """Construct the 16-qubit graph state adjacency used for DAGI validation."""
    x_matrix = np.array(
        [
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ]
    )
    z_matrix = np.array(
        [
            [0, 1, 0, 0, 1, 1],
            [1, 0, 1, 0, 0, 1],
            [0, 1, 0, 1, 0, 1],
            [0, 0, 1, 0, 1, 1],
            [1, 0, 0, 1, 0, 1],
            [1, 1, 1, 1, 1, 0],
        ]
    )
    position = [[4, 6], [10, 12], [16, 18], [22, 0]]

    gx, gz, vph, _ = f.StateAM(x_matrix, z_matrix, position, 200, figBM=False, layout=False)

    # Prepare the specific graph from the manuscript by applying Hadamards on
    # qubits 2, 6, 10, and 14 before stabilizer reduction.
    gates = [2, 6, 10, 14]
    gx_gates, gz_gates, vph_gates = LC_explore.HalamardGates(gx.copy(), gz.copy(), vph.copy(), gates)
    gx_gates, gz_gates, vph_gates = f.Triangular(gx_gates, gz_gates, vph_gates)
    gx_gates, gz_gates, vph_gates = f.CleanMatrix(gx_gates, gz_gates, vph_gates)

    graph = f.GraphDraw(gx_gates, gz_gates)
    adjacency = np.array(nx.to_numpy_array(graph, dtype=int))

    bulk_nodes = [3, 7, 11, 15]
    boundary_nodes = [node for node in range(adjacency.shape[0]) if node not in bulk_nodes]

    bulk_target = 15
    fragment_set = sorted(graph.neighbors(bulk_target))

    return GraphStateData(
        adjacency=adjacency,
        bulk_nodes=bulk_nodes,
        boundary_nodes=boundary_nodes,
        bulk_target=bulk_target,
        fragment_set=fragment_set,
    )


def gf2_rank(matrix: np.ndarray) -> int:
    """Return the rank over GF(2) using Gaussian elimination."""
    m = matrix.copy() % 2
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        pivot_row = None
        for r in range(rank, rows):
            if m[r, col] == 1:
                pivot_row = r
                break
        if pivot_row is None:
            continue
        if pivot_row != rank:
            m[[rank, pivot_row]] = m[[pivot_row, rank]]
        for r in range(rows):
            if r != rank and m[r, col] == 1:
                m[r] = (m[r] + m[rank]) % 2
        rank += 1
        if rank == rows:
            break
    return rank


def entropy_from_adjacency(adjacency: np.ndarray, subset: Sequence[int]) -> int:
    """Compute S(subset) in bits for a graph state using the rank formula."""
    if not subset:
        return 0
    n = adjacency.shape[0]
    keep = set(subset)
    complement = [i for i in range(n) if i not in keep]
    if not complement:
        return 0
    submatrix = adjacency[np.ix_(subset, complement)] % 2
    return gf2_rank(submatrix)


def mutual_information(adjacency: np.ndarray, a: Sequence[int], b: Sequence[int]) -> float:
    a_entropy = entropy_from_adjacency(adjacency, a)
    b_entropy = entropy_from_adjacency(adjacency, b)
    joint_entropy = entropy_from_adjacency(adjacency, list(set(a) | set(b)))
    return a_entropy + b_entropy - joint_entropy


def mutual_information_table(adjacency: np.ndarray, bulk: int, fragments: List[int]) -> Dict[Tuple[int, ...], float]:
    result: Dict[Tuple[int, ...], float] = {}
    bulk_set = [bulk]
    for size in range(1, len(fragments) + 1):
        for combo in combinations(fragments, size):
            result[combo] = mutual_information(adjacency, bulk_set, combo)
    return result


def serialize_table(table: Mapping[Tuple[int, ...], float]) -> List[Dict[str, object]]:
    """Produce a JSON-friendly, length-sorted table of subset values."""
    return [
        {"subset": list(subset), "value": value}
        for subset, value in sorted(table.items(), key=lambda item: (len(item[0]), item[0]))
    ]


def mobius_inversion(mi_table: Mapping[Tuple[int, ...], float]) -> Dict[Tuple[int, ...], float]:
    """Recover Möbius f-values given I(bulk : subset)."""
    f_values: Dict[Tuple[int, ...], float] = {}
    for subset in sorted(mi_table.keys(), key=len):
        total = mi_table[subset]
        for smaller, value in f_values.items():
            if set(smaller).issubset(subset):
                total -= value
        f_values[subset] = total
    return f_values


def aggregate_fk(f_values: Mapping[Tuple[int, ...], float]) -> Dict[int, float]:
    fk: Dict[int, float] = {}
    for subset, value in f_values.items():
        fk[len(subset)] = fk.get(len(subset), 0.0) + value
    return fk


def running_synergy_ratio(fk: Mapping[int, float], total_info: float) -> List[Tuple[int, float]]:
    ratios = []
    cumulative = 0.0
    for k in sorted(fk):
        if k >= 3:
            cumulative += fk[k]
        ratios.append((k, 0.0 if total_info == 0 else cumulative / total_info))
    return ratios


def analyze_fragments(
    adjacency: np.ndarray, bulk: int, fragments: List[int], label: str
) -> Dict[str, object]:
    """Compute mutual information, Möbius terms, and synergy ratios for a fragment set."""

    ordered_fragments = sorted(fragments)
    mi_table = mutual_information_table(adjacency, bulk, ordered_fragments)
    f_values = mobius_inversion(mi_table)
    fk = aggregate_fk(f_values)
    total_info = mi_table[tuple(ordered_fragments)]
    ratios = running_synergy_ratio(fk, total_info)

    return {
        "label": label,
        "bulk_target": bulk,
        "fragment_set": ordered_fragments,
        "mutual_information": serialize_table(mi_table),
        "f_values": serialize_table(f_values),
        "fk": fk,
        "synergy_ratio": ratios,
        "total_information": total_info,
    }


def run_experiment() -> Dict[str, object]:
    graph_data = build_graph_state()

    recovery_wedge = analyze_fragments(
        graph_data.adjacency,
        graph_data.bulk_target,
        graph_data.fragment_set,
        label="recovery_wedge",
    )

    outside_wedge = analyze_fragments(
        graph_data.adjacency,
        graph_data.bulk_target,
        [3, 4, 5, 6, 7],
        label="outside_wedge_control",
    )

    fragment_growth_sets = [
        ([0, 1, 2], "growth_3_qubits_012"),
        ([0, 2, 14], "growth_3_qubits_0_2_14"),
        ([0, 1, 2, 12], "growth_4_qubits_0_1_2_12"),
    ]

    fragment_growth = [
        analyze_fragments(graph_data.adjacency, graph_data.bulk_target, subset, label)
        for subset, label in fragment_growth_sets
    ]

    experiments = [recovery_wedge, outside_wedge, *fragment_growth]

    return {
        "graph": {
            "bulk_target": graph_data.bulk_target,
            "bulk_nodes": graph_data.bulk_nodes,
            "boundary_nodes": graph_data.boundary_nodes,
        },
        "experiments": experiments,
    }


if __name__ == "__main__":
    results = run_experiment()
    print("Graph bulk target:", results["graph"]["bulk_target"])
    for experiment in results["experiments"]:
        print("---", experiment["label"], "---")
        print("Fragments:", experiment["fragment_set"])
        print("Total I(bulk:fragments):", experiment["total_information"])
        print("f_k contributions:", experiment["fk"])
        print("Running R_{>=3}:", experiment["synergy_ratio"])
    with open("dagi_results.json", "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")
    print("Saved dagi_results.json")
