"""Decide positive one-in-three SAT when given an equivalent graphic parity system.

A clause is a triple of distinct nonnegative variable indices. An augmented
binary row stores its right-hand side in bit 0 and variable i in bit i+1.
"""

from collections import Counter
from itertools import combinations

import networkx as nx


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def canonical_basis(rows):
    """Reduced row basis of a GF(2) space, using the highest set bit as pivot."""
    pivots = {}
    for row in rows:
        for pivot in sorted(pivots, reverse=True):
            if row & (1 << pivot):
                row ^= pivots[pivot]
        if row:
            pivot = row.bit_length() - 1
            for index, other in pivots.items():
                if other & (1 << pivot):
                    pivots[index] = other ^ row
            pivots[pivot] = row
    return tuple(pivots[p] for p in sorted(pivots, reverse=True))


def equations(clauses):
    result = []
    for clause in clauses:
        if len(clause) != 3 or len(set(clause)) != 3 or any(i < 0 for i in clause):
            raise ValueError("expected three distinct nonnegative variable indices")
        result.append(1 | sum(1 << (i + 1) for i in clause))
    return tuple(result)


def columns(rows):
    used = 0
    for row in rows:
        used |= row >> 1
    return {
        i: sum(1 << j for j, row in enumerate(rows) if row & (1 << (i + 1)))
        for i in bits(used)
    }


def graphic_witness(clauses, rows):
    """Return a set of true variables, or None for UNSAT. Reject an invalid basis."""
    clauses = tuple(clauses)
    original = equations(clauses)
    rows = tuple(rows)
    if any(row < 0 for row in rows) or canonical_basis(rows) != canonical_basis(original):
        raise ValueError("augmented equations must have the same row space")
    incidence = columns(rows)
    if any(col.bit_count() > 2 for col in incidence.values()):
        raise ValueError("basis is not graphic")

    weights = Counter(i for clause in clauses for i in clause)
    root = len(rows)
    graph = nx.Graph()
    graph.add_nodes_from(range(root + 1))
    for i, col in incidence.items():
        ends = list(bits(col))
        u, v = ends if len(ends) == 2 else (ends[0], root)
        # Parallel edges have the same syndrome. An optimum needs only the cheapest one.
        if not graph.has_edge(u, v) or weights[i] < graph[u][v]["weight"]:
            graph.add_edge(u, v, weight=weights[i], variable=i)

    terminals = {j for j, row in enumerate(rows) if row & 1}
    if len(terminals) % 2:
        terminals.add(root)
    selected = set()
    for component in nx.connected_components(graph):
        ends = sorted(terminals & component)
        if len(ends) % 2:
            return None
        metric = nx.Graph()
        for u, v in combinations(ends, 2):
            cost = nx.shortest_path_length(graph, u, v, weight="weight")
            metric.add_edge(u, v, weight=cost)
        matching = nx.min_weight_matching(metric, weight="weight")
        if 2 * len(matching) != len(ends):
            raise RuntimeError("terminal pairing is incomplete")
        for u, v in matching:
            path = nx.shortest_path(graph, u, v, weight="weight")
            selected.symmetric_difference_update(
                graph[a][b]["variable"] for a, b in zip(path, path[1:])
            )

    if sum(weights[i] for i in selected) != len(clauses):
        return None
    if any(sum(i in selected for i in clause) != 1 for clause in clauses):
        raise RuntimeError("matching reconstruction produced an invalid witness")
    return selected


def small_graphic_basis(rows):
    """Exhaustive diagnostic for coefficient rank <= 4, with an augmented RHS."""
    rows = canonical_basis(rows)
    if 1 in rows:
        raise ValueError("parity equations are inconsistent")
    rank = len(rows)
    if rank > 4:
        raise ValueError("exhaustive basis search is limited to rank four")
    sums = [0]
    for row in rows:
        sums += [s ^ row for s in sums]
    for basis in combinations(sums[1:], rank):
        if canonical_basis(basis) == rows and all(
            col.bit_count() <= 2 for col in columns(basis).values()
        ):
            return basis
    return None
