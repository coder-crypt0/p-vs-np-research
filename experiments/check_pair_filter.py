"""Compose exact determinant counts across separators of two variables."""

from collections import Counter
from itertools import combinations, cycle, product
from random import Random

import networkx as nx

from check_determinant_filter import (
    chord_wires, compile_formula, determinant, direct_count, filtered_count, pair_block,
)


ASSIGNMENTS = tuple(product((0, 1), repeat=2))


def primal_graph(scopes):
    graph = nx.Graph()
    for scope in scopes:
        graph.add_nodes_from(scope)
        graph.add_edges_from(combinations(set(scope), 2))
    return graph


def incidence_rank(scopes):
    variables = set()
    forest = nx.utils.UnionFind()
    for scope in scopes:
        variables.update(scope)
        forest.union(*scope)
    components = sum(1 for _ in forest.to_sets())
    wires = sum(map(len, scopes))
    return wires - len(scopes) - len(variables) + components


def find_split(scopes, max_side_rank):
    graph = primal_graph(scopes)
    for boundary in combinations(sorted(graph), 2):
        interior = graph.subgraph(set(graph).difference(boundary))
        components = list(nx.connected_components(interior))
        if len(components) < 2:
            continue
        for side in components:
            selected = [i for i, scope in enumerate(scopes) if not side.isdisjoint(scope)]
            if incidence_rank([scopes[i] for i in selected]) <= max_side_rank:
                return boundary, set(selected)
    return None


def weighted_count(clauses, pairs, pins=None):
    weights = {v: (1 - bit, bit) for v, bit in (pins or {}).items()}
    entries, ends = compile_formula(clauses, weights=weights, pair_weights=pairs)
    return filtered_count(entries, chord_wires(ends))


def pair_filtered_count(clauses, max_side_rank=2):
    clauses = list(clauses)
    if any(not 1 <= len(c) <= 3 or 0 in c for c in clauses):
        raise ValueError("expected one to three nonzero signed literals per clause")
    pairs = {}
    while True:
        scopes = [tuple(map(abs, clause)) for clause in clauses] + list(pairs)
        split = find_split(scopes, max_side_rank)
        if split is None:
            return weighted_count(clauses, pairs)
        boundary, selected = split
        local_clauses = [clause for i, clause in enumerate(clauses) if i in selected]
        local_pairs = {pair: values for i, (pair, values) in enumerate(pairs.items(), len(clauses))
                       if i in selected}
        table = [weighted_count(local_clauses, local_pairs, dict(zip(boundary, bits)))
                 for bits in ASSIGNMENTS]
        if not any(table):
            return 0
        clauses = [clause for i, clause in enumerate(clauses) if i not in selected]
        for pair in local_pairs:
            del pairs[pair]
        pairs[boundary] = tuple(a * b for a, b in zip(
            pairs.get(boundary, (1, 1, 1, 1)), table,
        ))


def check_pair_blocks():
    for values in product((-2, 0, 1, 3), repeat=4):
        rows = pair_block(values)
        for bits, expected in zip(ASSIGNMENTS, values):
            kept = [i for i, bit in enumerate((*bits, 0)) if bit == 0]
            matrix = [[rows[i].get(j, 0) for j in kept] for i in kept]
            if determinant(matrix) != expected:
                raise AssertionError(f"wrong principal minor for {values}, {bits}")
        if weighted_count((), {(1, 2): values}) != sum(values):
            raise AssertionError(f"wrong weighted determinant for {values}")

    pairs = {(1, 2): (1, 2, 3, 4), (2, 3): (4, 3, 2, 1), (1, 3): (0, 1, 0, 2)}
    clauses = ((1, -2, 3),)
    entries, ends = compile_formula(clauses, pair_weights=pairs)
    if max(map(len, entries)) > 3 or max(
        Counter(col for row in entries for col in row).values()
    ) > 3:
        raise AssertionError("binary factors broke the sparsity bound")
    scopes = [(1, 2, 3), *pairs]
    if len(chord_wires(ends)) != incidence_rank(scopes):
        raise AssertionError("binary factors changed the incidence cycle rank")
    expected = sum(pairs[1, 2][2 * a + b] * pairs[2, 3][2 * b + c]
                   * pairs[1, 3][2 * a + c]
                   for a, b, c in product((0, 1), repeat=3) if a or not b or c)
    if filtered_count(entries, chord_wires(ends)) != expected:
        raise AssertionError("wrong count for overlapping binary factors")


def check_small():
    rng = Random(84992)
    literals = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    examples = [(), ((1,), (-1,)), ((1, -1, 2), (2, 3, -4), (3, 4, 5)),
                ((1, 2, 3), (-1, -2, 3), (2, 3, 4), (-3, 4, 5), (-4, 5, 6))]
    for _ in range(64):
        examples.append(tuple(tuple(rng.choices(literals, k=rng.randrange(1, 4)))
                              for _ in range(rng.randrange(1, 7))))
    for clauses in examples:
        if pair_filtered_count(clauses) != direct_count(clauses):
            raise AssertionError(f"wrong separator count for {clauses}")
    return len(examples)


def check_many_cycles():
    lobes = 48
    for conflicting in (False, True):
        signs = ASSIGNMENTS if conflicting else ((1, 1),)
        arms = [(x, 2 * a - 1, 2 * (2 * b - 1))
                for x, (a, b) in zip(range(3, lobes + 3), cycle(signs))]
        clauses = [clause for x, y, z in arms for clause in ((x, y, z), (-x, y, z))]
        graph = primal_graph([tuple(map(abs, clause)) for clause in clauses])
        _, ends = compile_formula(clauses)
        if list(nx.articulation_points(graph)) or list(nx.bridges(nx.MultiGraph(ends))):
            raise AssertionError("the example should have no articulation variable or wire bridge")
        if len(chord_wires(ends)) != 3 * lobes - 1:
            raise AssertionError("wrong cycle rank for the separator example")
        if pair_filtered_count(clauses) != (0 if conflicting else 3 * 2 ** lobes):
            raise AssertionError("wrong count for the separator example")


def check_three_port_example():
    clauses = ((1, 2, 3), (-5, 1), (-5, 2), (5, -1, -2), (-5, -3, -4))
    table = [direct_count(clauses, dict(zip((1, 2, 3), ((1 - b, b) for b in bits))))
             for bits in product((0, 1), repeat=3)]
    if table != [0, 2, 2, 2, 2, 2, 2, 1]:
        raise AssertionError("wrong conditional table for the three-port obstruction")


if __name__ == "__main__":
    check_pair_blocks()
    checked = check_small()
    check_many_cycles()
    check_three_port_example()
    print(f"Two-variable separators matched direct counting on {checked} formulas.")
    print("Two 143-cycle examples with no articulation variable were counted exactly.")
    print("Verified the SAT conditional counts for the three-port obstruction.")
