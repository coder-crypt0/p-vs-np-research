"""Condition on articulation variables and evaluate each block by its own determinant."""

from collections import defaultdict
from itertools import combinations
from math import prod
from random import Random

import networkx as nx

from check_determinant_filter import (
    bridge_filtered_count, chord_wires, compile_formula, direct_count, filtered_count,
)


def block_filtered_count(clauses):
    primal = nx.Graph()
    scopes = [set(map(abs, clause)) for clause in clauses]
    unary = defaultdict(lambda: [1, 1])
    for clause, scope in zip(clauses, scopes):
        if not 1 <= len(clause) <= 3 or 0 in clause:
            raise ValueError("expected one to three nonzero signed literals per clause")
        primal.add_nodes_from(scope)
        primal.add_edges_from(combinations(scope, 2))
        if len(scope) == 1:
            (variable,) = scope
            for bit in (0, 1):
                unary[variable][bit] *= any(bit == (literal > 0) for literal in clause)

    forest = nx.Graph()
    forest.add_nodes_from(primal)
    edge_block = {}
    for number, edges in enumerate(nx.biconnected_component_edges(primal)):
        node = ("block", number)
        members = {v for edge in edges for v in edge}
        forest.add_edges_from((v, node) for v in sorted(members))
        edge_block.update((frozenset(edge), node) for edge in edges)

    block_clauses = defaultdict(list)
    for clause, scope in zip(clauses, scopes):
        if len(scope) >= 2:
            a, b = sorted(scope)[:2]
            block_clauses[edge_block[frozenset((a, b))]].append(clause)

    counts = []
    for tree in nx.connected_components(forest):
        root = min(node for node in tree if isinstance(node, int))
        parents = nx.dfs_predecessors(forest, root)
        messages = {}
        for node in reversed(list(nx.dfs_preorder_nodes(forest, root))):
            children = {child: messages[child] for child in forest[node]
                        if parents.get(child) == node}
            if isinstance(node, int):
                messages[node] = tuple(unary[node][bit] * prod(
                    values[bit] for values in children.values()
                ) for bit in (0, 1))
            else:
                boundary = parents[node]
                conditional = []
                for pin in ((1, 0), (0, 1)):
                    entries, ends = compile_formula(
                        block_clauses[node], weights=children | {boundary: pin},
                    )
                    conditional.append(filtered_count(entries, chord_wires(ends)))
                messages[node] = tuple(conditional)
        counts.append(sum(messages[root]))
    return prod(counts)


def check_small():
    rng = Random(59073)
    literals = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    examples = [(), ((1,), (-1,)), ((1, 1, -1), (2,), (-3, 4), (4, -4, 5))]
    for _ in range(80):
        examples.append(tuple(tuple(rng.choices(literals, k=rng.randrange(1, 4)))
                              for _ in range(rng.randrange(1, 8))))
    for clauses in examples:
        if block_filtered_count(clauses) != direct_count(clauses):
            raise AssertionError(f"wrong block count for {clauses}")
        weights = {i: (rng.randrange(4), rng.randrange(4)) for i in range(1, 6)}
        entries, ends = compile_formula(clauses, weights=weights)
        expected = direct_count(clauses, weights)
        if filtered_count(entries, chord_wires(ends)) != expected:
            raise AssertionError(f"wrong weighted chord count for {clauses}")
        if bridge_filtered_count(entries, ends) != expected:
            raise AssertionError(f"wrong weighted bridge count for {clauses}")
    return len(examples)


def check_many_cycles():
    blocks = 48
    hub = blocks + 1
    clauses = [clause for x in range(1, hub) for clause in ((x, hub), (-x, hub))]
    for conflicting in (False, True):
        if conflicting:
            clauses[-2:] = [(blocks, -hub), (-blocks, -hub)]
        _, ends = compile_formula(clauses)
        graph = nx.MultiGraph(ends)
        if list(nx.bridges(graph)) or len(chord_wires(ends)) != blocks:
            raise AssertionError("the example should be bridgeless with 48 cycles")
        count = block_filtered_count(clauses)
        if count != (0 if conflicting else 1 << blocks):
            raise AssertionError("wrong count on the bridgeless example")


if __name__ == "__main__":
    checked = check_small()
    check_many_cycles()
    print(f"Block and weighted determinant methods: {checked} formulas checked.")
    print("Two connected bridgeless examples with 48 cycles counted exactly.")
