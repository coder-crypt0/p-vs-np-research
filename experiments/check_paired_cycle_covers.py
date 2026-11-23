"""Check SAT as a sparse cycle-cover problem with paired crossing edges."""

from collections import Counter
from itertools import product
from random import Random

import networkx as nx

from check_determinant_filter import compile_formula, direct_count


def support_clause_block(width):
    """Use auxiliary vertices to cover at most width minus one surviving ports."""
    if width == 1:
        return [{}]
    if width == 2:
        return [{2: 1}, {2: 1}, {0: 1, 1: 1, 2: 1}]
    return [{3: 1}, {3: 1, 4: 1}, {4: 1}, {0: 1, 1: 1, 3: 1}, {1: 1, 2: 1, 4: 1}]


def perfect_cover(vertices, edges):
    left = {(0, i) for i in vertices}
    graph = nx.Graph()
    graph.add_nodes_from(left)
    graph.add_nodes_from((1, j) for j in vertices)
    graph.add_edges_from(((0, i), (1, j)) for i, j in edges
                         if i in vertices and j in vertices)
    matching = nx.algorithms.bipartite.maximum_matching(graph, top_nodes=left)
    if len(matching) != 2 * len(left):
        return None
    return {(i, matching[0, i][1]) for _, i in left}


def wire_edges(entries):
    internal, wires = set(), {}
    for i, row in enumerate(entries):
        for j, (_, wire) in row.items():
            if wire is None:
                internal.add((i, j))
            else:
                wires.setdefault(wire, set()).add((i, j))
    return internal, wires


def paired_cover(vertices, internal, wires, selected):
    crossing = {edge for wire in selected for edge in wires[wire]}
    remaining = vertices.difference(i for i, _ in crossing)
    local = perfect_cover(remaining, internal)
    return None if local is None else local | crossing


def is_paired(cover, wires):
    return all(len(cover & directions) != 1 for directions in wires.values())


def decompose(vertices, multiplicities, degree):
    """Decompose an integral regular bipartite multigraph into perfect covers."""
    remaining = +Counter(multiplicities)
    covers = []
    for _ in range(degree):
        cover = perfect_cover(vertices, remaining)
        if cover is None:
            raise AssertionError("a regular bipartite multigraph should have a perfect matching")
        remaining.subtract(cover)
        remaining = +remaining
        covers.append(cover)
    if remaining:
        raise AssertionError("the covers must exhaust the multigraph")
    return covers


def check_small_formulas():
    rng = Random(95612)
    clauses = [(), ((1,), (-1,)), ((1, -1, 2),), ((1, 1), (-1, -1))]
    literals = (-2, -1, 1, 2)
    for _ in range(32):
        clauses.append(tuple(tuple(rng.choices(literals, k=rng.randrange(1, 4)))
                             for _ in range(2)))
    for formula in clauses:
        entries, ends = compile_formula(formula, clause_block=support_clause_block)
        vertices = set(range(len(entries)))
        internal, wires = wire_edges(entries)
        exists = any(
            paired_cover(vertices, internal, wires,
                         [wire for wire, used in enumerate(bits) if used]) is not None
            for bits in product((0, 1), repeat=len(ends))
        )
        if exists != bool(direct_count(formula)):
            raise AssertionError(f"wrong support answer for {formula}")
    return len(clauses)


def check_balanced_unsatisfiable_instance():
    clauses = tuple(product((1, -1), (2, -2), (3, -3)))
    entries, _ = compile_formula(clauses, clause_block=support_clause_block)
    vertices = set(range(len(entries)))
    _, wires = wire_edges(entries)
    multiplicities = {}
    for i, row in enumerate(entries):
        for j, (_, wire) in row.items():
            clause_edge = i < 5 * len(clauses) and wire is None
            quarter = clause_edge and (i == j or i % 5 == 1 or j % 5 == 1)
            multiplicities[i, j] = 1 if quarter else 2
    if direct_count(clauses):
        raise AssertionError("the eight sign patterns should forbid every assignment")
    row_degrees = Counter()
    column_degrees = Counter()
    for (i, j), multiplicity in multiplicities.items():
        row_degrees[i] += multiplicity
        column_degrees[j] += multiplicity
    if row_degrees != column_degrees or row_degrees != {i: 4 for i in vertices}:
        raise AssertionError("the integer multiplicities must have row and column sum four")
    covers = decompose(vertices, multiplicities, 4)
    if any(is_paired(cover, wires) for cover in covers):
        raise AssertionError("an unpaired cover should occur in every term of this witness")
    if Counter(edge for cover in covers for edge in cover) != multiplicities:
        raise AssertionError("the average must have the given strictly positive marginals")
    if any(sum(edge in cover for cover in covers) != 2
           for directions in wires.values() for edge in directions):
        raise AssertionError("each direction of a wire must occur exactly twice")


def check_linear_cost_obstruction():
    clauses = ((1, 2, 3), (1, 2, 4))
    entries, _ = compile_formula(clauses, clause_block=support_clause_block)
    vertices = set(range(len(entries)))
    internal, wires = wire_edges(entries)
    edges = internal.union(*wires.values())
    wire_variables = [v for v, degree in Counter(lit for clause in clauses for lit in clause).items()
                      for _ in range(degree)]
    good = []
    for bits in ((1, 1, 1, 1), (0, 0, 1, 1), (1, 0, 0, 0)):
        selected = [wire for wire, v in enumerate(wire_variables) if bits[v - 1]]
        cover = paired_cover(vertices, internal, wires, selected)
        if cover is None or not is_paired(cover, wires):
            raise AssertionError("the three satisfying assignments must have paired covers")
        good.append(cover)
    if set.union(*good) != edges:
        raise AssertionError("the three paired covers must use every supported edge")

    # Use one directed crossing cycle through variables 1 and 2 and both clauses.
    cycle = (10, 0, 3, 1, 12, 13, 6, 8, 5, 11)
    bad = set(zip(cycle, cycle[1:] + cycle[:1]))
    bad |= perfect_cover(vertices.difference(cycle), edges)
    if not bad <= edges or is_paired(bad, wires):
        raise AssertionError("the cover must be supported and use at least one wire singly")
    if ({i for i, _ in bad} != vertices or {j for _, j in bad} != vertices
            or len(bad) != len(vertices)):
        raise AssertionError("the unpaired edges must form a cycle cover")
    aggregate = Counter(edge for cover in good for edge in cover)
    residue = aggregate.copy()
    residue.subtract(bad)
    other_covers = decompose(vertices, residue, 2)
    if Counter(edge for cover in [bad, *other_covers] for edge in cover) != aggregate:
        raise AssertionError("the supported edge multiplicities must agree exactly")


if __name__ == "__main__":
    checked = check_small_formulas()
    check_balanced_unsatisfiable_instance()
    check_linear_cost_obstruction()
    print(f"Paired cycle-cover existence matched direct SAT on {checked} formulas.")
    print("Decomposed an unsatisfiable formula's balanced relaxation into four unpaired covers.")
    print("Verified the three-cover certificate against a linear edge-cost filter.")
