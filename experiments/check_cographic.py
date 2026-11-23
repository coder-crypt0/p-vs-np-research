"""Check the reduction to exact-one clauses on every triangle of a cone."""

from collections import Counter, defaultdict
from itertools import combinations, product

import networkx as nx

from check_parity import check
from graphic_sat import canonical_basis, equations


def cone_reduction(clauses):
    """Build the cone, one occurrence vertex per input occurrence, and equality diamonds."""
    occurrences = defaultdict(list)
    base_triangles = []
    for index, clause in enumerate(clauses):
        ports = tuple(range(3 * index + 1, 3 * index + 4))
        base_triangles.append(ports)
        for variable, port in zip(clause, ports):
            occurrences[variable].append(port)

    diamonds = []
    next_vertex = 3 * len(clauses) + 1
    for ports in occurrences.values():
        for left, right in zip(ports, ports[1:]):
            middle = (next_vertex, next_vertex + 1)
            next_vertex += 2
            diamonds.append((left, right, *middle))
            base_triangles.extend(((left, *middle), (right, *middle)))

    graph = nx.Graph()
    graph.add_nodes_from(range(next_vertex))
    interior_edges = {
        edge for triangle in base_triangles for edge in combinations(sorted(triangle), 2)
    }
    graph.add_edges_from(interior_edges)
    graph.add_edges_from((0, vertex) for vertex in range(1, next_vertex))
    triangles = base_triangles + [(0, u, v) for u, v in sorted(interior_edges)]
    return graph, triangles, occurrences, diamonds


def exact_port_assignments(clause_triangles, occurrences, diamonds):
    """Decode every coloring of occurrence vertices that can extend to the cone."""
    for ports in product(*clause_triangles):
        selected = set(ports)
        if any((left in selected) != (right in selected) for left, right, _, _ in diamonds):
            continue
        yield frozenset(var for var, var_ports in occurrences.items()
                        if var_ports[0] in selected), selected


def check_formula(clauses):
    graph, triangles, occurrences, diamonds = cone_reduction(clauses)
    actual = set()
    for clique in nx.enumerate_all_cliques(graph):
        if len(clique) > 3:
            break
        if len(clique) == 3:
            actual.add(tuple(sorted(clique)))
    check(actual == {tuple(sorted(t)) for t in triangles}, "unexpected cone triangle")
    edge_number = {tuple(sorted(e)): i for i, e in enumerate(graph.edges)}
    exact_clauses = [tuple(edge_number[e] for e in combinations(sorted(t), 2)) for t in triangles]
    basis = canonical_basis(row >> 1 for row in equations(exact_clauses))
    check(len(basis) == len(edge_number) - len(graph) + 1, "cycle rank differs")
    weights = Counter(i for clause in exact_clauses for i in clause)
    check(max(weights.values(), default=0) <= 6, "occurrence bound under cubic input")

    expected = set()
    for values in product((0, 1), repeat=len(occurrences)):
        chosen = frozenset(v for v, bit in zip(occurrences, values) if bit)
        if all(len(chosen.intersection(row)) == 1 for row in clauses):
            expected.add(chosen)

    decoded = set()
    for original, selected in exact_port_assignments(
        triangles[:len(clauses)], occurrences, diamonds
    ):
        decoded.add(original)
        same_as_apex = selected | {0}
        for left, _, interior, _ in diamonds:
            if left not in selected:
                same_as_apex.add(interior)
        cut = {number for (u, v), number in edge_number.items()
               if (u in same_as_apex) != (v in same_as_apex)}
        check(all(sum(e in cut for e in row) == 2 for row in exact_clauses), "bad cut")
        check(sum(weights[e] for e in cut) == 2 * len(triangles), "wrong cut weight")
    check(decoded == expected, f"assignment projection failed for {clauses}")
    return bool(expected)


if __name__ == "__main__":
    triples = tuple(combinations(range(5), 3))
    sat = unsat = 0
    for size in range(len(triples) + 1):
        for clauses in combinations(triples, size):
            if max(Counter(v for row in clauses for v in row).values(), default=0) > 3:
                continue
            incidence = nx.Graph()
            incidence.add_edges_from((("clause", i), ("variable", v))
                                     for i, row in enumerate(clauses) for v in row)
            if nx.check_planarity(incidence)[0]:
                cone, _, _, _ = cone_reduction(clauses)
                cone.remove_node(0)
                check(nx.check_planarity(cone)[0], "planar incidence gave nonplanar base")
            if check_formula(clauses):
                sat += 1
            else:
                unsat += 1
    check(sat > 0 and unsat > 0, "missing either outcome")
    print(f"Cone reduction: {sat} satisfiable and {unsat} unsatisfiable inputs checked.")
