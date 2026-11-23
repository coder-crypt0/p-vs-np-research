"""Reduce 3-SAT counting to a constant coefficient of a sparse determinant."""

from collections import Counter, defaultdict
from itertools import pairwise, product
from math import prod
from random import Random

import networkx as nx


def cycle_block(size, clause):
    """A signed cycle, with a unit diagonal when it represents OR."""
    if size == 1 and clause:
        return [{}]
    rows = [{i: 1} if clause else {} for i in range(size)]
    cycle_sign = (-1) ** (size if clause else size - 1)
    for i, row in enumerate(rows):
        row[(i + 1) % size] = cycle_sign if i == size - 1 else 1
    return rows


def compile_formula(clauses, *, pad_literals=False, weights=None):
    """Compile width-one through width-three clauses into a signed determinant."""
    occurrences = defaultdict(list)
    for i, clause in enumerate(clauses):
        if not 1 <= len(clause) <= 3 or 0 in clause:
            raise ValueError("expected one to three nonzero signed literals per clause")
        for j, literal in enumerate(clause):
            occurrences[abs(literal)].append((i, j, literal < 0))

    entries = []
    wire_ends = []
    owners = []

    def allocate(block):
        start = len(entries)
        entries.extend({start + col: (weight, None) for col, weight in row.items()}
                       for row in block)
        owners.extend([start] * len(block))
        return start

    def connect(left, right, selected_weight=1):
        wire = len(wire_ends)
        wire_ends.append((owners[left], owners[right]))
        entries[left][right] = (-selected_weight, wire)
        entries[right][left] = (1, wire)

    clause_starts = [allocate(cycle_block(len(clause), True)) for clause in clauses]
    for variable, ports in occurrences.items():
        absent, present = (weights or {}).get(variable, (1, 1))
        block = cycle_block(len(ports), False)
        block[0][1 % len(ports)] *= absent
        start = allocate(block)
        for offset, (clause, port, negate) in enumerate(ports):
            left = start + offset
            right = clause_starts[clause] + port
            selected_weight = present if offset == 0 else 1
            if negate or pad_literals:
                block = [{2: 1}, {2: 1}, {0: -1, 1: -1}] if negate else [
                    {1: -1}, {0: 1}, {2: 1}
                ]
                relay = allocate(block)
                connect(left, relay, selected_weight)
                left = relay + 1
                selected_weight = 1
            connect(left, right, selected_weight)
    return entries, wire_ends


def determinant(matrix):
    """Fraction-free elimination over integers."""
    size = len(matrix)
    if size == 0:
        return 1
    sign, previous = 1, 1
    for i in range(size - 1):
        pivot_row = next((r for r in range(i, size) if matrix[r][i]), None)
        if pivot_row is None:
            return 0
        if pivot_row != i:
            matrix[pivot_row], matrix[i] = matrix[i], matrix[pivot_row]
            sign = -sign
        pivot = matrix[i][i]
        for j in range(i + 1, size):
            row = matrix[j]
            for k in range(i + 1, size):
                numerator = row[k] * pivot - row[i] * matrix[i][k]
                row[k], remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("fraction-free division failed")
            row[i] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def chord_wires(wire_ends):
    """Select the wires outside a spanning forest of the gadget multigraph."""
    graph = nx.utils.UnionFind()
    chords = []
    for wire, (left, right) in enumerate(wire_ends):
        if graph[left] == graph[right]:
            chords.append(wire)
        else:
            graph.union(left, right)
    return chords


def signed_evaluation(entries, signs):
    matrix = [[0] * len(entries) for _ in entries]
    for i, row in enumerate(entries):
        for j, (weight, wire) in row.items():
            # The two Laurent exponents have identical values at signs.
            matrix[i][j] = weight * signs.get(wire, 1)
    return determinant(matrix)


def filtered_count(entries, wires):
    """Average determinant evaluations on signs for the specified wires."""
    total = sum(signed_evaluation(entries, dict(zip(wires, values)))
                for values in product((-1, 1), repeat=len(wires)))
    count, remainder = divmod(total, 2 ** len(wires))
    if remainder:
        raise ArithmeticError("the coefficient must be integral")
    return count


def bridge_filtered_count(entries, wire_ends):
    """Filter each component after deleting bridges, then compose its two boundary counts."""
    if not wire_ends:
        return filtered_count(entries, ())
    graph = nx.MultiGraph()
    for wire, (u, v) in enumerate(wire_ends):
        graph.add_edge(u, v, key=wire)
    bridges = [(u, v, next(iter(graph[u][v]))) for u, v in nx.bridges(graph)]
    graph.remove_edges_from(bridges)
    component = {vertex: i for i, vertices in enumerate(nx.connected_components(graph))
                 for vertex in vertices}
    forest = nx.Graph()
    forest.add_nodes_from(component.values())
    indices = [[] for _ in forest]
    owner = []
    for start, stop in pairwise([*sorted(component), len(entries)]):
        indices[component[start]].extend(range(start, stop))
        owner.extend([start] * (stop - start))

    wire_ports = {wire: (row, col, -weight) for row, edges in enumerate(entries)
                  for col, (weight, wire) in edges.items()
                  if wire is not None and owner[row] == wire_ends[wire][0]}
    bridge_ports = {}
    bridge_weights = {}
    for _, _, wire in bridges:
        u, v = wire_ends[wire]
        a, b = component[u], component[v]
        left, right, selected_weight = wire_ports[wire]
        bridge_ports[a, b] = left
        bridge_ports[b, a] = right
        bridge_weights[a, b] = bridge_weights[b, a] = selected_weight
        forest.add_edge(a, b)

    internal_wires = [[] for _ in forest]
    for u, _, wire in graph.edges(keys=True):
        internal_wires[component[u]].append(wire)

    def count_part(part, deleted, weights, chords):
        index = {old: new for new, old in enumerate(i for i in indices[part] if i != deleted)}
        local = []
        for old, row in index.items():
            keep, remove = weights.get(old, (1, 0))
            entries_here = {index[col]: (keep * weight, wire)
                            for col, (weight, wire) in entries[old].items() if col in index}
            diagonal = entries_here.get(row, (0, None))[0]
            if remove:
                entries_here[row] = (diagonal + remove, None)
            local.append(entries_here)
        return filtered_count(local, chords)

    root_counts = []
    for tree in nx.connected_components(forest):
        root = min(tree)
        parent = nx.dfs_predecessors(forest, root)
        messages = {}
        for part in reversed(list(nx.dfs_preorder_nodes(forest, root))):
            weights = {
                bridge_ports[part, child]: (
                    messages[child][0], messages[child][1] * bridge_weights[part, child],
                )
                for child in forest[part] if parent.get(child) == part
            }
            wires = internal_wires[part]
            chords = [wires[i] for i in chord_wires([wire_ends[wire] for wire in wires])]
            keep = count_part(part, None, weights, chords)
            if part == root:
                root_counts.append(keep)
            else:
                remove = count_part(part, bridge_ports[part, parent[part]], weights, chords)
                messages[part] = (keep, remove)
    return prod(root_counts)


def direct_count(clauses, weights=None):
    variables = sorted({abs(lit) for clause in clauses for lit in clause})
    total = 0
    for bits in product((0, 1), repeat=len(variables)):
        assignment = dict(zip(variables, bits))
        if all(any(assignment[abs(lit)] == int(lit > 0) for lit in clause)
               for clause in clauses):
            total += prod((weights or {}).get(v, (1, 1))[assignment[v]] for v in variables)
    return total


def check_formula(clauses):
    entries, ends = compile_formula(clauses)
    if max((len(row) for row in entries), default=0) > 3 or max(
        Counter(col for row in entries for col in row).values(), default=0
    ) > 3:
        raise AssertionError("reduction broke the sparsity bound")
    expected = direct_count(clauses)
    if filtered_count(entries, chord_wires(ends)) != expected:
        raise AssertionError(f"chord filtering failed on {clauses}")
    if bridge_filtered_count(entries, ends) != expected:
        raise AssertionError(f"bridge filtering failed on {clauses}")
    if len(ends) <= 7 and filtered_count(entries, range(len(ends))) != expected:
        raise AssertionError(f"full filtering failed on {clauses}")


if __name__ == "__main__":
    examples = [(), ((1, 2, 3),), ((1, 1, 1), (-1, -1, -1)),
                ((1, -1, 2), (1, -2, 1)), ((1, 2, 3), (1, 2, 3)),
                ((-3, -4, 1), (-2, -3, 1), (1, 2, 3)),
                ((-1,), (1, 2), (-2, 3)), ((1,), (-1,))]
    random = Random(91971)
    literals = [-4, -3, -2, -1, 1, 2, 3, 4]
    for _ in range(70):
        examples.append(tuple(tuple(random.choices(literals, k=3))
                              for _ in range(random.randrange(1, 5))))
    for formula in examples:
        check_formula(formula)
    print(f"Sparse determinant filter: {len(examples)} formulas checked against enumeration.")
