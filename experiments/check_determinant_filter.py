"""Reduce 3-SAT counting to a constant coefficient of a sparse determinant."""

from collections import Counter, defaultdict
from itertools import product
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


def compile_formula(clauses, *, pad_literals=False):
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

    def connect(left, right):
        wire = len(wire_ends)
        wire_ends.append((owners[left], owners[right]))
        entries[left][right] = (-1, wire)
        entries[right][left] = (1, wire)

    clause_starts = [allocate(cycle_block(len(clause), True)) for clause in clauses]
    for ports in occurrences.values():
        start = allocate(cycle_block(len(ports), False))
        for offset, (clause, port, negate) in enumerate(ports):
            left = start + offset
            right = clause_starts[clause] + port
            if negate or pad_literals:
                block = [{2: 1}, {2: 1}, {0: -1, 1: -1}] if negate else [
                    {1: -1}, {0: 1}, {2: 1}
                ]
                relay = allocate(block)
                connect(left, relay)
                left = relay + 1
            connect(left, right)
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


def direct_count(clauses):
    variables = sorted({abs(lit) for clause in clauses for lit in clause})
    total = 0
    for bits in product((0, 1), repeat=len(variables)):
        assignment = dict(zip(variables, bits))
        total += all(any(assignment[abs(lit)] == int(lit > 0) for lit in clause)
                     for clause in clauses)
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
