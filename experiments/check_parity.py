"""Reproduce the finite checks in notes/parity-matching.md."""

from itertools import combinations
from math import comb

from graphic_sat import (
    bits, canonical_basis, columns, equations, graphic_witness, small_graphic_basis,
)


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def exact_and_parity_costs(clauses, variables):
    masks = [row >> 1 for row in equations(clauses)]
    best = None
    satisfiable = False
    for assignment in range(1 << variables):
        counts = [(clause & assignment).bit_count() for clause in masks]
        if all(count == 1 for count in counts):
            satisfiable = True
        if all(count % 2 for count in counts):
            cost = sum(counts)
            best = cost if best is None else min(best, cost)
    return satisfiable, best


def check_small_instances():
    triples = tuple(combinations(range(5), 3))
    matching_checks = 0
    examples = {}
    for size in range(len(triples) + 1):
        for clauses in combinations(triples, size):
            sat, optimum = exact_and_parity_costs(clauses, 5)
            check(sat == (optimum == size), f"cost criterion failed on {clauses}")
            rows = equations(clauses)
            raw_graphic = all(col.bit_count() <= 2 for col in columns(rows).values())
            if raw_graphic:
                basis = rows
            elif len(canonical_basis(rows)) <= 4:
                basis = small_graphic_basis(rows)
            else:
                continue
            if basis is None:
                continue
            witness = graphic_witness(clauses, basis)
            check((witness is not None) == sat, f"matching decision failed on {clauses}")
            if witness is not None:
                check(all(len(witness.intersection(c)) == 1 for c in clauses), "bad witness")
            if not raw_graphic:
                examples.setdefault(sat, clauses)
            matching_checks += 1
    check(set(examples) == {False, True}, "missing a hidden-graphic SAT or UNSAT case")
    print(f"Five variables: {2**len(triples)} formulas, {matching_checks} matching comparisons.")


def fano_lines():
    return tuple((i, (i + 1) % 7, (i + 3) % 7) for i in range(7))


def check_fano():
    lines = fano_lines()
    rows = equations(lines)
    cols = tuple(columns(rows).values())
    check(len(canonical_basis(cols)) == 4, "Fano incidence rank")
    check(canonical_basis(rows[:4]) == canonical_basis(rows), "wrong row basis")
    check(all(cols) and len(set(cols)) == 7, "Fano loop or parallel pair")
    check(all(a ^ b ^ c for a, b, c in combinations(cols, 3)), "Fano triangle")
    check(small_graphic_basis(rows) is None, "unexpected graphic basis")
    check(exact_and_parity_costs(lines, 7) == (False, 9), "wrong Fano optimum")
    print(f"Fano obstruction: rank 4, seven columns, no triangle, no graphic basis "
          f"among {comb(15, 4)} candidate row sets. Minimum cost 9 for 7 clauses.")


def check_integer_relaxation():
    clauses = (
        (0, 6, 9), (0, 7, 8), (0, 10, 11), (1, 6, 7), (1, 8, 11), (1, 9, 10),
        (2, 3, 7), (2, 4, 11), (2, 5, 10), (3, 4, 6), (3, 5, 9), (4, 5, 8),
    )
    signed = [-1] * 2 + [0] * 4 + [1] * 6
    check(all(sum(signed[i] for i in c) == 1 for c in clauses), "bad integer witness")
    check(len({tuple(sorted(p)) for c in clauses for p in combinations(c, 2)}) == 36,
          "two clauses overlap in more than one variable")
    check(all(col.bit_count() == 3 for col in columns(equations(clauses)).values()),
          "instance is not cubic")
    check(exact_and_parity_costs(clauses, 12) == (False, 24), "wrong parity optimum")
    print("12-variable cubic obstruction: signed integer witness; "
          "all 4,096 Boolean assignments checked; minimum parity cost 24 for budget 12.")


def planted_ring(blocks):
    clauses = []
    planted = set()
    for block in range(blocks):
        start = 17 * block
        for j, (a, b, c) in enumerate(fano_lines()):
            link = start + 7 + j
            if j == 0:
                connector = 17 * ((block - 1) % blocks) + 16
            elif j == 6:
                connector = start + 16
            elif j <= 3:
                connector = start + 14
            else:
                connector = start + 15
            clauses.extend(((start + a, start + b, link), (start + c, connector, link)))
            planted.add(link)
    return clauses, planted


def check_ring(blocks):
    clauses, planted = planted_ring(blocks)
    check(all(len(planted.intersection(c)) == 1 for c in clauses), "plant is invalid")
    rows = equations(clauses)
    incidence = columns(rows)
    check(len(incidence) == 17 * blocks and len(rows) == 14 * blocks, "wrong size")
    check(all(2 <= col.bit_count() <= 3 for col in incidence.values()), "wrong degree")
    # Traverse clause vertices, with variables as hyperedges.
    seen = {0}
    stack = [0]
    while stack:
        clause = stack.pop()
        for variable in clauses[clause]:
            for other in bits(incidence[variable]):
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
    check(len(seen) == len(rows), "ring disconnected")
    for block in range(blocks):
        originals = sum(1 << (17 * block + i) for i in range(7))
        for j, line in enumerate(fano_lines()):
            left = rows[14 * block + 2 * j] >> 1
            right = rows[14 * block + 2 * j + 1] >> 1
            check(incidence[17 * block + 7 + j] == 3 << (14 * block + 2 * j),
                  "contraction column is not the two-row link")
            recovered = ((left ^ right) & originals) >> (17 * block)
            check(recovered == sum(1 << i for i in line), "incorrect minor")


if __name__ == "__main__":
    check_small_instances()
    check_fano()
    check_integer_relaxation()
    for count in range(1, 13):
        check_ring(count)
    print("Planted connected family: witnesses, degrees and disjoint minors checked for k=1..12.")
