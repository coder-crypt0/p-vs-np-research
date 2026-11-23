"""Check the point-mass family and the scope of its evaluation obstruction."""

from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

from check_determinant_filter import (
    bridge_filtered_count, chord_wires, compile_formula, determinant, direct_count,
    signed_evaluation,
)


def make_formula(target, unsatisfiable=False):
    hub = len(target) + 1
    clauses = [(i, -i, hub) if sign == -1 else (-i, i, hub)
               for i, sign in enumerate(target, 1)]
    clauses.extend((-i,) for i in range(1, hub))
    clauses.extend(((-hub,), (hub if unsatisfiable else -hub,)))
    return clauses


def check(test, message):
    if not test:
        raise AssertionError(message)


def check_determinant_signs():
    evaluations = 0
    rng = Random(21907)
    for rank in range(1, 4):
        sign_vectors = tuple(product((-1, 1), repeat=rank))
        zero_formula = make_formula(sign_vectors[0], unsatisfiable=True)
        zero_entries, zero_ends = compile_formula(zero_formula, pad_literals=True)
        chords = chord_wires(zero_ends)
        check(direct_count(zero_formula) == 0 and len(chords) == rank, "bad zero instance")
        for signs in sign_vectors:
            check(signed_evaluation(zero_entries, dict(zip(chords, signs))) == 0,
                  "unsatisfiable family had a nonzero sign evaluation")
            evaluations += 1

        for target in sign_vectors:
            formula = make_formula(target)
            entries, ends = compile_formula(formula, pad_literals=True)
            check(direct_count(formula) == 1 and ends == zero_ends, "bad hidden spike input")
            check(len(entries) == len(zero_entries), "matrices have different dimensions")
            for signs in sign_vectors:
                expected = (1 << rank) if signs == target else 0
                check(signed_evaluation(entries, dict(zip(chords, signs))) == expected,
                      f"wrong evaluation for {target} at {signs}")
                evaluations += 1

            # Changing a vertex sign only moves signs between incident wires.
            vertices = {vertex for edge in ends for vertex in edge}
            gauge = {v: rng.choice((-1, 1)) for v in sorted(vertices)}
            signs = dict(zip(chords, target))
            twisted = {i: signs.get(i, 1) * gauge[u] * gauge[v]
                       for i, (u, v) in enumerate(ends)}
            check(signed_evaluation(entries, twisted) == 1 << rank, "sign gauge changed det")
    return evaluations


def scaled_evaluation(entries, chord):
    """Evaluate at two on one chord, then clear a factor two from each row."""
    matrix = [[0] * len(entries) for _ in entries]
    for i, row in enumerate(entries):
        for j, (weight, wire) in row.items():
            scale = 2
            if wire == chord:
                scale = 4 if weight == -1 else 1
            matrix[i][j] = scale * weight
    return Fraction(determinant(matrix), 1 << len(entries))


def one_chord_coefficients(clauses):
    entries, ends = compile_formula(clauses, pad_literals=True)
    (chord,) = chord_wires(ends)
    positive = signed_evaluation(entries, {chord: 1})
    negative = signed_evaluation(entries, {chord: -1})
    constant = Fraction(positive + negative, 2)
    total = Fraction(positive - negative, 2)
    at_two = scaled_evaluation(entries, chord)
    forward = (at_two - constant - total / 2) * Fraction(2, 3)
    return constant, forward, total - forward


def scalar_compression_matrix(size):
    """Return a degree-three matrix whose determinant has a zero constant coefficient."""
    entries = [{} for _ in range(4 * size)]
    for i in range(size):
        u, v, a, b = range(4 * i, 4 * i + 4)
        following = 4 * ((i + 1) % size)
        cycle_weight = -2 if i == size - 1 and size % 2 else 2
        entries[u].update({v: (-1, i), a: (1, None), following: (cycle_weight, None)})
        entries[v].update({u: (1, i), b: (-1, None), v: (1, None)})
        entries[a].update({a: (1, None), v: (1, None)})
        entries[b].update({b: (1, None), u: (1, None)})
    return entries


def evaluate_at_grid(entries, values):
    matrix = [[0] * len(entries) for _ in entries]
    for i, row in enumerate(entries):
        for j, (weight, variable) in row.items():
            factor = 2
            if variable is not None:
                z = values[variable]
                factor = (4 if weight == -1 else 1) if z == 2 else 2 * z
            matrix[i][j] = factor * weight
    return Fraction(determinant(matrix), 1 << len(entries))


def recovered_constant(weights):
    coefficients = Counter({0: 1})
    for a in weights:
        following = Counter()
        for exponent, value in coefficients.items():
            following[exponent] += 2 * value
            following[exponent + a] += value
            following[exponent - a] += value
        coefficients = following
    return coefficients[0] - (1 << len(weights))


def check_scalar_compression():
    for size in range(2, 5):
        entries = scalar_compression_matrix(size)
        for values in product((1, -1, 2), repeat=size):
            expected = Fraction(1)
            for value in values:
                expected *= 2 + value + Fraction(1, value)
            check(evaluate_at_grid(entries, values) == expected - (1 << size),
                  "incorrect scalar-compression determinant")
    check(recovered_constant((1, 2, 3)) > 0, "collision should change the zero answer")
    check(recovered_constant((1, 2, 4)) == 0, "separated subset weights should be exact")


if __name__ == "__main__":
    evaluations = check_determinant_signs()
    check(one_chord_coefficients([(1, 1, 2)]) == (3, 1, -2), "wrong two-way leakage")
    check(one_chord_coefficients([(1, 1, 2), (-2,), (-1,)]) == (0, 0, -1),
          "wrong unsatisfiable Laurent polynomial")
    for unsatisfiable in (False, True):
        formula = make_formula((-1, 1) * 24, unsatisfiable)
        entries, ends = compile_formula(formula, pad_literals=True)
        check(len(chord_wires(ends)) == 48, "wrong incidence cycle rank")
        check(bridge_filtered_count(entries, ends) == int(not unsatisfiable),
              "bridge filtering failed on the large oracle family")
    check_scalar_compression()
    print(f"Sign-only filter: {evaluations} exact oracle values and two Laurent examples checked.")
    print("Bridge filtering: both 48-cycle oracle instances solved exactly.")
    print("Scalar compression: 117 exact grid values and the zero-coefficient collision checked.")
