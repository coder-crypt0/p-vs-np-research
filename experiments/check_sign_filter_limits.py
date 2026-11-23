"""Check the point-mass family and the scope of its evaluation obstruction."""

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
    print(f"Sign-only filter: {evaluations} exact oracle values and two Laurent examples checked.")
    print("Bridge filtering: both 48-cycle oracle instances solved exactly.")
