"""Finite checks for weighted matching connectors and edge-dependent SAT bases."""

from itertools import combinations, product
from math import inf
from random import Random


def matching_values(vertices, edges):
    """Return the nonnegative matching sum and the minimum signed additive cost."""
    if not vertices:
        return 1, 0
    first, *rest = vertices
    total, minimum = 0, inf
    for i, second in enumerate(rest):
        weight = edges.get((first, second))
        if weight is None:
            continue
        subtotal, subminimum = matching_values(rest[:i] + rest[i + 1:], edges)
        total += abs(weight) * subtotal
        minimum = min(minimum, weight + subminimum)
    return total, minimum


def check_matching_connectors():
    rng = Random(82071)
    partitions = ((0, 15), (3, 12), (5, 10), (9, 6))
    checked = 0
    for size in (4, 6, 8):
        for bipartite in (False, True):
            for _ in range(80):
                edges = {(i, j): rng.choice((-3, -2, -1, 1, 2, 3))
                         for i, j in combinations(range(size), 2)
                         if (not bipartite or i % 2 != j % 2) and rng.random() < 0.65}
                metrics = {
                    mask: matching_values(
                        [i for i in range(size) if not mask & (1 << i)], edges,
                    )
                    for partition in partitions for mask in partition
                }
                costs = [metrics[a][1] + metrics[b][1] for a, b in partitions]
                if costs.count(min(costs)) < 2:
                    raise AssertionError("a minimum complementary cost sum needs a tie")
                products = [metrics[a][0] * metrics[b][0] for a, b in partitions]
                if 2 * max(products) > sum(products):
                    raise AssertionError("nonnegative matching sums violate path exchange")
                checked += 1
    return checked


def parity_operators(prime):
    """Conjugates of diag(1,-1) in odd characteristic or diag(0,1) in characteristic two."""
    trace, determinant = (1, 0) if prime == 2 else (0, prime - 1)
    return [(a, b, c, d) for a, b, c, d in product(range(prime), repeat=4)
            if (a + d) % prime == trace and (a * d - b * c) % prime == determinant]


def apply_axis(tensor, matrix, axis, prime):
    bit = 1 << axis
    return tuple(
        (matrix[2 * bool(i & bit)] * tensor[i & ~bit]
         + matrix[2 * bool(i & bit) + 1] * tensor[i | bit]) % prime
        for i in range(8)
    )


def transpose(matrix):
    a, b, c, d = matrix
    return a, c, b, d


def compatible_first_operators(tensor, operators, prime):
    compatible = set()
    negative = tuple(-v % prime for v in tensor)
    for first in operators:
        on_first = apply_axis(tensor, first, 0, prime)
        for second in operators:
            on_second = apply_axis(tensor if prime == 2 else on_first, second, 1, prime)
            for third in operators:
                on_third = apply_axis(tensor if prime == 2 else on_second, third, 2, prime)
                if prime == 2:
                    summed = tuple((a + b + c) % prime
                                   for a, b, c in zip(on_first, on_second, on_third))
                    matches = summed in ((0,) * 8, tensor)
                else:
                    matches = on_third in (tensor, negative)
                if matches:
                    compatible.add(first)
                    break
            if first in compatible:
                break
    return compatible


def check_edge_bases():
    equality = tuple(int(i in (0, 7)) for i in range(8))
    for prime in (2, 3, 5):
        operators = parity_operators(prime)
        variable = compatible_first_operators(equality, operators, prime)
        for forbidden in (0, 1, 6, 7):
            clause = tuple(int(i != forbidden) for i in range(8))
            allowed = compatible_first_operators(clause, operators, prime)
            if prime != 2 and (not variable or not allowed):
                raise AssertionError("each local tensor separately admits a parity basis")
            if any(transpose(matrix) in allowed for matrix in variable):
                raise AssertionError(f"a shared parity basis exists over F_{prime}")
        if prime == 2 and variable:
            raise AssertionError("ternary equality alone has no parity basis in characteristic two")


if __name__ == "__main__":
    checked = check_matching_connectors()
    check_edge_bases()
    print(f"Checked the weighted exchange inequalities on {checked} four-terminal graphs.")
    print("Checked the edge-basis obstruction for ternary SAT factors over F_2, F_3 and F_5.")
