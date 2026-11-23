"""Permanent of a block matrix with weighted rank-one channels, modulo z**cutoff."""

from collections import defaultdict
from functools import cache
from math import factorial, prod


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def multiply(left, right):
    cutoff = len(left)
    result = [0] * cutoff
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right[:cutoff - i]):
                result[i + j] += a * b
    return tuple(result)


def permanent(matrix, cutoff):
    """Subset dynamic program over truncated integer polynomials."""
    zero = (0,) * cutoff
    states = {0: (1,) + zero[1:]}
    for row in matrix:
        following = defaultdict(lambda: zero)
        for used, poly in states.items():
            for j, value in enumerate(row):
                bit = 1 << j
                if used & bit == 0 and any(value):
                    following[used | bit] = add(following[used | bit], multiply(poly, value))
        states = following
    return states.get((1 << len(matrix)) - 1, zero)


def side_labels(size, valuations, cutoff):
    """Use -1 for an index retained in its own block."""
    stack = [((), 0)]
    while stack:
        labels, weight = stack.pop()
        if len(labels) == size:
            yield labels
            continue
        stack.append((labels + (-1,), weight))
        for i, valuation in enumerate(valuations):
            if weight + valuation < cutoff:
                stack.append((labels + (i,), weight + valuation))


def sides(indices, vectors, valuations, cutoff):
    one = (1,) + (0,) * (cutoff - 1)
    for labels in side_labels(len(indices), valuations, cutoff):
        counts = [0] * len(valuations)
        kept = []
        value = one
        for index, label in zip(indices, labels):
            if label == -1:
                kept.append(index)
            else:
                counts[label] += 1
                value = multiply(value, vectors[label][index])
        if any(value):
            yield tuple(counts), tuple(kept), value


def block_rank_jet(base, partition, channels, cutoff):
    """Return per(base + sum(z**d * u * v.T)) modulo z**cutoff.

    The supplied partition lists pairs of row and column index sets. The base
    is supported only in those rectangles. Each polynomial is a tuple of
    cutoff integer coefficients. Channels are triples (nonnegative d, u, v).
    """
    zero = (0,) * cutoff
    valuations, left, right = zip(*channels) if channels else ((), (), ())
    empty = (0,) * len(valuations)

    @cache
    def weight(counts):
        return sum(v * count for v, count in zip(valuations, counts))

    totals = {(empty, empty): (1,) + zero[1:]}
    for rows, columns in partition:
        @cache
        def base_minor(kept_rows, kept_columns):
            return permanent([[base[i][j] for j in kept_columns] for i in kept_rows], cutoff)

        local = defaultdict(lambda: zero)
        row_options = list(sides(rows, left, valuations, cutoff))
        column_options = list(sides(columns, right, valuations, cutoff))
        for alpha, kept_rows, u in row_options:
            for beta, kept_columns, v in column_options:
                if len(kept_rows) != len(kept_columns):
                    continue
                value = multiply(multiply(u, v), base_minor(kept_rows, kept_columns))
                local[alpha, beta] = add(local[alpha, beta], value)

        extended = defaultdict(lambda: zero)
        for (alpha, beta), value in totals.items():
            for (a, b), local_value in local.items():
                new_alpha = add(alpha, a)
                new_beta = add(beta, b)
                if max(weight(new_alpha), weight(new_beta)) < cutoff:
                    key = (new_alpha, new_beta)
                    extended[key] = add(extended[key], multiply(value, local_value))
        totals = extended

    result = zero
    for (alpha, beta), value in totals.items():
        if alpha == beta:
            degree = weight(alpha)
            pairings = prod(factorial(count) for count in alpha)
            shifted = (0,) * degree + tuple(pairings * c for c in value[:cutoff - degree])
            result = add(result, shifted)
    return result
