"""Cross-check truncated permanent computations and Fourier derivative ranks."""

from fractions import Fraction
from itertools import permutations
from math import factorial
from random import Random

from permanent_jets import add, block_rank_jet, multiply


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def check_jets():
    rng = Random(24017)
    for trial in range(96):
        cutoff = rng.randrange(1, 7)
        size = rng.randrange(1, 7)
        orders = [trial % 4, cutoff - 1, cutoff, 0][:trial % 5]
        zero = (0,) * cutoff

        def coefficient_vector():
            return tuple(rng.randrange(-3, 4) for _ in range(cutoff))

        split_rows = rng.randrange(size + 1)
        split_columns = rng.randrange(size + 1)
        partition = [(range(split_rows), range(split_columns)),
                     (range(split_rows, size), range(split_columns, size))]
        base = [[coefficient_vector() if (i < split_rows) == (j < split_columns) else zero
                 for j in range(size)] for i in range(size)]
        channels = [(order, [coefficient_vector() for _ in range(size)],
                     [coefficient_vector() for _ in range(size)]) for order in orders]
        dense = [row[:] for row in base]
        for order, u, v in channels:
            if order < cutoff:
                for i in range(size):
                    for j in range(size):
                        update = multiply(u[i], v[j])[:cutoff - order]
                        dense[i][j] = add(dense[i][j], (0,) * order + update)
        expected = zero
        for permutation in permutations(range(size)):
            term = (1,) + zero[1:]
            for i, j in enumerate(permutation):
                term = multiply(term, dense[i][j])
            expected = add(expected, term)
        check(block_rank_jet(base, partition, channels, cutoff) == expected,
              f"permanent mismatch on trial {trial}")


def parity_polynomial_derivative(value, start, degree):
    """Derivative of the polynomial interpolating (-1)**j on degree+1 integers."""
    delta = value - start
    falling, slope, result = 1, 0, 0
    for k in range(1, degree + 1):
        falling, slope = falling * (delta - k + 1), slope * (delta - k + 1) + falling
        sign = -1 if start % 2 else 1
        result += Fraction(sign * (-2)**k, factorial(k)) * slope
    return result


def walsh(values):
    values = list(values)
    step = 1
    while step < len(values):
        for start in range(0, len(values), 2 * step):
            for i in range(start, start + step):
                left, right = values[i], values[i + step]
                values[i], values[i + step] = left + right, left - right
        step *= 2
    return values


def rational_rank(rows):
    pivots = {}
    for source in rows:
        row = list(source)
        for col, value in enumerate(row):
            if not value:
                continue
            if col not in pivots:
                pivots[col] = [x / value for x in row]
                break
            for j in range(col, len(row)):
                row[j] -= value * pivots[col][j]
    return len(pivots)


def derivative_kernel(width, shift, sample, start, degree):
    count = 1 << (width * len(shift))
    mask = (1 << width) - 1
    multipliers = []
    nonzero = 0
    for frequency in range(count):
        bits = [((frequency >> (width * j)) & mask & gamma).bit_count() % 2
                for j, gamma in enumerate(shift)]
        value = sum((1 - Fraction(sample, j)) * bit for j, bit in enumerate(bits, 1))
        slope = parity_polynomial_derivative(value, start, degree)
        nonzero += bool(slope)
        multipliers.append(-Fraction(bits[sample - 1], sample) * slope)

    # Fourier inversion gives the true derivative entry for each row XOR column.
    return [value / count for value in walsh(multipliers)], nonzero


def check_derivative(width, shift, sample, start, degree):
    kernel, nonzero = derivative_kernel(width, shift, sample, start, degree)
    count = len(kernel)
    gamma = shift[sample - 1]
    separating_bit = (gamma & -gamma) << (width * (sample - 1))
    rows = [i for i in range(count) if i & separating_bit == 0]
    columns = [j for j in range(count) if j & separating_bit]
    rank = rational_rank([kernel[i ^ j] for j in columns] for i in rows)
    check(rank == nonzero // 2, f"off-block rank mismatch for shift {shift} at {sample}")
    return rank, len(rows)


def check_cycle_deformation():
    width, length, sample = 2, 2, 2
    count = 1 << (width * length)
    coordinate_bit = 1 << (width * (sample - 1))
    shifts = [g for g in range(count) if g & coordinate_bit and g & 3][:3]
    colors = sorted(range(count), key=lambda color: (color & ~coordinate_bit,
                                                      bool(color & coordinate_bit)))
    row_colors = colors * len(shifts)
    column_colors = [a ^ gamma for gamma in shifts for a in colors]
    size = len(row_colors)
    edges = [(u, v) for u in range(size) for v in (u, (u + 1) % size)]
    pairs = {(row_colors[u], column_colors[v]) for u, v in edges}
    check(len(pairs) == len(edges), "the cycle coloring repeats a pair")

    by_shift = {}
    for u, v in edges:
        gamma = row_colors[u] ^ column_colors[v]
        by_shift.setdefault(gamma, []).append((u, v))
    selected_rows = {u: i for i, u in enumerate(
        u for u in range(size) if not row_colors[u] & coordinate_bit
    )}
    selected_columns = {v: i for i, v in enumerate(
        v for v in range(size) if column_colors[v] & coordinate_bit
    )}
    matrix = [[Fraction(0)] * len(selected_columns) for _ in selected_rows]
    expected = 0
    for gamma, support in by_shift.items():
        blocks = tuple((gamma >> (width * i)) & 3 for i in range(length))
        kernel, nonzero = derivative_kernel(width, blocks, sample, 0, 2)
        if gamma & coordinate_bit:
            expected += nonzero // 2
        for u, _ in support:
            if u not in selected_rows:
                continue
            for _, v in support:
                if v in selected_columns:
                    value = kernel[row_colors[u] ^ column_colors[v]]
                    matrix[selected_rows[u]][selected_columns[v]] += value
    rank = rational_rank(matrix)
    check(rank == expected, "assembled off-block rank mismatch")
    return rank, len(selected_rows)


if __name__ == "__main__":
    check_jets()
    tests = [(1, (1,) * 6, i, 2, 3) for i in range(1, 7)]
    tests += [(1, (1,) * 6, 2, 0, 2), (2, (1, 3, 0), 2, 0, 2), (2, (2, 1, 3), 1, 0, 1)]
    ranks = [check_derivative(*test) for test in tests]
    check(any(0 < rank < limit for rank, limit in ranks), "missing a partial-rank case")
    check(any(rank == limit for rank, limit in ranks), "missing a full-rank case")
    print("Permanent jets: 96 exact comparisons against direct permutation sums.")
    print("Fourier derivative ranks (actual, off-block dimension):", ranks)
    print("Assembled cycle derivative rank (actual, off-block dimension):",
          check_cycle_deformation())
