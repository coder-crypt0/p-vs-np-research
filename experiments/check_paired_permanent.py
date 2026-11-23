"""Verify the paired-minor theorem, the sparse family, and the cost search."""

from itertools import permutations, product
from math import prod
from random import Random

from paired_permanent import low_valuation_witness, paired_parity
from permanent_jets import permanent


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def minor(matrix, ports, remove_rows, remove_columns):
    rows = {r for (r, _), bit in zip(ports, remove_rows) if bit}
    cols = {c for (_, c), bit in zip(ports, remove_columns) if bit}
    kept = [i for i in range(len(matrix)) if i not in rows]
    return sum(prod(matrix[i][j] for i, j in zip(kept, order))
               for order in permutations(i for i in range(len(matrix)) if i not in cols))


def check_gadget(matrix, ports):
    patterns = tuple(product((0, 1), repeat=len(ports)))
    diagonal = set()
    leaks = set()
    for s in patterns:
        for t in patterns:
            if sum(s) == sum(t) and minor(matrix, ports, s, t) % 2:
                if s == t:
                    diagonal.add(s)
                else:
                    leaks.add((s, t))

    cube, certificate = paired_parity(matrix, ports)
    if leaks:
        check(cube is None and certificate in leaks, "incorrect leaking minor certificate")
    else:
        check(certificate is None, "unexpected leaking minor certificate")
        check((set(product(*cube)) if cube is not None else set()) == diagonal,
              "incorrect diagonal support")


def check_all_small_gadgets():
    rng = Random(1967)
    for bits in product((0, 1), repeat=9):
        matrix = [bits[3 * j:3 * j + 3] for j in range(3)]
        for count in range(4):
            check_gadget(matrix, tuple((i, i) for i in range(count)))
    for _ in range(100):
        size = rng.randrange(1, 6)
        count = rng.randrange(size + 1)
        matrix = [[rng.randrange(-3, 4) for _ in range(size)] for _ in range(size)]
        ports = tuple(zip(rng.sample(range(size), count), rng.sample(range(size), count)))
        check_gadget(matrix, ports)

    cycle = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    ports = ((0, 0), (1, 1), (2, 2))
    check(paired_parity(cycle, ports)[1] is not None, "cycle has an odd cross-state")


def chain_matrix(length):
    equality = [[-1, -1, 1], [-1, -1, 1], [1, 1, 1]]
    blocks = [[[0, 1], [1, 0]]] + [equality] * length + [[[1]]]
    size = sum(map(len, blocks))
    starts = []
    matrix = [[0] * size for _ in range(size)]
    start = 0
    for block in blocks:
        starts.append(start)
        for i, row in enumerate(block):
            matrix[start + i][start:start + len(block)] = row
        start += len(block)
    for i, (left, right) in enumerate(zip(starts, starts[1:])):
        right_port = left if i == 0 else left + 1
        matrix[right_port][right] = matrix[right][right_port] = 1
    return matrix


def check_high_valuation_chain():
    gadget = [[-1, -1, 1], [-1, -1, 1], [1, 1, 1]]
    ports = ((0, 0), (1, 1))
    weights = {(0, 0): -2, (1, 1): 1}
    patterns = tuple(product((0, 1), repeat=2))
    for rows in patterns:
        for columns in patterns:
            if sum(rows) == sum(columns):
                expected = weights.get(rows, 0) if rows == columns else 0
                check(minor(gadget, ports, rows, columns) == expected, "bad equality minor")
    for length in range(1, 6):
        equality = {(0, 0): 1, (1, 1): 0}
        constraints = [((0,), {(0,): 0})]
        constraints += [((i, i + 1), equality) for i in range(length)]
        check(low_valuation_witness(constraints, length - 1) is None, "valuation lower bound")
        check(low_valuation_witness(constraints, length) == dict.fromkeys(range(length + 1), 0),
              "unique normalized chain assignment")
        matrix = [[(value,) for value in row] for row in chain_matrix(length)]
        check(permanent(matrix, 1) == ((-2) ** length,), "wrong composed permanent")


def check_valuation_search():
    rng = Random(73889)
    for _ in range(180):
        size = rng.randrange(1, 9)
        constraints = []
        for _ in range(rng.randrange(1, 11)):
            scope = tuple(rng.sample(range(size), rng.randrange(1, min(4, size) + 1)))
            pins = [rng.choice((None, 0, 1)) for _ in scope]
            has_odd = rng.randrange(4) != 0
            costs = {}
            for bits in product((0, 1), repeat=len(scope)):
                if has_odd and all(pin is None or pin == bit for pin, bit in zip(pins, bits)):
                    costs[bits] = 0
                elif rng.randrange(3) != 0:
                    costs[bits] = rng.randrange(1, 5)
            constraints.append((scope, costs))

        budget = rng.randrange(7)
        unaffordable = budget + 1
        expected = any(
            sum(costs.get(tuple(bits[i] for i in scope), unaffordable)
                for scope, costs in constraints) <= budget
            for bits in product((0, 1), repeat=size)
        )
        witness = low_valuation_witness(constraints, budget)
        check((witness is not None) == expected, "valuation decision disagrees with brute force")
        if witness is not None:
            cost = sum(costs.get(tuple(witness[i] for i in scope), unaffordable)
                       for scope, costs in constraints)
            check(cost <= budget, "returned valuation witness is not affordable")


if __name__ == "__main__":
    check_all_small_gadgets()
    check_valuation_search()
    check_high_valuation_chain()
    print("Paired minors: 2,148 matrix and port combinations checked.")
    print("Valuation search: 180 comparisons and five exact signed cycle-cover sums checked.")
