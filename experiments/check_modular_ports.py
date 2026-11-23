"""Three-port determinant factors over quadratic finite fields."""

from itertools import count, product
from math import isqrt
from random import Random

from check_determinant_filter import chord_wires, compile_formula, direct_count
from check_pair_filter import find_split, separator_filtered_count


ZERO = (0, 0)
ONE = (1, 0)


class QuadraticField:
    """Represent a + b*t by (a, b), with t*t a nonsquare modulo an odd prime."""

    def __init__(self, prime):
        self.prime = prime
        self.roots = {a * a % prime: a for a in range(prime)}
        self.nonsquare = next(a for a in range(2, prime) if a not in self.roots)
        self.inverse_nonsquare = pow(self.nonsquare, -1, prime)

    def lift(self, value):
        return (value % self.prime, 0) if isinstance(value, int) else value

    def add(self, left, right):
        a, b = left
        c, d = right
        return ((a + c) % self.prime, (b + d) % self.prime)

    def sub(self, left, right):
        a, b = left
        c, d = right
        return ((a - c) % self.prime, (b - d) % self.prime)

    def mul(self, left, right):
        a, b = left
        c, d = right
        return ((a * c + self.nonsquare * b * d) % self.prime,
                (a * d + b * c) % self.prime)

    def inv(self, value):
        a, b = value
        scale = pow((a * a - self.nonsquare * b * b) % self.prime, -1, self.prime)
        return (a * scale % self.prime, -b * scale % self.prime)

    def sqrt_base(self, value):
        value %= self.prime
        if value in self.roots:
            return (self.roots[value], 0)
        return (0, self.roots[value * self.inverse_nonsquare % self.prime])


def ternary_factor(scope, table, field):
    """Compile the table in 000 through 111 order, allowing negated port wires."""
    p = field.prime
    anchor = next((mask for mask, value in enumerate(table) if value % p), None)
    if anchor is None:
        return scope, [{}, {}, {}, {}]
    weight = table[anchor] % p
    inverse_weight = pow(weight, -1, p)
    minors = [table[anchor ^ kept] * inverse_weight % p for kept in range(8)]
    a, b, c = (minors[i] for i in (4, 2, 1))
    pair_products = (a * b - minors[6], b * c - minors[3], c * a - minors[5])
    ab, bc, ca = pair_products
    cycle_sum = minors[7] - a * b * c + a * bc + b * ca + c * ab
    root = field.sqrt_base(cycle_sum * cycle_sum - 4 * ab * bc * ca)
    half = field.lift(pow(2, -1, p))
    u = field.mul(half, field.add(field.lift(cycle_sum), root))
    v = field.sub(field.lift(cycle_sum), u)
    if u == v == ZERO:
        forward = [ONE if value % p else ZERO for value in pair_products]
        backward = [field.lift(value) for value in pair_products]
    else:
        cycle_product = u if u != ZERO else v
        forward = [cycle_product, ONE, ONE]
        backward = [field.mul(field.lift(ab), field.inv(cycle_product)),
                    field.lift(bc), field.lift(ca)]
    matrix = [{0: a, 1: forward[0], 2: backward[2]},
              {0: backward[0], 1: b, 2: forward[1]},
              {0: forward[2], 1: backward[1], 2: c}, {3: weight}]
    flips = anchor ^ 7
    ports = tuple(-v if flips & mask else v for v, mask in zip(scope, (4, 2, 1)))
    return ports, matrix


def modular_determinant(matrix, field):
    result = ONE
    size = len(matrix)
    for i in range(size):
        pivot = next((j for j in range(i, size) if matrix[j][i] != ZERO), None)
        if pivot is None:
            return ZERO
        if pivot != i:
            matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
            result = field.sub(ZERO, result)
        pivot_row = matrix[i]
        result = field.mul(result, pivot_row[i])
        inverse = field.inv(pivot_row[i])
        for row in matrix[i + 1:]:
            if row[i] == ZERO:
                continue
            ratio = field.mul(row[i], inverse)
            for j in range(i + 1, size):
                row[j] = field.sub(row[j], field.mul(ratio, pivot_row[j]))
    return result


def modular_count(clauses, tables, pins=None, *, field):
    weights = {v: (1 - bit, bit) for v, bit in (pins or {}).items()}
    factors = [ternary_factor(scope, table, field) for scope, table in tables.items()]
    entries, ends = compile_formula(clauses, weights=weights, port_factors=factors)
    chords = chord_wires(ends)
    rows = [{j: (field.lift(weight), wire) for j, (weight, wire) in row.items()}
            for row in entries]
    total = ZERO
    for signs in product((False, True), repeat=len(chords)):
        negative = {wire for wire, flip in zip(chords, signs) if flip}
        matrix = [[ZERO] * len(rows) for _ in rows]
        for i, row in enumerate(rows):
            for j, (weight, wire) in row.items():
                matrix[i][j] = field.sub(ZERO, weight) if wire in negative else weight
        total = field.add(total, modular_determinant(matrix, field))
    residue, extension_part = total
    if extension_part:
        raise ArithmeticError("the weighted count must lie in the prime field")
    return residue * pow(2, -len(chords), field.prime) % field.prime


def ternary_filtered_count(clauses, max_side_rank=2):
    clauses = tuple(clauses)
    limit = 1 << len({abs(lit) for clause in clauses for lit in clause})
    value, modulus = 0, 1
    for p in count(3, 2):
        if any(p % divisor == 0 for divisor in range(3, isqrt(p) + 1, 2)):
            continue
        field = QuadraticField(p)

        def count_mod_p(local_clauses, tables, pins=None):
            return modular_count(local_clauses, tables, pins, field=field)

        residue = separator_filtered_count(clauses, 3, max_side_rank, count_mod_p)
        value += modulus * ((residue - value) * pow(modulus, -1, p) % p)
        modulus *= p
        if modulus > limit:
            return value


def check_ternary_factors():
    rng = Random(84022)
    scope = (1, 2, 3)
    for p in (3, 5, 7, 11):
        field = QuadraticField(p)
        tables = [(0,) * 8, (0, 2, 2, 2, 2, 2, 2, 1),
                  *(tuple(int(mask == i) for mask in range(8)) for i in range(8))]
        tables.extend(tuple(rng.randrange(p) for _ in range(8)) for _ in range(48))
        for table in tables:
            for bits, expected in zip(product((0, 1), repeat=3), table):
                actual = modular_count((), {scope: table}, dict(zip(scope, bits)), field=field)
                if actual != expected % p:
                    raise AssertionError(f"wrong count modulo {p} for {table}, {bits}")


def check_ternary_separators():
    rng = Random(99036)
    for _ in range(24):
        clauses = tuple(tuple(rng.choice((-1, 1)) * v for v in scope) for scope in
                        ((1, 2, 4), (2, 3, 4), (1, 3, 5), (1, 2, 5), (3, 5, 6)))
        if ternary_filtered_count(clauses) != direct_count(clauses):
            raise AssertionError(f"wrong exact count for {clauses}")

    lobes = 32
    for conflicting in (False, True):
        clauses = []
        for i, x in enumerate(range(4, lobes + 4)):
            flips = i % 8 if conflicting else 0
            a, b, c = (v * (-1 if flips & (1 << (v - 1)) else 1) for v in (1, 2, 3))
            clauses.extend(((x, a, b), (-x, b, c)))
        _, ends = compile_formula(clauses)
        if len(chord_wires(ends)) != 3 * lobes - 2:
            raise AssertionError("wrong cycle rank for the three-separator example")
        if ternary_filtered_count(clauses) != (0 if conflicting else 5 * 2 ** lobes + 2):
            raise AssertionError("wrong exact count for the three-separator example")


def check_four_port_table():
    clauses = []
    table = [int(kept == 0 or kept.bit_count() == 2 or kept == 7) for kept in range(16)]
    for kept, allowed in enumerate(table):
        if allowed:
            continue
        a, b, c, d = (v if kept & (1 << (v - 1)) else -v for v in range(1, 5))
        t = 5 + kept
        clauses.extend(((-a, t), (-b, t), (a, b, -t), (t, c, d)))
    for kept, expected in enumerate(table):
        pins = {v: (0, 1) if kept & (1 << (v - 1)) == 0 else (1, 0) for v in range(1, 5)}
        if direct_count(clauses, pins) != expected:
            raise AssertionError("wrong SAT table for the four-port obstruction")


def check_separator_padding():
    variables = (1, 2, 3, 4)
    a, b, c, d = (5, 6, 7, 8)
    force_hubs = [literals for literals in product((a, -a), (b, -b), (c, -c))
                  if literals != (-a, -b, -c)]
    force_hubs.append((-a, -b, d))
    padding = [clause for x in variables for clause in ((x, a, b), (x, c, d))]
    for formula in (((1, 2, 3), (2, 3, 4)),
                    (*product((1, -1), (2, -2), (3, -3)), (2, 3, 4))):
        padded = [*formula, *force_hubs, *padding]
        if direct_count(padded) != direct_count(formula):
            raise AssertionError("the padding must extend every assignment uniquely")
        scopes = [tuple(map(abs, clause)) for clause in padded]
        if any(find_split(scopes, size, float("inf")) for size in (1, 2, 3)):
            raise AssertionError("a separator of at most three variables survived padding")
        _, ends = compile_formula(padded)
        if len(chord_wires(ends)) != 2 * len(formula) + 3 * len(variables) + 13:
            raise AssertionError("wrong cycle rank for the padded example")


if __name__ == "__main__":
    check_ternary_factors()
    check_ternary_separators()
    check_four_port_table()
    check_separator_padding()
    print("Checked 232 ternary tables over four prime fields through their quadratic extensions.")
    print("Matched direct counting on 24 formulas and exact counts on two 94-cycle examples.")
    print("Verified the Boolean SAT table for the four-port obstruction.")
    print("Verified count-preserving padding with no separator of at most three variables.")
