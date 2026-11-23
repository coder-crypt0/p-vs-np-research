"""Parity of paired deletion minors and search below a 2-adic cost budget."""

from graphic_sat import canonical_basis


def reduce_vector(vector, basis):
    for row in basis:
        pivot = 1 << (row.bit_length() - 1)
        if vector & pivot:
            vector ^= row
    return vector


def deletion_patterns(labels, ports):
    """Decode a basis of identity columns (tag 1) and matrix columns (tag 0)."""
    selected = set(labels)
    rows = tuple(int((i, 1) in selected) for i in range(ports))
    columns = tuple(int((i, 0) not in selected) for i in range(ports))
    return rows, columns


def paired_parity(matrix, ports):
    """Return (odd subcube, odd mismatched minor certificate).

    A port is a (row, column) pair. The empty support returns (None, None).
    A certificate is a pair of distinct port bit strings for row and column deletions.
    """
    size = len(matrix)
    paired_rows, paired_columns = zip(*ports) if ports else ((), ())
    if (any(len(row) != size for row in matrix)
            or len(set(paired_rows)) != len(ports)
            or len(set(paired_columns)) != len(ports)
            or any(not 0 <= index < size for index in (*paired_rows, *paired_columns))):
        raise ValueError("expected a square matrix and distinct valid row and column ports")

    matrix_columns = [sum((matrix[i][j] & 1) << i for i in range(size))
                      for j in range(size)]
    internal = [col for j, col in enumerate(matrix_columns) if j not in paired_columns]
    internal_basis = canonical_basis(internal)
    if len(internal_basis) != len(internal):
        return None, None
    quotient = {}
    for i, (row, col) in enumerate(ports):
        quotient[i, 0] = reduce_vector(matrix_columns[col], internal_basis)
        quotient[i, 1] = reduce_vector(1 << row, internal_basis)

    labels = []
    basis = ()
    for label, vector in quotient.items():
        if reduce_vector(vector, basis):
            labels.append(label)
            basis = canonical_basis((*basis, vector))
    if len(labels) != len(ports):
        return None, None

    rows, columns = deletion_patterns(labels, len(ports))
    if rows != columns:
        return None, (rows, columns)

    allowed = []
    for i, bit in enumerate(rows):
        alternate = quotient[i, 1 - bit]
        if alternate == 0:
            allowed.append(frozenset((bit,)))
        elif alternate == quotient[i, bit]:
            allowed.append(frozenset((0, 1)))
        else:
            for position, (j, _) in enumerate(labels):
                if j == i:
                    continue
                candidate = labels[:position] + [(i, 1 - bit)] + labels[position + 1:]
                if len(canonical_basis(quotient[label] for label in candidate)) == len(ports):
                    return None, deletion_patterns(candidate, len(ports))
            raise RuntimeError("a nonzero alternative must have a basis coordinate")
    return tuple(allowed), None


def low_valuation_witness(constraints, budget):
    """Find an assignment whose sum of local nonnegative costs is at most budget.

    Each constraint is (scope, {bit tuple: cost}). Omitted tuples are forbidden.
    The cost-zero tuples must be empty or a subcube.
    """
    prepared = []
    variables = set()
    for scope, costs in constraints:
        scope = tuple(scope)
        if len(set(scope)) != len(scope) or any(
            len(bits) != len(scope) or not set(bits) <= {0, 1} or cost < 0
            for bits, cost in costs.items()
        ):
            raise ValueError("invalid constraint scope or costs")
        variables.update(scope)
        odd_states = [bits for bits, cost in costs.items() if cost == 0]
        options = [set(bits[i] for bits in odd_states) for i in range(len(scope))]
        if odd_states:
            if len(odd_states) != 2 ** sum(len(bits) == 2 for bits in options):
                raise ValueError("zero-cost states must form a subcube")
            pins = {name: next(iter(bits)) for name, bits in zip(scope, options)
                    if len(bits) == 1}
        else:
            pins = None
        evens = [(dict(zip(scope, bits)), cost) for bits, cost in costs.items() if cost]
        prepared.append((pins, evens))

    def first_conflict(marked, fixed):
        values = dict(fixed)
        owners = {}
        for i, (pins, _) in enumerate(prepared):
            if i in marked:
                continue
            if pins is None:
                return (i,), values
            for name, bit in pins.items():
                if name in values and values[name] != bit:
                    return (i,) if name in fixed else (i, owners[name]), values
                values[name] = bit
                owners[name] = i
        return (), values

    states = [(frozenset(), {}, budget)]
    while states:
        marked, fixed, remaining = states.pop()
        if remaining < 0:
            continue
        conflict, values = first_conflict(marked, fixed)
        if not conflict:
            return {name: values.get(name, 0) for name in variables}
        for i in conflict:
            for assignment, cost in prepared[i][1]:
                if cost <= remaining and all(fixed.get(name, bit) == bit
                                             for name, bit in assignment.items()):
                    states.append((marked | {i}, fixed | assignment, remaining - cost))
    return None
