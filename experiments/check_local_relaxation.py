"""Reproduce the parity-to-exact-one gadget and its local distribution checks."""

from collections import Counter
from itertools import combinations, product
from math import gcd

import networkx as nx

from check_parity import check, fano_lines


def clause_witness(values):
    """Return the selector and slack values for an ordinary OR clause."""
    first_true = next((i for i, value in enumerate(values) if value), 0)
    choices = tuple(int(i == first_true) for i in range(3))
    slacks = tuple(value & (1 - choice) for value, choice in zip(values, choices))
    return choices + slacks


def translate(incident, charges):
    """Return exact-one clauses, local Boolean functions and a signed witness."""
    functions = []
    signed = []
    clauses = []

    def fresh(support, evaluate, integer_value):
        index = len(functions)
        functions.append((support, evaluate))
        signed.append(integer_value)
        return index

    def constant(value):
        return fresh((), lambda _: value, value)

    def zero_port():
        # The six remaining Fano lines force the three points on the deleted line true.
        deleted, *remaining = fano_lines()
        points = [constant(int(i in deleted)) for i in range(7)]
        clauses.extend(tuple(points[i] for i in line) for line in remaining)
        port, spare = constant(0), constant(0)
        clauses.append((points[deleted[0]], port, spare))
        return port

    literals = {}
    for edge in sorted({e for scope in incident for e in scope}):
        wire = []
        for _ in range(4):
            wire.append(fresh((edge,), lambda values: values[0], 0))
            wire.append(fresh((edge,), lambda values: 1 - values[0], 1))
        for index, endpoint in enumerate(wire):
            clauses.append((endpoint, wire[(index + 1) % len(wire)], zero_port()))
        literals[edge] = (wire[::2], wire[1::2])

    def make_or(scope, forbidden):
        def keep_output(index):
            def output(values):
                lit_values = tuple(x ^ y for x, y in zip(values, forbidden))
                return clause_witness(lit_values)[index]
            return output

        # At the signed assignment every original edge bit is zero.
        signed_values = (1, 0, 0, forbidden[0] - 1, forbidden[1], forbidden[2])
        aux = [fresh(scope, keep_output(i), value) for i, value in enumerate(signed_values)]
        clauses.append(tuple(aux[:3]))
        for i, edge in enumerate(scope):
            opposite = literals[edge][1 - forbidden[i]].pop()
            clauses.append((aux[i], opposite, aux[3 + i]))

    for scope, charge in zip(incident, charges):
        for forbidden in product((0, 1), repeat=3):
            if sum(forbidden) % 2 != charge:
                make_or(scope, forbidden)
    return clauses, functions, signed


def check_or():
    combinations_checked = 0
    for literals in product((0, 1), repeat=3):
        extensions = []
        for aux in product((0, 1), repeat=6):
            if sum(aux[:3]) == 1 and all(
                aux[i] + 1 - value + aux[i + 3] == 1 for i, value in enumerate(literals)
            ):
                extensions.append(aux)
        check(bool(extensions) == any(literals), "OR gadget Boolean projection")
        if extensions:
            check(clause_witness(literals) in extensions, "local witness function")
        signed = (1, 0, 0, literals[0] - 1, literals[1], literals[2])
        check(all(signed[i] + 1 - value + signed[i + 3] == 1
                  for i, value in enumerate(literals)), "signed witness")
        combinations_checked += 64
    return combinations_checked


def project(solutions, scope):
    counts = Counter(tuple(x[i] for i in scope) for x in solutions)
    factor = gcd(*counts.values())
    return {values: count // factor for values, count in counts.items()}


def check_k4():
    graph = nx.complete_graph(4)
    edges = tuple(sorted(graph.edges))
    incident = [tuple(i for i, edge in enumerate(edges) if v in edge) for v in graph]
    charges = (1, 0, 0, 0)
    clauses, functions, signed = translate(incident, charges)
    check(all(len(set(c)) == 3 for c in clauses), "repeated variable")
    pairs = [tuple(sorted(pair)) for c in clauses for pair in combinations(c, 2)]
    check(len(set(pairs)) == len(pairs), "clauses share two variables")
    check(all(sum(signed[i] for i in c) == 1 for c in clauses), "signed witness")
    check(max(Counter(i for c in clauses for i in c).values()) == 3, "occurrence bound")

    by_defect = {v: [] for v in graph}
    for x in product((0, 1), repeat=len(edges)):
        mismatch = [v for v, scope in enumerate(incident)
                    if sum(x[i] for i in scope) % 2 != charges[v]]
        check(bool(mismatch), "odd-charge source is unexpectedly satisfiable")
        if len(mismatch) == 1:
            by_defect[mismatch[0]].append(x)

    giant = {}
    for size in range(4):
        for scope in combinations(range(len(edges)), size):
            remainder = graph.copy()
            remainder.remove_edges_from(edges[i] for i in scope)
            component = max(nx.connected_components(remainder), key=len)
            check(len(component) > len(graph) / 2, "missing giant")
            giant[scope] = component
            law = project(by_defect[min(component)], scope)
            check(all(project(by_defect[v], scope) == law for v in component),
                  "defect choice changed the projected law")
    marginals = 0
    for scope, component in giant.items():
        for size in range(len(scope) + 1):
            for smaller in combinations(scope, size):
                check(project(by_defect[min(component)], smaller) ==
                      project(by_defect[min(giant[smaller])], smaller), "marginal mismatch")
                marginals += 1

    for clause in clauses:
        scope = tuple(sorted({edge for i in clause for edge in functions[i][0]}))
        for assignment in by_defect[min(giant[scope])]:
            values = (evaluate(tuple(assignment[i] for i in support))
                      for support, evaluate in (functions[i] for i in clause))
            check(sum(values) == 1, "clause not supported by its local distribution")

    return len(giant), marginals, len(clauses)


if __name__ == "__main__":
    gadget_cases = check_or()
    scopes, marginals, clauses = check_k4()
    print(f"OR gadget: {gadget_cases} Boolean cases checked and all signed witnesses checked.")
    print(f"K4 lift: {scopes} source scopes, {marginals} marginal identities, "
          f"{clauses} lifted clause distributions and one global signed witness checked.")
