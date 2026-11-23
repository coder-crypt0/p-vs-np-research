# Limits of querying the determinant at signs

The [sparse determinant reduction](determinant-filter.md) recovers the SAT count by averaging over (2^b) signs, where (b) is the incidence cycle rank. On a fixed connected gadget graph of cycle rank (b), there are uniquely satisfiable instances for which precisely one signing of the chords yields a nonzero determinant. The nonzero signing can be any of the (2^b) sign patterns. An unsatisfiable instance on the same gadget graph has zero at every signing. Thus even deciding satisfiability can require exponentially many evaluations for an oracle limited to sign queries.

## Same wiring for both polarities

In the determinant construction, subdivide every literal connection through a three-index block. For negative literals, use the negation block from the reduction. For positive ones, use this identity relation, with ports at the first two indices.

```text
0 -1  0
1  0  0
0  0  1
```

Its principal deletion determinants give one for 00 and 11 and zero for 01 and 10. Subdivision preserves cycle rank. A unit clause is a single zero diagonal entry, since deleting its only port gives the empty determinant one. This lets us add pins without adding cycles.

## A hidden nonzero signing

Take (b) variables (x_i) and one shared variable (y). For each (i), choose one ordered clause ((x_i,\neg x_i,y)) or ((\neg x_i,x_i,y)). Add the unit clauses (\neg x_i) for all (i) and two copies of the unit clause (\neg y). All variables are forced to zero, giving exactly one satisfying assignment. The wiring graph contains one four-wire cycle per (x_i), formed by the two occurrences of (x_i) and their intervening blocks. The cycles connect to the shared (y) gadget through bridges. Unit clauses also attach through bridges.

In any determinant term, unpaired crossings form a circulation, so the bridges are paired. Summing over bridge choices enforces all unit pins. Equality at the (y) gadget then leaves the (b) cycles independent. The (x_i) equality block has three ports in the order of its first occurrence, second occurrence and unit pin. With the pin set to zero, the only nonzero directed leak starts out of the first port and comes back into the second. The reverse orientation fails at the equality block, since deleting its second row and first column leaves a zero row. The cycle cover for the forward leak uses the first-to-second edge in the clause block and the two remaining edges in the (x_i) block. Its unused internal indices use unit diagonal entries.

With one positive and one negative relay, this permutation cycle has length ten and determinant sign −1. The four crossing edges have product +1 after factoring out the Laurent monomial. The first relay contributes −1 in either polarity. The second relay contributes −1 if negative and +1 if positive. Therefore the coefficient of the nonconstant monomial is −1 for the first clause order and +1 for the second. On chord signs (s_i\in\{-1,1\}), the local determinant after the pins is respectively (1-s_i) and (1+s_i). The product is (2^b) at the chosen sign pattern and zero elsewhere.

Replacing just one of the duplicate (\neg y) units by (y) yields the identically zero case. The two bridge pins disagree at the (y) equality block. With the fixed-size relay encoding this changes neither the wiring graph nor the dimension.

## Query bound and scope

Any assignment of signs to all wires is equivalent to a chord signing. Conjugating a gadget by a scalar sign leaves its internal block unchanged and multiplies every incident wire sign by that scalar. Starting from a root in each spanning tree, these conjugations make all forest signs positive. They preserve the determinant.

On the identically zero instance, a deterministic oracle algorithm sees only zeros. If it stops without testing a particular chord signing, the uniquely satisfiable instance with its nonzero value there gives the identical transcript. Every correct deterministic algorithm therefore needs (2^b) queries in the worst case. A randomized algorithm testing at most (q) signs has success probability at most (1/2+q/2^{b+1}) on a distribution giving half its mass to the zero instance and the other half to a uniform nonzero location. A success guarantee of (2/3) needs at least (2^b/3) queries.

These are oracle bounds. The algorithm may know the wiring, but must obtain the determinant values from sign evaluations and cannot read the blocks. The family has simple Boolean descriptions, and different complex or rational evaluations can reveal more than a sign query. Conversely, testing whether the Laurent determinant is identically zero does not suffice in general. The unsatisfiable formula ((x,x,y)\land\neg y\land\neg x) has Laurent determinant (-z^{-1}) in the one remaining chord variable. The single clause ((x,x,y)) has determinant (3+z-2z^{-1}), so in general leakage can occur in both directions and no choice of nonzero univariate weight makes the constant term the lowest Laurent term. Both finite identities follow by evaluating at (1,-1,2), since only exponents −1, 0 and 1 can appear. [Exact checks](../experiments/check_sign_filter_limits.py)
