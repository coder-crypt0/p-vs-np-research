# Limits of querying the determinant at signs

The [sparse determinant reduction](determinant-filter.md) recovers the SAT count by averaging over $2^b$ signs, where $b$ is the incidence cycle rank. On a fixed connected gadget graph of cycle rank $b$, there are uniquely satisfiable instances for which precisely one signing of the chords yields a nonzero determinant. The nonzero signing can be any of the $2^b$ sign patterns. An unsatisfiable instance on the same gadget graph has zero at every signing. Thus even deciding satisfiability can require exponentially many evaluations for an oracle limited to sign queries.

## Same wiring for both polarities

In the determinant construction, subdivide every literal connection through a three-index block. For negative literals, use the negation block from the reduction. For positive ones, use this identity relation, with ports at the first two indices.

```text
0 -1  0
1  0  0
0  0  1
```

Its principal deletion determinants give one for 00 and 11 and zero for 01 and 10. Subdivision preserves cycle rank. A unit clause is a single zero diagonal entry, since deleting its only port gives the empty determinant one. This lets us add pins without adding cycles.

## A hidden nonzero signing

Take $b$ variables $x_i$ and one shared variable $y$. For each $i$, choose one ordered clause $(x_i,\neg x_i,y)$ or $(\neg x_i,x_i,y)$. Add the unit clauses $\neg x_i$ for all $i$ and two copies of the unit clause $\neg y$. All variables are forced to zero, giving exactly one satisfying assignment. The wiring graph contains one four-wire cycle per $x_i$, formed by the two occurrences of $x_i$ and their intervening blocks. The cycles connect to the shared $y$ gadget through bridges. Unit clauses also attach through bridges.

In any determinant term, unpaired crossings form a circulation, so the bridges are paired. Summing over bridge choices enforces all unit pins. Equality at the $y$ gadget then leaves the $b$ cycles independent. The $x_i$ equality block has three ports in the order of its first occurrence, second occurrence and unit pin. With the pin set to zero, the only nonzero directed leak starts out of the first port and comes back into the second. The reverse orientation fails at the equality block, since deleting its second row and first column leaves a zero row. The cycle cover for the forward leak uses the first-to-second edge in the clause block and the two remaining edges in the $x_i$ block. Its unused internal indices use unit diagonal entries.

With one positive and one negative relay, this permutation cycle has length ten and determinant sign −1. The four crossing edges have product +1 after factoring out the Laurent monomial. The first relay contributes −1 in either polarity. The second relay contributes −1 if negative and +1 if positive. Therefore the coefficient of the nonconstant monomial is −1 for the first clause order and +1 for the second. On chord signs $s_i\in\{-1,1\}$, the local determinant after the pins is respectively $1-s_i$ and $1+s_i$. The product is $2^b$ at the chosen sign pattern and zero elsewhere.

Replacing just one of the duplicate $\neg y$ units by $y$ yields the identically zero case. The two bridge pins disagree at the $y$ equality block. With the fixed-size relay encoding this preserves the wiring graph and the dimension.

## Query bound and scope

Any assignment of signs to all wires is equivalent to a chord signing. Conjugating a gadget by a scalar sign leaves its internal block unchanged and multiplies every incident wire sign by that scalar. Starting from a root in each spanning tree, these conjugations make all forest signs positive. They preserve the determinant.

On the identically zero instance, a deterministic oracle algorithm sees only zeros. If it stops without testing a particular chord signing, the uniquely satisfiable instance with its nonzero value there gives the identical transcript. Every correct deterministic algorithm therefore needs $2^b$ queries in the worst case. A randomized algorithm testing at most $q$ signs has success probability at most $1/2+q/2^{b+1}$ on a distribution giving half its mass to the zero instance and the other half to a uniform nonzero location. A success guarantee of $2/3$ needs at least $2^b/3$ queries.

These oracle bounds apply when the algorithm knows the wiring and gets determinant values exclusively from sign evaluations. Inspecting the blocks allows the bridge algorithm described in the reduction note. Other complex or rational evaluations can reveal more than a sign query. Testing whether the Laurent determinant is identically zero also fails in general. The unsatisfiable formula $(x,x,y)\land\neg y\land\neg x$ has Laurent determinant $-z^{-1}$ in the one remaining chord variable. The single clause $(x,x,y)$ has determinant $3+z-2z^{-1}$. Its nonconstant terms go in both directions under any nonzero univariate weight, so the constant term lies between them. Both finite identities follow by evaluating at $1,-1,2$, since only exponents −1, 0 and 1 can appear. [Bridge algorithm](determinant-filter.md) · [Exact checks](../experiments/check_sign_filter_limits.py)

## Monomial substitution can change a zero constant coefficient

For every $b\ge 2$, construct a sparse determinant equal to

```text
∏ᵢ (2 + zᵢ + zᵢ⁻¹) − 2ᵇ.
```

Each factor comes from a four-index block with two main indices $u,v$ and two auxiliary indices $a,c$. Put $-z_i$ and $z_i^{-1}$ on the two opposite entries between $u,v$. Put one on the diagonals of $v,a,c$. Add paths $u\to a\to v$ with both weights +1 and $v\to c\to u$ with weights −1,+1. Eliminating the two auxiliary indices gives main entries $-z_i-1$ and $z_i^{-1}+1$, with diagonals 0 and 1. This gives the desired factor. Its principal minor after deleting $u$ has determinant one.

Connect the $u$ indices of the blocks in a directed cycle of weight magnitude two on every edge. Choose the sign of the last edge so the full crossing cycle has signed product $-2^b$. Only the all-used and all-unused crossing states balance the rows and columns at every block. Their contributions yield the displayed expression. Each row and column has at most three nonzeros, and the two appearances of every variable remain on opposite entries.

The original constant coefficient is zero. Under $z_i=t^{w_i}$, if two subsets of the $w_i$ have equal sum, their difference gives a nonzero exponent vector in $\{-1,0,1\}^b$ mapped to zero. All nonconstant coefficients are positive, so this makes the new constant coefficient strictly positive. Consequently separation requires $2^b$ distinct subset sums and $\sum_i |w_i|\ge 2^b-1$. Even an input-dependent monomial substitution cannot keep polynomial degree and preserve the zero answer for this family. The matrices here share the sparsity and reciprocal-pair structure of the SAT reduction. They are constructed as matrices, with no claim that they are outputs of the CNF compiler. [Exact checks](../experiments/check_sign_filter_limits.py)
