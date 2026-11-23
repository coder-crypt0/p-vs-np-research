# The parity of permanent gadgets

Computing the permanent modulo two is ordinary determinant computation. For any fixed k, the permanent modulo 2^k is also computable in polynomial time. If a SAT reduction preserved nonzeroness at small precision, that would offer a way around exact counting. Known local counting reductions can instead introduce a large power of two. For example, the sparse signed reduction of Dell et al. uses even constants at equality interfaces. The argument below explains an obstruction that applies to any such interface with paired ports. [Valiant](https://doi.org/10.1016/0304-3975(79)90044-6) · [Modular permanent algorithms](https://doi.org/10.4230/LIPIcs.MFCS.2021.36) · [Sparse reduction, Section 3](https://www.cs.bme.hu/~dmarx/papers/dell-tutte-talg.pdf)

## A local parity theorem

Take an n by n integer matrix M with k designated ports, each consisting of a row and a column. The designated rows are distinct, as are the designated columns. For port sets S and T with the same cardinality, let g(S,T) be the permanent after removing the rows of S and columns of T. The diagonal g(S,S) represents selecting the same ports for deletion on both sides. Assume g(S,T) is even whenever S differs from T. Then the S with g(S,S) odd form either the empty set or a subcube of the k bits. Each bit is pinned to zero, pinned to one, or free independently.

To prove this, reduce M modulo two and let eᵢ be the unit column of the row at port i. Given S,T, take all non-port columns of M, the port columns of M outside T, and the unit columns at ports in S. Expanding the determinant along the unit columns gives g(S,T) modulo two. If the non-port columns are dependent, every minor in question is even. Otherwise quotient their span out of the n-dimensional column space. Its dimension is k. The projected port columns of M and projected unit columns form k pairs.

A basis chosen from those pairs corresponds to an odd minor. Our assumption forces every such basis to choose exactly one column from each pair. Fix one basis, if one exists. Consider the alternative vector in pair i. If its basis expansion has a nonzero coefficient at a different pair j, replacing the basis vector at j by that alternative produces a basis with both columns from i and none from j. This would give an odd mismatched minor. Therefore the alternative in every pair lies in the span of its own selected vector. Over the binary field, it is either zero or the same vector. The selections are independent, as claimed.

This proof also gives a polynomial algorithm. Gaussian elimination either finds the subcube or gives the two deletion masks of an odd mismatched minor. It does not enumerate the minor table. The mismatch hypothesis is essential. A three-cycle permutation matrix gives odd diagonal values only for selecting no ports or all ports, but has odd mismatched minors. [Implementation](../experiments/paired_permanent.py)

## Consequence for local reductions

When all mismatched minors vanish as integers, such gadgets compose along paired ports. A selected wire uses both crossing edges. Once the crossing edges are fixed, the remaining sum factors into the gadget minors. Any mismatched wire pattern has a gadget with mismatched row and column deletions and contributes zero.

A uniformly weighted Boolean relation whose allowed masks do not form a subcube must have an even common weight. Equality on two bits, XOR on two bits, and exactly-one on three bits are examples. In a construction using b independent occurrences of uniformly weighted coupled gadgets, each allowed product is divisible by 2^b, before summing over assignments. This obstruction concerns interfaces with local cancellation of *all* mismatched minors. It leaves room for reductions in which mismatched states cancel only after composing gadgets, or which encode solutions in another quantity.

Unequal weights alone also provide no generic small precision guarantee. The integer matrix

```text
-1 -1  1
-1 -1  1
 1  1  1
```

with ports at the first two diagonal positions realizes equality with weight −2 on 00 and weight 1 on 11. All its other paired minors and all mismatched minors are zero. Connect b copies in a chain using two unit-weight crossing edges per wire, pin the left end to zero with a two-by-two swap matrix, and terminate the right end with the one-by-one matrix [1]. The only supported assignment is all zero. The signed permanent is (−2)^b, so its residues vanish at every precision through 2^b. This is a connected family of signed integer matrices with at most four nonzero entries per row or column.

## Search without cancellation under a valuation budget

Given explicit local paired signatures, give each nonzero state the exponent of two dividing its weight. Zero states are forbidden. The valuation of a product is the sum of its local costs. Each gadget's cost-zero states are a subcube or empty, so imposing cost zero on all undecided gadgets amounts to Boolean pins. A contradiction to those pins has one of three forms. One gadget has no zero-cost state, a gadget conflicts with an already fixed bit, or two gadgets pin the same bit oppositely.

Any affordable completion must use a positive-cost state at one of the at most two unresolved gadgets in that conflict. Branch on that gadget and on each affordable positive-cost state in its table, fix its bits, and subtract its cost. A branch without conflicts extends to an assignment because every remaining set of allowed zero-cost states is a subcube. With arity at most a, the time is polynomial in the input size times O((1+2^(a+1))^K) for budget K. For constant a and K = O(log n), this directly detects an affordable nonzero product in polynomial time. The calculation bypasses cancellation between products. It supplies no bound on the minimum valuation for unrestricted inputs, and the signed chain shows that the minimum can grow linearly even when only one assignment is supported.

The experiments check all binary three-by-three matrices under four choices of port count, 100 further signed matrices with independently positioned ports, 180 valuation-search instances against direct enumeration, and five exact permanents of the connected chain. [Checks](../experiments/check_paired_permanent.py)
