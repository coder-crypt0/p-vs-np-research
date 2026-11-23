# Can matching decide positive one-in-three SAT?

Let A be the clause-by-variable incidence matrix of a positive one-in-three instance with m clauses of three distinct variables. Let dᵢ be the number of occurrences of variable i. A Boolean vector x is a solution exactly when Ax = 1 over the integers. Gaussian elimination over the binary field instead permits a clause with all three variables true.

## An exact decision criterion

Minimize W(x) = Σᵢ dᵢxᵢ subject to Ax = 1 modulo two and x binary. For each feasible x, let b(x) count clauses in which all three variables are true. Counting chosen occurrences first by variable and then by clause gives

```text
W(x) = m + 2 b(x).
```

Thus the original instance is satisfiable exactly when the minimum is m. Every positive ternary input has at least one parity-feasible assignment, since the all-ones vector works. The weights are measured on the original clauses throughout.

The Lean theorem `sat_iff_parity_at_budget` proves the equivalent existence statement for arbitrary finite lists of clauses with shared Boolean variables. Its proof uses that every odd row count is at least one, with equality in the sum exactly when each count equals one. It needs no assumptions about complexity classes. [Lean proof](../proofs/ExactOne.lean)

## When the minimization is polynomial

Suppose binary row operations and removal of redundant equations give an equivalent system Bx = b in which every variable occurs in at most two rows. Transform the right-hand side with the rows. For each weight-two column draw an edge between its row vertices. Join a weight-one column to an extra root vertex. Zero columns can be set to zero. Edge weights are the original dᵢ. The terminals T are the rows with right-hand side one, together with the root when needed to make |T| even.

The feasible edge sets are exactly those with odd-degree vertex set T. A minimum-weight such set is a minimum T-join. Within each component there must be an even number of terminals. Pair them by a minimum-weight perfect matching in their shortest-path metric and take the symmetric difference of the corresponding paths. Any feasible edge set decomposes into paths pairing its terminals and cycles, so its cost is at least the metric matching optimum. Conversely the reconstructed paths yield a feasible edge set costing at most that optimum. Nonnegative weights suffice. This gives a polynomial decision algorithm given the graphic representation. Graphic recognition of binary matroids is also classical. [Edmonds and Johnson](https://doi.org/10.1007/BF01580113) · [Bixby and Cunningham](https://doi.org/10.1287/moor.5.3.321)

The implementation verifies equivalence of the augmented row spaces before solving and checks each returned SAT witness against the input. Its small-rank basis search is an exhaustive diagnostic capped at rank four. [Solver](../experiments/graphic_sat.py)

## A seven-clause obstruction

Number seven variables modulo seven and use the seven cyclic translates of the triple {0,1,3}. This is the incidence matrix of the Fano plane. Its binary column rank is four. One certificate uses the first four rows (indices 0 through 3) as coordinates. The remaining three rows are their binary sums, as the accompanying experiment verifies by row reduction. Every column has odd weight three, all columns are different, and no three sum to zero because the parity of their sum is odd.

If the columns had a graph-incidence representation, this rank-four graphic matroid could be represented by a connected graph on five vertices. Indeed, discard isolated vertices and identify one chosen vertex from each component into a single vertex. This preserves cycles and leaves rank plus one vertices. The graph would have seven edges, with no loop, parallel pair or triangle, because each would give a dependency of size at most three. A triangle-free simple graph on five vertices has at most six edges. For example, if the maximum degree is at most two there are at most five edges. At degree three its three neighbors are independent, leaving at most three further edges through the fifth vertex. At degree four the graph is a star. This contradiction survives arbitrary changes of row basis. The optimum on this particular instance is nine for a budget of seven.

## A linear number of exceptional columns

A natural repair is to delete a small set of exceptional columns, enumerate their values, and use a graphic parity system on the rest. A connected satisfiable family already forces at least k deletions on 17k variables and 14k positive triples. Each variable appears two or three times.

**Construction.** Take k blocks arranged in a ring. Each starts with its own seven Fano variables. Replace each ordered Fano line (a,b,c) with two clauses (a,b,t) and (c,s,t). Use a fresh link t for each line. In each block, the seven slots for s in cyclic line order use the preceding boundary connector, a local connector three times, a second local connector twice, and the following boundary connector. Neighboring blocks share their boundary connector. At k = 1 both boundary slots name the same variable. Per block this uses seven original variables, seven links, two local connectors and one boundary connector. Setting all links true and all other variables false gives an explicit witness. The Fano incidence graph is connected, splitting a line leaves its two parts joined by t, and the boundary connectors join successive blocks.

**Minor certificate.** In one block keep just its seven original columns and its seven links. Delete all other columns and discard zero rows. Contract each link column. Over the binary field this means adding the two incident rows, using one as the pivot row, then dropping the pivot row and column. The surviving seven rows are exactly the Fano incidence matrix. A minor of a graphic column matroid is graphic, since deletion and contraction correspond to deleting and contracting graph edges. The seven-clause obstruction therefore certifies that these fourteen columns are nongraphic even in isolation.

The k sets of fourteen columns are pairwise disjoint. Deleting fewer than k columns leaves one entire set intact, including its nongraphic minor. Consequently a universal logarithmic bound on the size of this deletion set is impossible even under the promises above. This argument concerns a fixed matrix before value-dependent simplification. The explicit satisfying assignment illustrates how easily a decision procedure can bypass that representation requirement.

## Can separate linear feasibility tests replace optimization?

A second candidate checks the real equations with box bounds and also checks the unrestricted integer equations. Every positive ternary instance passes the first test at xᵢ = 1/3. Passing the second test implies consistency modulo every modulus. We found a cubic, pairwise-linear instance that passes both and has no Boolean solution. Here each row gives the indices of a clause.

```text
0 6 9     0 7 8     0 10 11
1 6 7     1 8 11    1 9 10
2 3 7     2 4 11    2 5 10
3 4 6     3 5 9     4 5 8
```

An integer solution is (−1,−1,0,0,0,0,1,1,1,1,1,1). For a short Boolean impossibility proof, suppose x₀ = 1. Then all of x₆,…,x₁₁ vanish by the first three clauses. The last six clauses force every pair among x₂,…,x₅ to sum to one, which is impossible already on a triangle. If x₀ = 0, the first three clauses select three of x₆,…,x₁₁. In the last six clauses, each of x₂,…,x₅ occurs three times and each of x₆,…,x₁₁ once. Summing gives 3(x₂+x₃+x₄+x₅) + 3 = 6, so exactly one of the four is true. In the last six clauses the three pairs avoiding that true variable force three of the six high-index variables true. For each possible true variable 2, 3, 4, 5 these triples are respectively {6,8,9}, {8,10,11}, {7,9,10}, {6,7,11}. Each contains an entire forbidden pair from either the first three clauses or the next three. The next three also require x₁ = 0 because the six high-index variables include a true one. This exhausts both cases. The executable check independently tests all 4,096 assignments and finds minimum parity cost 24 against budget 12.

The same two relaxations fail on a family that also survives linearly many levels of local consistency. [Family and proof](local-relaxation.md)

## Verification and scope

- Exhausted all 1,024 distinct-clause formulas over five variables, including the empty formula. Compared direct Boolean search with parity minimization on every one. Compared the graph solver to direct search on 491 with a tested graphic representation, including satisfiable and unsatisfiable cases whose original columns exceed weight two.
- Enumerated all 1,365 choices of four distinct nonzero row vectors from the rank-four Fano row space, rejecting dependent choices and testing every basis. Checked rank, triangle exclusion, the optimum of nine, and the planted construction for k = 1 through 12. The minor proof above covers every positive integer k.
- Related work already gives strong lower bounds for a flat disjunction of affine systems, which would preserve all models. The decision algorithm here optimizes a single weighted parity system. Its frontier is the optimization on nongraphic matrices that arise from triples, especially when preprocessing can depend on the assigned values. We claim no literature priority for these deductions. [Cohen and Shinkar](https://eccc.weizmann.ac.il/report/2014/099/) · [Reproduction script](../experiments/check_parity.py)

The missing step for P = NP is a uniform, worst-case polynomial algorithm on **every** instance of this NP-complete language, including the nongraphic cases. The obstruction imposes no time bound on other algorithms. [Moore and Robson](https://arxiv.org/abs/math/0003039) · [Clay problem statement](https://www.claymath.org/millennium/p-vs-np/)
