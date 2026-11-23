# Local linear relaxations with a global integer check

Positive one-in-three SAT can be phrased as finding a Boolean vector satisfying a system of integer equations Ax = 1. The construction below gives unsatisfiable instances that pass both signed integer feasibility and a local linear relaxation. Each variable occurs at most three times, and any two clauses share at most one variable.

## Statement

Fix a connected simple cubic graph G with n vertices and edge expansion h. Every set U of at most n/2 vertices has at least h|U| edges to its complement. For any odd-total binary charges on its vertices, construct a positive one-in-three instance with 144n variables and 112n clauses.

1. The instance is unsatisfiable.
2. The raw equations have an unrestricted signed integer solution. Every integer choice of the original edge variables extends to one.
3. For every integer r ≥ 3 with r < hn/6, there is a consistent family of probability distributions on all subsets of at most r variables. Each clause is satisfied with probability one in the distribution on its three variables.

Property 3 is the scope-size formulation of the Sherali–Adams relaxation for Boolean CSPs. Constant-expansion cubic graphs exist for infinitely many n, giving instances feasible through linearly many levels under that convention. Signed integer feasibility also implies raw feasibility modulo every modulus. [Kolesnik and Wormald on expansion](https://arxiv.org/abs/1311.6555) · [Local consistency as a CSP reduction](https://arxiv.org/abs/2301.05084)

## Translate the parity constraints

Put one bit on each edge. At vertex v the XOR of its three bits must equal its charge. Summing all equations over the binary field gives zero on the left because each edge has two endpoints, and one on the right. To express the equation at one vertex as ordinary clauses, forbid each of the four bad triples with one OR clause.

For an OR with literal values ℓ₁, ℓ₂, ℓ₃, use six auxiliaries pᵢ and sᵢ. Require exactly one of p₁, p₂, p₃ and impose the three further exact-one equations

```text
pᵢ + (1 − ℓᵢ) + sᵢ = 1.
```

Over Boolean inputs and auxiliaries, the selected literal must be true. Conversely, select the first true literal and set sᵢ = ℓᵢ − pᵢ. The auxiliary values under this rule are functions of the three input edge bits. Over signed integers, choose p = (1,0,0) for every input and set sᵢ by the same difference.

Each edge needs four occurrences of its bit and four of its complement. Use eight copies in a cycle, alternating the bit and its complement. For each consecutive pair cᵢ,cⱼ, use a clause (cᵢ,cⱼ,z) with a fresh port z forced to zero.

To force zero, take the six Fano lines other than {0,1,3}. Every retained line contains two of the four points outside the deleted line. Those six pairs exhaust the edges of a complete graph on four points, so at most one outside point can be true. If an outside point were true, it would occur on one of the two retained lines through each missing-line point, forcing these missing-line points false. The other line through any one of them would then have no true points. Thus all outside points are false and all three missing-line points are true. Add the clause (0,z,w) with two fresh variables to force z and w to zero.

Every copy belongs to two cycle clauses and one OR gadget. Each port gadget uses nine fresh variables and seven internal clauses. The cycle clauses share at most one variable per pair and their ports are distinct. Fano lines intersect pairwise in one point, and each attachment meets its Fano lines only through point 0. Within an OR gadget, different clauses share at most a selector, and each data copy appears in one OR gadget only. This verifies the occurrence and overlap bounds. The construction uses 80 variables and 64 clauses per edge plus 24 variables and 16 clauses per vertex. A cubic graph has 3n/2 edges.

Every Boolean solution would induce consistent original edge bits and satisfy the contradictory parity system. Every signed integer edge assignment extends by alternating the given value and its integer complement around each cycle. The port variables take their Boolean solution and the signed OR extension completes the rest. [Lean checks for the gadgets](../proofs/LocalGadget.lean)

## Construct the local distributions

For a set S of fewer than hn/2 edges, the graph G − S has a unique component containing more than n/2 vertices. If all components had at most n/2 vertices, summing their expansion inequalities would give at least hn crossing edge incidences. Every crossing edge belongs to S and contributes two incidences, contradicting 2|S| < hn.

Choose a vertex v inside this large component and flip its charge. The resulting even-total parity system has a solution. Pair the vertices of charge one and take the symmetric difference of paths between paired vertices. Use the uniform distribution on all solutions to the corrected system and project it onto S. The projection is independent of the choice of v inside the large component. For two choices v and w, toggling a path in G − S between them bijects the corrected solution sets and leaves all S bits unchanged.

If S lies in S′, the large component for S′ lies inside the large component for S. Use the same defect vertex for both to see that the distribution for S is the S-marginal of the distribution for S′. When S contains the three edges of a vertex, that vertex is isolated in G − S. It cannot be the defect vertex, so its parity equation holds in the projected distribution.

Every output variable depends on at most three source bits. A scope U of at most r output variables therefore uses a union S of at most 3r source edges. Evaluate the local functions of the construction under the distribution on S to obtain the output distribution for U. Marginal agreement follows from the source distributions. Every clause of an OR gadget includes an auxiliary whose scope contains all three incident edges, so the parity equation holds and the first-true selector works. Port and cycle clauses hold directly for all values of their source edges.

The Python experiment checks 512 cases of the OR gadget. For the lift of the complete graph on four vertices, it checks 42 source scopes, 233 marginal identities, all 448 output clause distributions and a global signed witness. The asymptotic statement follows from the argument above. [Reproduction script](../experiments/check_local_relaxation.py)

## Consequence for the algorithm search

Combining the two feasibility tests with an AND operation still accepts this family. The translated formulas retain a route to a short refutation. Project each constant-size vertex gadget onto its three input bits, extract its XOR relation, and use global Gaussian elimination. This uses relations derived from Boolean feasibility across several clauses. A general solver based on that idea would need a polynomial rule for finding and combining enough relations on unrestricted inputs.
