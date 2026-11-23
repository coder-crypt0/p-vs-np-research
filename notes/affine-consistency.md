# Local rational and integer consistency together

The separate feasibility tests can be strengthened by keeping signed integer weights for all local assignments and allowing support only where a nonnegative rational relaxation permits it. The construction below passes this stronger test through a linear number of levels. It adapts Lichter and Pago's coprime-denominator argument to positive one-in-three SAT. [Limitations of Affine Integer Relaxations for Solving Constraint Satisfaction Problems](https://doi.org/10.4230/LIPIcs.ICALP.2025.166)

## Two contradictory systems and a switch

Use two copies of a connected cubic expander G with n vertices and edge expansion h. Orient each edge and use the outflow-minus-inflow equations at vertices, over the fields of order 2 and 3 respectively. In each copy assign charge 1 to one vertex and 0 to the rest. Both systems are impossible since all edge terms cancel in the sum of the equations. Add a Boolean switch which selects one of the two systems. Each value of the switch enforces its corresponding system and leaves the other unconstrained, so the disjunction is also impossible.

Encode the order-two edge values by bits and the order-three values by a triple of indicators constrained to contain exactly one true bit. For a vertex equation and each violating tuple of field values, make a clause that forbids that tuple unless the switch selects the other system. A forbidden value of a binary edge uses one literal. A forbidden value in the order-three case uses the negation of its indicator. These guarded clauses have four literals. There are a constant number per vertex. The one-hot indicator constraints apply in both modes.

Reduce an OR of four literals to two ORs of three by using (a ∨ b ∨ t) and (¬t ∨ c ∨ d), with the canonical choice t = ¬(a ∨ b). Use the [OR encoding](local-relaxation.md) for these ternary ORs. Complement pairs use an exact-one clause containing a forced zero, with the zero supplied by the same six-line Fano construction. All target clauses are positive exact-one constraints on three distinct variables. For each clause there is a deterministic extension from any satisfying source tuple for that clause. Every output variable is a function of the switch and at most three edge values. The switch may occur many times in this encoding. The number of output variables and clauses is O(n).

## Construct the two local laws

Fix the switch to one mode and the inactive system's edge values to zero. For any S of fewer than hn/2 queried active edges, choose a vertex in the unique component of G − S of size greater than n/2. Subtract one from its charge and solve the corrected system. Its charges sum to zero. A solution exists by routing each non-root vertex's charge to a root along a path, with signs determined by orientation.

The corrected solution set is a translate of the kernel of the incidence map, so its size is a power of the field order. Sample it uniformly and project to S. Moving the defect within the large component translates this set by a signed unit flow along a path avoiding S, so the local law is independent of the defect location within that component. The component for a larger queried set lies within the one for a smaller set. This gives marginal consistency. The three edges at any one vertex isolate that vertex in G − S, so its equation holds under the local law on those three edges.

For any r with 3r < hn/2, evaluate the encoding's deterministic functions on the source laws to get target local distributions μ₂ and μ₃ through width r. Each output clause is valid in its marginal. For a guarded clause in the selected system, its gadget auxiliaries depend on all three edges of its source vertex, and thus their law enforces that equation. The other system's guarded clauses are true under the switch. The one-hot, complement and constant constraints hold pointwise. Every entry of μ₂ has denominator dividing 2 raised to |E(G)|. The analogous denominator bound for μ₃ uses 3.

## Satisfy the signed integer test on allowed support

Write A = 2 raised to |E(G)| and B = 3 raised to |E(G)|. They are coprime, so choose integers u,v with uA + vB = 1. On every scope use the two tables

```text
rational probability = (μ₂ + μ₃) / 2
signed integer weight = uA μ₂ + vB μ₃
```

The first table is nonnegative and the second has integer entries. Both have total weight one and agree under marginalization. Every tuple to which the signed integer table assigns nonzero weight has positive probability in the first table, and every violating tuple has zero weight in both. The width-r local LP and the signed integer relaxation restricted to its feasible support therefore both accept. Taking signed first moments gives simultaneous feasibility of every affine integer equation valid on the retained local supports. Reducing the signed tables in any field gives the same result for affine equations over that field.

This closes the variant that strengthens local LP by adding an independent global consistency check on affine consequences from bounded scopes, including arbitrary fields obtained by reduction of these signed tables. A different variant conditions on specific tuples before asking for the signed integer tables. Fixing the switch pins the field, so the construction above supplies no lower bound for that variant. [Comparison with cohomological consistency](https://doi.org/10.4230/LIPIcs.ICALP.2025.166)
