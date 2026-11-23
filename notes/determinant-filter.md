# Exact SAT counts as a constant coefficient of a sparse determinant

The local parity obstruction for integer permanent gadgets applies when every mismatched port deletion cancels inside its gadget. Here we allow mismatches, record them using Laurent monomials, and discard them by extracting a coefficient. The result is a count-preserving reduction from 3-SAT to the constant coefficient of a determinant. Every row and column has at most three nonzero entries, and every nonconstant symbol appears in only two entries, once with exponent +1 and once with exponent −1. This is a coefficient problem on an efficiently evaluable determinant.

## Local determinants

Use one matrix index per port, with the row and column at that index treated as a pair. When a port is selected, delete its row and column. The following blocks have principal deletion determinants equal to the specified Boolean indicators.

- **Variable equality.** For k occurrences, use the matrix of a single directed k-cycle. Give its final edge weight (−1)^(k−1), so the full cycle contributes +1. Every proper nonempty subset of the cycle's indices has determinant zero, while the empty determinant is one. The supported port choices are therefore all zero and all one. A single occurrence uses [1].
- **Clause OR.** Use the identity plus a directed three-cycle whose edge weights multiply to −1. With all ports unselected the identity and the directed cycle contribute +1 and −1. Deleting any nonempty set breaks the directed cycle, leaving determinant one.
- **Literal negation.** Use the following matrix with ports at its first two indices. The first two rows are equal, so its full determinant is zero. Deleting either one of the ports gives determinant one, and deleting both leaves a zero.

```text
 0  0  1
 0  0  1
-1 -1  0
```

A source occurrence joins its variable and clause ports directly for a positive literal. For a negative literal insert a two-port negation block in the wire. These choices specify exactly the satisfying assignments to the occurring variables. All local determinants on supported masks are one.

## Connect the blocks and filter

Start with the block diagonal matrix. For a wire e between ports u and v, add the entries −zₑ at (u,v) and zₑ^−1 at (v,u). Each port receives only one wire. In every determinant term, the zₑ exponent is zero exactly when the term uses both crossing edges or neither. The product of a selected pair is −1, and its transposition sign is also −1. The net factor is +1, independently of all other choices. Once all paired crossing choices are fixed, the remaining determinant products are exactly the principal minors of the gadgets. Consequently the constant coefficient is the number of satisfying assignments. There is no power-of-two factor in that coefficient.

The number of rows and wires is linear in the number of literal occurrences. Each clause row or column has at most two local entries and one crossing entry. Equality and negation have at most two total per row or column. The entries that carry symbols have unit coefficients with signs. This supplies an exact signed-matrix representation for the zero decision we need.

## An exact algorithm in the cycle rank

For each determinant term, regard wires used in only one direction as directed edges between gadgets. Each square gadget loses the same number of rows and columns to crossing entries. Thus unpaired directions form a circulation. Choose a spanning forest of the gadget multigraph. If all chords have exponent zero, any remaining circulation would be supported on a forest. A leaf of a nonempty forest cannot balance a nonzero edge, so the circulation must vanish.

Set every forest wire variable to one and average the determinant over the sign assignments in {−1,1} to the chords. Each exponent is in {−1,0,1}, so the independent sign averages discard exactly the terms with a nonzero chord exponent. With cycle rank b the algorithm uses 2^b polynomial-time determinant evaluations. It is polynomial when b = O(log n). Inserting literal negations only subdivides wires, so b is the cycle rank of the source incidence multigraph. For m ternary clauses on n occurring variables in c connected components, b = 2m − n + c.

This shifts the hard step to extracting a particular coefficient. A single-variable substitution zₑ = t^{wₑ} cannot generically isolate that coefficient with small exponents if the method requires separating every nonzero circulation by its total weight. For k edge-disjoint cycles, each subset of the cycles has a total weight. If two subsets agree, their difference gives a nonzero signed circulation of weight zero. If all integer weights have magnitude at most W and the graph has m′ wires, there are at most 2m′W+1 possible totals. Thus a separation guarantee for every circulation needs 2m′W+1 ≥ 2^k on these graphs. This bound concerns the substitution guarantee. A coefficient algorithm may use the structure of the actual determinant and its cancellations.

The finite check compares the full sign filter on smaller formulas and the chord filter on all examples with direct enumeration. It covers empty, unsatisfiable, repeated-literal and mixed-sign inputs. [Implementation and checks](../experiments/check_determinant_filter.py)

For a fixed connected wiring graph, an exact sign-evaluation oracle can require all (2^b) chord signs even to distinguish a unique satisfying assignment from an unsatisfiable formula. [Query lower bound](sign-filter-limits.md)
