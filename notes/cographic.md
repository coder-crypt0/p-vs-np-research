# Exact-one SAT on graph cuts

The [graphic representation](parity-matching.md) makes the parity optimization a minimum T-join, which is polynomial-time solvable. In the dual representation, the same objective becomes a maximum cut. Here the exact-one decision question stays NP-complete even when the graph is supplied, deleting one universal vertex leaves a planar graph of maximum degree six, and the clauses are exactly its triangles.

## The cut formulation

Let H be a connected simple graph. Give every edge one Boolean variable and every triangle an exact-one clause on its three edges. Suppose the triangles span the binary cycle space. The parity system then has solutions of the form x = 1 + δ(S) over the binary field, where δ(S) is the edge cut between S and its complement. The all-ones vector satisfies every triangle equation, and the difference between any two parity solutions is orthogonal to the cycle space. Thus the nullspace is the cut space, and the triangle matrix represents the cographic matroid of H.

The vector x marks the edges whose endpoints are on the same side of the cut. A triangle has exactly one such edge if its vertices use both colors. If dₑ is the number of triangles containing edge e and there are m triangles, the original cost is 3m minus the cut weight for weights dₑ. The exact decision asks whether the maximum weight is 2m. Each triangle contributes at most two cut edges, so this is the exact threshold. The graph problem is vertex two-coloring without a monochromatic triangle. Its unrestricted form was already known to be NP-complete. The claim here is for a supplied graph of the form below. [Karpiński and Piecuch](https://arxiv.org/abs/1710.07132)

## Reduction with a planar base and a universal vertex

Start from positive one-in-three SAT with three distinct variables per clause and at most three occurrences per variable. NP-hardness already holds when the bipartite incidence graph is planar and cubic. [Moore and Robson](https://arxiv.org/abs/math/0003039)

1. Give each occurrence its own vertex. For each input clause fill in the triangle on its three occurrence vertices.
2. For each variable, link consecutive occurrence vertices p and q using a diamond. Add fresh vertices u,v and the two triangles (p,u,v) and (q,u,v). There is no edge pq. The exact-one conditions on the two vertex triples imply p = q. When both ports are selected, choose neither interior vertex. When both ports are unselected, choose either one interior vertex.
3. Add a universal vertex a to this base graph. Make one target clause for every triangle of the resulting graph H, using edge variables.

The only base triangles are the ones in steps 1 and 2. A triangle on occurrence vertices lies inside one input clause. A triangle containing a diamond interior vertex uses its other interior vertex and one of the two ports, since its two ports are nonadjacent. Every other triangle of H is (a,u,v) for a base edge uv.

Under a two-coloring, call a base vertex selected if its color agrees with a. The apex triangles forbid adjacent selected vertices. Each base triangle must use both colors, so it has exactly one selected vertex. The diamonds enforce occurrence equality, and the clause triangles then satisfy the original formula. Conversely, a satisfying assignment extends to selections in every diamond. Every base edge lies in a designated base triangle, so these selections also contain no adjacent pair. The resulting coloring splits every triangle of H.

The apex triangles form a fundamental cycle basis for the spanning tree consisting of the edges incident to a. The construction has O(m) size for m input clauses. An occurrence vertex has at most two neighbors from its clause and four from its diamonds. An interior vertex has degree three. A base edge appears in at most two base triangles and one apex triangle. Each spoke appears in at most six apex triangles. In particular every variable occurs in at most six target clauses.

For planarity of the base, take a planar drawing of the input incidence graph and replace each node by a disk and each edge by a disjoint corridor. Place one occurrence vertex at the variable end of each corridor. Inside a clause disk, draw the three sides of a triangle using two parallel arcs per incident corridor. Inside a variable disk of degree at most three, join consecutive occurrence vertices with diamonds arranged along the disk boundary. Any chosen order of three ports uses two of the three boundary intervals, so the two diamonds fit without intersecting. This gives a planar drawing after removing only a. Membership in NP follows by checking the edge assignment.

For comparison, a planar graph admits a polynomial-time maximum cut algorithm. The planar-plus-one-vertex hardness of unrestricted weighted Max-Cut was already known. This reduction enforces both the triangle occurrence weights and the precise upper-bound threshold 2m. [Barahona et al.](https://www.zib.de/userpage/groetschel/pubnew/paper/barahonagroetscheljuengeretal1988.pdf)

The finite checks compare source and recovered assignments, check the cycle-space rank and all graph triangles, and verify planarity of the base for tested planar source graphs. Lean checks the truth-table equivalences used for the triangles and diamonds. [Experiment](../experiments/check_cographic.py) · [Lean proof](../proofs/Cographic.lean)
