# Count across two-variable separators

Fixing two variables separates the counts on either side of a cut in the primal graph. For each of the four assignments, sum over the interior of one component. The resulting table is a binary factor. Multiplying it by the surviving factors and summing over the survivors gives the original count. This is variable elimination. The determinant-specific step is realizing an arbitrary resulting binary factor with three rows and columns and the same two wires as a binary clause. [Dechter on variable elimination](https://ics.uci.edu/~dechter/publications/r76A.pdf)

## A binary factor with no extra incidence cycles

For integer table entries $w_{ab}$ indexed by the bits of the two variables, use

$$
M=\begin{pmatrix}
0 & 0 & 1\\
-w_{00} & 0 & w_{10}\\
-w_{01} & -1 & w_{11}
\end{pmatrix}.
$$

The first two indices are ports. Delete the corresponding row and column when a port bit is $1$. The whole determinant is $w_{00}$ by expansion in its first row. Keeping the first port alone gives $w_{01}$, keeping the second alone gives $w_{10}$, and keeping neither gives $w_{11}$. Each row and column has at most three entries including the port wire. The third index requires no wires, so the gadget contributes one factor vertex of degree two to the incidence graph.

The algorithm searches the primal graph for a pair that leaves two or more components upon deletion. It selects one component of bounded incidence cycle rank, collects the factors that touch it, and evaluates the four boundary assignments using [unary pins](block-filter.md). Fixing the boundary separates those factors from the other components because every factor's scope forms a clique in the primal graph. A pin on a boundary variable absent from the local factors leaves the local count independent of that bit. Factors entirely on the boundary stay outside the selected region. If the replacement has the same scope as an existing binary factor, their tables multiply entry by entry.

Each step eliminates at least one variable. If the selected regions have incidence cycle ranks $b_1,\dots,b_t$ and the remaining core has rank $b_\ast$, the algorithm uses at most $4\sum_i 2^{b_i}+2^{b_\ast}$ determinant evaluations. Searching for eligible pairs takes polynomial time per step. The remaining core can still have large rank.

For a check, take $m$ pairs of clauses $(x_i\lor y\lor z)$ and $(\neg x_i\lor y\lor z)$ with distinct $x_i$. The primal graph has no articulation variable, and the incidence graph is bridgeless with rank $3m-1$. Removing $y,z$ leaves one component per $x_i$, each with local rank two. In the satisfiable version, each pair enforces $y\lor z$, giving $3\cdot 2^m$ assignments. Using the four sign patterns for $y,z$ across the pairs gives count zero with the same rank and separator. The checks cover both cases at $m=48$ and compare the algorithm with direct counting on smaller formulas. [Implementation and checks](../experiments/check_pair_filter.py)

## A limit at three ports over the reals

Consider a three-bit conditional table with $w_{000}=0$, $w_{111}=1$, and all six other entries equal to $2$. This table arises from a width-three CNF. Require $a\lor b\lor c$, enforce $t=a\land b$ with three clauses, and include $\neg t\lor\neg c\lor\neg d$. The interior bit $t$ is uniquely determined and $d$ is free except at $abc=111$, where it must be zero.

Suppose a real matrix realizes the table by the same principal-deletion rule with three ports and any number of auxiliary indices. Its auxiliary block has determinant $w_{111}=1$, so it is invertible. Eliminating it by a Schur complement leaves a real $3\times3$ matrix $S$ whose principal minors are exactly the table with the indices given by the surviving ports. The one-index minors give $S_{ii}=2$. The two-index minors give $S_{ij}S_{ji}=2$ for each distinct pair. Put $u=S_{12}S_{23}S_{31}$ and $v=S_{13}S_{32}S_{21}$. Then $uv=8$, while the full determinant gives

$$
0=\det S=8-3(2\cdot2)+u+v=-4+u+v.
$$

This would give $(u-v)^2=(u+v)^2-4uv=16-32=-16$. Thus no single real matrix allows arbitrary three-variable boundary tables through this port rule. The obstruction concerns one matrix with one principal deletion per boundary variable. The cycle rank of a general SAT instance remains uncontrolled.
