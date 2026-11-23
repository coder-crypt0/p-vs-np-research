# Modular determinant factors for three-variable separators

The real-matrix obstruction at three ports can be bypassed in the quadratic extension of an odd prime field. Every three-variable table over the prime field has a four-row determinant factor using one wire per variable. The construction fails for some four-variable tables over every field.

## Three ports in a quadratic extension

Let $p$ be an odd prime and let $w$ be the table of boundary counts modulo $p$. A table that vanishes modulo $p$ can be represented by a zero auxiliary row. Otherwise choose a nonzero entry $h$. Complement the boundary bits as needed so this entry corresponds to deleting all ports. The literal relays implement these complements and subdivide incidence edges. Index the resulting table by the subset $K$ of surviving ports and divide by $h$, giving values $m_K$ with $m_\varnothing=1$.

We seek a $3\times3$ matrix $S$ with these principal minors over $\mathbb F_{p^2}$. Put $S_{ii}=d_i=m_{\{i\}}$ and set the required pair products to $r_{ij}=d_i d_j-m_{\{i,j\}}$. The two oriented triple products must satisfy

$$
u+v=m_{\{1,2,3\}}-d_1d_2d_3+d_1r_{23}+d_2r_{31}+d_3r_{12},
\qquad uv=r_{12}r_{23}r_{31}.
$$

They are roots of a quadratic over $\mathbb F_p$. Every discriminant in $\mathbb F_p$ is a square in $\mathbb F_{p^2}$. In particular, for a fixed nonsquare $q$, write the extension as $\mathbb F_p[t]/(t^2-q)$. A nonsquare discriminant $\Delta$ has a square root $t\sqrt{\Delta/q}$ in that extension. All local factors for the same prime therefore use the same degree-two extension.

If one triple root is nonzero, call it $u$. Put $(S_{12},S_{23},S_{31})=(u,1,1)$ and $(S_{21},S_{32},S_{13})=(r_{12}/u,r_{23},r_{31})$. The pair products and both oriented triple products then have their required values. If both roots vanish, some $r_{ij}=0$. On the forward cycle put a $1$ for each nonzero pair product and $0$ for each zero pair product. On the reverse cycle put the pair products themselves. This realizes both zero triple products. A fourth row and column with a single diagonal entry $h$ restores the original scale.

Signs still filter the Laurent coefficient. Characteristic is odd, so averaging divides by a power of $2$. The resulting coefficient belongs to $\mathbb F_p$ because the paired-wire expansion gives a sum of products of the original counts modulo $p$. Iterating over odd primes and taking their product past $2^n$ recovers the exact count for $n$ variables by the Chinese remainder theorem. The prime computation and reconstruction cost polynomial time in $n$ in addition to the determinant evaluations. Each prime needs at most $8\sum_i 2^{b_i}+2^{b_\ast}$ evaluations, where $b_i$ and $b_\ast$ are the cycle ranks of the chosen regions and remaining core for that run.

In the example with $m$ interior variables, use $(x_i\lor a\lor b)$ and $(\neg x_i\lor b\lor c)$ for each $x_i$. The primal graph stays connected after deleting any two variables. Removing $a,b,c$ leaves one interior variable per component, with local incidence rank one. The full incidence rank is $3m-2$. Five hub assignments leave every $x_i$ free, two force every $x_i$, and the eighth is impossible, giving $5\cdot2^m+2$. Distributing all eight literal sign patterns over the hubs makes every hub assignment impossible without changing the unsigned graph. Both counts are checked at $m=32$. [Implementation and checks](../experiments/check_modular_ports.py)

## Four ports over any field

Take a Boolean table indexed by the set of surviving ports. Give value $1$ to the empty set, all six pairs, and the triple $\{1,2,3\}$. Give value $0$ to all other subsets. This is a width-three CNF boundary table. For each disallowed assignment, add a width-four clause forbidding that assignment, then split it using an auxiliary variable equal to the OR of its first two literals. The auxiliary variable is unique for each boundary assignment.

If a matrix with four ports and any number of auxiliary indices realizes the table, the auxiliary principal minor equals $1$. Its Schur complement $S$ on the four ports must have the specified principal minors. The singletons give $S_{ii}=0$ and the pairs give $S_{ij}S_{ji}=-1$.

Write $u_{ijk}=S_{ij}S_{jk}S_{ki}$ for $i<j<k$. Its product with the reverse triple is $-1$. The principal minor on the triple is the sum of the two triple products. Whenever that minor is zero, $u_{ijk}^2=1$. The relations on the pairs also give

$$
u_{123}u_{134}
=(S_{12}S_{23}S_{34}S_{41})(S_{31}S_{13})
=(S_{12}S_{23}S_{34}S_{41})(S_{24}S_{42})
=u_{124}u_{234}.
$$

The triples other than $\{1,2,3\}$ have zero minors, so squaring this identity forces $u_{123}^2=1$. Its reverse triple product is then $-u_{123}$, giving minor zero in every characteristic. The required minor is $1$, a contradiction. This gives a specific obstruction to the one-wire-per-variable construction, consistent with the general algebraic relations among four-port principal minors. [Lin and Sturmfels](https://arxiv.org/abs/0812.0601)

## Arbitrary SAT can be hidden behind four connected hubs

Take a width-three CNF with $n\ge4$ variables and $m$ clauses of three literals each. Introduce four hub variables $a,b,c,d$. The seven clauses on $a,b,c$ that exclude every assignment except $111$ force all three bits. Add $\neg a\lor\neg b\lor d$ to force $d=1$. For each original variable $x$, add $x\lor a\lor b$ and $x\lor c\lor d$. These two clauses hold under the uniquely forced hub assignment, so the extension preserves the exact count. The new clauses use distinct variables and contain no unit clause or tautology.

The primal graph contains the complete bipartite graph between the four hubs and the $n$ original variables. After deleting at most three vertices, each side still has a vertex, so the remaining primal graph stays connected. The incidence cycle rank is $2m+3n+13$. A procedure restricted to separators of at most three variables and then sign averaging on the core will enumerate exponentially many signs on this family. This is a limitation of that procedure on unsimplified inputs. The padding does not constrain other reductions or SAT algorithms.
