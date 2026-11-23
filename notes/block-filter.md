# Count across articulation variables

The [bridge method](determinant-filter.md) can evaluate a graph with many small cyclic regions separately. It needs a single edge between regions. This extension conditions on articulation variables, so it can also split bridgeless inputs. The decomposition uses standard graph structure. The determinant calculation supplies the weighted count inside each part. [Biconnected decomposition for constraint solving](https://ijcai.org/proceedings/2001/41)

## Put weights on variable values

For each variable $x$, take two integer weights $w_x(0),w_x(1)$. The desired count sums the product of the selected weights over satisfying assignments. The determinant compiler realizes these weights by changing only two entries per variable.

- If $x=0$, the local equality block is kept in full. Multiply one edge of its directed cycle by $w_x(0)$.
- If $x=1$, the equality block is entirely deleted. Multiply one of the two directed crossing entries of its first wire by $w_x(1)$. The paired transposition still has net positive sign, and this wire is selected exactly in this case. With a negation relay, use the wire on the variable side.

The constant coefficient therefore gives the weighted count. In particular, weights $(1,0)$ and $(0,1)$ pin a variable. Zero weights still make sense because the algebra sums over all paired states and the unsupported one contributes zero. In bridge composition each selected bridge's factor is included once, when its child contribution is applied at the parent port.

## Put complete clauses in the blocks

Make the primal graph on variables and join two distinct variables if they occur in the same clause. Every clause's variables form a clique, so all edges of a clause with at least two distinct variables lie in one biconnected block. Assign that clause to the block. Treat clauses with just one distinct variable as unary weights. This also accounts for repeated literals and tautologies.

Form the block-cut forest with a node for each block and each variable. A variable is joined to the blocks that contain it. Non-articulation variables are leaves or isolated nodes. Root each tree at a variable. A variable sends a message with two entries by multiplying its unary weights with its child-block messages separately for each bit. A block sends a message back by counting its assigned clauses twice using the child-variable messages as weights, once for each value of the parent variable. Sum the two entries at each root and multiply across connected components.

Correctness follows by induction on this forest. Fixing a variable's bit separates the assignments in its child blocks. Each block count sums over its own variables, with the incoming messages accounting for assignments in their disjoint subtrees. No clause is split across blocks and each unary factor is counted once.

If $b_B$ is the incidence cycle rank of the clauses in block $B$, this uses at most $2\sum_B 2^{b_B}$ determinant evaluations before any further bridge decomposition internal to the block. A single large biconnected block can still have linear cycle rank.

## A bridgeless check

For $i=1,\dots,b$, use the clauses $(x_i\lor y)$ and $(\neg x_i\lor y)$. The incidence graph has no bridges and cycle rank $b$. Its cycles meet at the articulation variable $y$. The formula requires $y=1$ and leaves all $x_i$ free, giving $2^b$ solutions. Complementing $y$ in just one pair makes the count zero without changing the graph. The checks cover both versions at $b=48$ and compare the weighted and block algorithms with direct enumeration on smaller mixed-width formulas. [Implementation and checks](../experiments/check_block_filter.py)
