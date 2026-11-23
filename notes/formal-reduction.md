# A checked reduction from 3-SAT to the parity budget

The Lean compiler takes an arbitrary list of three-literal clauses, with either sign on each literal. Each source variable has two named copies, constrained to carry opposite bits. For each source clause, three selector variables choose one literal and three slack variables enforce that the selected literal is true. The complement rows share one zero variable, forced by three additional exact-one clauses. This is a semantic formalization of a standard type of reduction used to establish the NP-completeness of one-in-three SAT. [Schaefer](https://doi.org/10.1145/800133.804350)

For an input formula with m clauses, Lean checks the following statements.

| Theorem | Checked statement |
| --- | --- |
| `compile_correct` | The 3-SAT formula is satisfiable exactly when the compiled positive one-in-three formula is satisfiable. |
| `compile_length` | The compiler emits exactly 3 + 7m rows. |
| `compile_distinct` | Each emitted row contains three distinct variable names, even if an input clause repeats a literal or a variable. |
| `compile_parity_budget` | The 3-SAT formula is satisfiable exactly when the compiled rows have an odd-parity assignment with total true occurrences 3 + 7m. |

A satisfying 3-SAT assignment extends by selecting the first true literal in each clause and filling its slack variables. An exact-one assignment recovers the source assignment from the positive literal copies. Whole clauses serve as keys for their auxiliary variables, so repeated source clauses can share selectors and slacks. The proof allows repeated output rows, which preserves satisfiability. [Lean implementation](../proofs/ReduceThreeSAT.lean)

The reduction reads each encoded source clause and emits a constant number of tagged copies of its literal names and clause key. This gives linear output length with self-delimiting variable names and fixed tags. A deterministic polynomial-time decision algorithm for reaching the budget on these inputs would decide 3-SAT in polynomial time. The checked file proves the semantic equivalence and clause count. A full Lean proof of P = NP would also need the algorithm with a machine-level time bound, the compiled bit encoding, and a composition proof under the target complexity-class definitions. Existing Lean Cook–Levin developments provide reusable interfaces, including `threeSAT_InP_implies_PeqNP` in a separate TM0 formalization. DeepMind's pinned statement uses Mathlib's TM2 polynomial-time definition. [Cook–Levin interface](https://github.com/starpos/lean-np-complete/tree/3817bf9a8292062863a24110cf648ae6c05a3429) · [DeepMind class definitions](https://github.com/google-deepmind/formal-conjectures/blob/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjecturesForMathlib/Computability/Complexity.lean)
