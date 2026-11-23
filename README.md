# Independent P versus NP research

We investigate deterministic algorithms for NP-complete decision problems. The first experiment uses positive one-in-three SAT. Every input clause names three distinct Boolean variables and requires exactly one of them to be true. This language remains NP-complete on cubic planar instances. [Moore and Robson](https://arxiv.org/abs/math/0003039)

## Formal target

DeepMind's Lean statement is `ComplexityTheory.P_ne_NP : P ≠ NP`. At revision `2207c54d` it is tagged open and its proof body uses `sorry`. The imported definitions model decisions on Boolean lists, deterministic polynomial time with encoded Turing machines, and NP with polynomially bounded witnesses. A positive resolution would prove the opposite equality. We use the definitions as a reference and keep our own lemmas independent of the conjecture declaration. [Pinned statement](https://github.com/google-deepmind/formal-conjectures/blob/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjectures/Millennium/PvsNP.lean) · [Class definitions](https://github.com/google-deepmind/formal-conjectures/blob/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjecturesForMathlib/Computability/Complexity.lean)

## Research record

- [Public fleet snapshot and overlap](notes/fleet.md)
- [One-in-three SAT, parity, and graph matching](notes/parity-matching.md)

## Reproduce

From the repository root, run the following with Python 3.10+ and Lean 4.33.1 on your path. The Python checks use the standard library. The Lean file imports only the bundled standard library.

```sh
python experiments/check_parity.py
lean -DwarningAsError=true proofs/ExactOne.lean
```

The Python experiment checks finite instances and the explicit obstruction. The Lean file proves the parity-cost criterion for arbitrarily many clauses. The graph-minor argument is written in the research note.
