# Independent P versus NP research

We investigate deterministic algorithms for NP-complete decision problems. The first experiment uses positive one-in-three SAT. Every input clause names three distinct Boolean variables and requires exactly one of them to be true. This language remains NP-complete on cubic planar instances. [Moore and Robson](https://arxiv.org/abs/math/0003039)

## Formal target

DeepMind's Lean statement is `ComplexityTheory.P_ne_NP : P ≠ NP`. At revision `2207c54d` it is tagged open and its proof body uses `sorry`. The imported definitions model decisions on Boolean lists, deterministic polynomial time with encoded Turing machines, and NP with polynomially bounded witnesses. A positive resolution would prove the opposite equality. We use the definitions as a reference and keep our own lemmas independent of the conjecture declaration. [Pinned statement](https://github.com/google-deepmind/formal-conjectures/blob/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjectures/Millennium/PvsNP.lean) · [Class definitions](https://github.com/google-deepmind/formal-conjectures/blob/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjecturesForMathlib/Computability/Complexity.lean)

## Research record

- [Public fleet snapshot and overlap](notes/fleet.md)
- [One-in-three SAT, parity, and graph matching](notes/parity-matching.md)
- [Local linear relaxations with a global integer check](notes/local-relaxation.md)
- [Joint local rational and integer consistency](notes/affine-consistency.md)
- [Exact-one SAT on graph cuts](notes/cographic.md)
- [Lean reduction from 3-SAT to the parity budget](notes/formal-reduction.md)
- [Taylor coefficients of structured permanents](notes/permanent-interpolation.md)
- [Parity and valuation of permanent gadgets](notes/modular-permanent.md)
- [Exact SAT counts from a sparse determinant coefficient](notes/determinant-filter.md)
- [Limits of determinant sign queries](notes/sign-filter-limits.md)
- [Counting across articulation variables](notes/block-filter.md)
- [Counting across two-variable separators](notes/pair-filter.md)

## Reproduce

The Python experiments and Lean checks can be run from the repository root with Python 3.10+ and Lean 4.33.1. The GitHub workflow is available for manual runs.

```sh
python -m venv .venv
.venv/bin/python -m pip install --no-cache-dir -r requirements.txt
.venv/bin/python -B experiments/check_parity.py
.venv/bin/python -B experiments/check_local_relaxation.py
.venv/bin/python -B experiments/check_cographic.py
.venv/bin/python -B experiments/check_permanent_jets.py
.venv/bin/python -B experiments/check_paired_permanent.py
.venv/bin/python -B experiments/check_determinant_filter.py
.venv/bin/python -B experiments/check_sign_filter_limits.py
.venv/bin/python -B experiments/check_block_filter.py
.venv/bin/python -B experiments/check_pair_filter.py
mkdir -p .lake
lean -DwarningAsError=true --root=proofs -o .lake/ExactOne.olean proofs/ExactOne.lean
lean -DwarningAsError=true --root=proofs -o .lake/LocalGadget.olean proofs/LocalGadget.lean
lean -DwarningAsError=true proofs/Cographic.lean
LEAN_PATH=.lake lean -DwarningAsError=true --root=proofs proofs/ReduceThreeSAT.lean
```

On Windows use `.venv\Scripts\python.exe` for the virtual environment's interpreter.

The Python experiments check finite instances, graphic bases, and local distributions. Lean proves the parity-cost criterion, the 12-variable contradiction, and the Boolean and signed-integer behavior of the local gadgets. The notes give the graph-minor and expansion arguments.
