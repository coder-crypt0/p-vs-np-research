# Computing Taylor coefficients of structured permanents

The October 2026 SAT algorithms change the baseline for this project. OpenAI's manuscript gives deterministic subexponential exact counting for fixed-width CNF. Its Lean endpoint checks a slightly weaker subexponential bound. Anthropic's independent tensor-contraction result also gives a subexponential bound. For our positive one-in-three inputs, replace each triple by its OR clause and the three binary clauses forbidding pairs of true variables. This preserves the count with no new variables. A worst-case polynomial algorithm is still the target. [OpenAI manuscript](https://github.com/openai/math/blob/7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/paper.pdf) · [OpenAI Lean endpoint](https://github.com/openai/math/blob/7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a/lean/OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean) · [Anthropic manuscript](https://hexagonmath.org/2610.00577)

Hatami's exposition of the sparse-permanent route describes a polynomial matrix A(t) whose value at each sample point is block diagonal plus low rank. Hermite interpolation recovers the desired permanent from Taylor coefficients at the sample points. The generic derivative algorithm enumerates row and column deletions. We test a way to replace this enumeration by dynamic programming and determine what rank information it needs. [Sparse-permanent exposition, Sections 2–5](https://www.cs.mcgill.ca/~hatami/expositions/ETH.pdf)

## Sum all coefficients together when the perturbation is structured

Work modulo z^h. Suppose

```text
A(z) = B(z) + ∑ᵢ z^{dᵢ} uᵢ(z) vᵢ(z)ᵀ.
```

Here B(z) is supported on disjoint row and column blocks, and dᵢ is a nonnegative integer order for each outer product. Extend Hatami's block-and-rank dynamic program over truncated polynomials. Label each row and column either as retained in its block or as assigned to an outer-product channel. In each block, multiply the selected vector entries by the permanent of the minor on the retained rows and columns. Pairing rows and columns in a channel contributes a factorial of the common channel count. Convolve the local tables over all blocks using the left and right counts as keys.

Let α and β be the accumulated count vectors and put W(α) = ∑ᵢ dᵢ αᵢ. Discard every state with W(α) ≥ h or W(β) ≥ h. Counts only increase, so such a state cannot contribute to degree less than h once the vectors match. At the end keep α = β and multiply each term by z^{W(α)} and the channel factorials. Block entries and vector entries can have arbitrary higher-order terms, since their products are also truncated. This computes all coefficients up to degree h − 1 at once.

If there are r₀ channels of order zero and r₊ of positive order, and at most K uses of positive-order channels fit under h, the number of possible count vectors is at most (n+1)^{r₀} times the binomial coefficient choosing K from r₊+K. On each side of a block of size at most s there are at most (r₀+1)^s times the sum of binom(s,j)r₊^j over j ≤ K admissible labelings. Minor permanents use a subset dynamic program of size polynomial in s times 2^s. This yields polynomial time for constant r₀ and K, s = O(log n), and polynomially bounded h and r₊. Over integer coefficients the bit lengths also grow polynomially in these parameters and the input bit lengths. [Implementation and independent finite checks](../experiments/check_permanent_jets.py)

## Why the first derivative obstructs a direct substitution

In Hatami's translation-matrix deformation let G = (F₂^k)^ℓ, let γ be the translation, and let p be the nonconstant approximating polynomial of degree d. At sample t = i, the first derivative of the Fourier multiplier at frequency ξ is

```text
−(1/i) zᵢ p′(S)
zⱼ = ξⱼ · γⱼ in {0,1}
S = ∑_{j≠i} (1 − i/j) zⱼ.
```

The time-independent low-rank correction has derivative zero. Assume γᵢ ≠ 0. Partition the i-th coordinate's symbols into pairs {a,a+γᵢ} and choose one symbol of each pair for the rows and the other for the columns. The off-block submatrix between these sets vanishes between different pairs. The Fourier transform of zᵢ on that coordinate is supported at 0 and γᵢ, with value −1/2 at γᵢ in the normalized convention. For each pair the remaining matrix is convolution on (F₂^k)^{ℓ−1}, diagonalized by its characters with eigenvalues p′(S)/(2i). Thus the off-block rank is exactly

```text
|G|/2 × Pr[p′(S) ≠ 0].
```

In particular, when m other coordinates have γⱼ ≠ 0, their independent fair bits give m nonzero coefficients with distinct signed values. At least half have the same sign and therefore distinct absolute values. Condition on the rest and apply the distinct-coefficient Littlewood–Offord bound. Every atom has probability O(m^(−3/2)), so the rank is at least |G|/2 times 1 − O((d−1)/m^(3/2)). For shifts with m proportional to ℓ and the parameters in the exposition, this approaches half the whole dimension. Subtracting any matrix supported in the original coordinate blocks leaves this submatrix unchanged. Even if the channel vectors vary with z, the first coefficient of a rank-one polynomial term of order zero has rank at most two, and a term of order one contributes rank at most one. Hence one needs a linear number of channels active through first order for this building block under the same partition. [Translation deformation, Lemma 3.2](https://www.cs.mcgill.ca/~hatami/expositions/ETH.pdf) · [Sárközy and Szemerédi](https://doi.org/10.4064/aa-11-2-205-208)

The exact-rational experiment compares the rank formula with direct elimination on the Fourier-inverted matrices, including partial-rank and full-rank cases. The lower bound concerns an individual translation and its sampled block partition. For a sparse matrix assembled from restricted translations, rank could drop under restrictions or cancellation. A faster derivative algorithm would need to use that additional structure, a different partition, or a different deformation. [Rank experiment](../experiments/check_permanent_jets.py)
