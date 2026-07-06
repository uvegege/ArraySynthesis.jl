# Mathematical Background

## Array Factor

For an $N$-element array with positions $\mathbf{p}_n$ (normalized to wavelengths) and complex
excitations $w_n$, the array factor in direction $\hat{r}$ is:

```math
AF(\hat{r}) = \sum_{n=1}^{N} w_n \, e^{j 2\pi \hat{r} \cdot \mathbf{p}_n}
```

Discretizing over $P$ directions and stacking into a matrix:

```math
\mathbf{AF} = A\,\mathbf{w}, \qquad A_{pn} = e^{j 2\pi \hat{r}_p \cdot \mathbf{p}_n}
```

The synthesis problem is to find $\mathbf{w}$ such that $\mathbf{AF}$ satisfies a set of constraints
(beam gain, sidelobe level, nulls), possibly optimizing some objective over the feasible set.

## Symmetric Arrays and Real Array Factor

A centrosymmetric array has positions in pairs $\pm\mathbf{p}_n$ (plus optionally the origin).
Imposing conjugate symmetric weights $w_{-n} = w_n^*$, each pair contributes:

```math
w_n e^{j2\pi\hat{r}\cdot\mathbf{p}_n} + w_n^* e^{-j2\pi\hat{r}\cdot\mathbf{p}_n}
= 2\,\text{Re}\!\left(w_n \, e^{j2\pi\hat{r}\cdot\mathbf{p}_n}\right)
```

so the array factor becomes:

```math
AF(\hat{r}) = \sum_n \Bigl[
  2\cos(2\pi\hat{r}\cdot\mathbf{p}_n)\,\text{Re}(w_n)
 -2\sin(2\pi\hat{r}\cdot\mathbf{p}_n)\,\text{Im}(w_n)
\Bigr]
```

This is **always real**. Advantages:

- Only $\lceil N/2 \rceil$ complex weights are free optimization variables.
- The modulus bound $|AF| \leq t$ collapses to $-t \leq AF \leq t$, exact and linear.

`SymmetricArray` stores only the representative half of the array;
`ConjugateSymmetricWeights` exploits this structure automatically.

## Excitation Types

The `AbstractExcitation` argument controls what the optimization variables are and
how the array factor is assembled.

| Type | Array | Variables | AF |
|---|---|---|---|
| `ComplexWeights` | `ArrayGeometry` | $w_n \in \mathbb{C}$, free | Complex |
| `RealWeights` | `ArrayGeometry` | $w_n \in \mathbb{R}$ | Complex |
| `ConjugateSymmetricWeights` | `SymmetricArray` | $\lceil N/2 \rceil$ complex weights | **Real** |
| `ProgressivePhaseAmplitude` | any | $a_n \in \mathbb{R}$, $\;w_n = a_n e^{j\boldsymbol{\beta}\cdot\mathbf{p}_n}$ | Complex (Real on `SymmetricArray`) |

**`ProgressivePhaseAmplitude`** factors out a fixed phase gradient toward a steering
direction $\boldsymbol{\beta}$, so the optimizer only handles real amplitudes $a_n$.
The steering matrix is evaluated in shifted coordinates $\hat{r} - \boldsymbol{\beta}$,
which makes the constraint at broadside exact and linear.
On a `SymmetricArray`, the AF is again real for the same reason as above.

## Modulus Constraints per Formulation

The fundamental sidelobe constraint is $|AF(\hat{r}_p)| \leq t_p$.

**SOCP — exact:**

```math
|AF(\hat{r}_p)| \leq t_p
```

**LP, QP, MILP — polyhedral (linear constraint):**
The unit disk is inner-approximated by a regular $F$-faced polygon:

```math
\mathrm{Re}\!\left(e^{j\theta_\ell}\, AF(\hat{r}_p)\right) \leq t_p \cos(\pi/F),
\qquad \theta_\ell = \frac{2\pi\ell}{F},\quad \ell = 0,\ldots,F-1
```

This constraint is linear in all three formulations; in QP the quadratic term appears
only in the objective, not in the constraints.
The approximation error is $1 - \cos(\pi/F)$ (~7.6% with the default $F = 8$);
increasing `polygon_faces` reduces it.

**Real AF** (from `ConjugateSymmetricWeights` or `ProgressivePhaseAmplitude` on a
`SymmetricArray`): the imaginary part is zero by construction, so the bound is exact
in all formulations:

```math
-t_p \leq \mathrm{Re}(AF(\hat{r}_p)) \leq t_p
```
