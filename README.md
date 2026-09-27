# sage-kahler

`sage-kahler` is an experimental SageMath package for symbolic computations in
local complex and Kähler geometry.

The initial version provides local holomorphic charts, formal conjugation,
and sparse complex differential forms backed by Sage's symbolic ring:

```python
from sage_kahler.all import *

X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()
z_bar, w_bar = X.conjugate_coordinates()

print(X)
print(bar(z))
print(abs2(z))

partial(z * bar(z))
dbar(z * bar(z))
X.dz(0).wedge(X.dzbar(0))

# ddbar factors its resulting coefficients automatically.
ddbar(abs(z))
```

Before differentiation, `partial`, `dbar`, and `ddbar` rewrite coordinate-dependent
`abs(f)` as `sqrt(f * bar(f))`, including in form coefficients. This makes the
dependence on the formal conjugate coordinates explicit; differentiating Sage's
`abs(z)` directly with respect to the independent symbol `z_bar` would give zero.
For example, `ddbar(abs(z))` equals `dz ∧ dbar(z) / (4*|z|)` for `z != 0`.
These are classical derivatives on the smooth locus; `abs(z)` is not
differentiable at the origin. Use `sqrt(abs2(z))` to write the same input
explicitly in formal coordinates.

Radial expressions are displayed using ordinary norm notation. Simplification
conditions remain available through `.conditions()` but are hidden in output:

```python
a = var("a")
result = ddbar(abs2(z) ** (a / 2))
result
# (1/4*a^2*|z|^(a - 2)) dz ∧ dbar(z)

result.raw()       # the unsimplified formal result
result.conditions()  # ('z != 0',)
```

Formal expressions can also be realized in real coordinates, which is useful
for checking Wirtinger calculations:

```python
X.real_coordinates()
# ((x1, y1), (x2, y2))

X.to_real(abs2(z))
# x1^2 + y1^2
```

The package treats `z`, `w`, `z_bar`, and `w_bar` as independent symbolic
variables for differentiation. Its conventions follow Huybrechts:

```text
d = partial + dbar
partial^2 = dbar^2 = 0
partial*dbar = -dbar*partial
```

Charts on the same coordinate names are interchangeable, so re-running
`X = ComplexChart(2, names=("z", "w"))` in a notebook keeps earlier forms
usable. A chart with the same names in a different order is a different chart,
and the most recently created one is used when a chart is inferred.

Declare real parameters with `var("beta", domain="real")`. Otherwise, `bar`
treats them as complex and produces `conjugate(beta)`.

`form1 == form2` is `True` only when every coefficient difference is proven to
vanish, by exact cancellation of coordinate square roots or by
`simplify_full`. On the chart, `sqrt(F * bar(F))` is taken to be `|F|`, so
for example `|z| |w| = |z w|`. A `False` result means no proof was found; a
difference that is clearly nonzero at a sample point is reported immediately.

## Development

Activate your Sage environment (for example, `conda activate sage` for a Conda
installation), then install the repository into that environment's Python:

```bash
python -m pip install -e .
```

Run the tests with Sage's Python:

```bash
python -m pytest tests
```

The `python` command must refer to the environment containing SageMath. Recent
Sage launchers may not support the older `sage -python` and `sage -pip` flags.

## Trying the package in Jupyter

From the repository directory, with the Sage environment active, run:

```bash
python -m jupyterlab
```

Create a notebook using the **SageMath** kernel. In its first cell, run:

```python
from sage_kahler.all import *

X = ComplexChart(1, names=("z",))
(z,) = X.coordinates()
ddbar(abs2(z))
# dz ∧ dbar(z)
```

Here `ddbar(f)` means `partial(dbar(f))`, with no extra factor of `I`.
Try the following in separate cells:

```python
ddbar(sqrt(abs2(z)))
# (1/4/|z|) dz ∧ dbar(z)
```

```python
alpha = var("alpha")
result = ddbar(abs2(z) ** (alpha / 2))
result
# (1/4*alpha^2*|z|^(alpha - 2)) dz ∧ dbar(z)
```

Use `result.terms()` to inspect the symbolic coefficients, `result.conditions()`
for the simplification conditions, and `result.raw()` for the result before
radial simplification. You can also check the sign convention:

```python
f = abs2(z) ** 2
partial(dbar(f)) == -dbar(partial(f))
# True
```


## LaTeX display in Sage notebooks

Use `%display latex` and leave the form as the last expression in a cell:

```python
%display latex
X = ComplexChart(1, names=("z",))
(z,) = X.coordinates()
ddbar(abs(z))
```

This renders as $\frac{1}{4|z|}\,dz\wedge d\bar{z}$ without a validity
message. `latex(result)` returns the LaTeX source; `print(result)` uses plain
text. Simplification conditions are still accessible with `result.conditions()`.


Forms use conjugate and norm notation throughout their displayed coefficients,
including `partial`, `dbar`, and `.raw()` results. Internal coefficients remain
formal expressions in independent symbols; inspect them with `.terms()`.

For a standalone scalar expression, use the chart's display wrapper:

```python
%display latex
X.display(bar(z) / sqrt(abs2(z)))
# Renders as conjugate(z) / |z| using mathematical notation.
```

`X.display(f)` is a zero-form whose coefficient is the original `f`. Native
Sage scalar expressions retain Sage's own norm formatting unless wrapped this
way. Conjugate coordinates have LaTeX names such as `\bar{z}`, while their
internal names and ordinary Sage text representation remain `z_bar`.


## Metric coefficients and determinant

For a Kähler potential `f`, `ddbar(f)` has coefficient matrix
$g_{i\bar j} = \partial^2 f / (\partial z_i\,\partial\bar z_j)$.
Use Sage's `det` directly:

```python
from sage.all import det

X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()
f = 2*abs(z)**2 + 3*abs(w)**2 + z*bar(w) + w*bar(z)
metric = ddbar(f)
metric.coefficient_matrix()  # [[2, 1], [1, 3]]
det(metric)                  # 5
# Equivalently: metric.det() or metric.determinant()
```

Rows correspond to `dz_i` and columns to `dbar(z_j)`, in chart coordinate
order. This is the complex Hessian, with no extra factor of `I`; the package
does not check whether the potential is real or the matrix is positive definite.
Only `(1,1)` forms (and zero) support this operation.

The determinant is a Sage symbolic expression, and the coefficient matrix is
a Sage symbolic matrix. Both carry chart information and display norms and
conjugates automatically, with no `X.display` call required.
Any domain restrictions on `metric` still apply to its determinant; recorded
restrictions are available through `metric.conditions()`.


Matrix and determinant output uses norm notation automatically:

```python
metric = ddbar(abs(z)**4 + abs(w)**2)
G = metric.coefficient_matrix()  # diagonal entries 4*|z|^2, 1
h = det(metric)                 # 4*|z|^2
h.diff(bar(z))                  # 4*z
G.inverse()
ddbar(log(h))
```

These are native Sage subclasses, so matrix arithmetic, indexing, scalar
arithmetic, substitutions, and differentiation continue to use independent
formal coordinates. Both plain text and `%display latex` use chart notation.
Use `G.raw()` or `h.raw()` to retrieve ordinary Sage objects. Modifying the raw
matrix does not modify `G` or the original form. `display=False` also retrieves
ordinary Sage results; previous explicit `display=True` calls remain supported.

Common result methods retain chart formatting. External Sage functions such
as `log(h)` may return ordinary Sage expressions: the calculation is unchanged,
but their output uses Sage's printer. `X.display(...)` remains available for
those standalone expressions. No global Sage printer is modified.


## Radial potentials

For potentials `F(s)` with `s = sum(abs(z_i)**2)`, `ddbar` recognizes the total
squared radius and differentiates the one-variable profile before substituting
coordinates back. It uses
$g_{i\bar j} = F'(s)\delta_{ij} + F''(s)\bar z_i z_j$, avoiding an expanded
multivariable Hessian followed by costly general simplification.

```python
X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()
s = abs(z)**2 + abs(w)**2
u = sqrt(s**2 + 1) + log(s / (sqrt(s**2 + 1) + 1))
metric = ddbar(u)
det(metric)  # 1
```

This example is defined for `s > 0`. The recorded condition excludes the joint
origin, not the individual coordinate axes. Inputs that cannot be recognized
as a function of the total squared radius use the ordinary Dolbeault path.

Radial coefficients are combined and factored before the total squared radius
is substituted back. Form and matrix displays preserve grouped powers such as
`(|z|^2 + |w|^2)^3`, rather than expanding them into degree-six monomials.


## Hermitian metrics

`HermitianMetric` follows Székelyhidi's *Introduction to Extremal Metrics*
(and Tian). A metric, extended complex-bilinearly, has components
$g_{j\bar k} = g(\partial_j, \partial_{\bar k})$, is Hermitian when
$\overline{g_{j\bar k}} = g_{k\bar j}$, and has form
$\omega = \sqrt{-1}\, g_{j\bar k}\, dz^j \wedge d\bar z^k$. A potential gives
$\omega = \sqrt{-1}\partial\bar\partial\varphi$, so its components are exactly
`ddbar(phi).coefficient_matrix()`:

```python
from sage.all import log
from sage_kahler.all import *

X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()
g = HermitianMetric.from_potential(log(1 + abs2(z) + abs2(w)))  # Fubini-Study

g.matrix()                 # g_{j kbar}: row j, column kbar
g.inverse()                # g^{i lbar}, with g^{i lbar} g_{k lbar} = delta^i_k
g.is_hermitian(), g.is_kahler()          # (True, True)
g.ricci_form() == 3 * g.fundamental_form()  # True: Ric = (n + 1) omega
g.scalar_curvature()       # 6 = n (n + 1)
```

Metrics can also be built from components, `HermitianMetric(X, [[1, 0],
[0, 1 + abs2(z)]])`, or from a form, `HermitianMetric.from_form(omega)`,
which divides the coefficients of `omega` by `I`. Neither checks the
Hermitian condition or positivity; use `is_hermitian()` and
`is_positive_definite_at({z: ..., w: ...})`, a numerical check at a point.

`inverse()` is the transpose of the ordinary matrix inverse of `matrix()`,
since the contraction is over the antiholomorphic index. `is_kahler()` tests
$\partial_i g_{j\bar k} = \partial_j g_{i\bar k}$. Also available are
`trace(alpha)` $= g^{j\bar k}\alpha_{j\bar k}$ for
$\alpha = \sqrt{-1}\alpha_{j\bar k}\,dz^j\wedge d\bar z^k$,
`laplacian(f)` $= g^{k\bar l}\partial_k\partial_{\bar l} f$ (half the
Riemannian Laplacian), `ricci_matrix()`
$R_{i\bar j} = -\partial_i\partial_{\bar j}\log\det g$, `ricci_form()` and
`scalar_curvature()`. For a non-Kähler metric these are the (first)
Chern–Ricci curvature and its trace.

`g(X, Y)` evaluates the metric on complex tangent vectors, and
`g.grad10(f)` $= g^{j\bar k}\partial_{\bar k} f\,\partial_j$ is the
$(1,0)$ gradient; `f` is a holomorphy potential when it is holomorphic.


## Complex tangent vectors

`VectorField` (module `vectors`) represents sections of
$T^{\mathbb C}M = T^{1,0}M \oplus T^{0,1}M$ in the frame
$\partial/\partial z^j, \partial/\partial\bar z^j$, dual to
$dz^j, d\bar z^j$:

```python
X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()
v = z * X.d_dz(0) + bar(w) * X.d_dzbar(1)

v(abs2(z))            # directional derivative: z*bar(z)
bar(v)                # formal conjugate, swapping the types
v.holomorphic_part(), v.vector_type()
v.contract(wedge(X.dz(0), X.dzbar(1)))   # interior product
(z * X.d_dz(0)).bracket(z**2 * X.d_dz(1))  # Lie bracket: 2*z^2 d/dw
(z * X.d_dz(0)).is_holomorphic()           # True
g(X.d_dz(0), X.d_dzbar(1)) == g[0, 1]      # True
```

Real tangent vectors are those with `bar(v) == v`, for example
$\partial/\partial x = \partial/\partial z + \partial/\partial\bar z$.


## Chern connection and curvature

`g.chern_connection()` (module `connection`) returns the Chern connection:
the connection compatible with `g` with
$\nabla_{\bar j}\,\partial_k = 0$. It is the Levi-Civita connection
exactly when `g` is Kähler (Székelyhidi, Remark 1.23).

```python
C = g.chern_connection()
C.christoffel(i, j, k)     # Gamma^i_{jk} = g^{i lbar} d_j g_{k lbar}
C.torsion(i, j, k)         # Gamma^i_{jk} - Gamma^i_{kj}; C.is_torsion_free()
C.covariant_derivative(v, T)   # nabla_v T for vector fields, forms, functions
C.curvature(i, j, k, l)    # R_{i jbar k lbar}
C.curvature_endomorphism(i, j, k, l)   # R_i^j_{k lbar} = -d_lbar Gamma^j_{ki}
C.curvature_operator(X, Y, Z)          # R(X, Y) Z from covariant derivatives
C.first_chern_ricci()      # g^{i jbar} R_{i jbar k lbar} = g.ricci_matrix()
C.second_chern_ricci()     # g^{k lbar} R_{i jbar k lbar}
```

Here $R_{i\bar j k\bar l} = -\partial_k\partial_{\bar l} g_{i\bar j} +
g^{p\bar q}\,\partial_k g_{i\bar q}\,\partial_{\bar l} g_{p\bar j}$ and
$(\nabla_k\nabla_{\bar l} - \nabla_{\bar l}\nabla_k)\partial_i =
R_i{}^j{}_{k\bar l}\,\partial_j$. For Fubini–Study,
$R_{i\bar j k\bar l} = g_{i\bar j}g_{k\bar l} + g_{i\bar l}g_{k\bar j}$.
The two Ricci contractions agree for Kähler metrics and differ in general.
Components are computed lazily and cached; `christoffel_symbols()` and
`curvature_tensor()` return dictionaries of the nonzero entries.
