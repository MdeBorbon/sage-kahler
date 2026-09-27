# sage-kahler

An experimental SageMath package for symbolic computations in local complex
and Kähler geometry: holomorphic charts, differential forms and Dolbeault
operators, Hermitian metrics, complex vector fields, and the Chern connection
with its curvature. Conventions follow Székelyhidi, *An Introduction to
Extremal Kähler Metrics*, and Tian.

## Installation

With the Sage environment active (for example `conda activate sage`), install
the package in editable mode:

```bash
python -m pip install -e .
```

Run the tests:

```bash
python -m pytest tests
```

Example notebooks are in `notebooks/`. Open them with the SageMath kernel:

```bash
python -m jupyterlab
```

The examples below assume a Sage session, such as a notebook.

## Charts and forms

```python
from sage_kahler.all import *

X = ComplexChart(2, names=("z", "w"))
z, w = X.coordinates()

ddbar(abs(z)**2)           # dz ∧ dbar(z)
ddbar(abs(z))              # (1/4/|z|) dz ∧ dbar(z)
X.dz(0).wedge(X.dzbar(1))  # dz ∧ dbar(w)
X.to_real(abs(z)**2)       # x1^2 + y1^2
```

- `z` and `bar(z)` are independent symbols. Before differentiating, `abs(f)`
  is rewritten as `sqrt(f*bar(f))`.
- `ddbar(f)` is `partial(dbar(f))`, with no factor of `I`. Signs follow
  `d = partial + dbar` and `partial dbar = -dbar partial`.
- Results are simplified away from singularities. `.conditions()` lists the
  assumptions, such as `z != 0`, and `.raw()` gives the unsimplified form.
- `form1 == form2` is `True` only when the difference is proven to vanish.
- `%display latex` shows results in mathematical notation.
- Declare real parameters with `var("beta", domain="real")`.

## Kähler potentials

`ddbar(f).coefficient_matrix()` is the complex Hessian
$\partial^2 f/\partial z_i\partial\bar z_j$:

```python
f = 2*abs(z)**2 + 3*abs(w)**2 + z*bar(w) + w*bar(z)
ddbar(f).coefficient_matrix()   # [[2, 1], [1, 3]]
det(ddbar(f))                   # 5
```

Potentials of $s = |z|^2 + |w|^2$ are differentiated as one-variable
functions, which keeps results compact. For the Eguchi–Hanson potential:

```python
s = abs(z)**2 + abs(w)**2
det(ddbar(sqrt(s**2 + 1) + log(s / (sqrt(s**2 + 1) + 1))))   # 1
```

## Hermitian metrics

`HermitianMetric` has components $g_{j\bar k} = g(\partial_j, \partial_{\bar k})$
and form $\omega = \sqrt{-1}\,g_{j\bar k}\,dz^j\wedge d\bar z^k$. A potential
gives $\omega = \sqrt{-1}\partial\bar\partial\varphi$.

```python
g = HermitianMetric.from_potential(log(1 + s))   # Fubini–Study
g.is_hermitian(), g.is_kahler()                  # (True, True)
g.ricci_form() == 3 * g.fundamental_form()       # True
g.scalar_curvature()                             # 6
```

Metrics can also be given by components, `HermitianMetric(X, [[1, 0], [0, 1 + abs(z)**2]])`,
or by a form, `HermitianMetric.from_form(omega)`. Other methods:

- `matrix()`, `inverse()` (with $g^{i\bar l}g_{k\bar l} = \delta^i_k$), `det()`
- `is_positive_definite_at({z: ..., w: ...})`, a numerical check at a point
- `trace(alpha)`, `laplacian(f)` $= g^{k\bar l}\partial_k\partial_{\bar l}f$
- `ricci_matrix()` $= -\partial_i\partial_{\bar j}\log\det g$, `ricci_form()`, `scalar_curvature()`
- `g(X, Y)` on vector fields, and `grad10(f)` $= g^{j\bar k}\partial_{\bar k}f\,\partial_j$

## Vector fields

`VectorField` represents complex vector fields in the frame
$\partial/\partial z^j, \partial/\partial\bar z^j$:

```python
v = z * X.d_dz(0) + bar(w) * X.d_dzbar(1)
v(abs(z)**2)                                # |z|^2
bar(v)                                      # (w) ∂/∂w + (bar(z)) ∂/∂bar(z)
(z * X.d_dz(0)).bracket(z**2 * X.d_dz(1))   # (2*z^2) ∂/∂w
(z * X.d_dz(0)).is_holomorphic()            # True
g(X.d_dz(0), X.d_dzbar(1))                  # equals g[0, 1]
g.grad10(abs(z)**2 / (1 + s))               # (z) ∂/∂z
```

Also available: `holomorphic_part()`, `vector_type()`, `is_real()`, and
`contract(form)` for the interior product.

## Chern connection and curvature

`g.chern_connection()` is the connection compatible with `g` with
$\nabla_{\bar j}\partial_k = 0$. It is the Levi-Civita connection exactly when
`g` is Kähler.

```python
C = g.chern_connection()
C.christoffel(0, 0, 0)     # -2*bar(z)/(|z|^2 + |w|^2 + 1)
C.curvature(0, 0, 0, 0)    # 2*(|w|^2 + 1)^2/(|z|^2 + |w|^2 + 1)^4
C.is_torsion_free()        # True
```

- `christoffel(i, j, k)` $= \Gamma^i_{jk} = g^{i\bar l}\partial_j g_{k\bar l}$, and `torsion(i, j, k)`
- `covariant_derivative(X, T)` for vector fields, forms and functions
- `curvature(i, j, k, l)` $= R_{i\bar j k\bar l}$, `curvature_endomorphism(i, j, k, l)`, `curvature_operator(X, Y, Z)`
- `first_chern_ricci()` and `second_chern_ricci()`, which agree for Kähler metrics

## License

GPL-2.0-or-later, the same license as SageMath.
