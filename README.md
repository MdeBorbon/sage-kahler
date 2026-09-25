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
ddbar(sqrt(abs2(z)))
```

Radial expressions are displayed using ordinary norm notation, together with
the domain condition used during simplification:

```python
a = var("a")
result = ddbar(abs2(z) ** (a / 2))
result
# (1/4*a^2*|z|^(a - 2)) dz ∧ dbar(z)  [valid where z != 0]

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
# (1/4/|z|) dz ∧ dbar(z)  [valid where z != 0]
```

```python
alpha = var("alpha")
result = ddbar(abs2(z) ** (alpha / 2))
result
# (1/4*alpha^2*|z|^(alpha - 2)) dz ∧ dbar(z)  [valid where z != 0]
```

Use `result.terms()` to inspect the symbolic coefficients, `result.conditions()`
for the simplification conditions, and `result.raw()` for the result before
radial simplification. You can also check the sign convention:

```python
f = abs2(z) ** 2
partial(dbar(f)) == -dbar(partial(f))
# True
```
