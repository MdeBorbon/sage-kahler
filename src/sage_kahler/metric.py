"""Hermitian metrics on a complex chart, following Székelyhidi's conventions.

A Hermitian metric ``g``, extended complex-bilinearly to complex tangent
vectors, is determined by its components

    g_{j kbar} = g(d/dz^j, d/dzbar^k),

with ``g(d/dz^j, d/dz^k) = g(d/dzbar^j, d/dzbar^k) = 0``. The Hermitian
condition is ``bar(g_{j kbar}) = g_{k jbar}`` and the associated form is

    omega = I * sum g_{j kbar} dz^j ∧ dzbar^k.

For a potential ``phi``, ``omega = I * partial(dbar(phi))``, so
``g_{j kbar} = d^2 phi / (dz^j dzbar^k)``: the coefficient matrix of
``ddbar(phi)``. The inverse ``g^{j kbar}`` satisfies
``g^{i lbar} g_{k lbar} = delta^i_k`` and the Laplacian is
``g^{k lbar} d_k d_lbar``, one half of the Riemannian Laplacian.
"""

# ****************************************************************************
#       Copyright (C) 2026 Martin de Borbon <martdeborbon@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import CDF, I, SR, log, matrix

from .conjugation import bar
from .forms import (
    DifferentialForm,
    _as_form,
    _formalize_absolute_values,
    _is_zero,
    _simplify_scalar,
    ddbar,
)
from .results import ChartExpression, ChartMatrix


class HermitianMetric:
    """A Hermitian metric on a :class:`ComplexChart`, given by ``g_{j kbar}``.

    Rows of the component matrix are holomorphic indices ``j`` and columns
    are antiholomorphic indices ``kbar``, in chart coordinate order. The
    constructor does not check the Hermitian condition or positivity; use
    ``is_hermitian()`` and ``is_positive_definite_at()``.
    """

    def __init__(self, chart, components, *, conditions=()):
        dimension = chart.dimension()
        components = matrix(SR, components)
        if components.dimensions() != (dimension, dimension):
            raise ValueError(
                f"expected a {dimension} x {dimension} matrix of components"
            )
        self._chart = chart
        # Expose conjugate dependence hidden in abs(), as the Dolbeault
        # operators do, so that derivatives of the components are correct.
        self._components = components.apply_map(
            lambda entry: _formalize_absolute_values(entry, chart)
        )
        self._conditions = tuple(dict.fromkeys(conditions))
        self._inverse = None
        self._determinant = None
        self._ricci = None

    @classmethod
    def from_form(cls, omega):
        """Return the metric whose form is ``omega = I g_{j kbar} dz^j ∧ dzbar^k``."""
        if not isinstance(omega, DifferentialForm):
            raise TypeError("from_form requires a differential form")
        components = omega.coefficient_matrix(display=False) / I
        return cls(omega.chart, components, conditions=omega.conditions())

    @classmethod
    def from_potential(cls, potential, chart=None):
        """Return the Kähler metric ``omega = I * partial(dbar(potential))``."""
        form = ddbar(_as_form(potential, chart))
        return cls(
            form.chart,
            form.coefficient_matrix(display=False),
            conditions=form.conditions(),
        )

    @property
    def chart(self):
        return self._chart

    def conditions(self):
        """Return domain restrictions inherited from the defining form."""
        return self._conditions

    def matrix(self):
        """Return the components ``g_{j kbar}``, row ``j`` and column ``kbar``."""
        return ChartMatrix.from_matrix(self._components, self._chart)

    def __getitem__(self, key):
        """Return ``g_{j kbar}`` for ``metric[j, k]``."""
        j, k = key
        return ChartExpression(SR, self._components[j, k], self._chart)

    def __call__(self, first, second):
        """Evaluate ``g(X, Y)`` on complex tangent vectors, bilinearly.

        With ``g = g_{j kbar} (dz^j ⊗ dzbar^k + dzbar^k ⊗ dz^j)`` this is
        ``g_{j kbar} (X^j Y^kbar + X^kbar Y^j)``; in particular
        ``g(d/dz^j, d/dzbar^k) = g_{j kbar}``.
        """
        from .vectors import VectorField

        for vector in (first, second):
            if not isinstance(vector, VectorField):
                raise TypeError("the metric is evaluated on vector fields")
            if vector.chart != self._chart:
                raise ValueError("the vector field belongs to a different chart")
        dimension = self._chart.dimension()
        total = sum(
            (
                self._components[j, k] * (
                    first[j] * second[dimension + k]
                    + first[dimension + k] * second[j]
                )
                for j in range(dimension)
                for k in range(dimension)
            ),
            SR.zero(),
        )
        return ChartExpression(SR, _simplify_scalar(total, self._chart), self._chart)

    def grad10(self, function):
        """Return ``grad^{1,0} f = g^{j kbar} (d f / dzbar^k) d/dz^j``.

        This is Székelyhidi's ``(1, 0)`` gradient; ``f`` is a holomorphy
        potential when the result is a holomorphic vector field.
        """
        from .vectors import VectorField, _raw

        self.inverse()
        function = _formalize_absolute_values(_raw(function), self._chart)
        dimension = self._chart.dimension()
        derivatives = [function.diff(zbar) for zbar in self._chart.conjugate_coordinates()]
        return VectorField(
            self._chart,
            {
                j: _simplify_scalar(
                    sum(self._inverse[j, k] * derivatives[k] for k in range(dimension)),
                    self._chart,
                )
                for j in range(dimension)
            },
            conditions=self._conditions,
        )

    def chern_connection(self):
        """Return the Chern connection of this metric (see ``connection``)."""
        from .connection import ChernConnection

        return ChernConnection(self)

    def inverse(self):
        """Return the inverse ``g^{i lbar}``, with ``g^{i lbar} g_{k lbar} = delta^i_k``.

        Row ``i`` and column ``lbar``. As a matrix this is the transpose of
        the ordinary inverse of ``matrix()``.
        """
        if self._inverse is None:
            inverse = self._components.inverse().transpose()
            self._inverse = inverse.apply_map(
                lambda entry: _simplify_scalar(entry, self._chart)
            )
        return ChartMatrix.from_matrix(self._inverse, self._chart)

    def determinant(self):
        """Return ``det(g_{j kbar})``."""
        if self._determinant is None:
            self._determinant = _simplify_scalar(
                self._components.determinant(), self._chart
            )
        return ChartExpression(SR, self._determinant, self._chart)

    det = determinant

    def fundamental_form(self):
        """Return ``omega = I * sum g_{j kbar} dz^j ∧ dzbar^k``."""
        dimension = self._chart.dimension()
        return DifferentialForm(
            self._chart,
            {
                (j, dimension + k): I * self._components[j, k]
                for j in range(dimension)
                for k in range(dimension)
            },
            conditions=self._conditions,
        )

    kahler_form = fundamental_form

    def is_hermitian(self):
        """Return True when ``bar(g_{j kbar}) = g_{k jbar}`` is proven.

        Like form equality, this treats undeclared symbolic parameters as
        real: Sage's simplifier reduces ``conjugate(a) - a`` to zero. Only the
        dependence on the chart coordinates is genuinely checked.
        """
        dimension = self._chart.dimension()
        return all(
            _is_zero(
                bar(self._components[j, k]) - self._components[k, j],
                self._chart,
            )
            for j in range(dimension)
            for k in range(j, dimension)
        )

    def is_kahler(self):
        """Return True when ``d_i g_{j kbar} = d_j g_{i kbar}`` is proven.

        For a Hermitian metric this is equivalent to ``d(omega) = 0``.
        """
        coordinates = self._chart.coordinates()
        dimension = self._chart.dimension()
        return all(
            _is_zero(
                self._components[j, k].diff(coordinates[i])
                - self._components[i, k].diff(coordinates[j]),
                self._chart,
            )
            for i in range(dimension)
            for j in range(i + 1, dimension)
            for k in range(dimension)
        )

    def is_positive_definite_at(self, point):
        """Numerically test positivity at ``point``, a dict ``{z_j: value}``.

        Conjugate coordinates are set to the complex conjugates of the given
        values; parameters must also be given numerical values.
        """
        substitutions = {}
        for coordinate, conjugate in zip(
            self._chart.coordinates(), self._chart.conjugate_coordinates()
        ):
            if coordinate not in point:
                raise ValueError(f"no value given for {coordinate}")
            value = CDF(point[coordinate])
            substitutions[coordinate] = value
            substitutions[conjugate] = value.conjugate()
        for variable, value in point.items():
            substitutions.setdefault(variable, value)
        numeric = matrix(
            CDF, self._components.subs(substitutions).apply_map(lambda e: CDF(e))
        )
        tolerance = 1e-10 * max(1, max(abs(entry) for entry in numeric.list()))
        if (numeric - numeric.conjugate_transpose()).norm() > tolerance:
            return False
        return all(
            eigenvalue.real() > tolerance
            for eigenvalue in numeric.eigenvalues()
        )

    def trace(self, alpha):
        """Return ``tr_omega alpha = g^{j kbar} alpha_{j kbar}``.

        Here ``alpha = I * alpha_{j kbar} dz^j ∧ dzbar^k`` is a (1, 1) form, as
        in Székelyhidi's Lemma 4.6, so ``trace(omega) = n``.
        """
        components = alpha.coefficient_matrix(display=False) / I
        return ChartExpression(SR, self._trace_matrix(components), self._chart)

    def laplacian(self, function):
        """Return ``g^{k lbar} d_k d_lbar function``."""
        form = _as_form(function, self._chart)
        if set(form.terms()) - {()}:
            raise ValueError("laplacian requires a scalar function")
        hessian = ddbar(form).coefficient_matrix(display=False)
        return ChartExpression(SR, self._trace_matrix(hessian), self._chart)

    def ricci_matrix(self):
        """Return ``R_{i jbar} = -d_i d_jbar log det(g)``.

        For a Kähler metric this is the Ricci curvature; for a general
        Hermitian metric it is the (first) Chern–Ricci curvature.
        """
        return ChartMatrix.from_matrix(self._ricci_components(), self._chart)

    def ricci_form(self):
        """Return ``Ric(omega) = I R_{i jbar} dz^i ∧ dzbar^j = -I ddbar(log det g)``."""
        form = -I * ddbar(log(self.determinant().raw()), chart=self._chart)
        return DifferentialForm(
            self._chart,
            form.terms(),
            conditions=self._conditions + form.conditions(),
        )

    def scalar_curvature(self):
        """Return ``S = g^{i jbar} R_{i jbar}``.

        With these conventions, the Fubini–Study metric on CP^n has
        ``S = n (n + 1)``; this is half of the Riemannian scalar curvature.
        """
        return ChartExpression(
            SR, self._trace_matrix(self._ricci_components()), self._chart
        )

    def _ricci_components(self):
        if self._ricci is None:
            form = ddbar(log(self.determinant().raw()), chart=self._chart)
            self._ricci = -form.coefficient_matrix(display=False)
        return self._ricci

    def _trace_matrix(self, components):
        self.inverse()
        dimension = self._chart.dimension()
        total = sum(
            self._inverse[j, k] * components[j, k]
            for j in range(dimension)
            for k in range(dimension)
        )
        return _simplify_scalar(SR(total), self._chart)

    def __repr__(self):
        return f"Hermitian metric with form {self.fundamental_form()!r}"

    def _latex_(self):
        from sage.all import latex

        return latex(self.fundamental_form())
