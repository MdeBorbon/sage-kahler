"""Complex tangent vector fields on a :class:`ComplexChart`.

Following Székelyhidi, vector fields are sections of the complexified tangent
bundle ``T^C M = T^{1,0} M ⊕ T^{0,1} M``, written in the coordinate frame

    d/dz^1, ..., d/dz^n, d/dzbar^1, ..., d/dzbar^n,

which is dual to ``dz^1, ..., dz^n, dzbar^1, ..., dzbar^n``. Combined indices
match those of differential forms: ``0 .. n-1`` are holomorphic and
``n .. 2n-1`` antiholomorphic. As with forms, coordinates and their formal
conjugates are independent symbols, and ``bar`` conjugates a vector field.
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

from sage.all import SR
from sage.symbolic.operators import add_vararg

from .forms import (
    DifferentialForm,
    _format_coefficient,
    _formalize_absolute_values,
    _is_zero,
)
from .conjugation import bar


def _raw(value):
    raw = getattr(value, "raw", None)
    return SR(raw() if raw is not None else value)


class VectorField:
    """A sparse complex vector field ``X = X^a d/dx^a`` on a chart."""

    def __init__(self, chart, components=None, *, conditions=()):
        self._chart = chart
        dimension = 2 * chart.dimension()
        if components is None:
            components = {}
        elif not isinstance(components, dict):
            components = list(components)
            if len(components) != dimension:
                raise ValueError(f"expected {dimension} components")
            components = dict(enumerate(components))
        combined = {}
        for index, coefficient in components.items():
            if not 0 <= index < dimension:
                raise IndexError("vector index is outside the chart")
            combined[index] = combined.get(index, SR.zero()) + _raw(coefficient)
        self._components = {
            index: coefficient
            for index, coefficient in combined.items()
            if not _is_zero(coefficient, chart)
        }
        self._conditions = tuple(dict.fromkeys(conditions))

    @property
    def chart(self):
        return self._chart

    def components(self):
        """Return a copy of the sparse index-to-coefficient dictionary."""
        return dict(self._components)

    def __getitem__(self, index):
        """Return the component ``X^a`` for a combined index ``a``."""
        if not 0 <= index < 2 * self._chart.dimension():
            raise IndexError("vector index is outside the chart")
        return self._components.get(index, SR.zero())

    def holomorphic_components(self):
        """Return ``(X^1, ..., X^n)``, the components along ``d/dz^j``."""
        return tuple(self[j] for j in range(self._chart.dimension()))

    def antiholomorphic_components(self):
        """Return ``(X^1bar, ..., X^nbar)``, the components along ``d/dzbar^j``."""
        dimension = self._chart.dimension()
        return tuple(self[dimension + j] for j in range(dimension))

    def conditions(self):
        return self._conditions

    def _restricted(self, holomorphic):
        dimension = self._chart.dimension()
        return VectorField(
            self._chart,
            {
                index: coefficient
                for index, coefficient in self._components.items()
                if (index < dimension) == holomorphic
            },
            conditions=self._conditions,
        )

    def holomorphic_part(self):
        """Return the ``(1, 0)`` part of the vector field."""
        return self._restricted(True)

    def antiholomorphic_part(self):
        """Return the ``(0, 1)`` part of the vector field."""
        return self._restricted(False)

    def vector_type(self):
        """Return ``(1, 0)``, ``(0, 1)``, or ``None`` for zero or mixed fields."""
        dimension = self._chart.dimension()
        kinds = {index < dimension for index in self._components}
        if kinds == {True}:
            return (1, 0)
        if kinds == {False}:
            return (0, 1)
        return None

    def __call__(self, function):
        """Return the derivative ``X(f) = X^a d f / dx^a`` of a function."""
        from .results import ChartExpression

        return ChartExpression(SR, self._apply(function), self._chart)

    def _apply(self, function):
        coordinates = self._chart.coordinates() + self._chart.conjugate_coordinates()
        function = _formalize_absolute_values(_raw(function), self._chart)
        return sum(
            (
                coefficient * function.diff(coordinates[index])
                for index, coefficient in self._components.items()
            ),
            SR.zero(),
        )

    def bracket(self, other):
        """Return the Lie bracket ``[X, Y]^a = X(Y^a) - Y(X^a)``."""
        other = self._coerce(other)
        indices = set(self._components).union(other._components)
        return VectorField(
            self._chart,
            {index: self._apply(other[index]) - other._apply(self[index])
             for index in indices},
            conditions=self._conditions + other._conditions,
        )

    def is_holomorphic(self):
        """Return True for a holomorphic vector field.

        That is a ``(1, 0)`` field whose components satisfy
        ``d v^i / dzbar^k = 0``; the zero field counts as holomorphic.
        """
        if self._components and self.vector_type() != (1, 0):
            return False
        return all(
            _is_zero(
                _formalize_absolute_values(coefficient, self._chart).diff(conjugate),
                self._chart,
            )
            for coefficient in self._components.values()
            for conjugate in self._chart.conjugate_coordinates()
        )

    def is_real(self):
        """Return True when ``bar(X) == X``, i.e. X is a real tangent vector."""
        return self == bar(self)

    def contract(self, form):
        """Return the interior product ``i_X form``, with ``dx^a(d/dx^b) = delta``.

        For a one-form ``alpha`` this is the zero-form ``alpha(X)``.
        """
        if not isinstance(form, DifferentialForm):
            raise TypeError("contract requires a differential form")
        if form.chart != self._chart:
            raise ValueError("the vector field and form belong to different charts")
        terms = {}
        for basis, coefficient in form.terms().items():
            for position, index in enumerate(basis):
                component = self._components.get(index)
                if component is None:
                    continue
                remaining = basis[:position] + basis[position + 1:]
                sign = -1 if position % 2 else 1
                terms[remaining] = terms.get(remaining, SR.zero()) + (
                    sign * component * coefficient
                )
        return DifferentialForm(
            self._chart,
            terms,
            conditions=self._conditions + form.conditions(),
        )

    interior_product = contract

    def _coerce(self, other):
        if not isinstance(other, VectorField):
            raise TypeError("expected a vector field")
        if other.chart != self._chart:
            raise ValueError("vector fields must belong to the same chart")
        return other

    def _formal_conjugate_(self):
        dimension = self._chart.dimension()
        return VectorField(
            self._chart,
            {
                (index + dimension if index < dimension else index - dimension):
                    bar(coefficient)
                for index, coefficient in self._components.items()
            },
            conditions=self._conditions,
        )

    def __add__(self, other):
        if isinstance(other, int) and other == 0:
            return self
        other = self._coerce(other)
        components = self.components()
        for index, coefficient in other._components.items():
            components[index] = components.get(index, SR.zero()) + coefficient
        return VectorField(
            self._chart, components, conditions=self._conditions + other._conditions
        )

    def __radd__(self, other):
        # Support sum(...) of vector fields.
        return self + other

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __mul__(self, scalar):
        if isinstance(scalar, (VectorField, DifferentialForm)):
            return NotImplemented
        scalar = _raw(scalar)
        return VectorField(
            self._chart,
            {index: coefficient * scalar for index, coefficient in self._components.items()},
            conditions=self._conditions,
        )

    def __rmul__(self, scalar):
        return self * scalar

    def __eq__(self, other):
        if isinstance(other, int) and other == 0:
            other = VectorField(self._chart)
        if not isinstance(other, VectorField) or other.chart != self._chart:
            return False
        indices = set(self._components).union(other._components)
        return all(
            _is_zero(self[index] - other[index], self._chart) for index in indices
        )

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return None

    def __repr__(self):
        if not self._components:
            return "0"
        pieces = []
        for index in sorted(self._components):
            coefficient = self._components[index]
            label = self._chart._vector_label(index)
            text = _format_coefficient(coefficient, self._chart)
            if (coefficient - 1).is_trivial_zero():
                pieces.append(label)
            elif (coefficient + 1).is_trivial_zero():
                pieces.append(f"-{label}")
            else:
                pieces.append(f"({text}) {label}")
        return " + ".join(pieces).replace("+ -", "- ")

    def _latex_(self):
        if not self._components:
            return "0"
        pieces = []
        for index in sorted(self._components):
            coefficient = self._components[index]
            label = self._chart._vector_latex(index)
            text = _format_coefficient(coefficient, self._chart, latex_mode=True)
            if (coefficient - 1).is_trivial_zero():
                pieces.append(label)
            elif (coefficient + 1).is_trivial_zero():
                pieces.append(f"-{label}")
            else:
                if coefficient.operator() == add_vararg:
                    text = rf"\left({text}\right)"
                pieces.append(rf"{text}\,{label}")
        return " + ".join(pieces).replace("+ -", "- ")


def basis_vector(chart, index):
    """Return the coordinate vector field with combined index ``index``."""
    return VectorField(chart, {index: SR.one()})
