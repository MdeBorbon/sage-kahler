"""The Chern connection and curvature of a Hermitian metric.

Conventions follow Székelyhidi, §1.4 and §1.6. The Chern connection on
``T^{1,0}`` is the unique connection compatible with the metric such that
``nabla_{jbar} d/dz^k = 0``; its Christoffel symbols are

    nabla_j d/dz^k = Gamma^i_{jk} d/dz^i,   Gamma^i_{jk} = g^{i lbar} d_j g_{k lbar}.

It acts on ``T^{0,1}`` by conjugation, and on forms by duality, so that
``nabla_i dz^k = -Gamma^k_{ij} dz^j``. For a Kähler metric it is the
Levi-Civita connection (Remark 1.23); in general its torsion
``T^i_{jk} = Gamma^i_{jk} - Gamma^i_{kj}`` vanishes exactly when the metric is
Kähler. The curvature is

    (nabla_k nabla_lbar - nabla_lbar nabla_k) d/dz^i = R_i^j_{k lbar} d/dz^j,
    R_i^j_{k lbar} = -d_lbar Gamma^j_{ki},
    R_{i jbar k lbar} = g_{p jbar} R_i^p_{k lbar}
        = -d_k d_lbar g_{i jbar} + g^{p qbar} (d_k g_{i qbar}) (d_lbar g_{p jbar}).

In ``R_{i jbar k lbar}`` the pair ``(i, jbar)`` is the endomorphism part and
``(k, lbar)`` the two-form part.
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

from sage.all import SR, matrix

from .conjugation import bar
from .forms import DifferentialForm, _is_zero, _simplify_scalar
from .results import ChartExpression, ChartMatrix
from .vectors import VectorField


def _sorted_basis(sequence):
    """Return ``(sign, sorted basis)``, or ``(0, ())`` for a repeated index."""
    if len(set(sequence)) != len(sequence):
        return 0, ()
    inversions = sum(
        sequence[a] > sequence[b]
        for a in range(len(sequence))
        for b in range(a + 1, len(sequence))
    )
    return (-1 if inversions % 2 else 1), tuple(sorted(sequence))


class ChernConnection:
    """The Chern connection of a :class:`HermitianMetric`."""

    def __init__(self, metric):
        self._metric = metric
        self._chart = metric.chart
        self._gamma = None
        self._gamma_bar = None
        self._curvature = {}
        self._metric_derivatives = None

    @property
    def metric(self):
        return self._metric

    @property
    def chart(self):
        return self._chart

    def _wrap(self, value):
        return ChartExpression(SR, value, self._chart)

    def _simplify(self, value):
        return _simplify_scalar(SR(value), self._chart)

    def _derivatives(self):
        """Return ``(d_k g, d_kbar g)`` as lists of matrices indexed by k."""
        if self._metric_derivatives is None:
            components = self._metric._components
            self._metric_derivatives = tuple(
                [
                    components.apply_map(lambda entry, v=variable: entry.diff(v))
                    for variable in variables
                ]
                for variables in (
                    self._chart.coordinates(),
                    self._chart.conjugate_coordinates(),
                )
            )
        return self._metric_derivatives

    def _christoffel_table(self):
        if self._gamma is None:
            inverse = self._metric.inverse().raw()
            holomorphic, _ = self._derivatives()
            dimension = self._chart.dimension()
            self._gamma = {
                (i, j, k): self._simplify(sum(
                    inverse[i, l] * holomorphic[j][k, l] for l in range(dimension)
                ))
                for i in range(dimension)
                for j in range(dimension)
                for k in range(dimension)
            }
            self._gamma_bar = {
                key: bar(value) for key, value in self._gamma.items()
            }
        return self._gamma

    def christoffel(self, i, j, k):
        """Return ``Gamma^i_{jk} = g^{i lbar} d_j g_{k lbar}``."""
        return self._wrap(self._christoffel_table()[i, j, k])

    def christoffel_symbols(self):
        """Return the nonzero symbols as ``{(i, j, k): Gamma^i_{jk}}``."""
        return {
            key: self._wrap(value)
            for key, value in self._christoffel_table().items()
            if not _is_zero(value, self._chart)
        }

    def torsion(self, i, j, k):
        """Return ``T^i_{jk} = Gamma^i_{jk} - Gamma^i_{kj}``."""
        table = self._christoffel_table()
        return self._wrap(self._simplify(table[i, j, k] - table[i, k, j]))

    def is_torsion_free(self):
        """Return True when all torsion components are proven to vanish.

        For a Hermitian metric this holds exactly when the metric is Kähler.
        """
        table = self._christoffel_table()
        dimension = self._chart.dimension()
        return all(
            _is_zero(table[i, j, k] - table[i, k, j], self._chart)
            for i in range(dimension)
            for j in range(dimension)
            for k in range(j + 1, dimension)
        )

    def covariant_derivative(self, direction, value):
        """Return ``nabla_X value`` for a vector field, form, or function.

        ``direction`` is the vector field ``X``. Functions are differentiated
        as ``X(f)``.
        """
        if not isinstance(direction, VectorField):
            raise TypeError("the direction must be a vector field")
        if direction.chart != self._chart:
            raise ValueError("the direction belongs to a different chart")
        if isinstance(value, VectorField):
            return self._derivative_of_vector(direction, value)
        if isinstance(value, DifferentialForm):
            return self._derivative_of_form(direction, value)
        return self._wrap(self._simplify(direction._apply(value)))

    def _derivative_of_vector(self, direction, vector):
        gamma = self._christoffel_table()
        dimension = self._chart.dimension()
        components = {}
        for offset, table in ((0, gamma), (dimension, self._gamma_bar)):
            for i in range(dimension):
                value = direction._apply(vector[offset + i]) + sum(
                    direction[offset + j] * table[i, j, k] * vector[offset + k]
                    for j in range(dimension)
                    for k in range(dimension)
                )
                components[offset + i] = self._simplify(value)
        return VectorField(
            self._chart,
            components,
            conditions=direction.conditions() + vector.conditions(),
        )

    def _derivative_of_basis_form(self, direction, index):
        """Return ``nabla_X dx^index`` as ``{index: coefficient}``."""
        self._christoffel_table()
        dimension = self._chart.dimension()
        offset, table = (
            (0, self._gamma) if index < dimension else (dimension, self._gamma_bar)
        )
        k = index - offset
        return {
            offset + p: -sum(
                direction[offset + j] * table[k, j, p] for j in range(dimension)
            )
            for p in range(dimension)
        }

    def _derivative_of_form(self, direction, form):
        if form.chart != self._chart:
            raise ValueError("the form belongs to a different chart")
        terms = {}
        for basis, coefficient in form.terms().items():
            terms[basis] = terms.get(basis, SR.zero()) + direction._apply(coefficient)
            for position, index in enumerate(basis):
                for new_index, factor in self._derivative_of_basis_form(
                    direction, index
                ).items():
                    sign, new_basis = _sorted_basis(
                        basis[:position] + (new_index,) + basis[position + 1:]
                    )
                    if sign:
                        terms[new_basis] = terms.get(new_basis, SR.zero()) + (
                            sign * coefficient * factor
                        )
        return DifferentialForm(
            self._chart,
            {basis: self._simplify(value) for basis, value in terms.items()},
            conditions=direction.conditions() + form.conditions(),
        )

    def curvature(self, i, j, k, l):
        """Return ``R_{i jbar k lbar}``."""
        key = (i, j, k, l)
        if key not in self._curvature:
            components = self._metric._components
            inverse = self._metric.inverse().raw()
            holomorphic, antiholomorphic = self._derivatives()
            dimension = self._chart.dimension()
            second = components[i, j].diff(self._chart.coordinates()[k]).diff(
                self._chart.conjugate_coordinates()[l]
            )
            self._curvature[key] = self._simplify(-second + sum(
                inverse[p, q] * holomorphic[k][i, q] * antiholomorphic[l][p, j]
                for p in range(dimension)
                for q in range(dimension)
            ))
        return self._wrap(self._curvature[key])

    def curvature_endomorphism(self, i, j, k, l):
        """Return ``R_i^j_{k lbar} = -d_lbar Gamma^j_{ki}``."""
        gamma = self._christoffel_table()[j, k, i]
        conjugate = self._chart.conjugate_coordinates()[l]
        return self._wrap(self._simplify(-gamma.diff(conjugate)))

    def curvature_tensor(self):
        """Return the nonzero ``{(i, j, k, l): R_{i jbar k lbar}}``."""
        dimension = self._chart.dimension()
        result = {}
        for i in range(dimension):
            for j in range(dimension):
                for k in range(dimension):
                    for l in range(dimension):
                        value = self.curvature(i, j, k, l)
                        if not _is_zero(value.raw(), self._chart):
                            result[i, j, k, l] = value
        return result

    def curvature_operator(self, first, second, vector):
        """Return ``R(X, Y) Z = nabla_X nabla_Y Z - nabla_Y nabla_X Z - nabla_[X,Y] Z``."""
        nabla = self.covariant_derivative
        return (
            nabla(first, nabla(second, vector))
            - nabla(second, nabla(first, vector))
            - nabla(first.bracket(second), vector)
        )

    def first_chern_ricci(self):
        """Return ``g^{i jbar} R_{i jbar k lbar}``, indexed by ``(k, lbar)``.

        This is the trace of the curvature endomorphism, and equals
        ``-d_k d_lbar log det g`` (Lemma 1.18), i.e. ``metric.ricci_matrix()``.
        """
        return self._contracted(endomorphism=True)

    def second_chern_ricci(self):
        """Return ``g^{k lbar} R_{i jbar k lbar}``, indexed by ``(i, jbar)``.

        Székelyhidi's definition of Ricci curvature. For a Kähler metric the
        two Chern–Ricci curvatures agree.
        """
        return self._contracted(endomorphism=False)

    def _contracted(self, endomorphism):
        inverse = self._metric.inverse().raw()
        dimension = self._chart.dimension()

        def entry(a, b):
            if endomorphism:
                terms = (
                    inverse[i, j] * self.curvature(i, j, a, b).raw()
                    for i in range(dimension)
                    for j in range(dimension)
                )
            else:
                terms = (
                    inverse[k, l] * self.curvature(a, b, k, l).raw()
                    for k in range(dimension)
                    for l in range(dimension)
                )
            return self._simplify(sum(terms))

        return ChartMatrix.from_matrix(
            matrix(SR, dimension, dimension, entry), self._chart
        )

    def __repr__(self):
        return f"Chern connection of the {self._metric!r}"
