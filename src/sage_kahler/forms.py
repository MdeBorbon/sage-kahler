"""Sparse complex differential forms and Dolbeault operators."""

from sage.all import SR

from .conjugation import bar, chart_for_expression


def _is_zero(expression):
    return bool(SR(expression).simplify_full() == 0)


def _wedge_basis(left, right):
    if set(left).intersection(right):
        return 0, ()
    inversions = sum(
        left_index > right_index
        for left_index in left
        for right_index in right
    )
    sign = -1 if inversions % 2 else 1
    return sign, tuple(sorted(left + right))


class DifferentialForm:
    """A sparse complex differential form on a :class:`ComplexChart`."""

    def __init__(
        self,
        chart,
        terms=None,
        *,
        conditions=(),
        radial_indices=(),
        raw_terms=None,
    ):
        self._chart = chart
        combined = {}
        for basis, coefficient in (terms or {}).items():
            basis = tuple(basis)
            if tuple(sorted(basis)) != basis or len(set(basis)) != len(basis):
                raise ValueError("basis indices must be distinct and sorted")
            if any(index < 0 or index >= 2 * chart.dimension() for index in basis):
                raise IndexError("basis index is outside the chart")
            combined[basis] = combined.get(basis, SR.zero()) + SR(coefficient)
        self._terms = {
            basis: coefficient
            for basis, coefficient in combined.items()
            if not _is_zero(coefficient)
        }
        self._conditions = tuple(dict.fromkeys(conditions))
        self._radial_indices = tuple(sorted(set(radial_indices)))
        self._raw_terms = None if raw_terms is None else dict(raw_terms)

    @property
    def chart(self):
        return self._chart

    def terms(self):
        """Return a copy of the sparse basis-to-coefficient dictionary."""
        return dict(self._terms)

    def conditions(self):
        """Return conditions used by assumption-aware simplification."""
        return self._conditions

    def raw(self):
        """Return the form before assumption-aware radial simplification."""
        return DifferentialForm(
            self._chart,
            self._terms if self._raw_terms is None else self._raw_terms,
        )

    def map_coefficients(self, function):
        """Return a form obtained by applying ``function`` to each coefficient."""
        return DifferentialForm(
            self._chart,
            {
                basis: function(coefficient)
                for basis, coefficient in self._terms.items()
            },
            conditions=self._conditions,
            radial_indices=self._radial_indices,
        )

    def factor(self):
        """Return the form with every symbolic coefficient factored by Sage."""
        return self.map_coefficients(lambda coefficient: coefficient.factor())

    def simplify_radial(self):
        """Simplify coefficients using ``z*z_bar = |z|^2`` away from zeros."""
        terms = {}
        radial_indices = set(self._radial_indices)
        nonzero_indices = set()
        for basis, coefficient in self._terms.items():
            simplified, used_indices, coefficient_nonzero_indices = (
                _simplify_radial_coefficient(coefficient, self._chart)
            )
            terms[basis] = simplified
            radial_indices.update(used_indices)
            nonzero_indices.update(coefficient_nonzero_indices)

        conditions = list(self._conditions)
        for index in sorted(nonzero_indices):
            condition = f"{self._chart.coordinates()[index]} != 0"
            if condition not in conditions:
                conditions.append(condition)

        return DifferentialForm(
            self._chart,
            terms,
            conditions=conditions,
            radial_indices=radial_indices,
            raw_terms=self.raw().terms(),
        )

    def degree(self):
        """Return the total degree, or ``None`` for an inhomogeneous form."""
        degrees = {len(basis) for basis in self._terms}
        if not degrees:
            return 0
        return degrees.pop() if len(degrees) == 1 else None

    def bidegrees(self):
        """Return all bidegrees occurring in the form."""
        dimension = self._chart.dimension()
        return {
            (
                sum(index < dimension for index in basis),
                sum(index >= dimension for index in basis),
            )
            for basis in self._terms
        }

    def bidegree(self):
        """Return the unique bidegree, or ``None`` for zero/mixed-type forms."""
        bidegrees = self.bidegrees()
        return next(iter(bidegrees)) if len(bidegrees) == 1 else None

    def wedge(self, other):
        """Return the exterior product with ``other``."""
        other = self._coerce(other)
        terms = {}
        for left_basis, left_coefficient in self._terms.items():
            for right_basis, right_coefficient in other._terms.items():
                sign, basis = _wedge_basis(left_basis, right_basis)
                if sign:
                    terms[basis] = terms.get(basis, SR.zero()) + (
                        sign * left_coefficient * right_coefficient
                    )
        return DifferentialForm(
            self._chart,
            terms,
            conditions=self._conditions + other._conditions,
            radial_indices=self._radial_indices + other._radial_indices,
        )

    def _coerce(self, other):
        if isinstance(other, DifferentialForm):
            if other.chart is not self._chart:
                raise ValueError("differential forms must belong to the same chart")
            return other
        return DifferentialForm(self._chart, {(): SR(other)})

    def _formal_conjugate_(self):
        dimension = self._chart.dimension()
        terms = {}
        for basis, coefficient in self._terms.items():
            conjugate_basis = tuple(
                index + dimension if index < dimension else index - dimension
                for index in basis
            )
            inversions = sum(
                conjugate_basis[i] > conjugate_basis[j]
                for i in range(len(conjugate_basis))
                for j in range(i + 1, len(conjugate_basis))
            )
            sign = -1 if inversions % 2 else 1
            canonical_basis = tuple(sorted(conjugate_basis))
            terms[canonical_basis] = terms.get(canonical_basis, SR.zero()) + (
                sign * bar(coefficient)
            )
        return DifferentialForm(
            self._chart,
            terms,
            conditions=self._conditions,
            radial_indices=self._radial_indices,
        )

    def __add__(self, other):
        other = self._coerce(other)
        terms = self.terms()
        for basis, coefficient in other._terms.items():
            terms[basis] = terms.get(basis, SR.zero()) + coefficient
        return DifferentialForm(
            self._chart,
            terms,
            conditions=self._conditions + other._conditions,
            radial_indices=self._radial_indices + other._radial_indices,
        )

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return DifferentialForm(
            self._chart,
            {basis: -coefficient for basis, coefficient in self._terms.items()},
            conditions=self._conditions,
            radial_indices=self._radial_indices,
        )

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, scalar):
        if isinstance(scalar, DifferentialForm):
            return NotImplemented
        return DifferentialForm(
            self._chart,
            {
                basis: coefficient * SR(scalar)
                for basis, coefficient in self._terms.items()
            },
            conditions=self._conditions,
            radial_indices=self._radial_indices,
        )

    def __rmul__(self, scalar):
        return self * scalar

    def __eq__(self, other):
        try:
            other = self._coerce(other)
        except (TypeError, ValueError):
            return False
        bases = set(self._terms).union(other._terms)
        return all(
            _is_zero(
                self._terms.get(basis, SR.zero())
                - other._terms.get(basis, SR.zero())
            )
            for basis in bases
        )

    def __repr__(self):
        if not self._terms:
            return "0"
        pieces = []
        for basis in sorted(self._terms, key=lambda item: (len(item), item)):
            coefficient = self._terms[basis]
            coefficient_text = _format_radial_coefficient(
                coefficient, self._chart, self._radial_indices
            )
            basis_text = " ∧ ".join(self._chart._basis_label(index) for index in basis)
            if not basis_text:
                pieces.append(coefficient_text)
            elif coefficient == 1:
                pieces.append(basis_text)
            elif coefficient == -1:
                pieces.append(f"-{basis_text}")
            else:
                pieces.append(f"({coefficient_text}) {basis_text}")
        result = " + ".join(pieces).replace("+ -", "- ")
        if self._conditions:
            result += f"  [valid where {', '.join(self._conditions)}]"
        return result


def _is_radial_in(expression, coordinate, conjugate_coordinate):
    variables = set(expression.variables())
    if coordinate not in variables or conjugate_coordinate not in variables:
        return False
    circle_derivative = (
        coordinate * expression.diff(coordinate)
        - conjugate_coordinate * expression.diff(conjugate_coordinate)
    )
    return _is_zero(circle_derivative)


def _simplify_radial_coefficient(expression, chart):
    simplified = SR(expression)
    used_indices = []
    nonzero_indices = []
    for index, (coordinate, conjugate_coordinate) in enumerate(
        zip(chart.coordinates(), chart.conjugate_coordinates())
    ):
        if not _is_radial_in(simplified, coordinate, conjugate_coordinate):
            continue
        radial_variable = SR.var(f"zz_sage_kahler_rho_{index}")
        candidate = simplified.subs(
            {coordinate: radial_variable / conjugate_coordinate}
        ).simplify()
        if (
            coordinate in candidate.variables()
            or conjugate_coordinate in candidate.variables()
        ):
            continue
        try:
            is_polynomial_at_zero = candidate.is_polynomial(radial_variable)
        except (AttributeError, TypeError, ValueError):
            is_polynomial_at_zero = False
        simplified = candidate.subs(
            {radial_variable: coordinate * conjugate_coordinate}
        )
        used_indices.append(index)
        if not is_polynomial_at_zero:
            nonzero_indices.append(index)
    return simplified.factor(), used_indices, nonzero_indices


def _format_radial_coefficient(expression, chart, radial_indices):
    display_expression = SR(expression)
    replacements = {}
    for index in radial_indices:
        coordinate = chart.coordinates()[index]
        conjugate_coordinate = chart.conjugate_coordinates()[index]
        display_name = f"zz_sage_kahler_norm_{index}"
        radius = SR.symbol(display_name, domain="positive")
        candidate = display_expression.subs(
            {coordinate: radius**2 / conjugate_coordinate}
        ).simplify()
        if (
            coordinate in candidate.variables()
            or conjugate_coordinate in candidate.variables()
        ):
            continue
        display_expression = candidate
        replacements[display_name] = f"|{coordinate}|"
    text = str(display_expression.factor())
    for variable_name, display_name in replacements.items():
        text = text.replace(variable_name, display_name)
    return text


def basis_form(chart, index):
    """Construct a coordinate basis one-form by its combined basis index."""
    return DifferentialForm(chart, {(index,): SR.one()})


def wedge(left, right):
    """Return the exterior product of two differential forms."""
    if not isinstance(left, DifferentialForm):
        if not isinstance(right, DifferentialForm):
            raise TypeError("at least one argument must be a differential form")
        left = right._coerce(left)
    return left.wedge(right)


def _as_form(value, chart):
    if isinstance(value, DifferentialForm):
        if chart is not None and value.chart is not chart:
            raise ValueError("the supplied chart does not match the differential form")
        return value
    if chart is None:
        chart = chart_for_expression(value)
    return DifferentialForm(chart, {(): SR(value)})


def _differentiate(value, kind, chart=None):
    form = _as_form(value, chart)
    chart = form.chart
    dimension = chart.dimension()
    coordinates = (
        chart.coordinates() if kind == "partial" else chart.conjugate_coordinates()
    )
    offset = 0 if kind == "partial" else dimension
    terms = {}
    for basis, coefficient in form._terms.items():
        for index, coordinate in enumerate(coordinates):
            derivative = coefficient.diff(coordinate)
            if _is_zero(derivative):
                continue
            sign, new_basis = _wedge_basis((offset + index,), basis)
            if sign:
                terms[new_basis] = terms.get(new_basis, SR.zero()) + sign * derivative
    return DifferentialForm(
        chart,
        terms,
        conditions=form._conditions,
        radial_indices=form._radial_indices,
    )


def partial(value, chart=None):
    """Apply the Dolbeault operator ``partial`` to an expression or form."""
    return _differentiate(value, "partial", chart=chart)


def dbar(value, chart=None):
    """Apply the Dolbeault operator ``dbar`` to an expression or form."""
    return _differentiate(value, "dbar", chart=chart)


def d(value, chart=None):
    """Apply the exterior differential ``d = partial + dbar``."""
    return partial(value, chart=chart) + dbar(value, chart=chart)


def ddbar(value, chart=None):
    """Apply ``partial(dbar(value))`` with readable radial simplification."""
    return partial(dbar(value, chart=chart)).simplify_radial()
