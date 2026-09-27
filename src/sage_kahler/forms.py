"""Sparse complex differential forms and Dolbeault operators."""

import operator
import re
from copy import copy

from sage.all import QQ, SR, ZZ, PolynomialRing, QuadraticField, latex, matrix, sqrt
from sage.functions.other import abs_symbolic
from sage.symbolic.operators import add_vararg, mul_vararg

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

    def coefficient_matrix(self, *, display=None):
        """Return the matrix of a (1, 1) form in ``dz_i wedge dbar(z_j)``.

        Rows are holomorphic indices and columns are antiholomorphic indices,
        in chart coordinate order. The zero form gives the zero matrix.
        For ``ddbar(f)`` this is the complex Hessian of ``f``, with no factor
        of ``I``. No Hermitian or positivity assumption is imposed.
        The native Sage matrix subclass formats norms automatically. Use
        ``raw()`` (or ``display=False``) for an ordinary Sage matrix.
        """
        if self.bidegrees() - {(1, 1)}:
            raise ValueError("coefficient_matrix requires a (1, 1) form")
        dimension = self._chart.dimension()
        result = matrix(SR, dimension, dimension, lambda i, j: self._terms.get(
            (i, dimension + j), SR.zero()
        ))
        if display is False:
            return result
        if display is True:
            return self._chart.display(result)
        from .results import ChartMatrix

        return ChartMatrix.from_matrix(result, self._chart)

    def determinant(self, *, display=None):
        """Return the symbolic determinant of the (1, 1) coefficient matrix.

        The scalar result is defined wherever the original form is defined;
        consult the form's ``conditions()`` for recorded domain restrictions.
        Norm notation is automatic; ``raw()`` retrieves the formal expression.
        """
        from .results import ChartExpression

        result = _simplify_scalar(
            self.coefficient_matrix(display=False).determinant(), self._chart
        )
        if display is False:
            return result
        if display is True:
            # Preserve the previous explicit presentation-wrapper API.
            return self._chart.display(result)
        return ChartExpression(SR, result, self._chart)

    def det(self, *, display=None):
        """Support Sage's ``det(form)`` for (1, 1) forms."""
        return self.determinant(display=display)

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
            coefficient_text = _format_coefficient(
                coefficient, self._chart
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
        return " + ".join(pieces).replace("+ -", "- ")

    def _latex_(self):
        """Render in Sage's ``%display latex`` mode, omitting conditions."""
        if not self._terms:
            return "0"
        pieces = []
        for basis in sorted(self._terms, key=lambda item: (len(item), item)):
            coefficient = self._terms[basis]
            coefficient_text = _format_coefficient(
                coefficient, self._chart, latex_mode=True
            )
            basis_text = r" \wedge ".join(
                self._chart._basis_latex(index) for index in basis
            )
            if not basis_text:
                pieces.append(coefficient_text)
            elif coefficient == 1:
                pieces.append(basis_text)
            elif coefficient == -1:
                pieces.append(f"-{basis_text}")
            else:
                # Group sums so the coefficient multiplies the entire basis.
                if coefficient.operator() == add_vararg:
                    coefficient_text = rf"\left({coefficient_text}\right)"
                pieces.append(rf"{coefficient_text}\,{basis_text}")
        return " + ".join(pieces).replace("+ -", "- ")


def _simplify_radial_coefficient(expression, chart):
    simplified = SR(expression)
    used_indices = []
    nonzero_indices = []
    for index, (coordinate, conjugate_coordinate) in enumerate(
        zip(chart.coordinates(), chart.conjugate_coordinates())
    ):
        variables = simplified.variables()
        if coordinate not in variables or conjugate_coordinate not in variables:
            continue
        # Test radial dependence by elimination, not by differentiating the
        # already differentiated coefficient and invoking simplify_full().
        radial_variable = SR.temp_var()
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


def _tree_size(value):
    return 1 + sum(_tree_size(operand) for operand in value.operands())


def _compact_polynomial_sums(expression):
    """Factor polynomial subexpressions only when their trees become smaller.

    This recovers expanded perfect powers in denominators without expanding
    already compact radial sums or running a general radical simplifier.
    """
    def compact(value):
        operands = value.operands()
        if operands:
            value = value.operator()(*(compact(operand) for operand in operands))
        if value.operator() == add_vararg and all(
            value.is_polynomial(variable) for variable in value.variables()
        ):
            candidate = value.factor()
            if _tree_size(candidate) < _tree_size(value):
                return candidate
        return value

    return compact(expression)


def _half_integer_exponent(value):
    """Return ``2*e`` when ``value`` is ``base^e`` with odd ``2*e``, else None."""
    if value.operator() is not operator.pow:
        return None
    try:
        doubled = ZZ(2 * QQ(value.operands()[1]))
    except (TypeError, ValueError):
        return None
    return doubled if doubled % 2 else None


def _cancel_radicals(expression, chart):
    """Cancel square roots of coordinate expressions using ``t^2 = A``.

    Each ``A^(k/2)`` becomes ``t^k`` for a new symbol ``t``, nested radicands
    first. The resulting rational function is reduced modulo the relations,
    its denominator is rationalized one radical at a time, and common factors
    are cancelled. Returns None unless the expression is rational in the
    coordinates, parameters and radicals.
    """
    coordinates = set(chart.coordinates() + chart.conjugate_coordinates())
    radicals = []  # (symbol, radicand in earlier symbols), innermost first

    def extract(value):
        if not coordinates.intersection(value.variables()):
            return value
        operands = value.operands()
        if not operands:
            return value
        doubled = _half_integer_exponent(value)
        if doubled is not None:
            radicand = extract(operands[0])
            for symbol, known in radicals:
                if (known - radicand).expand().is_zero():
                    return symbol**doubled
            symbol = SR.temp_var()
            radicals.append((symbol, radicand))
            return symbol**doubled
        return value.operator()(*(extract(operand) for operand in operands))

    rational = extract(SR(expression))
    if not radicals:
        return None
    variables = sorted(
        set(rational.variables()).union(
            *(radicand.variables() for _, radicand in radicals)
        ),
        key=str,
    )

    def convert(value, generators, field, base_ring):
        if value in generators:
            return generators[value]
        operands = value.operands()
        if not operands:
            return field(base_ring(value))
        converted = [convert(o, generators, field, base_ring) for o in operands]
        if value.operator() == add_vararg:
            return sum(converted, field.zero())
        if value.operator() == mul_vararg:
            result = field.one()
            for operand in converted:
                result *= operand
            return result
        if value.operator() is operator.pow:
            return converted[0] ** ZZ(operands[1])
        raise TypeError(f"not rational: {value}")

    for base_ring in (QQ, QuadraticField(-1, "I")):
        ring = PolynomialRing(base_ring, [f"v{i}" for i in range(len(variables))])
        field = ring.fraction_field()
        generators = dict(zip(variables, map(field, ring.gens())))
        try:
            value = convert(rational, generators, field, base_ring)
            relations = [
                (
                    ring(generators[symbol]),
                    ring(convert(radicand, generators, field, base_ring)),
                )
                for symbol, radicand in radicals
            ]
        except (TypeError, ValueError):
            continue
        break
    else:
        return None

    def reduce_relations(polynomial):
        # Outermost radicands contain inner symbols, so reduce outside in.
        for symbol, radicand in reversed(relations):
            polynomial = sum(
                (
                    polynomial.coefficient({symbol: n})
                    * radicand ** (n // 2) * symbol ** (n % 2)
                    for n in range(polynomial.degree(symbol) + 1)
                ),
                ring.zero(),
            )
        return polynomial

    numerator = reduce_relations(value.numerator())
    denominator = reduce_relations(value.denominator())
    if numerator.is_zero():
        return SR.zero()
    for symbol, radicand in reversed(relations):
        constant = denominator.coefficient({symbol: 0})
        linear = denominator.coefficient({symbol: 1})
        if linear.is_zero():
            continue
        conjugate = constant - linear * symbol
        numerator = reduce_relations(numerator * conjugate)
        denominator = reduce_relations(denominator * conjugate)
    common = numerator.gcd(denominator)
    numerator, denominator = numerator // common, denominator // common
    # Prefer 1/sqrt(A) to sqrt(A)/A.
    for symbol, radicand in relations:
        if (
            numerator.coefficient({symbol: 0}).is_zero()
            and (denominator % radicand).is_zero()
        ):
            numerator = numerator // symbol
            denominator = symbol * (denominator // radicand)

    values = {}
    for symbol, radicand in radicals:
        radicand = radicand.subs(values)
        factored = radicand.factor()
        if _tree_size(factored) < _tree_size(radicand):
            radicand = factored
        values[symbol] = sqrt(radicand)
    sr_values = [values.get(variable, variable) for variable in variables]

    def to_sr(polynomial):
        if polynomial.is_constant():
            return SR(polynomial.constant_coefficient())
        factorization = polynomial.factor()
        result = SR(ring(factorization.unit()).constant_coefficient())
        for factor, exponent in factorization:
            result *= SR(ring(factor)(*sr_values)) ** exponent
        return result

    return to_sr(numerator) / to_sr(denominator)


def _simplify_scalar(expression, chart):
    """Simplify a scalar, cancelling coordinate radicals when possible."""
    cancelled = _cancel_radicals(expression, chart)
    if cancelled is not None:
        return cancelled
    return SR(expression).simplify_full()


def _conjugate_norm(radicand, chart):
    """Write ``radicand`` as ``c * prod (F_i * bar(F_i))^e_i`` with ``c > 0``.

    Returns ``(c, [(F_i, e_i), ...])`` or None. Each ``F_i`` is the factor of
    its conjugate pair with fewer conjugate coordinates.
    """
    conjugates = set(chart.conjugate_coordinates())
    try:
        factor_list = radicand.factor_list()
    except (AttributeError, TypeError, ValueError):
        return None
    constant = QQ.one()
    remaining = []
    for factor, exponent in factor_list:
        if factor.variables():
            remaining.append((factor, exponent))
            continue
        try:
            constant *= QQ(factor) ** exponent
        except (TypeError, ValueError):
            return None
    norms = []
    while remaining:
        factor, exponent = remaining.pop(0)
        conjugate = bar(factor)
        for index, (other, other_exponent) in enumerate(remaining):
            if other_exponent != exponent:
                continue
            if (other - conjugate).expand().is_zero():
                sign = 1
            elif (other + conjugate).expand().is_zero():
                sign = -1
            else:
                continue
            del remaining[index]
            constant *= sign**exponent
            norms.append((
                min(
                    (factor, other),
                    key=lambda f: len(conjugates.intersection(f.variables())),
                ),
                exponent,
            ))
            break
        else:
            return None
    if not norms or constant <= 0:
        return None
    return constant, norms


def _format_coefficient(expression, chart, *, latex_mode=False):
    """Format a display-only copy, including norms inside mixed coefficients.

    Temporary positive radii permit simplification on the complex chart without
    replacing the independent coordinates in the stored symbolic expression.
    The radial metadata is deliberately not required for display. Square
    roots of ``c * F * bar(F)`` likewise display as ``sqrt(c) * |F|``.
    """
    pairs = tuple(zip(chart.coordinates(), chart.conjugate_coordinates()))
    radii = tuple(SR.temp_var(domain="positive") for _ in pairs)
    norms = []  # (F, positive token displayed as |F|)

    def norm_token(factor):
        for known, token in norms:
            if (known - factor).expand().is_zero():
                return token
        token = SR.temp_var(domain="positive")
        norms.append((factor, token))
        return token

    def rewrite(value):
        operands = value.operands()
        if operands:
            value = value.operator()(*(rewrite(operand) for operand in operands))
        for (coordinate, conjugate_coordinate), radius in zip(pairs, radii):
            variables = value.variables()
            if coordinate not in variables or conjugate_coordinate not in variables:
                continue
            candidate = value.subs(
                {coordinate: radius**2 / conjugate_coordinate}
            ).simplify()
            if (
                coordinate not in candidate.variables()
                and conjugate_coordinate not in candidate.variables()
            ):
                value = candidate
        doubled = _half_integer_exponent(value)
        if doubled is not None:
            norm = _conjugate_norm(value.operands()[0], chart)
            if norm is not None:
                constant, factors = norm
                value = sqrt(constant) ** doubled
                for factor, exponent in factors:
                    value *= norm_token(factor) ** (exponent * doubled)
        return value

    display_expression = _compact_polynomial_sums(rewrite(SR(expression)))
    # Use unique tokens rather than replacing coordinate-name substrings.
    conjugate_tokens = tuple(SR.temp_var() for _ in pairs)
    display_expression = display_expression.subs({
        conjugate: token
        for (_, conjugate), token in zip(pairs, conjugate_tokens)
    })
    render = (lambda value: str(latex(value))) if latex_mode else str
    replacements = {}
    for (coordinate, _), radius, token in zip(pairs, radii, conjugate_tokens):
        if latex_mode:
            replacements[render(radius)] = rf"\left|{latex(coordinate)}\right|"
            replacements[render(token)] = rf"\bar{{{latex(coordinate)}}}"
        else:
            replacements[render(radius)] = f"|{coordinate}|"
            replacements[render(token)] = f"bar({coordinate})"
    for factor, token in norms:
        text = _format_coefficient(factor, chart, latex_mode=latex_mode)
        replacements[render(token)] = (
            rf"\left|{text}\right|" if latex_mode else f"|{text}|"
        )
    # Replace all tokens in one pass, so inserted coordinate names stay intact.
    pattern = "|".join(
        re.escape(token) for token in sorted(replacements, key=len, reverse=True)
    )
    return re.sub(
        pattern, lambda match: replacements[match.group()], render(display_expression)
    )


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


def _formalize_absolute_values(expression, chart):
    """Expose the conjugate-coordinate dependence hidden inside Sage's abs.

    On the complex chart, ``abs(f) = sqrt(f * bar(f))``. Rewrite before
    either Dolbeault derivative, including abs nested in other functions or
    form coefficients. Absolute values of coordinate-independent parameters
    are left alone. Derivatives are classical, on the smooth locus.
    """
    coordinates = set(chart.coordinates() + chart.conjugate_coordinates())

    def rewrite(value):
        if not coordinates.intersection(value.variables()):
            return value
        operands = value.operands()
        if not operands:
            return value
        if value.operator() == abs_symbolic:
            argument = rewrite(operands[0])
            return sqrt(argument * bar(argument))
        return value.operator()(*(rewrite(operand) for operand in operands))

    return rewrite(expression)


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
        coefficient = _formalize_absolute_values(coefficient, chart)
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


def _simplify_univariate_derivative(expression):
    """Cancel radical denominators before expanding back into coordinates."""
    from sympy import cancel, expand, radsimp

    try:
        # Restrict rationalization to one radical term. This is algebraic
        # cancellation, without the trigonometric passes of simplify_full.
        symbolic = expression._sympy_()
        return SR(cancel(expand(radsimp(symbolic, max_terms=1))))
    except (TypeError, ValueError, NotImplementedError):
        return expression.simplify_rational()


def _radial_hessian(form):
    """Compute the complex Hessian of F(sum |z_i|^2) in one variable.

    Structural recognition is intentionally conservative: if the substitution
    leaves any chart coordinates, the ordinary Dolbeault path handles it.
    """
    chart = form.chart
    if chart.dimension() < 2 or set(form._terms) != {()}:
        return None
    coordinates = chart.coordinates()
    conjugates = chart.conjugate_coordinates()
    radius_squared = sum(z * zb for z, zb in zip(coordinates, conjugates))
    rho = SR.temp_var(domain="positive")
    expression = _formalize_absolute_values(form._terms[()], chart)
    remainder = radius_squared - coordinates[0] * conjugates[0]
    profile = expression.subs({coordinates[0]: (rho - remainder) / conjugates[0]})
    if set(profile.variables()).intersection(coordinates + conjugates):
        return None
    first_raw = profile.diff(rho)
    first = _simplify_univariate_derivative(first_raw).factor()
    second_raw = first_raw.diff(rho)
    second = _simplify_univariate_derivative(first.diff(rho)).factor()
    conditions = list(form._conditions)
    # Record radial singularities without excluding individual coordinate
    # hyperplanes: the total squared radius, not every z_i, must be nonzero.
    def singular_at_origin(derivative):
        try:
            return derivative.denominator().subs({rho: 0}).is_trivial_zero()
        except (ValueError, ZeroDivisionError):
            return True

    if any(singular_at_origin(derivative) for derivative in (first, second)):
        condition = " + ".join(f"|{z}|^2" for z in coordinates) + " != 0"
        if condition not in conditions:
            conditions.append(condition)
    first_raw, second_raw = (
        derivative.subs({rho: radius_squared})
        for derivative in (first_raw, second_raw)
    )
    terms = {}
    raw_terms = {}
    for i, zb in enumerate(conjugates):
        for j, z in enumerate(coordinates):
            basis = (i, chart.dimension() + j)
            # Combine diagonal terms while rho is still a single variable;
            # factoring after substitution expands powers of the total norm.
            coefficient = (second * zb * z + (first if i == j else 0)).factor()
            terms[basis] = coefficient.subs({rho: radius_squared})
            raw_terms[basis] = second_raw * zb * z + (first_raw if i == j else 0)
    return DifferentialForm(chart, terms, conditions=conditions, raw_terms=raw_terms)


def ddbar(value, chart=None):
    """Apply partial(dbar(value)), using a compact Hessian for radial scalars."""
    form = _as_form(value, chart)
    radial = _radial_hessian(form)
    if radial is not None:
        return radial
    return partial(dbar(form)).simplify_radial()



class _MatrixDisplay:
    """A chart-aware matrix presentation, preserving a copy of its entries."""

    def __init__(self, chart, value):
        self._chart = chart
        self._matrix = copy(value)

    def matrix(self):
        """Return a Sage matrix copy for further computation."""
        return copy(self._matrix)

    def _latex_(self):
        rows = [
            " & ".join(
                _format_coefficient(entry, self._chart, latex_mode=True)
                for entry in row
            )
            for row in self._matrix.rows()
        ]
        columns = "r" * self._matrix.ncols()
        return (
            rf"\left(\begin{{array}}{{{columns}}}"
            + r" \\ ".join(rows)
            + r"\end{array}\right)"
        )

    def __repr__(self):
        rows = [
            [_format_coefficient(entry, self._chart) for entry in row]
            for row in self._matrix.rows()
        ]
        widths = [max((len(row[j]) for row in rows), default=0)
                  for j in range(self._matrix.ncols())]
        return "\n".join(
            "[" + "  ".join(entry.rjust(width) for entry, width in zip(row, widths)) + "]"
            for row in rows
        ) or "[]"
