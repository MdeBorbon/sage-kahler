import pytest
from sage.all import SR, det, latex, log, sqrt

from sage_kahler.all import ComplexChart, abs2, bar, ddbar, partial, dbar
from sage_kahler.forms import _formalize_absolute_values


def potential(s):
    return sqrt(s**2 + 1) + log(s / (sqrt(s**2 + 1) + 1))


def test_nested_absolute_values_are_fully_rewritten():
    chart = ComplexChart(2, names=("nested_radial_z", "nested_radial_w"))
    z, w = chart.coordinates()
    expression = potential(abs(z)**2 + abs(w)**2)
    expected = potential(abs2(z) + abs2(w))
    assert _formalize_absolute_values(expression, chart) == expected


@pytest.mark.parametrize("use_abs", [False, True])
def test_radical_log_potential_uses_compact_radial_hessian(monkeypatch, use_abs):
    import sage_kahler.forms as forms

    chart = ComplexChart(2, names=("radial_z", "radial_w"))
    z, w = chart.coordinates()
    s = abs2(z) + abs2(w)
    input_s = abs(z)**2 + abs(w)**2 if use_abs else s

    def no_expanded_path(*args, **kwds):
        raise AssertionError("radial potentials must avoid expanded differentiation")

    monkeypatch.setattr(forms, "_differentiate", no_expanded_path)
    result = ddbar(potential(input_s))
    coefficients = result.coefficient_matrix(display=False)
    first = sqrt(s**2 + 1) / s
    second = -1 / (s**2 * sqrt(s**2 + 1))
    for i, zi in enumerate((z, w)):
        for j, zj in enumerate((z, w)):
            expected = second * bar(zi) * zj + (first if i == j else 0)
            assert (coefficients[i, j] - expected).simplify_full() == 0
    assert det(result) == 1
    assert result.conditions() == ("|radial_z|^2 + |radial_w|^2 != 0",)
    # Individual coordinate hyperplanes are allowed; only the origin is singular.
    axis = {z: 0, bar(z): 0, w: 1, bar(w): 1}
    assert (coefficients[0, 0].subs(axis) - sqrt(2)).simplify_full() == 0
    assert (coefficients[1, 1].subs(axis) - 1/sqrt(2)).simplify_full() == 0
    rendered = str(latex(result))
    assert "sqrt" in rendered
    assert "symbol" not in rendered
    # The raw terms retain the same Hessian before one-variable cancellation.
    for basis, coefficient in result.raw().terms().items():
        i, combined_j = basis
        assert (coefficient.subs(axis) - coefficients[i, combined_j - 2].subs(axis)).simplify_full() == 0


def test_radial_polynomials_and_nonradial_fallback():
    chart = ComplexChart(3, names=("poly_z", "poly_w", "poly_v"))
    z, w, v = chart.coordinates()
    s = sum(abs2(t) for t in (z, w, v))
    result = ddbar(s**3)
    assert result.conditions() == ()
    assert result == partial(dbar(s**3))
    nonradial = z*bar(w)**2 + abs2(v)**2
    assert ddbar(nonradial) == partial(dbar(nonradial))


def test_radial_output_keeps_powers_of_total_squared_radius_grouped():
    chart = ComplexChart(2, names=("compact_z", "compact_w"))
    z, w = chart.coordinates()
    result = ddbar(potential(abs(z)**2 + abs(w)**2))
    original = result.terms()
    # All display entry points must keep the cube of the sum, not expand it
    # into degree-six monomials. No undisclosed dummy radius may be printed.
    for value in (result, result.coefficient_matrix(), result.coefficient_matrix()[0, 0]):
        rendered = str(latex(value))
        assert r"\right)}^{3}" in rendered
        assert r"\right)}^{2} + 1" in rendered
        assert "^{6}" not in rendered
        assert "symbol" not in rendered
        assert "^6" not in str(value)
    assert result.terms() == original
    assert det(result) == 1



def test_determinant_display_recovers_expanded_polynomial_powers():
    chart = ComplexChart(2, names=("factor_z", "factor_w"))
    z, w = chart.coordinates()
    s = abs2(z) + abs2(w)
    form = ddbar(log(1 + s))
    determinant = det(form)
    original = determinant.raw()
    expanded = 1 / ((1 + s)**3).expand()
    for value in (
        determinant, det(form.coefficient_matrix()),
        chart.display(expanded), form.det(display=True),
    ):
        rendered = str(latex(value))
        assert r"\right)}^{3}" in rendered
        assert "^{6}" not in rendered
        assert "^6" not in str(value)
        assert "symbol" not in rendered
    assert determinant.raw() == original
    assert (original - 1/(1+s)**3).simplify_full() == 0


def test_determinant_cancels_nested_radicals_of_non_radial_norms():
    chart = ComplexChart(2, names=("conic_z", "conic_w"))
    z, w = chart.coordinates()
    phi = sqrt(abs(z)**2 + abs(w)**2 + 2*abs(z*w - 1) + 2)
    metric = ddbar(phi)
    expected = 1 / (16 * sqrt(abs2(z*w - 1)))
    for determinant in (det(metric), det(metric.coefficient_matrix())):
        assert bool(determinant.raw() == expected)
        assert str(determinant) == "1/16/|conic_w*conic_z - 1|"
        assert str(latex(determinant)) == (
            rf"\frac{{1}}{{16 \, \left|{latex(w*z - 1)}\right|}}"
        )
    # Hessian coefficients display sqrt(F * bar(F)) as |F| too.
    assert "sqrt" not in str(metric)
    assert "|conic_w*conic_z - 1|" in str(metric)


def test_radical_cancellation_agrees_with_the_unsimplified_determinant():
    from sage.all import CC, I

    chart = ComplexChart(2, names=("check_z", "check_w"))
    z, w = chart.coordinates()
    s = abs2(z) + abs2(w)
    form = ddbar(s + I*z*bar(w)**2 - I*bar(z)*w**2 + sqrt(1 + s + abs2(z - w)))
    raw = form.coefficient_matrix(display=False).determinant()
    point = {
        z: 1/2 + I/3, bar(z): 1/2 - I/3,
        w: -1/5 + 2*I/7, bar(w): -1/5 - 2*I/7,
    }
    assert abs(CC(det(form).raw().subs(point)) - CC(raw.subs(point))) < 1e-12
