import pytest

from sage.all import I, SR, det, exp, latex, matrix, sqrt, var

from sage_kahler.all import ComplexChart, abs2, bar, d, dbar, ddbar, partial, wedge


def test_coordinate_basis_forms_have_the_expected_types():
    chart = ComplexChart(2, names=("a", "b"))

    assert chart.dz(0).bidegree() == (1, 0)
    assert chart.dzbar(1).bidegree() == (0, 1)
    assert wedge(chart.dz(0), chart.dzbar(1)).bidegree() == (1, 1)


def test_wedge_product_is_antisymmetric():
    chart = ComplexChart(2, names=("c", "e"))
    dc = chart.dz(0)
    de = chart.dz(1)

    assert wedge(dc, dc) == 0
    assert wedge(dc, de) == -wedge(de, dc)


def test_partial_and_dbar_differentiate_independent_variables():
    chart = ComplexChart(1, names=("p",))
    (p,) = chart.coordinates()

    assert partial(p * bar(p)) == bar(p) * chart.dz(0)
    assert dbar(p * bar(p)) == p * chart.dzbar(0)
    assert d(p * bar(p)) == partial(p * bar(p)) + dbar(p * bar(p))


def test_dolbeault_identities_on_functions_and_forms():
    chart = ComplexChart(2, names=("r", "s"))
    r, s = chart.coordinates()
    expression = r * s * bar(r) ** 2 + bar(s)

    assert partial(partial(expression)) == 0
    assert dbar(dbar(expression)) == 0
    assert partial(dbar(expression)) == -dbar(partial(expression))
    assert d(d(expression)) == 0


def test_ddbar_and_form_conjugation():
    chart = ComplexChart(1, names=("q",))
    (q,) = chart.coordinates()

    assert ddbar(q * bar(q)) == wedge(chart.dz(0), chart.dzbar(0))
    assert bar(chart.dz(0)) == chart.dzbar(0)
    assert bar(wedge(chart.dz(0), chart.dzbar(0))) == -wedge(
        chart.dz(0), chart.dzbar(0)
    )


def test_factor_simplifies_form_coefficients_and_ddbar_factors_automatically():
    chart = ComplexChart(1, names=("t",))
    (t,) = chart.coordinates()
    expected = (1 / (4 * sqrt(abs2(t)))) * wedge(
        chart.dz(0), chart.dzbar(0)
    )

    unfactored = partial(dbar(sqrt(abs2(t))))

    assert unfactored.factor() == expected
    assert ddbar(sqrt(abs2(t))) == expected
    assert repr(ddbar(sqrt(abs2(t)))) == (
        "(1/4/|t|) dt ∧ dbar(t)"
    )


def test_symbolic_radial_power_has_readable_output_and_raw_form():
    chart = ComplexChart(1, names=("u",))
    (u,) = chart.coordinates()
    exponent = var("alpha")

    result = ddbar(abs2(u) ** (exponent / 2))

    assert result.conditions() == ("u != 0",)
    assert "|u|^(alpha - 2)" in repr(result)
    assert "u_bar" not in repr(result.raw())
    assert "u_bar" in str(result.raw().terms()[(0, 1)])
    assert result == (
        exponent**2 / 4 * abs2(u) ** (exponent / 2 - 1)
    ) * wedge(chart.dz(0), chart.dzbar(0))


def test_smooth_radial_polynomial_does_not_exclude_the_origin():
    chart = ComplexChart(1, names=("m",))
    (m,) = chart.coordinates()

    result = ddbar(abs2(m) ** 2)

    assert result.conditions() == ()
    assert repr(result) == "(4*|m|^2) dm ∧ dbar(m)"


def test_real_coordinate_realization():
    chart = ComplexChart(1, names=("v",), real_names=("xv", "yv"))
    (v,) = chart.coordinates()
    ((xv, yv),) = chart.real_coordinates()

    assert chart.to_real(v) == xv + I * yv
    assert chart.to_real(bar(v)) == xv - I * yv
    assert chart.to_real(abs2(v)) == xv**2 + yv**2


def test_wirtinger_derivatives_agree_with_real_coordinate_derivatives():
    chart = ComplexChart(1, names=("h",), real_names=("xh", "yh"))
    (h,) = chart.coordinates()
    ((xh, yh),) = chart.real_coordinates()
    expression = (1 + abs2(h)) ** 3
    real_expression = chart.to_real(expression)

    partial_coefficient = partial(expression).terms()[(0,)]
    dbar_coefficient = dbar(expression).terms()[(1,)]
    ddbar_coefficient = ddbar(expression).terms()[(0, 1)]

    expected_partial = (real_expression.diff(xh) - I * real_expression.diff(yh)) / 2
    expected_dbar = (real_expression.diff(xh) + I * real_expression.diff(yh)) / 2
    expected_ddbar = (
        real_expression.diff(xh, 2) + real_expression.diff(yh, 2)
    ) / 4

    assert SR(chart.to_real(partial_coefficient) - expected_partial).simplify_full() == 0
    assert SR(chart.to_real(dbar_coefficient) - expected_dbar).simplify_full() == 0
    assert SR(chart.to_real(ddbar_coefficient) - expected_ddbar).simplify_full() == 0


def test_absolute_value_has_both_wirtinger_derivatives():
    chart = ComplexChart(1, names=("norm_z",))
    (z,) = chart.coordinates()
    radius = sqrt(abs2(z))

    assert partial(abs(z)) == (bar(z) / (2 * radius)) * chart.dz(0)
    assert dbar(abs(z)) == (z / (2 * radius)) * chart.dzbar(0)
    result = ddbar(abs(z))
    assert result == (1 / (4 * radius)) * wedge(chart.dz(0), chart.dzbar(0))
    assert result.conditions() == ("norm_z != 0",)
    assert "1/4/|norm_z|" in repr(result)
    assert partial(dbar(abs(z))) == -dbar(partial(abs(z)))


def test_absolute_value_powers_and_nested_expressions():
    chart = ComplexChart(1, names=("nested_z",))
    (z,) = chart.coordinates()
    alpha = var("norm_alpha")

    assert ddbar(abs(z)**2) == wedge(chart.dz(0), chart.dzbar(0))
    assert ddbar(abs(z)**2).conditions() == ()
    assert ddbar(abs(z)**alpha) == ddbar(abs2(z)**(alpha / 2))
    assert ddbar(exp(abs(z))) == ddbar(exp(sqrt(abs2(z))))
    assert ddbar(abs(z + I)) == (
        1 / (4 * sqrt(abs2(z + I)))
    ) * wedge(chart.dz(0), chart.dzbar(0))
    assert dbar(abs(abs(z) + z)) == dbar(sqrt(abs2(sqrt(abs2(z)) + z)))


def test_absolute_values_in_multivariable_form_coefficients():
    chart = ComplexChart(2, names=("norm_u", "norm_v"))
    u, v = chart.coordinates()
    expression = abs(u + I*v)
    formal = sqrt(abs2(u + I*v))

    assert dbar(expression * chart.dz(0)) == dbar(formal * chart.dz(0))
    assert partial(expression * chart.dzbar(1)) == partial(formal * chart.dzbar(1))
    assert d(d(expression)) == 0


def test_display_normalizes_mixed_coefficients_without_changing_them():
    chart = ComplexChart(1, names=("display_z",), bar_names=("internal_conj",))
    (z,) = chart.coordinates()
    form = partial(abs(z))
    before = form.terms()
    for value in (form, form.raw(), form.factor()):
        text = str(latex(value))
        assert r"\bar{\mathit{display}_{z}}" in text
        assert r"\left|\mathit{display}_{z}\right|" in text
        assert "internal" not in text
        assert "sqrt" not in text
        assert "internal_conj" not in repr(value)
    assert form.terms() == before
    assert form == (bar(z) / (2 * sqrt(abs2(z)))) * chart.dz(0)
    assert partial(dbar(abs(z))) == -dbar(partial(abs(z)))


def test_scalar_display_handles_nested_norms_and_multiple_coordinates():
    chart = ComplexChart(2, names=("z_7", "z_8"))
    z, w = chart.coordinates()
    expression = bar(z) / sqrt(abs2(z)) + exp(sqrt(abs2(w)))
    displayed = chart.display(expression)
    text = str(latex(displayed))
    assert r"\bar{z_{7}}" in text
    assert r"\left|z_{7}\right|" in text
    assert r"\left|z_{8}\right|" in text
    assert "sqrt" not in text
    assert displayed.terms()[()] == expression
    assert str(latex(bar(z))) == r"{\bar{z_{7}}}"
    assert str(bar(z)) == "z_7_bar"
    # A mixed product is not a squared norm.
    mixed = str(latex(chart.display(sqrt(z * bar(w)))))
    assert r"\sqrt" in mixed
    assert r"\left|" not in mixed



def test_metric_coefficient_matrix_and_sage_determinant():
    chart = ComplexChart(2, names=("metric_z", "metric_w"))
    z, w = chart.coordinates()
    potential = 2*abs(z)**2 + 3*abs(w)**2 + I*z*bar(w) - I*w*bar(z)
    form = ddbar(potential)

    assert form.coefficient_matrix() == matrix(SR, [[2, I], [-I, 3]])
    assert det(form) == 5
    assert form.det() == form.determinant() == 5
    # Extracting a matrix must not expose mutable form storage.
    extracted = form.coefficient_matrix()
    extracted[0, 0] = 99
    assert det(form) == 5


def test_determinant_of_variable_and_degenerate_metrics():
    chart = ComplexChart(2, names=("hessian_z", "hessian_w"))
    z, w = chart.coordinates()
    form = ddbar(abs(z)**4 + abs(w)**2)
    assert (det(form) - 4*abs2(z)).simplify_full() == 0
    assert det(ddbar(abs(z)**2)) == 0
    zero = ddbar(0, chart=chart)
    assert zero.coefficient_matrix() == matrix(SR, 2, 2, 0)
    assert det(zero) == 0

    line = ComplexChart(1, names=("line_z",))
    (t,) = line.coordinates()
    radial = ddbar(abs(t))
    assert (det(radial) - 1/(4*sqrt(abs2(t)))).simplify_full() == 0
    assert radial.conditions() == ("line_z != 0",)


def test_coefficient_matrix_rejects_other_form_types():
    chart = ComplexChart(2, names=("type_z", "type_w"))
    forms = (
        chart.dz(0),
        wedge(chart.dz(0), chart.dz(1)),
        wedge(chart.dzbar(0), chart.dzbar(1)),
        wedge(chart.dz(0), chart.dzbar(0)) + 1,
    )
    for form in forms:
        with pytest.raises(ValueError, match=r"\(1, 1\)"):
            det(form)



def test_metric_matrix_and_determinant_norm_display():
    chart = ComplexChart(2, names=("display_metric_z", "display_metric_w"))
    z, w = chart.coordinates()
    metric = ddbar(abs(z)**4 + abs(w)**2)
    coefficients = metric.coefficient_matrix()
    original = matrix(coefficients)
    presentation = chart.display(coefficients)
    assert "|display_metric_z|^2" in repr(presentation)
    assert rf"\left|{latex(z)}\right|^{{2}}" in str(latex(presentation))
    assert presentation.matrix() == original
    assert coefficients == original
    assert repr(metric.coefficient_matrix(display=True)) == repr(presentation)
    assert "|display_metric_z|^2" in repr(metric.det(display=True))
    assert metric.det(display=True).terms()[()] == det(metric)
    assert str(latex(chart.display(det(metric)))) == str(latex(metric.det(display=True)))
    # Matrix computations still use the independent formal coordinates.
    assert coefficients[0, 0].diff(bar(z)) == 4*z
    assert det(metric).diff(bar(z)) == 4*z
    # Mixed off-diagonal entries retain their distinct conjugate indices.
    mixed = chart.display(matrix(SR, [[z*bar(w), sqrt(abs2(z))]]))
    rendered = str(latex(mixed))
    assert rf"\bar{{{latex(w)}}}" in rendered
    assert "sqrt" not in rendered
    assert mixed.matrix().ncols() == 2


def test_rational_coefficients_are_printed_without_parentheses():
    chart = ComplexChart(2, names=("print_z", "print_w"))
    z, w = chart.coordinates()

    assert repr(2 * chart.dz(0)) == "2 dprint_z"
    assert repr(-chart.dz(0) * (SR(1) / 2) + chart.dz(1)) == "-1/2 dprint_z + dprint_w"
    assert repr((1 + 2*I) * chart.dz(0)) == "(2*I + 1) dprint_z"
    assert repr(z * chart.dz(1)) == "(print_z) dprint_w"


def test_products_of_norms_cancel_in_coefficients():
    chart = ComplexChart(2, names=("prod_z", "prod_w"))
    z, w = chart.coordinates()
    # |z| |w| = |z w| on the chart, although the formal radicals differ.
    vanishing = sqrt(abs2(z)) * sqrt(abs2(w)) - sqrt(abs2(z * w))

    assert chart.display(vanishing).terms() == {}
    assert chart.display(vanishing) == 0
    assert sqrt(abs2(z)) * chart.dz(0) != sqrt(abs2(w)) * chart.dz(0)


def test_display_and_simplification_do_not_leak_temporary_variables():
    from sage.symbolic.assumptions import assumptions

    chart = ComplexChart(2, names=("leak_z", "leak_w"))
    z, w = chart.coordinates()
    metric = ddbar(sqrt(abs(z)**2 + abs(w)**2 + 2*abs(z*w - 1) + 2))
    before = (len(assumptions()), len(SR.symbols))
    for _ in range(3):
        repr(metric), latex(metric), repr(det(metric))
        ddbar(abs(z)**2 + abs(w)**2)

    assert (len(assumptions()), len(SR.symbols)) == before
