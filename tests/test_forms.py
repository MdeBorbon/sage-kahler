from sage.all import I, SR, sqrt, var

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
        "(1/4/|t|) dt ∧ dbar(t)  [valid where t != 0]"
    )


def test_symbolic_radial_power_has_readable_output_and_raw_form():
    chart = ComplexChart(1, names=("u",))
    (u,) = chart.coordinates()
    exponent = var("alpha")

    result = ddbar(abs2(u) ** (exponent / 2))

    assert result.conditions() == ("u != 0",)
    assert "|u|^(alpha - 2)" in repr(result)
    assert "u*u_bar" in repr(result.raw())
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
