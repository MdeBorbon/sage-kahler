from sage.all import I

from sage_kahler.all import ComplexChart, abs2, bar


def test_chart_creates_formal_conjugate_coordinates():
    chart = ComplexChart(2, names=("z", "w"))
    z, w = chart.coordinates()
    z_bar, w_bar = chart.conjugate_coordinates()

    assert tuple(map(str, (z_bar, w_bar))) == ("z_bar", "w_bar")
    assert bar(z) == z_bar
    assert bar(w) == w_bar


def test_bar_is_an_involution_and_conjugates_expressions():
    chart = ComplexChart(2, names=("u", "v"))
    u, v = chart.coordinates()
    u_bar, v_bar = chart.conjugate_coordinates()
    expression = (2 + I) * u**2 + v_bar

    assert bar(expression) == (2 - I) * u_bar**2 + v
    assert bar(bar(expression)) == expression


def test_abs2_uses_formal_conjugation():
    chart = ComplexChart(1, names=("xi",))
    (xi,) = chart.coordinates()
    (xi_bar,) = chart.conjugate_coordinates()

    assert abs2(xi) == xi * xi_bar
