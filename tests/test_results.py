from sage.all import SR, det, latex, log, matrix
from sage.symbolic.expression import Expression
from sage.structure.element import Matrix

from sage_kahler.all import ComplexChart, abs2, bar, ddbar
from sage_kahler.results import ChartExpression, ChartMatrix


def metric_example():
    chart = ComplexChart(2, names=("auto_z", "auto_w"))
    z, w = chart.coordinates()
    return chart, z, w, ddbar(abs(z)**4 + abs(w)**2)


def test_native_metric_results_display_automatically():
    chart, z, w, form = metric_example()
    coefficients = form.coefficient_matrix()
    scalar = det(form)
    assert isinstance(coefficients, Matrix)
    assert isinstance(scalar, Expression)
    norm = rf"\left|{latex(z)}\right|^{{2}}"
    for value in (coefficients, scalar, det(coefficients), coefficients[0, 0]):
        assert norm in str(latex(value))
        assert "|auto_z|^2" in str(value)
        assert "auto_z_bar" not in str(value)
    assert type(scalar.raw()) is Expression
    assert not isinstance(coefficients.raw(), ChartMatrix)
    assert scalar.raw() == 4*abs2(z)
    assert coefficients.raw() == matrix(SR, [[4*abs2(z), 0], [0, 1]])
    assert form.coefficient_matrix(display=False) == coefficients.raw()
    assert type(form.det(display=False)) is Expression
    coefficients.raw()[0, 0] = 123
    assert coefficients[0, 0] == 4*abs2(z)


def test_scalar_computations_preserve_formal_derivatives_and_common_display():
    chart, z, w, form = metric_example()
    scalar = det(form)
    raw = scalar.raw()
    pairs = (
        (scalar + scalar, raw + raw), (scalar + 1, raw + 1),
        (1 + scalar, 1 + raw), (scalar - 1, raw - 1),
        (1 - scalar, 1 - raw), (scalar * scalar, raw * raw),
        (2 * scalar, 2 * raw), (scalar / 2, raw / 2),
        (1 / scalar, 1 / raw), (scalar**2, raw**2), (-scalar, -raw),
    )
    for actual, expected in pairs:
        assert isinstance(actual, ChartExpression)
        assert (actual.raw() - expected).simplify_full() == 0
        assert "|auto_z|" in str(actual)
    assert scalar.diff(bar(z)) == 4*z
    assert scalar.subs({z: 2, bar(z): 2}) == 16
    assert scalar.factor().raw() == raw.factor()
    assert scalar.expand().raw() == raw.expand()
    assert scalar.simplify_full().raw() == raw.simplify_full()
    assert log(scalar) == log(raw)
    assert ddbar(log(scalar)) == ddbar(log(raw))
    assert ddbar(scalar) == ddbar(raw)
    assert scalar * chart.dz(0) == raw * chart.dz(0)
    assert SR(scalar).parent() is SR


def test_matrix_operations_and_entries_remain_computational():
    chart, z, w, form = metric_example()
    value = form.coefficient_matrix()
    raw = value.raw()
    for result, expected in (
        (value + value, raw + raw),
        (value * value, raw * raw),
        (value.transpose(), raw.transpose()),
        (value.inverse(), raw.inverse()),
        (value.subs({z: w, bar(z): bar(w)}), raw.subs({z: w, bar(z): bar(w)})),
    ):
        assert result == expected
        assert "|auto_" in str(result)
    assert det(value) == det(form)
    assert value[0, 0].diff(bar(z)) == 4*z
    value[0, 0] = 7
    assert det(value) == 7
    assert det(form) == 4*abs2(z)
    assert det(ddbar(0, chart=chart)) == 0
