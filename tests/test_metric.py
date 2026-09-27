import pytest

from sage.all import I, SR, identity_matrix, log, matrix, var

from sage_kahler.all import (
    ComplexChart,
    HermitianMetric,
    abs2,
    bar,
    d,
    ddbar,
    wedge,
)


def test_components_follow_the_potential_hessian():
    chart = ComplexChart(2, names=("metric_a", "metric_b"))
    a, b = chart.coordinates()
    potential = 2 * abs2(a) + 3 * abs2(b) + a * bar(b) + b * bar(a)
    metric = HermitianMetric.from_potential(potential)

    assert metric.matrix().raw() == matrix(SR, [[2, 1], [1, 3]])
    assert metric[0, 1] == 1
    assert metric.det() == 5
    assert metric.fundamental_form() == I * ddbar(potential)


def test_from_form_removes_the_factor_of_i():
    chart = ComplexChart(2, names=("form_a", "form_b"))
    omega = I * (
        wedge(chart.dz(0), chart.dzbar(0)) + 2 * wedge(chart.dz(1), chart.dzbar(1))
    )
    metric = HermitianMetric.from_form(omega)

    assert metric.matrix().raw() == matrix(SR, [[1, 0], [0, 2]])
    assert metric.fundamental_form() == omega


def test_inverse_contracts_on_the_antiholomorphic_index():
    chart = ComplexChart(2, names=("inverse_a", "inverse_b"))
    a, b = chart.coordinates()
    metric = HermitianMetric(chart, [[1, a], [bar(a), 2]])
    g = metric.matrix().raw()
    inverse = metric.inverse().raw()

    # g^{i lbar} g_{k lbar} = delta^i_k
    product = inverse * g.transpose()
    assert (product - identity_matrix(2)).simplify_full() == 0
    assert metric.is_hermitian()


def test_fubini_study_is_kahler_einstein():
    chart = ComplexChart(2, names=("fs_z", "fs_w"))
    z, w = chart.coordinates()
    metric = HermitianMetric.from_potential(log(1 + abs2(z) + abs2(w)))

    assert metric.is_hermitian()
    assert metric.is_kahler()
    assert metric.ricci_form() == 3 * metric.fundamental_form()
    assert metric.scalar_curvature() == 6
    assert metric.trace(metric.fundamental_form()) == 2
    assert metric.is_positive_definite_at({z: 0.3 + 0.2 * I, w: -1})


def test_flat_metric_has_zero_curvature_and_half_laplacian():
    chart = ComplexChart(2, names=("flat_z", "flat_w"))
    z, w = chart.coordinates()
    metric = HermitianMetric.from_potential(abs2(z) + abs2(w))

    assert metric.ricci_form() == 0
    assert metric.scalar_curvature() == 0
    # Half the Euclidean Laplacian of x^2 + y^2 + u^2 + v^2.
    assert metric.laplacian(abs2(z) + abs2(w)) == 2


def test_non_kahler_hermitian_metric():
    chart = ComplexChart(2, names=("herm_z", "herm_w"))
    z, w = chart.coordinates()
    metric = HermitianMetric(chart, [[1, 0], [0, 1 + abs2(z)]])

    assert metric.is_hermitian()
    assert not metric.is_kahler()
    assert d(metric.fundamental_form()) != 0
    expected = -1 / (1 + abs2(z)) ** 2
    assert (metric.ricci_matrix().raw()[0, 0] - expected).simplify_full() == 0
    assert (metric.scalar_curvature().raw() - expected).simplify_full() == 0
    assert metric.laplacian(abs2(w)) == 1 / (1 + abs2(z))


def test_hermitian_and_positivity_checks_reject_bad_matrices():
    chart = ComplexChart(2, names=("bad_z", "bad_w"))
    z, w = chart.coordinates()

    assert not HermitianMetric(chart, [[1, z], [0, 1]]).is_hermitian()
    negative = HermitianMetric(chart, [[1, 0], [0, -1]])
    assert negative.is_hermitian()
    assert not negative.is_positive_definite_at({z: 0, w: 0})


def test_real_parameters_are_hermitian():
    chart = ComplexChart(1, names=("param_z",))
    (z,) = chart.coordinates()
    real = var("metric_real_parameter", domain="real")

    assert HermitianMetric(chart, [[real]]).is_hermitian()
    assert not HermitianMetric(chart, [[real * I * z]]).is_hermitian()


def test_absolute_values_in_components_are_differentiated_correctly():
    chart = ComplexChart(2, names=("abs_z", "abs_w"))
    z, w = chart.coordinates()

    assert not HermitianMetric(chart, [[1, 0], [0, 1 + abs(z) ** 2]]).is_kahler()


def test_rejects_wrong_matrix_size():
    chart = ComplexChart(2, names=("size_z", "size_w"))
    with pytest.raises(ValueError):
        HermitianMetric(chart, [[1]])
