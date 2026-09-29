import pytest

from sage.all import log

from sage_kahler.all import ComplexChart, HermitianMetric, abs2, bar, wedge


@pytest.fixture(scope="module")
def fubini_study():
    chart = ComplexChart(2, names=("conn_fs_z", "conn_fs_w"))
    z, w = chart.coordinates()
    metric = HermitianMetric.from_potential(log(1 + abs2(z) + abs2(w)))
    return chart, metric, metric.chern_connection()


@pytest.fixture(scope="module")
def non_kahler():
    chart = ComplexChart(2, names=("conn_h_z", "conn_h_w"))
    z, _ = chart.coordinates()
    metric = HermitianMetric(chart, [[1, 0], [0, 1 + abs2(z)]])
    return chart, metric, metric.chern_connection()


def test_fubini_study_christoffel_symbols(fubini_study):
    chart, _, connection = fubini_study
    coordinates = chart.coordinates()
    radius = 1 + sum(abs2(z) for z in coordinates)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                expected = -(
                    (i == j) * bar(coordinates[k]) + (i == k) * bar(coordinates[j])
                ) / radius
                gamma = connection.christoffel(i + 1, j + 1, k + 1)
                assert (gamma - expected).raw().simplify_full() == 0
    assert connection.is_torsion_free()
    symbols = connection.christoffel_symbols()
    # Gamma^1_{22} vanishes; the keys use the chart's indices, from 1.
    assert (1, 2, 2) not in symbols and (2, 2, 2) in symbols
    assert symbols[1, 1, 2] == connection.christoffel(1, 1, 2)
    assert (1, 1, 1, 1) in connection.curvature_tensor()
    with pytest.raises(IndexError):
        connection.christoffel(0, 1, 1)


def test_fubini_study_has_constant_holomorphic_sectional_curvature(fubini_study):
    _, metric, connection = fubini_study
    g = metric.matrix().raw()
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    expected = g[i, j] * g[k, l] + g[i, l] * g[k, j]
                    curvature = connection.curvature(i + 1, j + 1, k + 1, l + 1)
                    difference = curvature.raw() - expected
                    assert difference.simplify_full() == 0


def test_kahler_ricci_contractions_agree(fubini_study):
    _, metric, connection = fubini_study
    ricci = metric.ricci_matrix().raw()
    assert (connection.first_chern_ricci().raw() - ricci).simplify_full() == 0
    assert (connection.second_chern_ricci().raw() - ricci).simplify_full() == 0


def test_lowered_curvature_matches_endomorphism(non_kahler):
    _, metric, connection = non_kahler
    g = metric.matrix().raw()
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    lowered = sum(
                        g[p, j] * connection.curvature_endomorphism(
                            i + 1, p + 1, k + 1, l + 1
                        ).raw()
                        for p in range(2)
                    )
                    curvature = connection.curvature(i + 1, j + 1, k + 1, l + 1)
                    difference = curvature.raw() - lowered
                    assert difference.simplify_full() == 0


def test_curvature_is_the_commutator_of_covariant_derivatives(non_kahler):
    chart, _, connection = non_kahler
    for i in chart.irange():
        for k in chart.irange():
            for l in chart.irange():
                expected = sum(
                    connection.curvature_endomorphism(i, j, k, l) * chart.d_dz(j)
                    for j in chart.irange()
                )
                result = connection.curvature_operator(
                    chart.d_dz(k), chart.d_dzbar(l), chart.d_dz(i)
                )
                assert result == expected


def test_non_kahler_metric_has_torsion_and_distinct_riccis(non_kahler):
    chart, metric, connection = non_kahler
    z, _ = chart.coordinates()

    assert not connection.is_torsion_free()
    assert (connection.torsion(2, 1, 2) - bar(z) / (1 + abs2(z))).raw().simplify_full() == 0
    first = connection.first_chern_ricci().raw()
    second = connection.second_chern_ricci().raw()
    assert (first - metric.ricci_matrix().raw()).simplify_full() == 0
    assert (second[1, 1] + 1 / (1 + abs2(z))).simplify_full() == 0
    assert (first - second).simplify_full() != 0


def test_connection_is_metric_compatible_and_of_type_1_0(non_kahler):
    chart, metric, connection = non_kahler
    z, w = chart.coordinates()
    nabla = connection.covariant_derivative
    first = z * chart.d_dz(2) + chart.d_dzbar(1)
    second = w * chart.d_dz(1) + bar(z) * chart.d_dz(2)
    third = bar(w) * chart.d_dzbar(2) + chart.d_dzbar(1)

    for direction in (chart.d_dz(1), chart.d_dzbar(1), first):
        assert direction(metric(second, third)) == (
            metric(nabla(direction, second), third)
            + metric(second, nabla(direction, third))
        )
    for k in chart.irange():
        for l in chart.irange():
            assert nabla(chart.d_dzbar(l), chart.d_dz(k)) == 0
            assert nabla(chart.d_dz(l), chart.d_dzbar(k)) == 0


def test_torsion_tensor_formula(non_kahler):
    chart, _, connection = non_kahler
    nabla = connection.covariant_derivative
    first = chart.d_dz(1)
    second = chart.d_dz(2)
    torsion = nabla(first, second) - nabla(second, first) - first.bracket(second)
    expected = sum(connection.torsion(i, 1, 2) * chart.d_dz(i) for i in chart.irange())
    assert torsion == expected


def test_covariant_derivative_of_forms_is_dual(non_kahler):
    chart, _, connection = non_kahler
    z, w = chart.coordinates()
    nabla = connection.covariant_derivative
    direction = chart.d_dz(1) + z * chart.d_dzbar(2)
    vector = w * chart.d_dz(2) + chart.d_dzbar(2)
    form = z * chart.dz(2) + bar(w) * chart.dzbar(2)

    # X(alpha(Y)) = (nabla_X alpha)(Y) + alpha(nabla_X Y)
    pairing = vector.contract(form).terms().get((), 0)
    left = direction(pairing)
    right = (
        vector.contract(nabla(direction, form)).terms().get((), 0)
        + nabla(direction, vector).contract(form).terms().get((), 0)
    )
    assert (left - right).raw().simplify_full() == 0
    assert nabla(chart.d_dz(1), chart.dz(2)) == -connection.christoffel(2, 1, 2) * chart.dz(2)

    # Leibniz rule on a two-form.
    two_form = wedge(chart.dz(1), chart.dz(2))
    assert nabla(direction, two_form) == (
        wedge(nabla(direction, chart.dz(1)), chart.dz(2))
        + wedge(chart.dz(1), nabla(direction, chart.dz(2)))
    )
