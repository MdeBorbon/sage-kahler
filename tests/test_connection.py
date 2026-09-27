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
                assert (connection.christoffel(i, j, k) - expected).raw().simplify_full() == 0
    assert connection.is_torsion_free()


def test_fubini_study_has_constant_holomorphic_sectional_curvature(fubini_study):
    _, metric, connection = fubini_study
    g = metric.matrix().raw()
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    expected = g[i, j] * g[k, l] + g[i, l] * g[k, j]
                    difference = connection.curvature(i, j, k, l).raw() - expected
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
                        g[p, j] * connection.curvature_endomorphism(i, p, k, l).raw()
                        for p in range(2)
                    )
                    difference = connection.curvature(i, j, k, l).raw() - lowered
                    assert difference.simplify_full() == 0


def test_curvature_is_the_commutator_of_covariant_derivatives(non_kahler):
    chart, _, connection = non_kahler
    for i in range(2):
        for k in range(2):
            for l in range(2):
                expected = sum(
                    connection.curvature_endomorphism(i, j, k, l) * chart.d_dz(j)
                    for j in range(2)
                )
                result = connection.curvature_operator(
                    chart.d_dz(k), chart.d_dzbar(l), chart.d_dz(i)
                )
                assert result == expected


def test_non_kahler_metric_has_torsion_and_distinct_riccis(non_kahler):
    chart, metric, connection = non_kahler
    z, _ = chart.coordinates()

    assert not connection.is_torsion_free()
    assert (connection.torsion(1, 0, 1) - bar(z) / (1 + abs2(z))).raw().simplify_full() == 0
    first = connection.first_chern_ricci().raw()
    second = connection.second_chern_ricci().raw()
    assert (first - metric.ricci_matrix().raw()).simplify_full() == 0
    assert (second[1, 1] + 1 / (1 + abs2(z))).simplify_full() == 0
    assert (first - second).simplify_full() != 0


def test_connection_is_metric_compatible_and_of_type_1_0(non_kahler):
    chart, metric, connection = non_kahler
    z, w = chart.coordinates()
    nabla = connection.covariant_derivative
    first = z * chart.d_dz(1) + chart.d_dzbar(0)
    second = w * chart.d_dz(0) + bar(z) * chart.d_dz(1)
    third = bar(w) * chart.d_dzbar(1) + chart.d_dzbar(0)

    for direction in (chart.d_dz(0), chart.d_dzbar(0), first):
        assert direction(metric(second, third)) == (
            metric(nabla(direction, second), third)
            + metric(second, nabla(direction, third))
        )
    for k in range(2):
        for l in range(2):
            assert nabla(chart.d_dzbar(l), chart.d_dz(k)) == 0
            assert nabla(chart.d_dz(l), chart.d_dzbar(k)) == 0


def test_torsion_tensor_formula(non_kahler):
    chart, _, connection = non_kahler
    nabla = connection.covariant_derivative
    first = chart.d_dz(0)
    second = chart.d_dz(1)
    torsion = nabla(first, second) - nabla(second, first) - first.bracket(second)
    expected = sum(connection.torsion(i, 0, 1) * chart.d_dz(i) for i in range(2))
    assert torsion == expected


def test_covariant_derivative_of_forms_is_dual(non_kahler):
    chart, _, connection = non_kahler
    z, w = chart.coordinates()
    nabla = connection.covariant_derivative
    direction = chart.d_dz(0) + z * chart.d_dzbar(1)
    vector = w * chart.d_dz(1) + chart.d_dzbar(1)
    form = z * chart.dz(1) + bar(w) * chart.dzbar(1)

    # X(alpha(Y)) = (nabla_X alpha)(Y) + alpha(nabla_X Y)
    pairing = vector.contract(form).terms().get((), 0)
    left = direction(pairing)
    right = (
        vector.contract(nabla(direction, form)).terms().get((), 0)
        + nabla(direction, vector).contract(form).terms().get((), 0)
    )
    assert (left - right).raw().simplify_full() == 0
    assert nabla(chart.d_dz(0), chart.dz(1)) == -connection.christoffel(1, 0, 1) * chart.dz(1)

    # Leibniz rule on a two-form.
    two_form = wedge(chart.dz(0), chart.dz(1))
    assert nabla(direction, two_form) == (
        wedge(nabla(direction, chart.dz(0)), chart.dz(1))
        + wedge(chart.dz(0), nabla(direction, chart.dz(1)))
    )
