import pytest

from sage.all import I, SR, latex, log

from sage_kahler.all import (
    ComplexChart,
    HermitianMetric,
    VectorField,
    abs2,
    bar,
    wedge,
)


def test_coordinate_vectors_are_dual_to_coordinate_forms():
    chart = ComplexChart(2, names=("dual_z", "dual_w"))
    for a, vector in enumerate(
        [chart.d_dz(0), chart.d_dz(1), chart.d_dzbar(0), chart.d_dzbar(1)]
    ):
        for b, form in enumerate(
            [chart.dz(0), chart.dz(1), chart.dzbar(0), chart.dzbar(1)]
        ):
            assert vector.contract(form) == (1 if a == b else 0)


def test_types_parts_and_conjugation():
    chart = ComplexChart(2, names=("type_z", "type_w"))
    z, w = chart.coordinates()
    vector = 2 * z * chart.d_dz(0) + bar(w) * chart.d_dzbar(1)

    assert chart.d_dz(0).vector_type() == (1, 0)
    assert chart.d_dzbar(1).vector_type() == (0, 1)
    assert vector.vector_type() is None
    assert vector.holomorphic_part() == 2 * z * chart.d_dz(0)
    assert vector.antiholomorphic_part() == bar(w) * chart.d_dzbar(1)
    assert bar(vector) == 2 * bar(z) * chart.d_dzbar(0) + w * chart.d_dz(1)
    assert not vector.is_real()
    assert (vector + bar(vector)).is_real()
    assert vector.holomorphic_components() == (2 * z, 0)


def test_vector_fields_differentiate_functions():
    chart = ComplexChart(1, names=("derivative_z",))
    (z,) = chart.coordinates()

    assert chart.d_dz(0)(abs2(z) ** 2) == 2 * z * bar(z) ** 2
    assert chart.d_dzbar(0)(z**3) == 0
    # abs() is rewritten through the formal conjugate, as for Dolbeault operators.
    assert (chart.d_dzbar(0)(abs(z) ** 2) - z).simplify_full() == 0


def test_real_coordinate_vectors():
    chart = ComplexChart(1, names=("real_z",))
    d_dx = chart.d_dz(0) + chart.d_dzbar(0)
    d_dy = I * (chart.d_dz(0) - chart.d_dzbar(0))

    assert d_dx.is_real()
    assert d_dy.is_real()


def test_interior_product_signs():
    chart = ComplexChart(2, names=("contract_z", "contract_w"))
    z, w = chart.coordinates()
    form = wedge(chart.dz(0), chart.dzbar(1))
    vector = 2 * z * chart.d_dz(0) + bar(w) * chart.d_dzbar(1)

    assert vector.contract(form) == 2 * z * chart.dzbar(1) - bar(w) * chart.dz(0)
    assert vector.contract(vector.contract(form)) == 0


def test_lie_bracket_and_holomorphy():
    chart = ComplexChart(2, names=("bracket_z", "bracket_w"))
    z, w = chart.coordinates()
    first = z * chart.d_dz(0)
    second = z**2 * chart.d_dz(1)

    assert first.bracket(second) == 2 * z**2 * chart.d_dz(1)
    assert first.bracket(second) == -second.bracket(first)
    assert chart.d_dz(0).bracket(chart.d_dzbar(0)) == 0
    assert first.is_holomorphic()
    assert not (abs2(z) * chart.d_dz(0)).is_holomorphic()
    assert not chart.d_dzbar(0).is_holomorphic()
    assert VectorField(chart).is_holomorphic()


def test_metric_evaluates_on_complex_vectors():
    chart = ComplexChart(2, names=("eval_z", "eval_w"))
    z, w = chart.coordinates()
    metric = HermitianMetric.from_potential(log(1 + abs2(z) + abs2(w)))

    for j in range(2):
        for k in range(2):
            assert metric(chart.d_dz(j), chart.d_dzbar(k)) == metric[j, k]
            assert metric(chart.d_dzbar(k), chart.d_dz(j)) == metric[j, k]
            assert metric(chart.d_dz(j), chart.d_dz(k)) == 0
            assert metric(chart.d_dzbar(j), chart.d_dzbar(k)) == 0

    real = chart.d_dz(0) + chart.d_dzbar(0)
    assert metric(real, real) == 2 * metric[0, 0]


def test_grad10_of_holomorphy_potentials():
    chart = ComplexChart(2, names=("grad_z", "grad_w"))
    z, w = chart.coordinates()
    flat = HermitianMetric.from_potential(abs2(z) + abs2(w))
    fubini_study = HermitianMetric.from_potential(log(1 + abs2(z) + abs2(w)))

    assert flat.grad10(abs2(z)) == z * chart.d_dz(0)
    gradient = fubini_study.grad10(abs2(z) / (1 + abs2(z) + abs2(w)))
    assert gradient == z * chart.d_dz(0)
    assert gradient.is_holomorphic()
    assert not flat.grad10(abs2(z) ** 2).is_holomorphic()


def test_display():
    chart = ComplexChart(1, names=("u",))
    (u,) = chart.coordinates()
    vector = chart.d_dz(0) - bar(u) * chart.d_dzbar(0)

    assert repr(vector) == "∂/∂u + (-bar(u)) ∂/∂bar(u)"
    assert latex(vector) == (
        r"\frac{\partial}{\partial u} - \bar{u}\,\frac{\partial}{\partial \bar{u}}"
    )
    assert repr(VectorField(chart)) == "0"


def test_rejects_mismatched_input():
    chart = ComplexChart(1, names=("bad_vector_z",))
    other = ComplexChart(1, names=("other_vector_z",))
    with pytest.raises(ValueError):
        VectorField(chart, [1])
    with pytest.raises(ValueError):
        chart.d_dz(0) + other.d_dz(0)
    with pytest.raises(IndexError):
        chart.d_dz(1)
