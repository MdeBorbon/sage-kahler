import pytest
from sage.all import SR

from sage_kahler.all import ComplexChart


def test_chart_with_explicit_names():
    chart = ComplexChart(2, names=("z", "w"))

    assert chart.dimension() == 2
    assert tuple(map(str, chart.coordinates())) == ("z", "w")
    assert repr(chart) == "Complex chart of dimension 2 with coordinates (z, w)"


def test_chart_with_default_names():
    chart = ComplexChart(2)

    assert tuple(map(str, chart.coordinates())) == ("z1", "z2")


def test_coordinates_are_sage_symbolic_expressions():
    z, w = ComplexChart(2, names=("z", "w")).coordinates()

    assert str(z + w) in {"w + z", "z + w"}
    assert z.parent() is SR
    assert w.parent() is SR


@pytest.mark.parametrize("dimension", [0, -1, 1.5])
def test_dimension_must_be_a_positive_integer(dimension):
    with pytest.raises(ValueError, match="positive integer"):
        ComplexChart(dimension)


def test_number_of_names_must_match_dimension():
    with pytest.raises(ValueError, match="expected 2 coordinate names"):
        ComplexChart(2, names=("z",))


def test_default_real_names_do_not_collide_with_complex_names():
    chart = ComplexChart(1, names=("x",))

    assert tuple(map(str, chart.real_coordinates()[0])) == ("re_x", "im_x")
