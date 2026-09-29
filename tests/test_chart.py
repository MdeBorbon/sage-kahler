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


def test_recreated_chart_accepts_forms_from_the_earlier_chart():
    from sage_kahler.all import abs2, ddbar

    first = ComplexChart(2, names=("rerun_z", "rerun_w"))
    form = ddbar(abs2(first.coordinates()[0]))
    # Re-running a notebook cell creates a new chart on the same symbols.
    second = ComplexChart(2, names=("rerun_z", "rerun_w"))
    z, w = second.coordinates()

    assert second == first
    assert hash(second) == hash(first)
    assert form + ddbar(abs2(w)) == (
        second.dz(1).wedge(second.dzbar(1)) + second.dz(2).wedge(second.dzbar(2))
    )
    # Coordinate order defines the basis, so a reordered chart differs.
    assert second != ComplexChart(2, names=("rerun_w", "rerun_z"))


def test_indices_start_at_one_by_default():
    chart = ComplexChart(2, names=("one_z", "one_w"))
    z, w = chart.coordinates()

    assert chart.start_index() == 1
    assert list(chart.irange()) == [1, 2]
    assert chart.d_dz(1)(z) == 1
    assert chart.dzbar(2)(chart.d_dzbar(2)) == 1
    for method in (chart.dz, chart.dzbar, chart.d_dz, chart.d_dzbar):
        with pytest.raises(IndexError):
            method(0)
        with pytest.raises(IndexError):
            method(3)
    assert tuple(map(str, ComplexChart(2).coordinates())) == ("z1", "z2")


def test_frames_and_coframes_follow_sage_manifolds():
    chart = ComplexChart(2, names=("frame_z", "frame_w"))
    frame, coframe = chart.frame(), chart.coframe()

    assert frame[1] == chart.d_dz(1)
    assert coframe[2] == chart.dz(2)
    assert chart.conjugate_frame()[:] == (chart.d_dzbar(1), chart.d_dzbar(2))
    assert chart.conjugate_coframe()[:] == (chart.dzbar(1), chart.dzbar(2))
    assert len(frame) == 2 and list(frame) == list(frame[:])
    assert frame[2:] == (chart.d_dz(2),)
    for i in chart.irange():
        for j in chart.irange():
            assert coframe[i](frame[j]) == (1 if i == j else 0)
    with pytest.raises(IndexError):
        frame[0]
    assert repr(frame) == "Coordinate frame (∂/∂frame_z, ∂/∂frame_w)"
    assert repr(coframe) == "Coordinate coframe (dframe_z, dframe_w)"


def test_start_index_zero_matches_the_sage_manifolds_default():
    chart = ComplexChart(2, names=("zero_z", "zero_w"), start_index=0)

    assert chart.frame()[0] == chart.d_dz(0)
    assert chart.dz(1).terms() == {(1,): 1}
    assert chart.dzbar(0).terms() == {(2,): 1}
    with pytest.raises(IndexError):
        chart.dz(2)
    assert chart != ComplexChart(2, names=("zero_z", "zero_w"))
    assert tuple(map(str, ComplexChart(2, start_index=0).coordinates())) == (
        "z0",
        "z1",
    )
