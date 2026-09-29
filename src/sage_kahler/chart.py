"""Local holomorphic coordinate charts."""

# ****************************************************************************
#       Copyright (C) 2026 Martin de Borbon <martdeborbon@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import I, SR, latex

from .conjugation import register_coordinate_pair


class ComplexChart:
    """A minimal complex chart with symbolic holomorphic coordinates.

    As in SageManifolds, indices run from ``start_index`` (here 1 by
    default): ``dz(1)``, ``frame()[1]``, ``metric[1, 1]``. Where a form or
    vector field is indexed by a combined index ``a``, the coordinates are
    ordered ``z^1, ..., z^n, zbar^1, ..., zbar^n``, so with ``start_index=1``
    the index of ``zbar^j`` is ``n + j``. Sage matrices, such as
    ``metric.matrix()``, keep Sage's own indexing from 0.
    """

    def __init__(
        self, dimension, names=None, bar_names=None, real_names=None, *,
        start_index=1,
    ):
        try:
            integer_dimension = int(dimension)
        except (TypeError, ValueError, OverflowError) as exc:
            raise TypeError("dimension must be a positive integer") from exc

        if integer_dimension != dimension or integer_dimension <= 0:
            raise ValueError("dimension must be a positive integer")
        try:
            integer_start = int(start_index)
        except (TypeError, ValueError, OverflowError) as exc:
            raise TypeError("start_index must be an integer") from exc
        if integer_start != start_index:
            raise TypeError("start_index must be an integer")
        labels = range(integer_start, integer_start + integer_dimension)

        if names is None:
            names = tuple(f"z{i}" for i in labels)
        elif isinstance(names, str):
            names = tuple(name.strip() for name in names.split(","))
        else:
            names = tuple(names)

        if len(names) != integer_dimension:
            raise ValueError(
                f"expected {integer_dimension} coordinate names, got {len(names)}"
            )
        if any(not isinstance(name, str) or not name.strip() for name in names):
            raise ValueError("coordinate names must be non-empty strings")
        if len(set(names)) != len(names):
            raise ValueError("coordinate names must be distinct")

        if bar_names is None:
            bar_names = tuple(f"{name}_bar" for name in names)
        elif isinstance(bar_names, str):
            bar_names = tuple(name.strip() for name in bar_names.split(","))
        else:
            bar_names = tuple(bar_names)
        if len(bar_names) != integer_dimension:
            raise ValueError(
                f"expected {integer_dimension} conjugate coordinate names, "
                f"got {len(bar_names)}"
            )
        if any(not isinstance(name, str) or not name.strip() for name in bar_names):
            raise ValueError("conjugate coordinate names must be non-empty strings")
        if len(set(bar_names)) != len(bar_names):
            raise ValueError("conjugate coordinate names must be distinct")
        if set(names).intersection(bar_names):
            raise ValueError("coordinate and conjugate-coordinate names must be disjoint")

        if real_names is None:
            if integer_dimension == 1:
                real_names = ("x", "y")
            else:
                real_names = tuple(
                    name
                    for index in labels
                    for name in (f"x{index}", f"y{index}")
                )
            if set(real_names).intersection(set(names).union(bar_names)):
                real_names = tuple(
                    name
                    for coordinate_name in names
                    for name in (
                        f"re_{coordinate_name}",
                        f"im_{coordinate_name}",
                    )
                )
        else:
            if isinstance(real_names, str):
                real_names = tuple(name.strip() for name in real_names.split(","))
            else:
                real_names = tuple(real_names)
        if len(real_names) != 2 * integer_dimension:
            raise ValueError(
                f"expected {2 * integer_dimension} real coordinate names, "
                f"got {len(real_names)}"
            )
        if any(not isinstance(name, str) or not name.strip() for name in real_names):
            raise ValueError("real coordinate names must be non-empty strings")
        if len(set(real_names)) != len(real_names):
            raise ValueError("real coordinate names must be distinct")
        if set(real_names).intersection(set(names).union(bar_names)):
            raise ValueError(
                "real, holomorphic, and conjugate-coordinate names must be disjoint"
            )

        self._dimension = integer_dimension
        self._start_index = integer_start
        self._coordinate_names = names
        self._coordinates = tuple(SR.var(name) for name in names)
        self._conjugate_coordinates = tuple(
            SR.var(name, latex_name=rf"\bar{{{latex(coordinate)}}}")
            for name, coordinate in zip(bar_names, self._coordinates)
        )
        real_variables = tuple(SR.var(name) for name in real_names)
        self._real_coordinates = tuple(
            (real_variables[2 * index], real_variables[2 * index + 1])
            for index in range(integer_dimension)
        )
        for coordinate, conjugate_coordinate in zip(
            self._coordinates, self._conjugate_coordinates
        ):
            register_coordinate_pair(self, coordinate, conjugate_coordinate)

    def dimension(self):
        """Return the complex dimension of the chart."""
        return self._dimension

    def start_index(self):
        """Return the first index, 1 by default, as in SageManifolds."""
        return self._start_index

    def irange(self):
        """Return the coordinate indices, as ``Manifold.irange`` does."""
        return range(self._start_index, self._start_index + self._dimension)

    def coordinates(self):
        """Return the chart's symbolic holomorphic coordinates."""
        return self._coordinates

    def conjugate_coordinates(self):
        """Return the formal conjugates, independent for differentiation."""
        return self._conjugate_coordinates

    def real_coordinates(self):
        """Return ``(x_i, y_i)`` pairs realizing ``z_i = x_i + I*y_i``."""
        return self._real_coordinates

    def to_real(self, expression):
        """Rewrite a formal expression in the chart's real coordinates.

        As for the Dolbeault operators, coordinate-dependent ``abs(f)`` is
        first rewritten as ``sqrt(f * bar(f))``.
        """
        from .forms import _formalize_absolute_values

        expression = _formalize_absolute_values(SR(expression), self)
        substitutions = {}
        for coordinate, conjugate_coordinate, (x, y) in zip(
            self._coordinates,
            self._conjugate_coordinates,
            self._real_coordinates,
        ):
            substitutions[coordinate] = x + I * y
            substitutions[conjugate_coordinate] = x - I * y
        return SR(expression).subs(substitutions).expand().simplify_full()

    def display(self, expression):
        """Display a scalar or matrix with chart-aware norms and conjugates.

        Scalars are wrapped as zero-forms; their original coefficient is in
        ``terms()[()]``. Matrix displays provide ``matrix()`` to retrieve a
        copy of the original Sage matrix. Neither changes the input.
        """
        from sage.structure.element import Matrix
        from .forms import DifferentialForm, _MatrixDisplay

        if isinstance(expression, Matrix):
            return _MatrixDisplay(self, expression)
        return DifferentialForm._internal(self, {(): SR(expression)})

    def dz(self, index):
        """Return the coordinate ``(1, 0)`` basis form ``dz^index``."""
        from .forms import basis_form

        return basis_form(self, self._position(index))

    def dzbar(self, index):
        """Return the coordinate ``(0, 1)`` basis form ``dzbar^index``."""
        from .forms import basis_form

        return basis_form(self, self._dimension + self._position(index))

    def d_dz(self, index):
        """Return the coordinate ``(1, 0)`` vector field ``d/dz^index``."""
        from .vectors import basis_vector

        return basis_vector(self, self._position(index))

    def d_dzbar(self, index):
        """Return the coordinate ``(0, 1)`` vector field ``d/dzbar^index``."""
        from .vectors import basis_vector

        return basis_vector(self, self._dimension + self._position(index))

    def frame(self):
        """Return the coordinate frame ``(d/dz^1, ..., d/dz^n)`` of ``T^{1,0}``.

        As ``Chart.frame()`` in SageManifolds, it is indexed from
        ``start_index`` and ``frame()[:]`` is the tuple of vector fields.
        """
        return CoordinateFrame(
            self, [self.d_dz(i) for i in self.irange()], "Coordinate frame"
        )

    def coframe(self):
        """Return the coordinate coframe ``(dz^1, ..., dz^n)``, dual to ``frame()``."""
        return CoordinateFrame(
            self, [self.dz(i) for i in self.irange()], "Coordinate coframe"
        )

    def conjugate_frame(self):
        """Return ``(d/dzbar^1, ..., d/dzbar^n)``, the frame of ``T^{0,1}``."""
        return CoordinateFrame(
            self,
            [self.d_dzbar(i) for i in self.irange()],
            "Conjugate coordinate frame",
        )

    def conjugate_coframe(self):
        """Return ``(dzbar^1, ..., dzbar^n)``, dual to ``conjugate_frame()``."""
        return CoordinateFrame(
            self,
            [self.dzbar(i) for i in self.irange()],
            "Conjugate coordinate coframe",
        )

    def partial(self, value):
        """Apply ``partial`` using this chart."""
        from .forms import partial

        return partial(value, chart=self)

    def dbar(self, value):
        """Apply ``dbar`` using this chart."""
        from .forms import dbar

        return dbar(value, chart=self)

    def _position(self, index):
        """Return the 0-based position of a coordinate index."""
        return _position(index, self._start_index, self._dimension)

    def _combined_position(self, index):
        """Return the 0-based position of a combined index, ``zbar`` after ``z``."""
        return _position(index, self._start_index, 2 * self._dimension)

    def _basis_label(self, combined_index):
        if combined_index < self._dimension:
            return f"d{self._coordinate_names[combined_index]}"
        return f"dbar({self._coordinate_names[combined_index - self._dimension]})"

    def _basis_latex(self, combined_index):
        coordinate = self._coordinates[combined_index % self._dimension]
        label = str(latex(coordinate))
        if combined_index >= self._dimension:
            label = rf"\bar{{{label}}}"
        return f"d{label}"

    def _vector_label(self, combined_index):
        if combined_index < self._dimension:
            return f"∂/∂{self._coordinate_names[combined_index]}"
        return f"∂/∂bar({self._coordinate_names[combined_index - self._dimension]})"

    def _vector_latex(self, combined_index):
        coordinate = self._coordinates[combined_index % self._dimension]
        label = str(latex(coordinate))
        if combined_index >= self._dimension:
            label = rf"\bar{{{label}}}"
        return rf"\frac{{\partial}}{{\partial {label}}}"

    def _key(self):
        return (
            tuple(map(str, self._coordinates + self._conjugate_coordinates)),
            self._start_index,
        )

    def __eq__(self, other):
        """Charts on the same symbolic coordinates are interchangeable.

        Re-running ``ComplexChart(2, names=("z", "w"))`` in a notebook gives
        a chart that accepts forms built on the earlier one. The start index
        is compared, since it changes the meaning of indices; real coordinate
        names only affect ``to_real`` and are not compared.
        """
        if not isinstance(other, ComplexChart):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self):
        return hash(self._key())

    def __repr__(self):
        coordinates = ", ".join(map(str, self._coordinates))
        return (
            f"Complex chart of dimension {self._dimension} "
            f"with coordinates ({coordinates})"
        )


def _position(index, start, length):
    """Return ``index - start``, checking that it lies in ``range(length)``."""
    try:
        integer_index = int(index)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError("index must be an integer") from exc
    if integer_index != index:
        raise TypeError("index must be an integer")
    if not start <= integer_index < start + length:
        raise IndexError(
            f"index out of range: {integer_index} not in [{start}, {start + length - 1}]"
        )
    return integer_index - start


class CoordinateFrame:
    """A coordinate frame or coframe, indexed from the chart's start index."""

    def __init__(self, chart, elements, description):
        self._chart = chart
        self._elements = tuple(elements)
        self._description = description

    @property
    def chart(self):
        return self._chart

    def __getitem__(self, key):
        """Return one element, or a tuple for a slice such as ``[:]``."""
        if isinstance(key, slice):
            start = self._chart.start_index()
            shift = lambda bound: None if bound is None else bound - start
            return self._elements[slice(shift(key.start), shift(key.stop), key.step)]
        return self._elements[self._chart._position(key)]

    def __iter__(self):
        return iter(self._elements)

    def __len__(self):
        return len(self._elements)

    def __repr__(self):
        return f"{self._description} ({', '.join(map(repr, self._elements))})"

    def _latex_(self):
        from sage.all import latex

        return (
            r"\left(" + ", ".join(str(latex(e)) for e in self._elements) + r"\right)"
        )
