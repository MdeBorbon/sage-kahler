"""Local holomorphic coordinate charts."""

from sage.all import I, SR, latex

from .conjugation import register_coordinate_pair


class ComplexChart:
    """A minimal complex chart with symbolic holomorphic coordinates."""

    def __init__(self, dimension, names=None, bar_names=None, real_names=None):
        try:
            integer_dimension = int(dimension)
        except (TypeError, ValueError, OverflowError) as exc:
            raise TypeError("dimension must be a positive integer") from exc

        if integer_dimension != dimension or integer_dimension <= 0:
            raise ValueError("dimension must be a positive integer")

        if names is None:
            names = tuple(f"z{i}" for i in range(1, integer_dimension + 1))
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
                    for index in range(1, integer_dimension + 1)
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
        """Rewrite a formal expression in the chart's real coordinates."""
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
        return DifferentialForm(self, {(): SR(expression)})

    def dz(self, index):
        """Return the coordinate ``(1, 0)`` basis form at ``index``."""
        from .forms import basis_form

        self._validate_index(index)
        return basis_form(self, index)

    def dzbar(self, index):
        """Return the coordinate ``(0, 1)`` basis form at ``index``."""
        from .forms import basis_form

        self._validate_index(index)
        return basis_form(self, self._dimension + index)

    def d_dz(self, index):
        """Return the coordinate ``(1, 0)`` vector field ``d/dz^index``."""
        from .vectors import basis_vector

        self._validate_index(index)
        return basis_vector(self, index)

    def d_dzbar(self, index):
        """Return the coordinate ``(0, 1)`` vector field ``d/dzbar^index``."""
        from .vectors import basis_vector

        self._validate_index(index)
        return basis_vector(self, self._dimension + index)

    def partial(self, value):
        """Apply ``partial`` using this chart."""
        from .forms import partial

        return partial(value, chart=self)

    def dbar(self, value):
        """Apply ``dbar`` using this chart."""
        from .forms import dbar

        return dbar(value, chart=self)

    def _validate_index(self, index):
        try:
            integer_index = int(index)
        except (TypeError, ValueError, OverflowError) as exc:
            raise TypeError("coordinate index must be an integer") from exc
        if integer_index != index or not 0 <= integer_index < self._dimension:
            raise IndexError(
                f"coordinate index must be between 0 and {self._dimension - 1}"
            )

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
        return tuple(map(str, self._coordinates + self._conjugate_coordinates))

    def __eq__(self, other):
        """Charts on the same symbolic coordinates are interchangeable.

        Re-running ``ComplexChart(2, names=("z", "w"))`` in a notebook gives
        a chart that accepts forms built on the earlier one. Real coordinate
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
