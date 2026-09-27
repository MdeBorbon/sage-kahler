"""Native Sage results with chart-aware presentation.

Computations use independent formal coordinates. Only the display is changed;
``raw()`` provides the native base-class object for interoperability.
"""

# ****************************************************************************
#       Copyright (C) 2026 Martin de Borbon <martdeborbon@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import SR, matrix
from sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
from sage.symbolic.expression import Expression


class ChartExpression(Expression):
    """A Sage symbolic expression carrying its chart for display."""

    def __init__(self, parent, value, chart=None):
        super().__init__(parent, value)
        self._chart = chart

    def raw(self):
        """Return an ordinary Sage expression with the same symbolic tree."""
        return Expression(SR, self)

    def _repr_(self):
        from .forms import _format_coefficient

        return _format_coefficient(self.raw(), self._chart)

    def _latex_(self):
        from .forms import _format_coefficient

        return _format_coefficient(self.raw(), self._chart, latex_mode=True)


    def _wrap(self, result):
        if isinstance(result, Expression):
            return ChartExpression(SR, result, self._chart)
        return result

    @staticmethod
    def _unwrap(value):
        return value.raw() if isinstance(value, ChartExpression) else value

    def __add__(self, other):
        return self._wrap(self.raw() + self._unwrap(other))

    def __radd__(self, other):
        return self._wrap(self._unwrap(other) + self.raw())

    def __sub__(self, other):
        return self._wrap(self.raw() - self._unwrap(other))

    def __rsub__(self, other):
        return self._wrap(self._unwrap(other) - self.raw())

    def __mul__(self, other):
        return self._wrap(self.raw() * self._unwrap(other))

    def __rmul__(self, other):
        return self._wrap(self._unwrap(other) * self.raw())

    def __truediv__(self, other):
        return self._wrap(self.raw() / self._unwrap(other))

    def __rtruediv__(self, other):
        return self._wrap(self._unwrap(other) / self.raw())

    def __pow__(self, other, modulo=None):
        return self._wrap(pow(self.raw(), self._unwrap(other), modulo))

    def __neg__(self):
        return self._wrap(-self.raw())

    def derivative(self, *args, **kwds):
        return self._wrap(self.raw().derivative(*args, **kwds))

    diff = derivative
    differentiate = derivative

    def subs(self, *args, **kwds):
        return self._wrap(self.raw().subs(*args, **kwds))

    substitute = subs

    def simplify_full(self, *args, **kwds):
        return self._wrap(self.raw().simplify_full(*args, **kwds))

    def factor(self, *args, **kwds):
        return self._wrap(self.raw().factor(*args, **kwds))

    def expand(self, *args, **kwds):
        return self._wrap(self.raw().expand(*args, **kwds))


class ChartMatrix(Matrix_symbolic_dense):
    """A native symbolic matrix with chart-aware entries and determinant."""

    _chart = None

    @classmethod
    def from_matrix(cls, value, chart):
        result = cls(value.parent(), value.list())
        result._chart = chart
        return result

    def raw(self):
        """Return an independent ordinary Sage matrix."""
        return matrix(SR, self.nrows(), self.ncols(), self.list())

    def _display_chart(self):
        # Sage matrix operations may allocate this subclass without copying
        # Python attributes. Infer the chart from the resulting entries then.
        if self._chart is not None:
            return self._chart
        from .conjugation import chart_for_expression

        variables = {
            variable for entry in self.list() for variable in entry.variables()
        }
        try:
            return chart_for_expression(sum(variables, SR.zero()))
        except ValueError:
            return None

    def _repr_(self):
        from .forms import _MatrixDisplay

        chart = self._display_chart()
        if chart is None:
            return repr(self.raw())
        return repr(_MatrixDisplay(chart, self.raw()))

    def __repr__(self):
        return self._repr_()

    def __str__(self):
        return self._repr_()

    def _latex_(self):
        from sage.all import latex
        from .forms import _MatrixDisplay

        chart = self._display_chart()
        if chart is None:
            return str(latex(self.raw()))
        return _MatrixDisplay(chart, self.raw())._latex_()

    def determinant(self, *args, **kwds):
        from .forms import _simplify_scalar

        result = self.raw().determinant(*args, **kwds)
        chart = self._display_chart()
        if chart is None:
            return result.simplify_full()
        return ChartExpression(SR, _simplify_scalar(result, chart), chart)

    det = determinant

    def transpose(self):
        return ChartMatrix.from_matrix(self.raw().transpose(), self._display_chart())

    def inverse(self, *args, **kwds):
        return ChartMatrix.from_matrix(
            self.raw().inverse(*args, **kwds), self._display_chart()
        )

    def subs(self, *args, **kwds):
        return ChartMatrix.from_matrix(
            self.raw().subs(*args, **kwds), self._display_chart()
        )

    substitute = subs

    def __getitem__(self, key):
        result = super().__getitem__(key)
        chart = self._display_chart()
        if chart is None:
            return result
        if isinstance(result, Matrix_symbolic_dense):
            return ChartMatrix.from_matrix(result, chart)
        if isinstance(result, Expression):
            return ChartExpression(SR, result, chart)
        return result
