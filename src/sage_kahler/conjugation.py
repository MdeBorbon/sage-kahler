"""Formal complex conjugation for local holomorphic coordinates."""

# ****************************************************************************
#       Copyright (C) 2026 Martin de Borbon <martdeborbon@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import SR


_CONJUGATES = {}
_VARIABLE_CHARTS = {}


def register_coordinate_pair(chart, coordinate, conjugate_coordinate):
    """Register a pair of formal conjugate variables belonging to ``chart``."""
    for left, right in (
        (coordinate, conjugate_coordinate),
        (conjugate_coordinate, coordinate),
    ):
        existing = _CONJUGATES.get(left)
        if existing is not None and existing != right:
            raise ValueError(f"the symbolic variable {left} is already registered")
        _CONJUGATES[left] = right
        _VARIABLE_CHARTS[left] = chart


def conjugate_variable(variable):
    """Return the registered formal conjugate of ``variable``, or None."""
    return _CONJUGATES.get(variable)


def chart_for_expression(expression):
    """Return the unique registered chart containing an expression's variables."""
    expression = SR(expression)
    charts = {
        _VARIABLE_CHARTS[variable]
        for variable in expression.variables()
        if variable in _VARIABLE_CHARTS
    }
    if not charts:
        raise ValueError(
            "cannot infer a chart from this expression; pass chart= explicitly"
        )
    if len(charts) != 1:
        raise ValueError("the expression contains coordinates from multiple charts")
    return charts.pop()


def bar(value):
    """Return the formal complex conjugate of an expression or differential form."""
    conjugate_method = getattr(value, "_formal_conjugate_", None)
    if conjugate_method is not None:
        return conjugate_method()

    expression = SR(value)
    registered = _CONJUGATES.get(expression)
    if registered is not None:
        return registered

    operands = expression.operands()
    if not operands:
        return expression.conjugate()

    operator = expression.operator()
    return operator(*(bar(operand) for operand in operands))


def abs2(expression):
    """Return ``expression * bar(expression)`` using formal conjugation."""
    expression = SR(expression)
    return expression * bar(expression)
