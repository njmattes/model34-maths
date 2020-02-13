#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from math import sqrt, atan
import sympy as sp
from sympy.geometry import Ellipse


class NeckEllipse(object):
    def __init__(self, width, depth, y_int, b=None):
        """Ellipse, either inner or outer, in section of neck.

        #TODO: Double check that this width arg is correct.
        #TODO: Change arg name to x_int.
        :param width: Width of ellipse---if outer, at underside of fretboard,
        if inner, at TKTKTK
        :type width: float
        :param depth: Depth from underside of fretboard to bottom apex
        of ellipse
        :type depth: float
        :param y_int: y-intercept where width crosses x-axis
        :type y_int: float
        #TODO: Is the minor axis always horizontal?
        #TODO: Look at nut section in CAD files.
        :param b: length of minor (horizontal) axis if known
        :type b: float
        """
        self.width = width
        self.depth = depth
        self.y_int = y_int
        self._a = None
        self._b = b
        self._k = None
        self._y = None

    @property
    def ellipse(self):
        """

        :return: Inner or outer elliptical profile.
        :rtype: sp.geometry.Ellipse
        """
        return self.get_ellipse()

    @property
    def a(self):
        """Major (vertical) axis.

        :return: Length of major axis, a, from center of ellipse.
        :rtype: float
        """
        if self._a is None:
            self._a = self.get_a()
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def b(self):
        """Minor (horizontal) axis.

        :return: Length of minor axis, b, from center of ellipse.
        :rtype:
        """
        if self._b is None:
            self._b = self.get_random_b()
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    @property
    def k(self):
        # TODO: Check position of origin.
        """Y-offset from center of ellipse to origin (top of fretboard)

        :return: y-offset (major axis - depth)
        :rtype: float
        """
        if self._k is None:
            self._k = self.a - self.depth
        return self._k

    @property
    def y(self):
        """Equation of ellipse expressed in terms of y.

        :return: Equation of ellipse expressed in terms of y.
        :rtype: #TODO: Wtf
        """
        if self._y is None:
            a, b, x, y, k = sp.symbols('a b x y k')
            eq = sp.solve(x ** 2 / b ** 2 + (y - k) ** 2 / a ** 2 - 1, y)[0]
            self._y = eq.subs(a, self.a).subs(k, self.k).subs(b, self.b)
        return self._y

    @property
    def hradius(self):
        """Unused.

        :return: Horizontal radius.
        :rtype: float
        """
        return self.ellipse.hradius

    @property
    def vradius(self):
        """Unused.

        :return: Vertical radius
        :rtype: float
        """
        return self.ellipse.vradius

    @property
    def center(self):
        """Unused.

        :return: Center of ellipse.
        :rtype: sp.geometry.Point
        """
        return self.ellipse.center

    @staticmethod
    def get_random_in_range(minimum, maximum):
        """Generate random number between minimum and maximum values. This
        function is inclusive of the minimum value, non-inclusive of the maximum.

        :param minimum: minimum inclusive value for random number
        :type minimum: float
        :param maximum: maximum non-inclusive value for random number
        :type maximum: float
        :return: Random number between [minimum, maximum)
        :rtype: float
        """
        return minimum + np.random.random() * (maximum - minimum)

    def get_a_from_b(self, _b, _x, _t):
        """From a known minor horizontal axis (b), and point in space (x, t),
        find the major vertical axis (b).

        :param _b: Minor horizontal axis
        :type _b: float
        :param _x: x-coordinate
        :type _x: float
        :param _t: Negative y-offset of major axis
        :type _t: float
        :return: Major vertical axis
        :rtype: float
        """
        a = sp.symbols('a')
        return max(sp.solve(
            _x ** 2 / _b ** 2 + (self.y_int - (a - _t)) ** 2 / a ** 2 - 1, a
        ))

    def get_random_b(self):
        """Generate random horizontal minor radius.

        :return: Random number greater than half the width of the neck,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(self.width, self.width * 2)

    def get_a(self):
        """Get value for major axis a at nut, given a value for minor axis b.

        :param _b: Minor (horizontal) axis
        :type _b: float
        :return:
        :rtype:
        """
        return self.get_a_from_b(self.b, self.width, self.depth)

    def get_ellipse(self):
        """Return Ellipse as sympy object.

        :return: sympy ellipse object.
        :rtype: sp.geometry.Ellipse
        """
        return Ellipse(
            center=(0, self.a - self.depth),
            vradius=self.a,
            hradius=self.b,
        )


