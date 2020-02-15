#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from math import sqrt, atan
import sympy as sp
from sympy.geometry import Ellipse


class NeckEllipse(object):
    def __init__(self, width, depth, y_int, a=None):
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
        :param a: length of minor (horizontal) axis if known
        :type a: float
        """
        self.width = width
        self.depth = depth
        self.y_int = y_int
        self._b = None
        self._a = a
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
    def b(self):
        """Major (vertical) axis.

        :return: Length of major axis, a, from center of ellipse.
        :rtype: float
        """
        if self._b is None:
            self._b = self.get_b()
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    @property
    def a(self):
        """Minor (horizontal) axis.

        :return: Length of minor axis, b, from center of ellipse.
        :rtype:
        """
        if self._a is None:
            self._a = self.get_random_a()
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def k(self):
        # TODO: Check position of origin.
        """Y-offset from center of ellipse to origin (top of fretboard)

        :return: y-offset (major axis - depth)
        :rtype: float
        """
        if self._k is None:
            self._k = self.b - self.depth
        return self._k

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

    def get_y(self, x):
        """

        :return:
        :rtype:
        """
        return self.k - self.b * sqrt(self.a ** 2 - x ** 2) / self.a

    def get_b_from_a(self, _a, _x, _t):
        """From a known horizontal axis (a), and point in space (x, t),
        find the vertical axis (b).

        :param _a: Horizontal axis
        :type _a: float
        :param _x: x-coordinate
        :type _x: float
        :param _t: Negative y-offset of major axis
        :type _t: float
        :return: Major vertical axis
        :rtype: float
        """
        a = sp.symbols('a')
        return max(sp.solve(
            _x ** 2 / _a ** 2 + (self.y_int - (a - _t)) ** 2 / a ** 2 - 1, a
        ))

    def get_random_a(self):
        """Generate random horizontal radius.

        :return: Random number greater than half the width of the neck,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(self.width, self.width * 2)

    def get_b(self):
        """Get value for major axis a at nut, given a value for minor axis b.

        :param _b: Minor (horizontal) axis
        :type _b: float
        :return:
        :rtype:
        """
        return self.get_b_from_a(self.a, self.width, self.depth)
