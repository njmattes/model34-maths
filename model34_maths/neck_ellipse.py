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
        self._n = None

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

    @property
    def n(self):
        if self._n is None:
            self._n = self.get_random_n()
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

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
        return self.k - self.b * sqrt(
            self.a ** self.n - x ** self.n) / self.a**(self.n / 2)

    def get_random_a(self):
        """Generate random horizontal radius.

        :return: Random number greater than half the width of the neck,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(self.width, self.width * 2)

    def get_random_n(self):
        """Generate random exponent.

        :return: Random number greater that 2 and less than 4,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(2, 4)

    def get_b(self):
        """Get value for major axis a at nut, given a value for minor axis b.

        :return:
        :rtype:
        """
        return (self.a * sqrt(1 / ((self.a - self.width) *
                                   (self.a + self.width))) *
                (self.a - self.depth - self.y_int))
