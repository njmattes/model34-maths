#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, atan
import numpy as np
from sympy.geometry import Ellipse
import sympy as sp


NUT_WIDTH = .8125
BODY_WIDTH = 1.125
NECK_LENGTH = 23.672
FACET_WIDTH = .62
NUT_DEPTH = .5
BODY_DEPTH = .75 + .06
FRETBOARD_THICKNESS = .125
RAIL_WIDTH = .125
NECK_THICKNESS = .063
SCALE_LENGTH = 34


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


class NeckSection(object):
    def __init__(self, neck_position, min_outer_b=None, max_outer_b=None,
                 min_inner_b=None, max_inner_b=None):
        """Section of neck made up of simple(r) shapes.

        #TODO: Finish docstring
        :param neck_position: position along neck, 0 <= x <= neck_length
        :type neck_position: float
        :param min_inner_b:
        :type min_inner_b: float
        :param max_inner_b:
        :type max_inner_b: float
        :param min_outer_b:
        :type min_outer_b: float
        :param max_outer_b:
        :type max_outer_b: float
        """
        self.neck_position = neck_position
        neck_position = neck_position / (SCALE_LENGTH / 2)
        self.width = NUT_WIDTH + neck_position * (BODY_WIDTH - NUT_WIDTH)
        self.depth = NUT_DEPTH + neck_position * (BODY_DEPTH - NUT_DEPTH)
        self.inner = NeckEllipse(
            self.width - RAIL_WIDTH, self.depth - NECK_THICKNESS,
            FRETBOARD_THICKNESS,
            b=None if min_inner_b is None else
            min_inner_b+(max_inner_b-min_inner_b)*neck_position)
        self.outer = NeckEllipse(
            self.width, self.depth, 0,
            b=None if min_outer_b is None else
            min_outer_b+(max_outer_b-min_outer_b)*neck_position)
        self._y_bar = None

    @property
    def y_bar(self):
        """y-bar. y component of centroid. (x component is assumed to be 0).

        :return: y-bar
        :rtype: float
        """
        if self._y_bar is None:
            self._y_bar = (
                    self.integral_yda_between_ellipses(self.inner.width) /
                    self.integral_da_between_ellipses(self.inner.width))
        return self._y_bar

    def integral_da_between_ellipses(self, _x):
        """Integral of dA in the denominator of x-bar and y-bar.

        y = k - a * sqrt(b**2 - x**2) / b
        dA = [y of inner ellipse] - [y of outer ellipse]

        Integral( (ki - ai * sqrt(bi**2 - x**2) / bi) -
                  (ko - ao * sqrt(bo**2 - x**2) / bo) )

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return (
                -self.inner.a * _x * sqrt(self.inner.b ** 2 - _x ** 2) /
                (2 * self.inner.b) + self.inner.a * self.inner.b *
                atan((_x * sqrt(self.inner.b ** 2 - _x ** 2)) /
                     (_x ** 2 - self.inner.b ** 2)) / 2 +
                self.outer.a * _x * sqrt(self.outer.b ** 2 - _x ** 2) /
                (2 * self.outer.b) - self.outer.a * self.outer.b *
                atan((_x * sqrt(self.outer.b ** 2 - _x ** 2)) /
                     (_x ** 2 - self.outer.b ** 2)) / 2 +
                _x * (self.inner.a - self.inner.depth) -
                _x * (self.outer.a - self.outer.depth)
        )

    def integral_xda_between_ellipses(self, _x):
        """Integral of x * dA in the numerator of x-bar.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return \
            ((2 * self.inner.a * (self.inner.b**2 - _x**2)**1.5) /
             self.inner.b -
             (2 * self.outer.a * (self.outer.b**2 - _x**2)**1.5) /
             self.outer.b +
             3 * _x**2 * (self.inner.k - self.outer.k)) / 6

    def integral_yda_between_ellipses(self, _x):
        """Integral of y * dA in the numerator of y-bar.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return \
            (self.inner.b**2 * self.outer.b**2 * _x * (
                self.inner.a**2 + self.inner.k**2 -
                self.outer.k**2 - self.outer.a**2) +
             _x**3 * (self.inner.b**2 * self.outer.a**2 -
                      self.outer.b**2 * self.inner.a**2) / 3 -
             self.inner.a * self.inner.b * self.outer.b**2 *
             self.inner.k * _x * sqrt(self.inner.b**2 - _x**2) -
             self.inner.a * self.inner.b**3 * self.outer.b**2 * self.inner.k *
             atan(_x / sqrt(self.inner.b**2 - _x**2)) +
             self.inner.b**2 * self.outer.b * self.outer.k * self.outer.a *
             _x * sqrt(self.outer.b**2 - _x**2) + self.inner.b**2 *
             self.outer.b**3 * self.outer.k * self.outer.a * atan(
                        _x / sqrt(self.outer.b**2 - _x**2))) / \
            (2 * self.inner.b**2 * self.outer.b**2)

    def second_moment_area(self, _x):
        """Integral of (y - y_bar)**2.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return \
            ((4 * self.y_bar**2 * self.inner.b**2 * self.outer.b**2) +
             (self.inner.a**2 * self.inner.b**2 * self.outer.b**2) +
             (-self.inner.a**2 * self.outer.b**2 * _x**2) +
             (self.inner.b**2 * self.outer.b**2 * self.inner.k**2) +
             (self.inner.b**2 * self.outer.b**2 * self.outer.k**2) +
             (self.inner.b**2 * self.outer.b**2 * self.outer.a**2) +
             (-self.inner.b**2 * self.outer.a**2 * _x**2) +
             (-4 * self.y_bar * self.inner.b**2 * self.outer.b**2 * self.inner.k) +
             (-4 * self.y_bar * self.inner.b**2 * self.outer.b**2 * self.outer.k) +
             (2 * self.inner.b**2 * self.outer.b**2 * self.inner.k * self.outer.k) +
             (4 * self.y_bar * self.inner.a * self.inner.b * self.outer.b**2 * sqrt(self.inner.b**2 - _x**2)) +
             (4 * self.y_bar * self.inner.b**2 * self.outer.b * self.outer.a * sqrt(self.outer.b**2 - _x**2)) +
             (-2 * self.inner.a * self.inner.b * self.outer.b**2 * self.inner.k * sqrt(self.inner.b**2 - _x**2)) +
             (-2 * self.inner.a * self.inner.b * self.outer.b**2 * self.outer.k * sqrt(self.inner.b**2 - _x**2)) +
             (-2 * self.inner.b**2 * self.outer.b * self.inner.k * self.outer.a * sqrt(self.outer.b**2 - _x**2)) +
             (-2 * self.inner.b**2 * self.outer.b * self.outer.k * self.outer.a * sqrt(self.outer.b**2 - _x**2)) +
             (2 * self.inner.a * self.inner.b * self.outer.b * self.outer.a * sqrt(self.inner.b**2 - _x**2) * sqrt(self.outer.b**2 - _x**2))
             )/(4 * self.inner.b**2 * self.outer.b**2)

    def second_moment_area2(self, _x):
        """Broken. Also stupid. Do not use.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return \
            ((-self.inner.a * _x * sqrt(self.inner.b**2 - _x**2) /
              2 * self.inner.b) + (self.inner.a * self.inner.b *
                                   atan(_x * sqrt(self.inner.b**2 - _x**2) /
                                        (_x**2 - self.inner.b**2))) / 2 -
             self.outer.a * _x * sqrt(self.outer.b**2 - _x**2) /
             (2 * self.outer.b) + (self.outer.b * self.outer.a * atan(
                        _x * sqrt(self.outer.b**2 - _x**2) /
                        (_x**2 - self.outer.b**2))) / 2 +
             self.inner.k * _x + self.outer.k * _x - 2 * _x * self.y_bar)


class Neck(object):
    def __init__(self, n):
        """Neck made up of n NeckSections.

        """
        nut = NeckSection(neck_position=0)
        body = NeckSection(neck_position=17)

        # Values from first prototype
        body.inner.b = 1.95
        body.outer.b = 2.8
        body.inner.b = .826
        body.outer.b = 1.66
        # End values from first prototype

        self.sections = [nut] + [
            NeckSection(neck_position=_x/10*17, min_inner_b=nut.inner.b,
                        max_inner_b=body.inner.b, min_outer_b=nut.outer.b,
                        max_outer_b=body.outer.b) for _x in np.arange(1, n-1)
        ] + [body]

        # sum of the load of all strings, Psi
        self.load = 36.5 + 42 + 51.3 + 42.8
        # Young's modulus of aluminum in Psi
        self.elastic_modulus = 10**6
        self.string_to_fretboard_distance = .125
        self._mean_load_centroid_distance = None

    @property
    def mean_load_centroid_distance(self):
        """Average distance of strings from y-centroid in each NeckSection.

        :return:
        :rtype: float
        """
        if self._mean_load_centroid_distance is None:
            self._mean_load_centroid_distance = (
                    np.array([_.y_bar for _ in self.sections]).mean() +
                    self.string_to_fretboard_distance)
        return self._mean_load_centroid_distance

    def curve_fit_area_moment(self):
        """Fit a polynomial to the area moments of inertia of each NeckSection.

        :return:
        :rtype:
        """
        moments = np.array([section.second_moment_area(section.inner.width)
                            for section in self.sections], dtype='float')
        xs = np.linspace(0, 1, 10, dtype='float')
        fit = np.polynomial.polynomial.Polynomial.fit(xs, moments, 2)
        return fit

    def curvature(self, _x):
        """Find deflection of neck at position _x.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        m = self.mean_load_centroid_distance * self.load
        e = self.elastic_modulus
        (a, b, c) = self.curve_fit_area_moment()
        return (m / e) * (
                (2 * atan(
                    (2 * a * _x + b) /
                    sqrt(4 * a * c - b**2))) /
                sqrt(4 * a * c - b**2))


if __name__ == '__main__':
    neck = Neck(n=10)
    print(neck.sections[0].inner.b, neck.sections[0].outer.b, )
    print(neck.sections[-1].inner.b, neck.sections[-1].outer.b, )
    for s in neck.sections:
        print(s.neck_position, neck.curvature(s.neck_position))
