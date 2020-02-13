#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, atan


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
             (-4 * self.y_bar * self.inner.b**2 * self.outer.b**2 *
              self.inner.k) +
             (-4 * self.y_bar * self.inner.b**2 * self.outer.b**2 *
              self.outer.k) +
             (2 * self.inner.b**2 * self.outer.b**2 * self.inner.k *
              self.outer.k) +
             (4 * self.y_bar * self.inner.a * self.inner.b * self.outer.b**2 *
              sqrt(self.inner.b**2 - _x**2)) +
             (4 * self.y_bar * self.inner.b**2 * self.outer.b * self.outer.a *
              sqrt(self.outer.b**2 - _x**2)) +
             (-2 * self.inner.a * self.inner.b * self.outer.b**2 *
              self.inner.k * sqrt(self.inner.b**2 - _x**2)) +
             (-2 * self.inner.a * self.inner.b * self.outer.b**2 *
              self.outer.k * sqrt(self.inner.b**2 - _x**2)) +
             (-2 * self.inner.b**2 * self.outer.b * self.inner.k *
              self.outer.a * sqrt(self.outer.b**2 - _x**2)) +
             (-2 * self.inner.b**2 * self.outer.b * self.outer.k *
              self.outer.a * sqrt(self.outer.b**2 - _x**2)) +
             (2 * self.inner.a * self.inner.b * self.outer.b * self.outer.a *
              sqrt(self.inner.b**2 - _x**2) * sqrt(self.outer.b**2 - _x**2))
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


