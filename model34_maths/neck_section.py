#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, atan
from model34_maths.neck_ellipse import NeckEllipse
from model34_maths.config import BODY_DEPTH, BODY_WIDTH
from model34_maths.config import NUT_DEPTH, NUT_WIDTH
from model34_maths.config import RAIL_WIDTH, SCALE_LENGTH
from model34_maths.config import NECK_THICKNESS, FRETBOARD_THICKNESS


class NeckSection(object):
    def __init__(self, neck_position, min_outer_a=None, max_outer_a=None,
                 min_inner_a=None, max_inner_a=None):
        """Section of neck made up of simple(r) shapes.

        #TODO: Finish docstring
        :param neck_position: position along neck, 0 <= x <= neck_length
        :type neck_position: float
        :param min_inner_a:
        :type min_inner_a: float
        :param max_inner_a:
        :type max_inner_a: float
        :param min_outer_a:
        :type min_outer_a: float
        :param max_outer_a:
        :type max_outer_a: float
        """
        self.neck_position = neck_position
        neck_position = neck_position / (SCALE_LENGTH / 2)
        self.width = NUT_WIDTH + neck_position * (BODY_WIDTH - NUT_WIDTH)
        self.depth = NUT_DEPTH + neck_position * (BODY_DEPTH - NUT_DEPTH)
        self.inner = NeckEllipse(
            self.width - RAIL_WIDTH, self.depth - NECK_THICKNESS,
            FRETBOARD_THICKNESS,
            a=None if min_inner_a is None else
            min_inner_a + (max_inner_a - min_inner_a) * neck_position)
        self.outer = NeckEllipse(
            self.width, self.depth, 0,
            a=None if min_outer_a is None else
            min_outer_a + (max_outer_a - min_outer_a) * neck_position)
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
                    (2*self.integral_da_between_ellipses(self.inner.width)))
        return self._y_bar

    def integral_da_between_ellipses(self, x):
        """Integral of dA in the denominator of x-bar and y-bar.

        y = k - a * sqrt(b**2 - x**2) / b
        dA = [y of inner ellipse] - [y of outer ellipse]

        Integral( (ki - ai * sqrt(bi**2 - x**2) / bi) -
                  (ko - ao * sqrt(bo**2 - x**2) / bo) )

        integrate (k - b * sqrt(a^2-x^2)/a) - (K - B * sqrt(A^2-x^2)/A) dx

        :param x:
        :type x:
        :return:
        :rtype:
        """
        a = self.inner.a
        A = self.outer.a
        b = self.inner.b
        B = self.outer.b
        k = self.inner.k
        K = self.outer.k
        eq = (
            -b * x * sqrt(a**2 - x**2) / (2 * a) +
            (a * b * atan(x * sqrt(a**2 - x**2) / (x**2 - a**2))) / 2 +
            B * x * sqrt(A**2 - x**2) / (2 * A) -
            (A * B * atan(x * sqrt(A**2 - x**2) / (x**2 - A**2))) / 2 +
            k * x - K * x
        )
        return eq

    def integral_yda_between_ellipses(self, x):
        """Integral of y * dA in the numerator of y-bar.

        :param _x:
        :type _x:
        :return:
        :rtype:

        integrate ((k - b * sqrt(a^2-x^2)/a) - ((k - b * sqrt(a^2-x^2)/a) - (K - B * sqrt(A^2-x^2)/A)) / 2) * ((k - b * sqrt(a^2-x^2)/a) - (K - B * sqrt(A^2-x^2)/A)) dx


        """
        a = self.inner.a
        A = self.outer.a
        b = self.inner.b
        B = self.outer.b
        k = self.inner.k
        K = self.outer.k
        eq = (-1 / (6 * a**2 * A**2)) * (
                -3 * a**2 * A**2 * b**2 * x + 3 * a * A**2 * b * k * x *
                sqrt(a**2 - x**2) + 3 * a**2 * A**2 * B**2 * x -
                3 * a**2 * A * B * K * x * sqrt(A**2 - x**2) -
                3 * a**2 * A**2 * k**2 * x +
                3 * a**2 * A**2 * K**2 * x -
                3 * a**2 * A**3 * B * K * atan(x / sqrt(A**2 - x**2)) -
                a**2 * B**2 * x**3 +
                3 * a**3 * A**2 * b * k * atan(x / sqrt(a**2 - x**2)) +
                A**2 * b**2 * x**3)
        return eq

    def second_moment_area(self, _x):
        """Integral of (y - y_bar)**2.

        :param _x:
        :type _x:
        :return:
        :rtype:
        """
        return \
            ((4 * self.y_bar ** 2 * self.inner.a ** 2 * self.outer.a ** 2) +
             (self.inner.b ** 2 * self.inner.a ** 2 * self.outer.a ** 2) +
             (-self.inner.b ** 2 * self.outer.a ** 2 * _x ** 2) +
             (self.inner.a ** 2 * self.outer.a ** 2 * self.inner.k ** 2) +
             (self.inner.a ** 2 * self.outer.a ** 2 * self.outer.k ** 2) +
             (self.inner.a ** 2 * self.outer.a ** 2 * self.outer.b ** 2) +
             (-self.inner.a ** 2 * self.outer.b ** 2 * _x ** 2) +
             (-4 * self.y_bar * self.inner.a ** 2 * self.outer.a ** 2 *
              self.inner.k) +
             (-4 * self.y_bar * self.inner.a ** 2 * self.outer.a ** 2 *
              self.outer.k) +
             (2 * self.inner.a ** 2 * self.outer.a ** 2 * self.inner.k *
              self.outer.k) +
             (4 * self.y_bar * self.inner.b * self.inner.a * self.outer.a ** 2 *
              sqrt(self.inner.a ** 2 - _x ** 2)) +
             (4 * self.y_bar * self.inner.a ** 2 * self.outer.a * self.outer.b *
              sqrt(self.outer.a ** 2 - _x ** 2)) +
             (-2 * self.inner.b * self.inner.a * self.outer.a ** 2 *
              self.inner.k * sqrt(self.inner.a ** 2 - _x ** 2)) +
             (-2 * self.inner.b * self.inner.a * self.outer.a ** 2 *
              self.outer.k * sqrt(self.inner.a ** 2 - _x ** 2)) +
             (-2 * self.inner.a ** 2 * self.outer.a * self.inner.k *
              self.outer.b * sqrt(self.outer.a ** 2 - _x ** 2)) +
             (-2 * self.inner.a ** 2 * self.outer.a * self.outer.k *
              self.outer.b * sqrt(self.outer.a ** 2 - _x ** 2)) +
             (2 * self.inner.b * self.inner.a * self.outer.a * self.outer.b *
              sqrt(self.inner.a ** 2 - _x ** 2) * sqrt(self.outer.a ** 2 - _x ** 2))
             )/(4 * self.inner.a ** 2 * self.outer.a ** 2)
