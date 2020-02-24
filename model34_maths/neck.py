#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, atan, log
import numpy as np
from numpy.polynomial.polynomial import polyfit
from model34_maths.neck_section import NeckSection
from model34_maths.neck_area import NeckArea
from model34_maths.neck_curve import NeckCurve
from model34_maths.config import SCALE_LENGTH, NeckConfig


class Neck(NeckConfig):
    def __init__(self, n):
        """Neck made up of n NeckSections.

        """

        minor_as = self.get_random_in_range(self.widths * 1.1,
                                            self.widths * 2)
        minor_bs = self.get_b(
            minor_as, self.widths, self.depths,
            np.array([self.rail_depth, 0, self.rail_depth, 0])
        )

        self.sections = [
            NeckSection([NeckArea(
                0,
                self.scale_value(i / n, *self.widths[0::2]),
                NeckCurve(
                    self.scale_value(i / n, *minor_as[0::2]),
                    self.scale_value(i / n, *minor_bs[0::2]),
                    self.scale_value(i / n, *minor_bs[0::2] -
                                             self.depths[0::2]),
                    self.n0),
                NeckCurve(
                    self.scale_value(i / n, *minor_as[1::2]),
                    self.scale_value(i / n, *minor_bs[1::2]),
                    self.scale_value(i / n, *minor_bs[1::2] - self.depths[1::2]),
                    self.n1),
            ), NeckArea(
                self.scale_value(i / n, *self.widths[0::2]),
                self.scale_value(i / n, *self.widths[1::2]),
                NeckCurve(0, 0, 0, 0),
                NeckCurve(
                    self.scale_value(i / n, *minor_as[1::2]),
                    self.scale_value(i / n, *minor_bs[1::2]),
                    self.scale_value(i / n, *minor_bs[1::2] - self.depths[1::2]),
                    self.n1), ), ])
            for i in reversed(range(n + 1))
        ]

        # sum of the load of all strings, Psi
        self.load = 36.5 + 42 + 51.3 + 42.8
        self.elastic_modulus = 10 * 10 ** 6
        self.string_to_fretboard_distance = .125
        self._mean_load_centroid_distance = None
        self._p1 = None
        self._p2 = None
        self._p3 = None
        self._deflection = None
        self.n = n

    @property
    def p1(self):
        if self._p1 is None:
            (self._p3, self._p2, self._p1) = self.curve_fit_area_moment()
        return self._p1

    @property
    def p2(self):
        if self._p2 is None:
            (self._p3, self._p2, self._p1) = self.curve_fit_area_moment()
        return self._p2

    @property
    def p3(self):
        if self._p3 is None:
            (self._p3, self._p2, self._p1) = self.curve_fit_area_moment()
        return self._p3

    @property
    def deflection(self):
        if self._deflection is None:
            self._deflection = [
                self.curvature(i, section)
                for i, section in enumerate(self.sections)
            ]
        return self._deflection

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
        moments = np.array([section.second_moment
                            for section in self.sections], dtype='float')
        xs = np.linspace(0, 1, self.n+1, dtype='float') * SCALE_LENGTH / 2
        fit = polyfit(xs, moments, 2)
        return list(fit)

    def curvature(self, x, section):
        """Find deflection of neck at position _x.

        :param x:
        :type x:
        :param section:
        :type section:
        :return:
        :rtype:
        """
        m = (self.string_to_fretboard_distance - section.y_bar) * self.load
        return (m / self.elastic_modulus) * (
            ((2 * self.p1 * x + self.p2) *
             atan((2 * self.p1 * x + self.p2) /
                  sqrt(4 * self.p1 * self.p3 - self.p2 ** 2))) /
            (self.p1 * sqrt(4 * self.p1 * self.p3 - self.p2 ** 2)) -
            log(x * (self.p1 * x + self.p2) + self.p3) / (2 * self.p1)
            # x
        )

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

    def get_random_n(self):
        """Generate random exponent.

        :return: Random number greater that 2 and less than 4,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(2, 3)

    def get_b(self, a, w, d, y):
        """Get value for major axis a at nut, given a value for minor axis b.

        :return:
        :rtype:
        """
        return a * np.sqrt(1 / ((a - w) * (a + w))) * (a - d - y)

    @staticmethod
    def scale_value(x, a, b):
        return a + (b - a) * x


if __name__ == '__main__':
    neck = Neck(10)
    # for _ in neck.sections:
    #     print('A:', _.area, 'y:', _.y_bar, 'I:', _.second_moment, )
    print(np.array(neck.deflection))
