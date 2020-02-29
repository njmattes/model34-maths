#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import dblquad
from numpy.polynomial.polynomial import polyfit
from model34_maths.neck_section import NeckSection
from model34_maths.neck_area import NeckArea
from model34_maths.neck_curve import NeckCurve
from model34_maths.config import NeckConfig


class Neck(NeckConfig):
    def __init__(self, a_min=1.2, a_max=2.5, n=10):
        """Neck made up of n NeckSections.

        """

        self.n = n
        self.minor_as = self.get_random_in_range(self.widths * a_max,
                                                 self.widths * a_min)
        self.minor_bs = None
        self.sections = None
        self._mean_load_centroid_distance = None
        self._p1 = None
        self._p2 = None
        self._p3 = None
        self._deflection = None
        self.set_model()

    def set_model(self):
        self.minor_bs = self.get_b(
            self.minor_as, self.widths, self.depths,
            np.array([self.rail_depth, 0, self.rail_depth, 0])
        )
        self.sections = [
            NeckSection(self.build_areas(i), i / self.n)
            for i in range(self.n + 1)
        ]

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
                self.curvature(section)
                for section in self.sections
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
        xs = np.linspace(0, 1, self.n+1, dtype='float') * self.scale / 2
        fit = polyfit(xs, moments, 2)
        return list(fit)

    def curvature(self, section):
        """Find deflection of neck at position _x.

        :param section:
        :type section:
        :return:
        :rtype:
        """
        def f(y, x):
            return 1 / (x ** 2 * self.p1 + x * self.p2 + self.p3)

        return (((self.string_to_fretboard_distance - section.y_bar) *
                 self.load / self.elastic_modulus) *
                dblquad(f, 0, self.scale / 2 * section.z,
                        lambda x: 0, lambda x: self.scale / 2 * section.z)[0])

    @staticmethod
    def scale_value(x, a, b):
        return a + (b - a) * x

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

    @staticmethod
    def get_b(a, w, d, y):
        """Get value for major axis a at nut, given a value for minor axis b.

        :return:
        :rtype:
        """
        return a * np.sqrt(1 / ((a - w) * (a + w))) * (a - d - y)

    def get_random_n(self):
        """Generate random exponent.

        :return: Random number greater that 2 and less than 4,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(2, 3)

    def build_areas(self, i):
        areas = list()

        outer_rail = self.scale_value(i / self.n, *self.widths[0::2])
        inner_rail = self.scale_value(i / self.n, *self.widths[1::2])

        outer_curve = NeckCurve(
            self.scale_value(i / self.n, *self.minor_as[0::2]),
            self.scale_value(i / self.n, *self.minor_bs[0::2]),
            self.scale_value(i / self.n,
                             *self.minor_bs[0::2] - self.depths[0::2]),
            self.n0)
        inner_curve = NeckCurve(
            self.scale_value(i / self.n, *self.minor_as[1::2]),
            self.scale_value(i / self.n, *self.minor_bs[1::2]),
            self.scale_value(i / self.n,
                             *self.minor_bs[1::2] - self.depths[1::2]),
            self.n1)

        edge_area = NeckArea(
            inner_rail, outer_rail, outer_curve, NeckCurve(1, 0, 0, 2),
        )
        rib_area = NeckArea(  # inner rail
            0, .25, inner_curve, NeckCurve(1, 0, 0, 2),
        )

        areas.append(edge_area)
        areas.append(rib_area)

        if outer_curve.y(0) < -.75 - 1e-8:
            break_x = outer_curve.x(-.75)
            print('break_x:', break_x)
            areas.append(NeckArea(
                0,
                break_x,
                NeckCurve(0, 0, -.75, 0),
                inner_curve))
            areas.append(NeckArea(
                break_x, inner_rail,
                outer_curve, inner_curve))
        else:
            areas.append(NeckArea(
                0, inner_rail,
                outer_curve, inner_curve))

        return areas


if __name__ == '__main__':
    neck = Neck()
    for _ in neck.sections:
        print('A:', _.area, 'y:', _.y_bar, 'I:', _.second_moment)
    print(np.array(neck.deflection))
