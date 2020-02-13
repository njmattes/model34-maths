#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, atan, log
import numpy as np
from model34_maths.neck_section import NeckSection
from model34_maths.config import SCALE_LENGTH


class Neck(object):
    def __init__(self, n):
        """Neck made up of n NeckSections.

        """
        nut = NeckSection(neck_position=0)
        body = NeckSection(neck_position=SCALE_LENGTH/2)

        # Values from first prototype
        body.inner.b = 1.95
        body.outer.b = 2.8
        nut.inner.b = .826
        nut.outer.b = 1.66
        # End values from first prototype

        self.sections = [nut] + [
            NeckSection(neck_position=_x/10*SCALE_LENGTH/2,
                        min_inner_b=nut.inner.b, max_inner_b=body.inner.b,
                        min_outer_b=nut.outer.b, max_outer_b=body.outer.b)
            for _x in np.arange(1, n-1)
        ] + [body]

        # sum of the load of all strings, Psi
        self.load = 36.5 + 42 + 51.3 + 42.8
        # Young's modulus of aluminum in Psi
        self.elastic_modulus = 10**6
        self.string_to_fretboard_distance = .125
        self._mean_load_centroid_distance = None
        self.n = n

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
        xs = np.linspace(0, 1, self.n, dtype='float') * SCALE_LENGTH / 2
        fit = np.polynomial.polynomial.Polynomial.fit(xs, moments, 2)
        return list(fit)

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
            (2 * (2 * a * _x + b) * atan((2 * a * _x + b) /
                                         sqrt(4 * a * c - b**2)) -
             sqrt(4 * a * c - b**2) * log(_x * (a * _x + b) + c)) /
            2 * a * sqrt(4 * a * c - b**2)
        )

    def curvature2(self, _x):
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
    print('inner widths:', [s.inner.width for s in neck.sections])
    print('outer widths:', [s.outer.width for s in neck.sections])
    print('inner bs:', neck.sections[0].inner.b, neck.sections[-1].inner.b, )
    print('outer bs:', neck.sections[0].outer.b, neck.sections[-1].outer.b, )
    print('I_xx:', [
        '{0:.4f}'.format(_.second_moment_area(_.inner.width)) for _ in neck.sections
    ])
    for s in neck.sections:
        print('{0:.2f}'.format(s.neck_position),
              '{0:.6f}'.format(neck.curvature(s.neck_position)))
