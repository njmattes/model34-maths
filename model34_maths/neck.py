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
        body.inner.a = 1.95
        body.outer.a = 2.8
        nut.inner.a = .826
        nut.outer.a = 1.66
        # End values from first prototype

        self.sections = [nut] + [
            NeckSection(neck_position=_x/10*SCALE_LENGTH/2,
                        min_inner_a=nut.inner.a, max_inner_a=body.inner.a,
                        min_outer_a=nut.outer.a, max_outer_a=body.outer.a)
            for _x in np.arange(1, n-1)
        ] + [body]

        # sum of the load of all strings, Psi
        self.load = 36.5 + 42 + 51.3 + 42.8
        # Young's modulus of aluminum in Psi
        self.elastic_modulus = 10**6
        self.string_to_fretboard_distance = .125
        self._mean_load_centroid_distance = None
        self.n = n
        self._p1 = None
        self._p2 = None
        self._p3 = None
        self._deflection = None

    @property
    def p1(self):
        if self._p1 is None:
            (self._p1, self._p2, self._p3) = self.curve_fit_area_moment()
        return self._p1

    @property
    def p2(self):
        if self._p2 is None:
            (self._p1, self._p2, self._p3) = self.curve_fit_area_moment()
        return self._p2

    @property
    def p3(self):
        if self._p3 is None:
            (self._p1, self._p2, self._p3) = self.curve_fit_area_moment()
        return self._p3

    @property
    def deflection(self):
        if self._deflection is None:
            self._deflection = [
                neck.curvature(s.neck_position) for s in neck.sections
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
        moments = np.array([section.second_moment_area
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
        return 2 * (m / e) * (
            (2 * (2 * self.p1 * _x + self.p2) *
             atan((2 * self.p1 * _x + self.p2) /
                  sqrt(4 * self.p1 * self.p3 - self.p2**2)) -
             sqrt(4 * self.p1 * self.p3 - self.p2**2) *
             log(_x * (self.p1 * _x + self.p2) + self.p3)) /
            (2 * self.p1 * sqrt(4 * self.p1 * self.p3 - self.p2**2))
        )


if __name__ == '__main__':
    neck = Neck(n=10)
    print('inner widths:', [s.inner.width for s in neck.sections])
    print('outer widths:', [s.outer.width for s in neck.sections])
    print('inner as:', [s.inner.a for s in neck.sections])
    print('inner bs:', [s.inner.b for s in neck.sections])
    print('outer As:', [s.outer.a for s in neck.sections])
    print('outer Bs:', [s.outer.b for s in neck.sections])
    print('inner ks:', [s.inner.k for s in neck.sections])
    print('outer Ks:', [s.outer.k for s in neck.sections])
    print('y-bars:', [s.y_bar for s in neck.sections], )
    print('I_xx:', [
        '{0:.4f}'.format(_.second_moment_area) for _ in neck.sections
    ])
    for s in neck.sections:
        print('{0:.2f}'.format(s.neck_position),
              '{0:.6f}'.format(neck.curvature(s.neck_position)))
