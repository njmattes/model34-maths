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
    def __init__(self, a_min=1.2, a_max=6, n=10):
        """Neck made up of n NeckSections.

        """

        self.n = n
        self.minor_bs = self.get_random_in_range(self.depths,
                                                 self.depths * a_max)
        self.minor_as = None
        self.sections = None
        self._p1 = None
        self._p2 = None
        self._p3 = None
        self._deflection = None
        self.set_model()

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
    def get_a(b, w, d, y, n):
        """Get value for major axis a at nut, given a value for minor axis b.

        :return:
        :rtype:
        """
        return w * (1 - np.abs((b - d + y) / b) ** n) ** (-1 / n)

    def get_random_n(self):
        """Generate random exponent.

        :return: Random number greater that 2 and less than 4,
        and less than twice that.
        :rtype: float
        """
        return self.get_random_in_range(2, 3)

    def get_inner_curve(self, i):
        return NeckCurve(
            self.scale_value(i / self.n, *self.minor_as[1::2]),
            self.scale_value(i / self.n, *self.minor_bs[1::2]),
            self.scale_value(i / self.n,
                             *self.minor_bs[1::2] - self.depths[1::2]),
            self.r_inner)

    def get_outer_curve(self, i):
        return NeckCurve(
            self.scale_value(i / self.n, *self.minor_as[0::2]),
            self.scale_value(i / self.n, *self.minor_bs[0::2]),
            self.scale_value(i / self.n,
                             *self.minor_bs[0::2] - self.depths[0::2]),
            self.r_outer)

    def get_inner_rail(self, i):
        """Calculate horizontal location of inner edge of the fretbaord rail.

        :param i: Position index along neck
        :type i: int
        :return: x-coordinate of inner rail edge
        :rtype: float
        """
        return self.scale_value(i / self.n, *self.widths[1::2])

    def get_outer_rail(self, i):
        """Calculate horizontal location of outer edge of the fretbaord rail.

        :param i: Position index along neck
        :type i: int
        :return: x-coordinate of outer rail edge
        :rtype: float
        """
        return self.scale_value(i / self.n, *self.widths[0::2])

    def build_areas(self, i):
        areas = list()

        # Area under the horizontal rail and above the curved face along
        # the outer edge of the neck.
        edge_area = NeckArea(
            self.get_inner_rail(i),   # x0: inside of rail
            self.get_outer_rail(i),   # x1: outside edge of rail
            self.get_outer_curve(i),  # y0: curved face of neck
            NeckCurve(1, 0, 0, 2),    # y1: top face of rail
        )

        areas.append(edge_area)

        # If the curved face at position i extends below 0.75 inches
        if self.get_outer_curve(i).y(0) < -.75 - 1e-8:
            # Get the horizontal point at which the curved face intersects
            # the face at 0.75 inches
            break_x = self.get_outer_curve(i).x(-.75)
            areas.append(NeckArea(
                0,                         # x0: middle of neck
                break_x,                   # x1: intersection
                NeckCurve(0, 0, -.75, 0),  # y0: -0.75 face
                self.get_inner_curve(i)))  # y1: inner curved face
            areas.append(NeckArea(
                break_x,                   # x0: intersection
                self.get_inner_rail(i),    # x1: inner rail edge
                self.get_outer_curve(i),   # y0: outer curved face
                self.get_inner_curve(i)))  # y1: inner curved face
        else:
            areas.append(NeckArea(
                0,                         # x0: middle of neck
                self.get_inner_rail(i),    # x1: inner rail edge
                self.get_outer_curve(i),   # y0: outer curved face
                self.get_inner_curve(i)))  # y1: inner curved face

        return areas

    def build_areas_with_solid_rib(self, i):
        areas = self.build_areas(i)
        rib_area = NeckArea(  # inner rail
            0, .125+self.inner_rib,
            self.get_inner_curve(i), NeckCurve(1, 0, 0, 2),
        )
        areas.append(rib_area)
        return areas

    def build_areas_with_open_rib(self, i):
        areas = self.build_areas(i)
        rib_area = NeckArea(  # inner rail
            .125,
            .125+self.inner_rib,
            self.get_inner_curve(i),
            NeckCurve(1, 0, 0, 2),
        )
        areas.append(rib_area)
        return areas

    def build_areas_with_medial_rib(self, i):
        areas = self.build_areas(i)
        rib_thickness = .06
        medial_rib_area = NeckArea(
            (self.get_inner_rail(i)-.25)/2+.25-rib_thickness/2,
            (self.get_inner_rail(i)-.25)/2+.25+rib_thickness/2,
            self.get_inner_curve(i), NeckCurve(1, 0, 0, 2),
        )
        areas.append(medial_rib_area)
        return areas

    def build_areas_with_longitudinal_ribs(self, i, n):
        areas = self.build_areas(i)
        rib_thickness = .06
        for j in range(n):
            areas.append(NeckArea(
                (self.get_inner_rail(i)-(.125+self.inner_rib))/n*(j+1)-rib_thickness/2,
                (self.get_inner_rail(i)-(.125+self.inner_rib))/n*(j+1)+rib_thickness/2,
                self.get_inner_curve(i), NeckCurve(1, 0, 0, 2),
            ))
        return areas

    def set_model(self, solid_rib=False,
                  medial_rib=False, longitudinal_ribs=False):
        self.depths = np.array([
            self.octave_depth,
            self.octave_depth - self.octave_shell,
            self.nut_depth,
            self.nut_depth - self.nut_shell,
        ])
        self.minor_as = self.get_a(
            self.minor_bs, self.widths, self.depths,
            np.array([self.rail_depth, 0, self.rail_depth, 0]),
            np.array([self.r_outer, self.r_inner, self.r_outer, self.r_inner]),
        )
        if solid_rib:
            self.sections = [
                NeckSection(self.build_areas_with_solid_rib(i), i / self.n)
                for i in range(self.n + 1)
            ]
        elif medial_rib:
            self.sections = [
                NeckSection(self.build_areas_with_medial_rib(i), i / self.n)
                for i in range(self.n + 1)
            ]
        elif longitudinal_ribs and longitudinal_ribs > 0:
            self.sections = [
                NeckSection(self.build_areas_with_longitudinal_ribs(i, longitudinal_ribs), i / self.n)
                for i in range(self.n + 1)
            ]
        else:
            self.sections = [
                NeckSection(self.build_areas_with_open_rib(i), i / self.n)
                for i in range(self.n + 1)
            ]

    def reset(self):
        self._p1 = None
        self._p2 = None
        self._p3 = None
        self._deflection = None


if __name__ == '__main__':
    neck = Neck()
    neck.reset()
    neck.minor_bs = np.array([3.63100775, 0.9375, 2.2802993, 0.625])
    neck.r_outer = 2.7
    neck.r_inner = 3.
    neck.rail_width = 0.16
    neck.nut_shell = 0.032
    neck.octave_shell = .0625 + 0.032
    neck.nut_depth = 0.55
    neck.octave_depth = 0.817
    neck.inner_rib = 0.25
    neck.set_model()

    for _ in neck.sections:
        print('A:', _.area, 'y:', _.y_bar, 'I:', _.second_moment)
    print('Asum:', sum([section.area for section in neck.sections]))
    print(np.array(neck.deflection))
    print(neck.depths)

    neck = Neck()
    neck.reset()
    neck.minor_bs = np.array([6.75, 0.9375, 4.875, 0.625 ])
    neck.r_outer = 3
    neck.r_inner = 3.
    neck.rail_width = 0.16
    neck.nut_shell = 0.032
    neck.octave_shell = .0625 + 0.032
    neck.nut_depth = 0.55
    neck.octave_depth = 0.817
    neck.inner_rib = 0.25
    neck.set_model()
    self.minor_as = self.get_a(
        self.minor_bs, self.widths, self.depths,
        np.array([self.rail_depth, 0, self.rail_depth, 0]),
        np.array([self.r_outer, self.r_inner, self.r_outer, self.r_inner]),
    )