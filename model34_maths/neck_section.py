#!/usr/bin/env python
# -*- coding: utf-8 -*-
from scipy.integrate import dblquad
from model34_maths.neck_area import NeckArea
from model34_maths.neck_curve import NeckCurve


class NeckSection(object):
    def __init__(self, areas, z):
        self.areas = areas
        self._area = None
        self._y_bar = None
        self._first_moment = None
        self._second_moment = None
        self.x_min = min([area.x_min for area in self.areas])
        self.x_max = max([area.x_max for area in self.areas])
        self.z = z

    @property
    def area(self):
        if self._area is None:
            return 2 * sum([_.area for _ in self.areas])
        return self._area

    @property
    def y_bar(self):
        if self._y_bar is None:
            self._y_bar = (
                    sum([_.area * _.y_bar for _ in self.areas]) /
                    sum([_.area for _ in self.areas]))
        return self._y_bar

    @property
    def first_moment(self):
        if self._first_moment is None:
            pass
        return self._first_moment

    @property
    def second_moment(self):
        if self._second_moment is None:
            self._second_moment = 2 * sum([
                self.get_second_moment_of_area(area) for area in self.areas
            ])
        return self._second_moment

    def get_second_moment_of_area(self, area):
        def f(y, x):
            return (y - self.y_bar) ** 2

        return dblquad(f, area.x_min, area.x_max,
                       lambda x: area.y_min(x),
                       lambda x: area.y_max(x))[0]


if __name__ == '__main__':
    y0 = NeckCurve(1, 1, 0, 2)
    y1 = NeckCurve(.75, .75, 0, 2)
    area0 = NeckArea(0, .75, y0, y1)
    y1 = NeckCurve(0, 0, 0, 0)
    area1 = NeckArea(.75, 1, y0, y1)
    section = NeckSection([area0, area1])
    print('A:', section.area)
    print('y_bar:', section.y_bar)
    print('I:', section.second_moment)
