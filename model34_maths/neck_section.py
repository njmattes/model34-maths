#!/usr/bin/env python
# -*- coding: utf-8 -*-
from model34_maths.neck_area import NeckArea
from model34_maths.neck_curve import NeckCurve


class NeckSection(object):
    def __init__(self, areas):
        self.areas = areas
        self._area = None
        self._y_bar = None
        self._first_moment = None
        self._second_moment = None
        self.x_min = min([area.x_min for area in self.areas])
        self.x_max = max([area.x_max for area in self.areas])

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
            self._second_moment = sum([
                2 * (_.second_moment + _.area * (_.y_bar - self.y_bar) ** 2)
                for _ in self.areas])
        return self._second_moment


if __name__ == '__main__':
    y0 = NeckCurve(1, 1, 0, 2)
    y1 = NeckCurve(.75, .75, 0, 2)
    area0 = NeckArea(0, .75, y0, y1)
    y1 = NeckCurve(0, 0, 0, 0)
    area1 = NeckArea(.75, 1, y0, y1)
    section = NeckSection([area0, area1])
    print(section.area)
    print(section.y_bar)
    print(section.second_moment)
