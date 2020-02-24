#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt
from model34_maths.decorators import finite_method
from scipy.integrate import quad


N = 3


class NeckArea(object):
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self._area = None
        self._first_moment = None
        self._second_moment = None
        self._y_bar = None

    @property
    def area(self):
        if self._area is None:
            self._area = self.integral_da()[0]
        return self._area

    @property
    def first_moment(self):
        if self._first_moment is None:
            def f(x):
                return self.y(x) ** 2 / 2
            self._first_moment = quad(f, self.x_min, self.x_max)[0]
        return self._first_moment

    @property
    def second_moment(self):
        if self._second_moment is None:
            def f(x):
                return (((self.y_max(x) - self.y_min(x)) ** 3) / 12 +
                        (self.y_max(x) - self.y_min(x)) *
                        ((self.y_max(x) + self.y_min(x)) / 2) ** 2)
            self._second_moment = quad(f, self.x_min, self.x_max)[0]
        return self._second_moment

    @property
    def y_bar(self):
        if self._y_bar is None:
            self._y_bar = self.integral_yda()[0] / self.integral_da()[0]
        return self._y_bar

    def y(self, x):
        return (self.y_max(x) + self.y_min(x)) / 2

    def integral_da(self):
        def f(x):
            return self.y_max(x) - self.y_min(x)
        return quad(f, self.x_min, self.x_max)

    def integral_yda(self):
        def f(x):
            return self.y(x) * (self.y_max(x) - self.y_min(x))
        return quad(f, self.x_min, self.x_max)


if __name__ == '__main__':
    from model34_maths.neck_curve import NeckCurve
    y0 = NeckCurve(1, 1, 0, 2)
    y1 = NeckCurve(.75, .75, 0, 2)
    area = NeckArea(0, .75, y0, y1)
    print(area.integral_y())
    print('area:', area.area)
    print('y_bar:', area.y_bar)
    print('1st:', area.first_moment)
    print('2nd:', area.second_moment)
