#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import quad
from mystic.penalty import quadratic_inequality
from mystic.constraints import as_constraint
from model34_maths.neck import Neck


x0 = np.array([.8125, .8125, 1.125, 1.125, 2, 2, 2, 2, .125, .05, .5, .8125, .125])
ub = np.array([4, 4, 6, 6, 4, 4, 4, 4, .25, .125, .625, 1, .25])
lb = np.array([.6875, .8125, 1, 1.125, 2, 2, 2, 2, .125, .05, .5, .8125, .125])
bounds = np.stack((lb, ub)).T
neck = Neck()


def nut_inner_b(x):
    return (x[10] - x[9]) - x[0]


def nut_outer_b(x):
    return x[10] - x[1]


def body_inner_b(x):
    return (x[11] - x[9]) - x[2]


def body_outer_b(x):
    return x[11] - x[3]


def inner_a(x):
    return ((1.125 - x[8]) * (1 - np.abs((x[2] - (x[11] - x[9])) / x[2]) **
                              x[6])) ** (-1 / x[6])


def inner_y(x, xx):
    print(inner_a(x))
    return (x[2] - (x[11] - x[9])) - x[2] * (
            1 - np.abs(xx / inner_a(x)) ** x[6]) ** (1 / x[6])


def outer_a(x):
    return (1.125 * (1 - np.abs((x[3] - x[11]) / x[3]) ** x[7])) ** (-1 / x[7])


def outer_y(x, xx):
    print(outer_a(x))
    return (x[3] - x[11]) - x[3] * (
            1 - np.abs(xx / outer_a(x)) ** x[7]) ** (1 / x[7])


def curves(x):
    set_neck(x)
    return outer_y(x, 1.125 / 2) - inner_y(x, 1.125 / 2)


def area(x):
    print('x', x)
    def f(xx):
        print(xx, inner_y(x, xx), outer_y(x, xx))
        return inner_y(x, xx) - outer_y(x, xx)
    return quad(f, 0, 1.125 - x[8])[0] #- 0.05281695966578498


@quadratic_inequality(area, k=1e4)
@quadratic_inequality(curves, k=1e4)
@quadratic_inequality(nut_inner_b, k=1e4)
@quadratic_inequality(nut_outer_b, k=1e4)
@quadratic_inequality(body_inner_b, k=1e4)
@quadratic_inequality(body_outer_b, k=1e4)
def penalty(x):
    return 0.0


solver = as_constraint(penalty)


def constraint(x):
    return x


def set_neck(x):
    neck.reset()
    neck.endpoint_bs = np.array(x)[0:4].reshape((2, 2))
    neck.exponents = np.array(x)[4:8].reshape((2, 2))
    neck.rail_width = x[8]
    neck.nut_shell = x[9]
    # self.neck.octave_shell = .0625 + x[9]
    neck.octave_shell = x[11] - .75 + x[9]
    neck.nut_depth = x[10]
    neck.octave_depth = x[11]
    neck.inner_rib = x[12]
    neck.set_model()


def objective(x):
    set_neck(x)
    return 10 * np.sum(
        np.arange(len(neck.deflection)) *
        np.array(neck.deflection) ** 2)


def solve1():
    from mystic.solvers import diffev2
    from mystic.math import almostEqual
    from mystic.monitors import Monitor, VerboseMonitor
    mon = VerboseMonitor(10)

    result = diffev2(objective, x0=x0, bounds=bounds,
                     penalty=penalty, #constraints=constraint,
                     npop=20, ftol=1e-8, gtol=200, disp=True,
                     full_output=True, cross=0.8, scale=0.9,
                     itermon=mon)

    print(result[0])


if __name__ == '__main__':
    print(area([3.99997929, 0.81250005, 5.99999573, 1.12500008, 2.00000162,
                3.99999996, 2.00000007, 3.9999999, 0.19250634, 0.12499998,
                0.625, 0.99999999, 0.25]))
    print(area([.625, 4, 1, 2,
                2, 2, 2, 2,
                0.8754, 0.0625,
                0.5, 0.8125, 0.125]))
    # solve1()
