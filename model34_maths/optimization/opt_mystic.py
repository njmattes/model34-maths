#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import quad
from mystic.penalty import quadratic_inequality, barrier_inequality
from mystic.constraints import as_constraint, and_
from mystic.symbolic import generate_constraint, generate_solvers, simplify
from model34_maths.neck import Neck


x0 = np.array([1.125, 1.125, .8125, .8125, 2, 2, 2, 2, .125, .05, .5, .8125, .125])
ub = np.array([4, 4, 4, 4, 4, 4, 4, 4, .25, .125, .625, 1, .25])
lb = np.array([1, 1.125, .6875, .8125, 2, 2, 2, 2, .125, .05, .5, .8125, .125])
bounds = np.stack((lb, ub)).T
n = Neck()

equations = """
(x11 - x9) - x0 <= 0
x11 - x1 <= 0
(x10 - x9) - x2 <= 0
x10 - x3 <= 0  
"""
eqn = simplify(equations, all=True)
constraints = generate_constraint(generate_solvers(eqn), join=and_)


def area(x, neck):
    set_neck(neck, x)
    target_area = 3.7
    return sum([section.area for section in neck.sections]) - target_area


def curves(x, neck):
    set_neck(neck, x)
    xs = np.linspace(0, neck.octave_width - neck.rail_width, 20)
    return np.max(
        (neck.get_outer_curve(0).y(xs) + .08) -
        neck.get_inner_curve(0).y(xs))


@quadratic_inequality(area, kwds=dict(neck=n))
@quadratic_inequality(curves, kwds=dict(neck=n))
def penalty(x):
    return 0.0


def set_neck(neck, x):
    neck.reset()
    neck.endpoint_bs = np.array(x)[0:4].reshape((2, 2))
    neck.exponents = np.array(x)[4:8].reshape((2, 2))
    neck.rail_width = x[8]
    neck.nut_shell = x[9]
    neck.octave_shell = x[11] - .75 + x[9]
    neck.nut_depth = x[10]
    neck.octave_depth = x[11]
    neck.inner_rib = x[12]
    neck.set_model()


def objective(x, neck):
    set_neck(neck, x)
    return 10 * np.sum(
        # np.arange(len(neck.deflection), 0, -1) *
        np.arange(len(neck.deflection)) *
        np.array(neck.deflection) ** 2)


def solve1():
    from mystic.solvers import diffev2
    from mystic.monitors import VerboseMonitor
    mon = VerboseMonitor(10)

    result = diffev2(objective, x0=x0, bounds=bounds,
                     penalty=penalty,
                     constraints=constraints,
                     args=(n, ),
                     npop=20, ftol=1e-8,
                     gtol=20,
                     disp=True,
                     full_output=True, cross=0.8, scale=0.9,
                     itermon=mon)

    print(result[0])


if __name__ == '__main__':
    # print(area([5.99999573, 1.12500008, 3.99997929, 0.81250005, 2.00000162,
    #             3.99999996, 2.00000007, 3.9999999, 0.19250634, 0.12499998,
    #             0.625, 0.99999999, 0.25]))
    # print(area([1, 2, .625, 4,
    #             2, 2, 2, 2,
    #             0.1875, 0.0625,
    #             0.5, 0.8125, 0.125]))
    solve1()
    print(n.deflection)
