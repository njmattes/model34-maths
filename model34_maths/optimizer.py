#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import fabs
import numpy as np
from scipy.optimize import minimize
from model34_maths.neck import Neck


X0 = np.array([
    [.5, 1.5, 0, 20, ],  # b octave inner
    [2, 1.5, 0, 20, ],  # b octave outer
    [.5, 1.5, 0, 20, ],       # b nut inner
    [2, 1.5, 0, 20, ],        # b nut outer
    [2, 2, 2, 3.5, ],               # r nut inner
    [1, 2, 2, 3.5, ],               # r nut outer
    [2, 2, 2, 3.2, ],               # r octave inner
    [1.5, 2, 2.5, 3.2, ],             # r octave outer
    [.125, .125, .125, .5, ],     # rail width
    [.03125, .03125, 0, .125, ],  # shell
    [.0625, .5, .5, .625, ],      # nut depth
    [.0625, .8125, .8125, 1, ],   # octave depth
    [.0625, .0625, .0625, .125, ],  # inner rib
])
x0 = np.random.ranf(X0.shape[0]) * X0[:, 0] + X0[:, 1]


class OptimalNeck(object):
    def __init__(self, neck):
        self.neck = neck
        self.x0 = x0
        # self.target_area = 5.0634014771373765  # Weight with full inner rib
        self.target_area = 3.68  # Weight with open inner rib
        # self.target_area = 3.68  # Weight with open inner rib
        # self.target_area = 3.93  # Weight with open inner rib??????
        # self.target_area = 8.542788824745898  # New curve fitting for sections
        self.target_deflection = .191

    def __call__(self, *args, **kwargs):
        self.minimize()

    @property
    def bounds(self):
        bounds = X0[:, 2:]
        bounds[0] = self.neck.widths[0, 0] * np.array([1, 4])
        bounds[1] = self.neck.widths[0, 1] * np.array([1, 6])
        bounds[2] = self.neck.widths[1, 0] * np.array([1, 4])
        bounds[3] = self.neck.widths[1, 1] * np.array([1, 6])
        return bounds

    @property
    def constraints(self):
        return [
            # Area less than target
            dict(type='ineq', fun=self.get_area),
            # inner faces inside outer faces
            dict(type='ineq', fun=self.reasonable_curves),
            # octave inner b > .75 thickness of neck
            dict(type='ineq', fun=lambda x: x[0] - .75),
            # octave outer b > octave depth
            dict(type='ineq', fun=lambda x: x[1] - x[11]),
            # nut inner b > nut depth - shell
            dict(type='ineq', fun=lambda x: x[2] - (x[10] - x[9])),
            # nut outer b > nut depth
            dict(type='ineq', fun=lambda x: x[3] - x[10]),
        ]

    def get_area(self, x):
        self.set_neck(x)
        scalar = .9
        total_area = sum([s.area for s in self.neck.sections])
        diff = self.target_area * scalar - total_area
        return diff

    def reasonable_curves(self, x):
        self.set_neck(x)
        reasonable_thickness = .1
        xs = np.linspace(0, self.neck.nut_width - x[8], 20)
        return np.min(
            self.neck.get_inner_curve(0).y(xs) -
            (self.neck.get_outer_curve(0).y(xs) + reasonable_thickness))

    def get_deflection(self, x):
        self.set_neck(x)
        return self.target_deflection * .8 - self.neck.deflection[-1]

    def set_neck(self, x):
        self.neck.reset()
        self.neck.endpoint_bs = x[0:4].reshape((2, 2))
        self.neck.exponents = x[4:8].reshape((2, 2))
        self.neck.rail_width = x[8]
        self.neck.nut_shell = x[9]
        # self.neck.octave_shell = .0625 + x[9]
        self.neck.octave_shell = x[11] - .75 + x[9]
        self.neck.nut_depth = x[10]
        self.neck.octave_depth = x[11]
        self.neck.inner_rib = x[12]
        self.neck.set_model()

    def objective(self, x):
        self.set_neck(x)
        return 10 * np.sum(
            np.arange(len(self.neck.deflection)) *
            np.array(self.neck.deflection) ** 2)
        # (sum([s.area for s in self.neck.sections]) / self.target_area) *
        # 100 * (self.neck.deflection[-1]) ** 2

    def minimize(self):
        return minimize(
            self.objective,
            self.x0,
            method='SLSQP',
            bounds=self.bounds,
            constraints=self.constraints,
            options={
                'disp': True,
                'iprint': 10,
                'maxiter': 1000,
                # 'eps': 2e-10,
                'ftol': 1e-4,
            })
