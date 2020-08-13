#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import fabs
import numpy as np
from scipy.optimize import minimize
from model34_maths.neck import Neck


class OptimalNeck(object):
    def __init__(self, neck):
        self.neck = neck
        self.x0 = np.random.ranf(11) * \
                  np.array([
                      1.5, 1.5, 1.5, 1.5, 1, 1, .125, .03125, .0625, .0625,
                      .125]) + \
                  np.array([
                      2, .5, 2, .5, 2, 2, .125, .03125, .5, .75 + .0625, .125])
        self.target_area = 5.0634014771373765  # Weight with full inner rib
        self.target_area = 3.93  # Weight with open inner rib
        self.target_deflection = .191

    def __call__(self, *args, **kwargs):
        self.minimize()

    @property
    def bounds(self):
        return np.array([
            [self.neck.widths[0], self.neck.widths[0] * 6],
            [self.neck.widths[1], self.neck.widths[1] * 2],
            [self.neck.widths[2], self.neck.widths[2] * 6],
            [self.neck.widths[3], self.neck.widths[3] * 2],
            [2, 3],  # Outer curvature
            [2, 4],  # Inner curvature
            [.125, .5],  # Rail
            [.0, .125],  # Shell thickness
            [.5, .625],  # Depth at nut
            [.75 + .0625, .75 + 4 * .0625],  # Depth at body
            [.0625, .125],  # Inner rib
            # [.0625, .1875],
        ])

    @property
    def constraints(self):
        return dict(type='ineq', fun=self.get_area)
            #      fun=self.get_deflection),

    def get_area(self, x):
        self.set_neck(x)
        scalar = .95
        diff = self.target_area * scalar - sum([s.area for s in self.neck.sections])
        return diff

    def get_deflection(self, x):
        self.set_neck(x)
        return self.target_deflection * .8 - self.neck.deflection[-1]

    def set_neck(self, x):
        self.neck.reset()
        self.neck.minor_bs[0] = x[0]
        self.neck.minor_bs[1] = x[1]
        self.neck.minor_bs[2] = x[2]
        self.neck.minor_bs[3] = x[3]
        self.neck.r_outer = x[4]
        self.neck.r_inner = x[5]
        self.neck.rail_width = x[6]
        self.neck.nut_shell = x[7]
        self.neck.octave_shell = .0625 + x[7]
        self.neck.nut_depth = x[8]
        self.neck.octave_depth = x[9]
        self.neck.inner_rib = x[10]
        # self.neck.inner_rib = x[8]
        self.neck.set_model(medial_rib=True)

    def objective(self, x):
        self.set_neck(x)
        return (
            # (sum([s.area for s in self.neck.sections]) / self.target_area) *
            # 100 * (self.neck.deflection[-1]) ** 2
            np.sum(
                np.arange(len(self.neck.deflection)) *
                np.array(self.neck.deflection) ** 2)
        )

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
