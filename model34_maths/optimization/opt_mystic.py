#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial.distance import cdist
from mystic.differential_evolution import DifferentialEvolutionSolver2
from mystic.penalty import quadratic_inequality, barrier_inequality
from mystic.constraints import and_
from mystic.symbolic import generate_constraint, generate_solvers, simplify
from mystic.monitors import VerboseMonitor, Monitor
from mystic.termination import ChangeOverGeneration
from model34_maths.third_party.intersection.intersect import intersection
from model34_maths.neck import Neck


def set_neck(neck, x):
    """Set parameters of Neck object.

    :param neck:
    :type neck:
    :param x:
    :type x:
    :return:
    :rtype:
    """
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


def get_solver(n=100):
    # Initial guess represents parameters of current neck
    x0 = np.array([1.125, 1.125, .8125, .8125, 2, 2, 2, 2,
                   .125, .05, .5, .8125, .125])
    # Upper bounds
    ub = np.array([4, 6, 4, 6, 4, 4, 4, 4,
                   .25, .125, .625, 1, .25])
    # Lower bounds
    lb = np.array([1, 1.125, .6875, .8125, 2, 2, 2, 2,
                   .125, .01, .5, .8125, .125])
    _neck = Neck(n=n)

    equations = """
    x0 >= (.75 - x9)
    x1 >= x11
    x2 >= (x10 - x9)    x3 >= x10  
    """
    # x9 >= .05
    # x9 <= .125
    # x10 >= .5
    # x10 <= .625
    # x11 >= .8125
    # x11 <= 1
    eqn = simplify(equations, all=True)
    constraints = generate_constraint(generate_solvers(eqn), join=and_)

    def area(x, neck):
        """Apply penalty if area of sections is larger than target area.

        :param x:
        :type x:
        :param neck:
        :type neck:
        :return:
        :rtype:
        """
        target_area = {100: 40.27, 50: 20.13, 20: 8.05, 10: 3.7}[n]
        pen = sum([section.area for section in neck.sections]) / target_area
        return pen - 1

    def intersect(x, neck):
        """

        :param x:
        :type x:
        :param neck:
        :type neck:
        :return: 2 if intersection, 0 if none
        :rtype:
        """
        t = np.linspace(neck.endpoint_ts[1, 0], np.pi / 2, 100)
        x1 = (neck.endpoint_as[1, 0] * np.cos(t) **
              (2 / neck.exponents[1, 0]))
        y1 = (-neck.endpoint_bs[1, 0] * np.sin(t) **
              (2 / neck.exponents[1, 0]) + neck.endpoint_ks[1, 0])
        x2 = (neck.endpoint_as[1, 1] * np.cos(t) **
              (2 / neck.exponents[1, 1]))
        y2 = (-neck.endpoint_bs[1, 1] * np.sin(t) **
              (2 / neck.exponents[1, 1]) + neck.endpoint_ks[1, 1])
        return max([len(xy) for xy in intersection(x1, y1, x2, y2)])

    def thickness(x, neck, i):
        xys = np.empty((4, 100))
        for j in range(2):
            t = np.linspace(neck.endpoint_ts[i, j], np.pi / 2, 100)
            idx = j * 2
            xys[idx, :] = (neck.endpoint_as[i, j] * np.cos(t) **
                           (2 / neck.exponents[i, j]))
            xys[idx + 1, :] = (-neck.endpoint_bs[i, j] * np.sin(t) **
                               (2 / neck.exponents[i, j]) +
                               neck.endpoint_ks[i, j])
        return np.min(cdist(xys[0:2, :].T, xys[2:4, :].T))

    def body_min_thickness(x, neck):
        return 1 - (thickness(x, neck, 0) / neck.nut_shell)

    def body_max_thickness(x, neck):
        return 1 - (ub[9] / thickness(x, neck, 0))

    def nut_min_thickness(x, neck):
        return 1 - (thickness(x, neck, 1) / neck.nut_shell)

    def nut_max_thickness(x, neck):
        return 1 - (ub[9] / thickness(x, neck, 1))

    # (2 * k) ** (5 ** n) * f(x) ** 2
    @quadratic_inequality(area, k=1e+3, kwds=dict(neck=_neck))
    @quadratic_inequality(body_min_thickness, k=1e+3, kwds=dict(neck=_neck))
    @quadratic_inequality(body_max_thickness, k=1e+3, kwds=dict(neck=_neck))
    @quadratic_inequality(nut_min_thickness, k=1e+3, kwds=dict(neck=_neck))
    @quadratic_inequality(nut_max_thickness, k=1e+3, kwds=dict(neck=_neck))
    @quadratic_inequality(intersect, k=5, kwds=dict(neck=_neck))  #
    def penalty(x):
        return 0.0

    def objective(x, neck):
        # Penalty functions are called after objective, so set neck model here.
        set_neck(neck, x)
        obj = np.sum(
            1 / np.arange(len(neck.deflection), 0, -1) *
            np.array(neck.deflection) ** 2 ** .5)
        return 10 * obj

    eval_monitor = Monitor()
    step_monitor = VerboseMonitor(10)
    n_pop = 10
    f_tol = 1e-8
    g_tol = 80
    n_x0 = len(x0)

    _solver = DifferentialEvolutionSolver2(n_x0, n_pop * n_x0)
    _solver.SetRandomInitialPoints(min=lb, max=ub)
    _solver.SetStrictRanges(min=lb, max=ub)
    _solver.SetPenalty(penalty)
    _solver.SetConstraints(constraints)
    _solver.SetEvaluationMonitor(eval_monitor)
    _solver.SetGenerationMonitor(step_monitor)
    _solver.SetObjective(objective, ExtraArgs=(_neck, ))
    _solver.enable_signal_handler()
    _solver.SetTermination(ChangeOverGeneration(f_tol, g_tol))

    return _solver


if __name__ == '__main__':
    from datetime import datetime
    cross_probability = 0.9
    scaling_factor = 0.8
    solver = get_solver()
    solver.Solve(disp=True,
                 CrossProbability=cross_probability,
                 ScalingFactor=scaling_factor)
    # solver.Step(disp=True,
    #             CrossProbability=cross_probability,
    #             ScalingFactor=scaling_factor)
    np.save('{}-solutions.npy'.format(
        datetime.now().strftime('%Y%m%d%H%M%S')),
        solver.solution_history)
    np.save('{}-energy.npy'.format(
        datetime.now().strftime('%Y%m%d%H%M%S')),
        solver.energy_history)
    print(', '.join([str(x) for x in solver.bestSolution]))
    print(solver.bestEnergy)
