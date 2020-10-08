#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from model34_maths.optimizer import OptimalNeck
from model34_maths.neck import Neck


def original_neck():
    original_neck = Neck()

    original_neck.endpoint_bs = np.array([[1, 2], [.625, 4], ])
    original_neck.inner_rib = .125
    original_neck.rail_depth = 0
    # Calculate area / weight with solid inner rib
    original_neck.set_model(solid_rib=True)
    print('A0[2016]: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
        original_neck.sections[0].area, original_neck.sections[-1].area,
        sum([section.area for section in original_neck.sections])))

    # Set target area for optimal neck
    original_neck.reset()

    # Calculate deflection with open inner rib (screw hole voids affect
    # deflection).
    # original_neck.inner_rib = .125
    original_neck.set_model()
    print('A0[2016]: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
        original_neck.sections[0].area, original_neck.sections[-1].area,
        sum([section.area for section in original_neck.sections])))
    return original_neck


def optimal_neck():
    new_neck = Neck()

    # Set new neck parameters to optimal parameters
    # FIXME::
    new_neck.endpoint_bs = optimal.x[0:4].reshape((2, 2))
    new_neck.exponents = optimal.x[4:8].reshape((2, 2))
    new_neck.rail_width = optimal.x[8]
    new_neck.nut_shell = optimal.x[9]
    new_neck.octave_shell = optimal.x[11] - .75 + optimal.x[9]
    # new_neck.octave_shell = .0625 + optimal.x[9]
    new_neck.nut_depth = optimal.x[10]
    new_neck.octave_depth = optimal.x[11]
    new_neck.inner_rib = optimal.x[12]

    new_neck.set_model(solid_rib=True)
    print('A0[2020]: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
        new_neck.sections[0].area, new_neck.sections[-1].area,
        sum([section.area for section in new_neck.sections])))

    # new_neck.set_model()
    # print('A0: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
    #     new_neck.sections[0].area, new_neck.sections[-1].area,
    #     sum([section.area for section in new_neck.sections])))
    return new_neck


if __name__ == '__main__':

    optimize_run = True

    if optimize_run:
        opt_neck = OptimalNeck(Neck())
        optimal = opt_neck.minimize()
    else:
        optimal = None

    original_neck = original_neck()
    deflection_a = original_neck.deflection[-1]
    print(deflection_a)

    if optimal.success:

        new_neck = optimal_neck()
        deflection_b = new_neck.deflection[-1]
        print('d[nut][2016]: {:.3f}, d[sum]: {:.3f}'.format(
            deflection_a, sum(original_neck.deflection)))
        print('d[nut][2020]: {:.3f}, d[sum]: {:.3f}'.format(
            deflection_b, sum(new_neck.deflection)))

        # Reduction in deflection
        print('∂d: {:.1f}%\n'.format(
            (deflection_a - deflection_b) * 100 / deflection_a))

        print('b: {}\n'
              'r: {}\n'
              'rail: {}\n'
              'shell: {}\n'
              'depth: {}\n'
              'rib: {}'.format(
                    optimal.x[0:4], optimal.x[4:8], optimal.x[8], optimal.x[9],
                    optimal.x[10:12], optimal.x[12]))

        # print(new_neck.aa)
        # print(new_neck.bs)
        # print(new_neck.depths)
        # print(new_neck.widths)
        print(original_neck.deflection)
        print(new_neck.deflection)

        import matplotlib.pyplot as plt
        for i in [0, len(new_neck.sections) - 1]:
            section = new_neck.sections[i]
            for area in section.areas:
                x = np.linspace(area.x_min, area.x_max, 100)
                plt.plot(x, [area.y_min(n) for n in x])
                plt.plot(x, [area.y_max(n) for n in x])
        plt.show()

    else:

        print('Optimization unsuccessful.')
