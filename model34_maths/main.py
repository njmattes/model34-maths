#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from model34_maths.optimizer import OptimalNeck
from model34_maths.neck import Neck


if __name__ == '__main__':
    opt_neck = OptimalNeck(Neck())
    optimal = opt_neck.minimize()

    if optimal.success:

        original_neck = Neck()
        new_neck = Neck()

        original_neck.minor_bs = np.array([2, 1, 4, .625])
        original_neck.inner_rib = .125
        original_neck.rail_depth = 0
        # Calculate area / weight with solid inner rib
        original_neck.set_model(solid_rib=True)
        print('A0: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
            original_neck.sections[0].area, original_neck.sections[-1].area,
            sum([section.area for section in original_neck.sections])))

        # Set target area for optimal neck
        new_neck.target_area = sum([
            section.area for section in original_neck.sections])
        original_neck.reset()

        # Calculate deflection with open inner rib (screw hole voids affect
        # deflection).
        original_neck.inner_rib = .125
        original_neck.set_model()
        print('A0: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
            original_neck.sections[0].area, original_neck.sections[-1].area,
            sum([section.area for section in original_neck.sections])))
        deflection_a = original_neck.deflection[-1]

        # Set new neck parameters to optimal parameters
        new_neck.minor_bs[0] = optimal.x[0]
        new_neck.minor_bs[1] = optimal.x[1]
        new_neck.minor_bs[2] = optimal.x[2]
        new_neck.minor_bs[3] = optimal.x[3]
        new_neck.r_outer = optimal.x[4]
        new_neck.r_inner = optimal.x[5]
        new_neck.rail_width = optimal.x[6]
        new_neck.nut_shell = optimal.x[7]
        new_neck.octave_shell = .0625 + optimal.x[7]
        new_neck.nut_depth = optimal.x[8]
        new_neck.octave_depth = optimal.x[9]
        new_neck.inner_rib = optimal.x[10]

        new_neck.set_model(solid_rib=True, medial_rib=True)
        print('A0: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
            new_neck.sections[0].area, new_neck.sections[-1].area,
            sum([section.area for section in new_neck.sections])))

        new_neck.set_model()
        print('A0: {:.3f}, A1 {:.3f}, Asum: {:.2f}'.format(
            new_neck.sections[0].area, new_neck.sections[-1].area,
            sum([section.area for section in new_neck.sections])))

        deflection_b = new_neck.deflection[-1]
        print('d[nut]: {:.3f}, d[sum]: {:.3f}'.format(
            deflection_a, sum(original_neck.deflection)))
        print('d[nut]: {:.3f}, d[sum]: {:.3f}'.format(
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
                    optimal.x[0:4], optimal.x[4:6], optimal.x[6], optimal.x[7],
                    optimal.x[8:10], optimal.x[10]))

        print(new_neck.minor_as)
        print(new_neck.minor_bs)
        print(new_neck.depths)
        print(new_neck.widths)
        print(original_neck.deflection)
        print(new_neck.deflection)

    else:

        print('Optimization unsuccessful.')
