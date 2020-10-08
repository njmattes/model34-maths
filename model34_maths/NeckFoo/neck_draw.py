#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from model34_maths.neck import Neck


def draw_areas(neck):
    for i in [0, len(neck.sections)-1]:
        section = neck.sections[i]
        for area in section.areas:
            # print(area.x_min, area.x_max)
            x = np.linspace(area.x_min, area.x_max, 100)
            plt.plot(x, [area.y_min(n) for n in x])
            plt.plot(x, [area.y_max(n) for n in x])
    #     x = np.linspace(0, neck.octave_width -
    #                     (neck.octave_width - neck.nut_width) *
    #                     i / (len(neck.sections) - 1), 50)
    #     plt.plot(x, curve.y(x))
    # plt.axvline(x=neck.octave_width)
    # plt.axvline(x=neck.nut_width)
    # plt.axhline(y=0)
    # plt.axhline(y=neck.rail_depth)
    # plt.axhline(y=-neck.nut_depth)
    # plt.axhline(y=-neck.octave_depth)
    plt.show()


def draw_sections(neck):
    for i in range(len(neck.sections)):
    # for i in range(1):
        curve = neck.get_outer_curve(i)
        x = np.linspace(0, neck.octave_width -
                        (neck.octave_width - neck.nut_width) *
                        i / (len(neck.sections) - 1), 50)
        plt.plot(x, curve.y(x))
    plt.axvline(x=neck.octave_width)
    plt.axvline(x=neck.nut_width)
    plt.axhline(y=0)
    plt.axhline(y=neck.rail_depth)
    plt.axhline(y=-neck.nut_depth)
    plt.axhline(y=-neck.octave_depth)
    plt.show()

    # for i in range(1):
    for i in range(len(neck.sections)):
        curve = neck.get_inner_curve(i)
        x = np.linspace(0, (neck.octave_width - neck.rail_width) +
                        ((neck.nut_width - neck.rail_width) -
                         (neck.octave_width - neck.rail_width)) *
                        i / (len(neck.sections) - 1), 50)
        plt.plot(x, curve.y(x))
    plt.axvline(x=neck.octave_width - neck.rail_width)
    plt.axvline(x=neck.nut_width - neck.rail_width)
    plt.axhline(y=0)
    plt.axhline(y=neck.rail_depth)
    plt.axhline(y=-neck.nut_depth + neck.nut_shell)
    plt.axhline(y=-neck.octave_depth + neck.octave_shell)
    plt.show()


def draw_original():
    original_neck = Neck()
    original_neck.endpoint_bs = np.array([[.625, 4], [1, 2]])
    original_neck.inner_rib = .125
    original_neck.rail_depth = 0
    # Calculate area / weight with solid inner rib
    original_neck.set_model(solid_rib=True)
    # draw_sections(original_neck)
    draw_areas(original_neck)

    from scipy.integrate import quad

    def f(x):
        # return (original_neck.get_inner_curve(0).y(x) +
        #         original_neck.get_outer_curve(0).y(x)) / 2 * \
        #        (original_neck.get_inner_curve(0).y(x) -
        #         original_neck.get_outer_curve(0).y(x))
        return (original_neck.get_inner_curve(0).y(x) -
                original_neck.get_outer_curve(0).y(x))

    # print(quad(f, 0, original_neck.octave_width - original_neck.rail_width))
    from model34_maths.neck_area import NeckArea
    area = NeckArea(
        0,
        original_neck.get_inner_rail(0),  # x1: inner rail edge
        original_neck.get_outer_curve(-1),  # y0: outer curved face
        original_neck.get_inner_curve(-1)
    )
    print(area.area)
    print(area.x_min)
    print(area.x_max)
    x = np.linspace(0, original_neck.get_inner_rail(-1), 50)
    plt.plot(x, original_neck.get_outer_curve(-1).y(x))
    plt.plot(x, original_neck.get_inner_curve(-1).y(x))
    plt.show()
    plt.plot(x, area.y_max(x))
    plt.plot(x, area.y_min(x))
    plt.show()


    # print(original_neck.sections[0].areas)



def draw_optimal():
    opt_neck = Neck()
    opt_neck.endpoint_bs = np.array([0.625, 2.7197315,
                                     1.55822206, 4.04447216]).reshape((2, 2))
    opt_neck.exponents = np.array([[3.82064514, 3.,
                                    3.68748353, 2.]]).reshape((2, 2))
    opt_neck.rail_width = 0.21883729056706025
    opt_neck.nut_shell = 0.0
    opt_neck.octave_shell = 0.010039654395593034
    opt_neck.nut_depth = 0.62481776
    opt_neck.octave_depth = 0.92106462
    opt_neck.inner_rib = 0.125
    # Calculate area / weight with solid inner rib
    opt_neck.set_model(solid_rib=True)
    # draw_sections(opt_neck)
    draw_areas(opt_neck)


if __name__ == '__main__':
    draw_original()
    # draw_optimal()
