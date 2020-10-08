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
    original_neck.endpoint_bs = np.array([[1, 2], [.625, 4], ])
    original_neck.inner_rib = .125
    original_neck.rail_depth = 0
    # Calculate area / weight with solid inner rib
    original_neck.set_model(solid_rib=True)
    # draw_sections(original_neck)
    draw_areas(original_neck)
    print(sum([s.area for s in original_neck.sections]))
    original_neck.set_model()
    print(original_neck.deflection)


def draw_optimal():
    opt_neck = Neck()
    opt_neck.endpoint_bs = np.array([1.55822206, 4.04447216,
                                     0.625, 2.7197315, ]).reshape((2, 2))
    opt_neck.exponents = np.array([[3.68748353, 2.,
                                    3.82064514, 3., ]]).reshape((2, 2))
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


def draw_any(x):
    neck = Neck()
    neck.endpoint_bs = x[0:4].reshape((2, 2))
    neck.exponents = x[4:8].reshape((2, 2))
    neck.rail_width = x[8]
    neck.nut_shell = x[9]
    neck.octave_shell = x[11] - .75 + x[9]
    neck.nut_depth = x[10]
    neck.octave_depth = x[11]
    neck.inner_rib = x[12]
    # Calculate area / weight with solid inner rib
    neck.set_model(solid_rib=True)
    # draw_sections(opt_neck)
    print(sum([s.area for s in neck.sections]))
    draw_areas(neck)
    neck.set_model()
    print(neck.deflection)


if __name__ == '__main__':
    # draw_original()
    # draw_optimal()
    draw_any(np.array([
        1.18552443,
        3.08404271,
        0.82819804,
        0.8125,
        3.82721312,
        2.,
        3.50399796,
        2.,
        0.13062725,
        0.05011796,
        0.5,
        0.93869175,
        0.14332651])
    )
    draw_original()
