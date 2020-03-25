#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np


class NeckConfig(object):
    # n = n
    nut_width = .8125
    nut_depth = .5
    nut_shell = .0625
    octave_width = 1.125
    octave_depth = .75 + .0625
    octave_shell = .0625 + .0625
    rail_width = .1875
    rail_depth = .0625
    scale = 34
    fretboard_thickness = .125
    inner_rib = .125

    load = 36.5 + 42 + 51.3 + 42.8
    elastic_modulus = 10 * 10 ** 6
    string_to_fretboard_distance = .125

    widths = np.array([
        octave_width,
        octave_width - rail_width,
        nut_width,
        nut_width - rail_width,
    ])

    depths = np.array([
        octave_depth,
        octave_depth - octave_shell,
        nut_depth,
        nut_depth - nut_shell,
    ])

    r_outer = 2
    r_inner = 2


if __name__ == '__main__':
    pass
