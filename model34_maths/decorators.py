#!/usr/bin/env python
# -*- coding: utf-8 -*-
import functools


def finite_method(n, a):
    def finite_method_decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            b = args[0].inner.width
            v0 = 0
            for k in range(n):
                _x = a + (b - a) / (2 * n) + (
                        b - a) / n * k
                v0 += f(*args, **kwargs)
            return ((b - a) / n) * v0
        return wrapper
    return finite_method_decorator
