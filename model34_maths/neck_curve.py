#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt


class NeckCurve(object):
    def __init__(self, a, b, k, n):
        self.a = a
        self.b = b
        self.k = k
        self.n = n

    def __call__(self, *args, **kwargs):
        return self.y(args[0])

    def y(self, x):
        """

        Note: If y = k, b = 0.

        :param x:
        :type x:
        :return:
        :rtype:
        """
        try:
            return self.k - self.b / self.a ** (self.n / 2) * sqrt(
                self.a ** self.n - x ** self.n
            )
        except ValueError:
            print('domain error')
            print('a:', self.a, 'b:', self.b, 'k:', self.k, 'n:', self.n, 'x:', x)
            return self.k - self.b / self.a ** (self.n / 2) * sqrt(
                self.a ** self.n - x ** self.n
            )
