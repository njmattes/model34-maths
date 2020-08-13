#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import sqrt, fabs


class NeckCurve(object):
    def __init__(self, a, b, k, n):
        """Hyperelliptical curve in the form


        |  x  |^n | (y - k) |^n
        | --- | + | ------- | = 1
        |  a  |   |    b    |

        :param a: Minor horizontal radius
        :type a: float
        :param b: Minor vertical radius
        :type b: float
        :param k: Vertical offset
        :type k: float
        :param n: Degree of curvature
        :type n: float
        """
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
        if self.b == 0:
            return self.k
        try:
            # print(self.a, self.b, self.k, self.n)
            return (self.k -
                    self.b * (1 - fabs(x / self.a) ** self.n) ** (1 / self.n))
        except ValueError:
            print('domain error, √a^n-x^n:',
                  self.a ** self.n - x ** self.n, 'k:', self.k, 'a:', self.a)
            return self.k

    def x(self, y):
        """

        :param y:
        :type y:
        :return:
        :rtype:
        """
        try:
            return (
                self.a *
                (1 - fabs((self.k - y) / self.b) ** self.n) ** (1 / self.n))
            # return (self.a * sqrt(self.b**2 - (self.k - y)**2) / self.b)
        except ValueError:
            print('domain error')
            print('a:', self.a, 'b:', self.b, 'k:', self.k,
                  'n:', self.n, 'y:', y)
            return self.k - self.b / self.a ** (self.n / 2) * sqrt(
                self.a ** self.n - y ** self.n
            )
