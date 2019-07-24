#
# femtofitter/pyfsi.py
#

from mpmath import hyp2f2
import numpy as np

HBAR_C = 0.19732697
HBAR_C_SQR = HBAR_C * HBAR_C
ETA = 1.0 / 388


def gamov(qinv):
    """
    Gamov factor
    """
    x = np.true_divide(2.0 * np.pi * HBAR_C * ETA, qinv,
                       where=qinv!=0,
                       out=np.zeros_like(qinv))
    return np.true_divide(x, np.exp(x) - 1.0, where=x!=0.0, out=x)


def calculate(qinv, R):
    """
    Coulomb factor for single (qinv, Rinv) pion pair
    """
    z = -4 * (qinv * R / HBAR_C) ** 2
    try:
        h = float(hyp2f2(.5, 1, 1.5, 1.5, z))
    except TypeError:
        h = np.fromiter(map(float, [hyp2f2(.5, 1, 1.5, 1.5, _z)
                                    for _z in z]), dtype=float)
    k = gamov(qinv) * (1 + 8 * ETA * R / np.sqrt(np.pi) * h)
    return k
