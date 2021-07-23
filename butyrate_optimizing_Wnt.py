import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fmin

import pickle


from butyrate_optimizing_Wnt_constants import *

def fitfunction(k):

    def diff(x,T):

        # Differential equations

        y = [0.0 for i in range(len(x))]

        # butyrate in intestine
        y[0] = FB1 - muB * x[0] - AB12 * x[0]

        # butyrate in blood
        y[1] = AB12 * x[0] - muB * x[1] - AB23 * x[1]

        # butyrate in bone
        y[2] = AB23 * x[1] - muB * x[2]

        # naive T cells in intestine
        y[3] = FN1_minus - b_minus * gamma * x[3] - muN * x[3] - b_plus * x[3] * x[0] + FN1_plus

        # naive T cells in Blood
        y[4] = FN2_minus - b_minus * gamma * x[4] - muN * x[4] - b_plus * x[4] * x[1] + FN2_plus

        # naive T cells in Bone
        y[5] = FN3_minus - b_minus * gamma * x[5] - muN * x[5] - b_plus * x[5] * x[2] + FN3_plus

        # regulatory T cells in intestine
        y[6] = b_minus * gamma * x[3] + b_plus * x[3] * x[0] - deltaT12 * x[6] - muT * x[6]

        # regulatory T cells in blood
        y[7] = deltaT12 * x[6] - muT * x[7] - deltaT23 * x[7] + b_minus * gamma * x[4] + b_plus * x[4] * x[1]

        # regulatory T cells in bone
        y[8] = deltaT23 * x[7] - muT * x[8] + b_plus * x[5] * x[2] + b_minus * gamma * x[5]

        # TGF-beta production and decay
        y[9] = muT_beta * (1 - x[9]) + VT_beta * (x[8] / x8 - 1) / (kT_beta - (x[8] / x8 - 1))

        # Wnt10b formation
        y[10] = muW * (1 - x[10]) + k * (x[9] - 1)

        return y

    T = np.arange(0.0, N, 0.001)
    x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
    result = np.array(odeint(diff, x, T))

    residual = (result[-1, 10] - 3.2)**2

    return residual


minimum = fmin(fitfunction, x0 = 1)

print(minimum[0])