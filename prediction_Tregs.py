import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pickle


from prediction_Tregs_constants import *



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

    return y

x = (x0,x1,x2,x3,x4,x5,x6,x7,x8)

T = np.arange(0.0, N, 0.001)

result = np.array(odeint(diff, x, T))
#pickle.dump(result, open('result1a.p', 'wb'))


# plot
plt.rcParams.update({'font.size': 25})
print(result[:,3])
plt.plot(T, (result[:,8]/(result[:,5] + result[:,8]))*100, linewidth=3)
plt.xlabel('Time (days)')
plt.ylabel('Bone Tregs (%)')
plt.show()

