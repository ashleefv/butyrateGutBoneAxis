from scipy.optimize import fsolve, least_squares
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random,pickle
import numpy as np



# Input values
# butyrate dose in the intestine
bI1 = 0.18

# Data for parameter estimation
# value of distribution
blood_fraction = 0.13       # fraction regulatory T cell in blood
bone_fraction = 0.42        # fraction regulatory T cell in bone

# after butyrate
blood_fraction1 = 0.16      # fraction regulatory T cell in blood after butyrate

# Delta butyrate for parameter estimation (rate constant parameter for with_butyrate)
bI = 0.18
bB = 0.29
bb = 0.29

# evaluating amount from percentage
# assume constant amount of CD4+ T cell in intestine, blood and bone = 1
x7 = blood_fraction/(1-blood_fraction)        # regulatory T cell in the Blood
x8 = bone_fraction/(1-bone_fraction)          # regulatory T cell in the Bone

x71 = blood_fraction1/(1-blood_fraction1)

# constant parameters
b_minus = 0.1718       # without butyrate naive T cell differentiation
muN = 0.02             # Half-life of naive T cells
muT = 0.064            # Half-life of T lymphocytes
muB = 166.3            # butyrate half life day-1

# evaluating migration rate, activity and intestine percentage of Tregs
# x = intestine content, y = migration rate, z = activity
def f(variable):
    x, y, z = variable
    first = b_minus*z - y*x - muT*x
    second = y * x - muT * x7 - y * x7  + b_minus*z
    third = y * x7 - muT * x8 + b_minus*z
    return (first, second, third)

solution =  fsolve(f, (0.1, 0.1, 0.1))

# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb , mb, m1b = variable
    first = fb - muB * bI - mb * bI
    second = mb * bI - muB * bB - m1b * bB
    third = m1b * bB - muB * bb
    return (first, second, third)

solution1 =  fsolve(f, (0.1,0.1,0.1))

# evaluating rate constant parameter for butyrate
def f(variable):
    xnew, b = variable
    first = b_minus*solution[2] - xnew*solution[1] - muT*xnew + b*bI
    second = xnew*solution[1] - muT * x71 - solution[1]*x71  + b_minus*solution[2] + b*bB
    return (first, second)

solution2 =  fsolve(f, (0.1,0.1))

# Updating butyrate dose
def f(variable):
    fb = variable
    first = fb - muB * bI1 -  solution1[1] * bI1
    return (first)

solution3 =  fsolve(f, (0.1))

# evaluated parameters
gamma = solution[2]          # activity
deltaT12 = solution[1]       # Migration of regulatory T cells from the intestine T to the blood
deltaT23 = solution[1]       # Migration of regulatory T cells from the blood T to the bone
FB1 = solution3[0]           # constant formation rate of butyrate
AB12 = solution1[1]          # Absorption of butyrate in blood
AB23 = solution1[2]          # Absorption of butyrate in bone
b_plus = solution2[1]        # rate constant parameter for butyrate

# Initial values
x0 = 0                        # butyrate in the intestine
x1 = 0                        # butyrate in the blood
x2 = 0                        # butyrate in the bone
x3 = 1                        # naive CD4+ T cell in the intestine
x4 = 1                        # naive CD4+ T cell in the Blood
x5 = 1                        # naive CD4+ T cell in the Bone
x6 = solution[0]*x3           # regulatory T cell in the intestine
x7 = x7*x4                    # regulatory T cell in the Blood
x8 = x8*x5                    # regulatory T cell in the Bone

# formation of naive CD4+ T cell without and with butyrate
# without butyrate
FN1_minus = b_minus*gamma*x3 + muN*x3
FN2_minus = b_minus*gamma*x4 + muN*x4
FN3_minus = b_minus*gamma*x5 + muN*x5

# check for influx of naive CD4+ T cells
bB1 = (AB12 * bI1)/(muB + AB23)

if solution3[0] != 0:
    FN1_plus = b_plus*x3*bI1
    FN2_plus = b_plus*x4*bB1
    FN3_plus = b_plus*x5*bB1
else:
    FN1_plus = 0
    FN2_plus = 0
    FN3_plus = 0


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
    y[3] =  FN1_minus - b_minus*gamma*x[3] - muN*x[3] - b_plus*x[3]*x[0] + FN1_plus

    # naive T cells in Blood
    y[4] = FN2_minus - b_minus * gamma * x[4] - muN * x[4] - b_plus * x[4] * x[1] + FN2_plus

    # naive T cells in Bone
    y[5] = FN3_minus - b_minus * gamma * x[5] - muN * x[5] - b_plus * x[5] * x[2] + FN3_plus

    # regulatory T cells in intestine
    y[6] = b_minus*gamma*x[3] + b_plus*x[3]*x[0]- deltaT12 * x[6] - muT * x[6]

    # regulatory T cells in blood
    y[7] = deltaT12 * x[6] - muT * x[7] - deltaT23 * x[7] + b_minus*gamma*x[4] + b_plus*x[4]*x[1]

    # regulatory T cells in bone
    y[8] = deltaT23 * x[7] - muT * x[8] + b_plus * x[5] * x[2] + b_minus*gamma*x[5]

    return y

N = 28        # Total number of days
x = (x0,x1,x2,x3,x4,x5,x6,x7,x8)

T = np.arange(0.0, N, 0.001)
resulth = np.array(odeint(diff, x, T))

con = resulth[-1,8]/(resulth[-1,5] + resulth[-1,8])         # homeostasis value of bone Tregs


# local sensitivity analysis

k = 8           # total number of parameters
store = []

for j in range(k):
    A = [1.0 for i in range(k)]
    A[j] = 1.1

    # Input values

    # butyrate dose in the intestine
    bI1 = 0.18

    # Data for parameter estimation
    # value of distribution
    blood_fraction = 0.13       # fraction regulatory T cell in blood
    bone_fraction = 0.42        # fraction regulatory T cell in bone

    # after butyrate
    blood_fraction1 = 0.16      # fraction regulatory T cell in blood after butyrate

    # Delta butyrate for parameter estimation (rate constant parameter for with_butyrate)
    bI = 0.18
    bB = 0.29
    bb = 0.29

    # evaluating amount from percentage
    # assume constant amount of CD4+ T cell in intestine, blood and bone = 1
    x7 = blood_fraction/(1-blood_fraction)        # regulatory T cell in the Blood
    x8 = bone_fraction/(1-bone_fraction)          # regulatory T cell in the Bone

    x71 = blood_fraction1/(1-blood_fraction1)


    # constant parameters
    b_minus = 0.1718       # without butyrate naive T cell differentiation
    muN = 0.02                 # Half-life of naive T cells
    muT = 0.064                # Half-life of T lymphocytes
    muB = 166.3            # butyrate half life day-1


    # evaluating migration rate, activity and intestine percentage of Tregs
    # x = intestine content, y = migration rate, z = activity
    def f(variable):
        x, y, z = variable
        first = b_minus*z - y*x - muT*x
        second = y * x - muT * x7 - y * x7  + b_minus*z
        third = y * x7 - muT * x8 + b_minus*z
        return (first, second, third)

    solution =  fsolve(f, (0.1, 0.1, 0.1))

    # evaluating constant formation and absorption rate for distribution of butyrate
    def f(variable):
        fb , mb, m1b = variable
        first = fb - muB * bI - mb * bI
        second = mb * bI - muB * bB - m1b * bB
        third = m1b * bB - muB * bb
        return (first, second, third)

    solution1 =  fsolve(f, (0.1,0.1,0.1))


    # evaluating rate constant parameter for butyrate
    def f(variable):
        xnew, b = variable
        first = b_minus*solution[2] - xnew*solution[1] - muT*xnew + b*bI
        second = xnew*solution[1] - muT * x71 - solution[1]*x71  + b_minus*solution[2] + b*bB
        return (first, second)

    solution2 =  fsolve(f, (0.1,0.1))


    # Updating butyrate dose
    def f(variable):
        fb = variable
        first = fb - muB * bI1 -  solution1[1] * bI1
        return (first)

    solution3 =  fsolve(f, (0.1))

    # evaluated parameters
    gamma = solution[2]*A[2]           # activity
    deltaT12 = solution[1]*A[6]       # Migration of regulatory T cells from the intestine T to the blood
    deltaT23 = solution[1]*A[7]           # Migration of regulatory T cells from the blood T to the bone
    FB1 = solution3[0]*A[3]              # constant formation rate of butyrate
    AB12 = solution1[1]*A[4]              # Absorption of butyrate in blood
    AB23 = solution1[2]*A[5]              # Absorption of butyrate in bone
    b_plus = solution2[1]*A[1]      # rate constant parameter for butyrate

    # Initial values
    x0 = 0                        # butyrate in the intestine
    x1 = 0                        # butyrate in the blood
    x2 = 0                        # butyrate in the bone
    x3 = 1                        # naive CD4+ T cell in the intestine
    x4 = 1                        # naive CD4+ T cell in the Blood
    x5 = 1                        # naive CD4+ T cell in the Bone
    x6 = solution[0]*x3              # regulatory T cell in the intestine
    x7 = x7*x4                       # regulatory T cell in the Blood
    x8 = x8*x5                       # regulatory T cell in the Bone

    b_minus = 0.1718 * A[0]
    # formation of naive CD4+ T cell without and with butyrate
    # without butyrate
    FN1_minus = b_minus*gamma*x3 + muN*x3
    FN2_minus = b_minus*gamma*x4 + muN*x4
    FN3_minus = b_minus*gamma*x5 + muN*x5

    # check for influx of naive CD4+ T cells
    bB1 = (AB12 * bI1)/(muB + AB23)

    if solution3[0] != 0:
        FN1_plus = b_plus*x3*bI1
        FN2_plus = b_plus*x4*bB1
        FN3_plus = b_plus*x5*bB1
    else:
        FN1_plus = 0
        FN2_plus = 0
        FN3_plus = 0


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
        y[3] =  FN1_minus - b_minus*gamma*x[3] - muN*x[3] - b_plus*x[3]*x[0] + FN1_plus

        # naive T cells in Blood
        y[4] = FN2_minus - b_minus * gamma * x[4] - muN * x[4] - b_plus * x[4] * x[1] + FN2_plus

        # naive T cells in Bone
        y[5] = FN3_minus - b_minus * gamma * x[5] - muN * x[5] - b_plus * x[5] * x[2] + FN3_plus

        # regulatory T cells in intestine
        y[6] = b_minus*gamma*x[3] + b_plus*x[3]*x[0]- deltaT12 * x[6] - muT * x[6]

        # regulatory T cells in blood
        y[7] = deltaT12 * x[6] - muT * x[7] - deltaT23 * x[7] + b_minus*gamma*x[4] + b_plus*x[4]*x[1]

        # regulatory T cells in bone
        y[8] = deltaT23 * x[7] - muT * x[8] + b_plus * x[5] * x[2] + b_minus*gamma*x[5]

        return y


    N = 28        # Total number of days
    x = (x0,x1,x2,x3,x4,x5,x6,x7,x8)


    T = np.arange(0.0, N, 0.001)
    result = np.array(odeint(diff, x, T))
    store.append(result[-1,8]/(result[-1,5] + result[-1,8]))


# for 10% increase
multiplier = 0.1

# for 10% decrease
#multiplier = -0.1

sensitivity =  np.subtract(store,con)/(con*multiplier)

# plot
# Set position of bar on X axis
r1 = np.arange(8)
plt.rcParams.update({'font.size': 15})
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.yscale('log')
matplotlib.pyplot.bar(r1, sensitivity, width=0.5, bottom=None, align='center', data=None)
matplotlib.pyplot.bar(4, abs(sensitivity[4]), width=0.5, bottom=None, align='center', data=None, color = 'r')
plt.xticks([r for r in range(8)], ['b-','b+','γ','FB1', 'AB12', 'AB23', 'δT12', 'δT23'])
plt.ylabel('Local sensitivity of bone Tregs')
plt.ylim([0.0001,1])
#plt.savefig("images/IECR/" + "10.png", dpi = 300, bbox_inches='tight')
plt.show()
