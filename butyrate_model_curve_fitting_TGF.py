import numpy as np
import math
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random, pickle
import pandas as pd
from scipy.optimize import fsolve, least_squares
from scipy.optimize import curve_fit


# Evaluating steady state time of Tregs

# butyrate dose in the intestine
bI1 = 0.18  # with antibiotic -0.11

# Data for parameter estimation
# value of distribution
blood_fraction = 0.13  # fraction regulatory T cell in blood
bone_fraction = 0.42  # fraction regulatory T cell in bone

# after butyrate
blood_fraction1 = 0.16  # fraction regulatory T cell in blood after butyrate

# Delta butyrate for parameter estimation (rate constant parameter for with_butyrate)
bI = 0.18
bB = 0.29
bb = 0.29

# evaluating amount from percentage
# assume constant amount of CD4+ T cell in intestine, blood and bone = 1
x7 = blood_fraction / (1 - blood_fraction)  # regulatory T cell in the Blood
x8 = bone_fraction / (1 - bone_fraction)  # regulatory T cell in the Bone

x71 = blood_fraction1 / (1 - blood_fraction1)

# constant parameters
b_minus = 0.1718  # without butyrate naive T cell differentiation
muN = 0.02  # Half-life of naive T cells
muT = 0.064  # Half-life of T lymphocytes
muT_beta = 2.0  # Half-life of TGF-beta
muB = 166.3  # butyrate half life day-1


# evaluating migration rate, activity and intestine percentage of Tregs
# x = intestine content, y = migration rate, z = activity
def f(variable):
    x, y, z = variable
    first = b_minus * z - y * x - muT * x
    second = y * x - muT * x7 - y * x7 + b_minus * z
    third = y * x7 - muT * x8 + b_minus * z
    return (first, second, third)


solution = fsolve(f, (0.1, 0.1, 0.1))

# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb, mb, m1b = variable
    first = fb - muB * bI - mb * bI
    second = mb * bI - muB * bB - m1b * bB
    third = m1b * bB - muB * bb
    return (first, second, third)


solution1 = fsolve(f, (0.1, 0.1, 0.1))


# evaluating rate constant parameter for butyrate
def f(variable):
    xnew, b = variable
    first = b_minus * solution[2] - xnew * solution[1] - muT * xnew + b * bI
    second = xnew * solution[1] - muT * x71 - solution[1] * x71 + b_minus * solution[2] + b * bB
    return (first, second)


solution2 = fsolve(f, (0.1, 0.1))


# Updating butyrate dose
def f(variable):
    fb = variable
    first = fb - muB * bI1 - solution1[1] * bI1
    return (first)

solution3 = fsolve(f, (0.1))

# evaluated parameters
gamma = solution[2]  # activity
deltaT12 = solution[1]  # Migration of regulatory T cells from the intestine T to the blood
deltaT23 = solution[1]  # Migration of regulatory T cells from the blood T to the bone
FB1 = solution3[0]  # constant formation rate of butyrate
AB12 = solution1[1]  # Absorption of butyrate in blood
AB23 = solution1[2]  # Absorption of butyrate in bone
b_plus = solution2[1]  # rate constant parameter for butyrate

# Initial values
x0 = 0  # butyrate in the intestine
x1 = 0  # butyrate in the blood
x2 = 0  # butyrate in the bone
x3 = 1  # naive CD4+ T cell in the intestine
x4 = 1  # naive CD4+ T cell in the Blood
x5 = 1  # naive CD4+ T cell in the Bone
x6 = solution[0] * x3  # regulatory T cell in the intestine
x7 = x7 * x4  # regulatory T cell in the Blood
x8 = x8 * x5  # regulatory T cell in the Bone


# calculation of steady state time
def diff(x, T, bI1,):
    # Updating butyrate dose
    def f(variable):
        fb = variable
        first = fb - muB * bI1 - solution1[1] * bI1
        return (first)

    solution3 = fsolve(f, (0.1))

    FB1 = solution3[0]  # constant formation rate of butyrate

    # formation of naive CD4+ T cell without and with butyrate
    # without butyrate
    FN1_minus = b_minus * gamma * x3 + muN * x3
    FN2_minus = b_minus * gamma * x4 + muN * x4
    FN3_minus = b_minus * gamma * x5 + muN * x5

    # with butyrate
    bB1 = (solution1[1] * bI1) / (muB + solution1[2])
    bb1 = bB1 * solution1[2] / muB
    # print(bB1,bb1)

    if solution3[0] != 0:
        FN1_plus = b_plus * x3 * bI1
        FN2_plus = b_plus * x4 * bB1
        FN3_plus = b_plus * x5 * bB1
    else:
        FN1_plus = 0
        FN2_plus = 0
        FN3_plus = 0

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


store = []
N = 280  # Total number of days
x = (x0, x1, x2, x3, x4, x5, x6, x7, x8)

T = np.arange(0.0, N, 0.001)
result = np.array(odeint(diff, x, T, args=(0.18,), atol=1e-22))
result1 = np.array(odeint(diff, x, T, args=(-0.11,), atol=1e-22))

count = int(N / 0.001)
count1 = 20
count2 = 0
count3 = 0
tp = 0
tn = 0
tss = 0

# check to see if the value of bone Tregs remain same for 20 iterations for positive butyrate dose to evaluate steady state time
for i in range(count - 1):
    for j in range(20):
        ss = result[(i + j + 1), 8] - result[(i + j), 8]
        if ss <= 1e-7:
            count2 = count2 + 1
    if count2 == count1:
        tp = i
        break
    count2 = 0

# check to see if the value of bone Tregs remain same for 20 iterations for negative butyrate dose to evaluate steady state time
for i in range(count - 1):
    for j in range(20):
        ss = result1[(i + j), 8] - result1[(i + j + 1), 8]
        if ss <= 1e-7:
            count2 = count2 + 1
    if count2 == count1:
        tn = i
        break
    count2 = 0

if tp > tn:
    tss = tp
else:
    tss = tn

print('Steady state time for bone Tregs is: ', tss* 0.001)


# Estimation of parameter VT_beta and KT_beta
VT_beta = 1             # initial guess of parameter
kT_beta = 1             # initial guess of parameter
rms = 1
rms1 = 1
break_check1 = 1
break_check2 = 0
break_check_counter = 0
break_check_counter1 = 20
counter = 0

while (break_check_counter != break_check_counter1):

    x_guess = 1 + (VT_beta * (result[tss, 8] / x8 - 1)) / (muT_beta * (kT_beta - (result[tss, 8] / x8 - 1)))
    x_guess1 = 1 + (VT_beta * (result1[tss, 8] / x8 - 1)) / (muT_beta * (kT_beta - (result1[tss, 8] / x8 - 1)))

    # curve fit
    y_data = np.zeros((2, 3))
    c_sol1 = np.zeros((2, 3))

    tspan = np.array([0, int(28), int(tss * 0.001)])
    C_data = [1, 3.4, x_guess]
    C_data1 = [1, 0.9, x_guess1]

    y_data[0] = np.array(C_data)
    y_data[1] = np.array(C_data1)
    C_data3 = np.ravel((y_data))


    def myode1(x, t, bI1, k1, k2):
        def f(variable):
            fb = variable
            first = fb - muB * bI1 - solution1[1] * bI1
            return (first)

        solution3 = fsolve(f, (0.1))

        FB1 = solution3[0]  # constant formation rate of butyrate

        # formation of naive CD4+ T cell without and with butyrate
        # without butyrate
        FN1_minus = b_minus * gamma * x3 + muN * x3
        FN2_minus = b_minus * gamma * x4 + muN * x4
        FN3_minus = b_minus * gamma * x5 + muN * x5

        # with butyrate
        bB1 = (solution1[1] * bI1) / (muB + solution1[2])
        bb1 = bB1 * solution1[2] / muB

        if solution3[0] != 0:
            FN1_plus = b_plus * x3 * bI1
            FN2_plus = b_plus * x4 * bB1
            FN3_plus = b_plus * x5 * bB1
        else:
            FN1_plus = 0
            FN2_plus = 0
            FN3_plus = 0

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

        # TGF-beta in bone
        y[9] = muT_beta * (1 - x[9]) + (k1 * (x[8] / x8 - 1)) / (k2 - (x[8] / x8 - 1))

        return y

    def fitfunc(t, k1, k2):

        x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, C_data[0])  # initial value
        C_sol = np.array(odeint(myode1, x, t, args=(0.18, k1, k2)))
        c_sol1[0] = C_sol[:, 9]
        x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, C_data1[0])  # initial value
        C_sol = np.array(odeint(myode1, x, t, args=(-0.11, k1, k2)))
        c_sol1[1] = C_sol[:, 9]

        return np.ravel(c_sol1)


    k_fit, kcov = curve_fit(fitfunc, tspan, C_data3, p0=([VT_beta,kT_beta]))

    VT_beta = k_fit[0]
    kT_beta = k_fit[1]
    x_fit = 1 + (VT_beta * (result[tss, 8] / x8 - 1)) / (muT_beta * (kT_beta - (result[tss, 8] / x8 - 1)))


    def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())

    break_check2 = break_check1
    fit = fitfunc(tspan, k_fit[0], k_fit[1])
    rms = rmse(fit, C_data3)
    break_check1 = rms

    if break_check1 == break_check2:
        break_check_counter +=1


    def fitfunc1(t, k1, k2):
        x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,C_data[0])     #initial value
        C_sol = odeint(myode1, x, t,args=(0.18, k1, k2))
        return C_sol[:,9]

    def fitfunc2(t, k1, k2):

        x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,C_data[0])     #initial value
        C_sol = odeint(myode1, x, t,args=(-0.11, k1, k2))
        return C_sol[:,9]

    tspan1 = np.arange(0.0, tss * 0.001, 0.001)
    fit1 = fitfunc1(tspan1, k_fit[0], k_fit[1])
    fit2 = fitfunc2(tspan1, k_fit[0], k_fit[1])

    plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
    plt.plot(tspan[0:2], C_data[0:2], 'ro', label='Data', markersize=12)
    plt.plot(tspan[0:2], C_data1[0:2], 'ro', markersize=12)
    plt.plot(tspan[2], C_data[2], 'ro', fillstyle = 'none', label='S.S.', markersize=12)
    plt.plot(tspan[2], C_data1[2], 'ro', fillstyle = 'none', markersize=12)
    plt.plot(tspan1, fit1, 'b-', label='Fitted f(t)', linewidth=3)
    plt.plot(tspan1, fit2, 'b-', linewidth=3)
    plt.legend(loc='upper left')
    plt.xlabel('Time (days)')
    plt.ylabel('TGF-β fold change')
    plt.title('RMSE= %f' %rms)
    plt.ylim([0,7.3])
    #plt.savefig('images/smb/' + str(counter) + '.png', dpi=300, bbox_inches='tight')
    counter += 1
    plt.show(block=False)
    plt.pause(1)
    plt.close()

print("parameters are: ", VT_beta, kT_beta)


