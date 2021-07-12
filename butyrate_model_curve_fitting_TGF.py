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

#result111 = np.array(pickle.load(open('result1.p', 'rb')))
#result222 = np.array(pickle.load(open('result2.p', 'rb')))
#T1 = np.array(pickle.load(open('time.p', 'rb')))


# butyrate dose in the intestine
bI1 = 0.18  # with antibiotic -0.11, -0.093, -0.13

# butyrate dose in the intestine
# bI1 = 0.18  # with antibiotic -0.11, -0.093, -0.13

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
without_but1 = 0.1718  # without butyrate naive T cell differentiation
Eta1 = 0.02  # Half-life of naive T cells
Eta2 = 0.064  # Half-life of T lymphocytes
Eta3 = 0.064  # Half-life of T lymphocytes
# q4 = 5.0e-4                 # TGF-beta production by Tregs in the target organ
q4 = 2.762  # TGF-beta production by Tregs in the target organ
Eta4 = 2.0  # Half-life of TGF-beta
but_half = 166.3  # butyrate half life day-1
# rhoT = 20.6
# rhoTV = 0.49888397
# rhoTk = 0.2539052
#rhoTV = 0.51664553
#rhoTk = 0.33229081
# rhoT1 = 39.
rhoW = 1.7


# evaluating migration rate, activity and intestine percentage of Tregs
# x = intestine content, y = migration rate, z = activity
def f(variable):
    x, y, z = variable
    first = without_but1 * z - y * x - Eta2 * x
    second = y * x - Eta2 * x7 - y * x7 + without_but1 * z
    third = y * x7 - Eta3 * x8 + without_but1 * z
    return (first, second, third)


solution = fsolve(f, (0.1, 0.1, 0.1))
# solution = least_squares(f, (0.1, 0.1, 0.1), bounds = ((0, 0, 0), (1, 1, 1)))
print(solution)


# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb, mb, m1b = variable
    first = fb - but_half * bI - mb * bI
    second = mb * bI - but_half * bB - m1b * bB
    third = m1b * bB - but_half * bb
    return (first, second, third)


solution1 = fsolve(f, (0.1, 0.1, 0.1))
print(solution1)


# evaluating rate constant parameter for butyrate
def f(variable):
    xnew, b = variable
    first = without_but1 * solution[2] - xnew * solution[1] - Eta2 * xnew + b * bI
    second = xnew * solution[1] - Eta2 * x71 - solution[1] * x71 + without_but1 * solution[2] + b * bB
    return (first, second)


solution2 = fsolve(f, (0.1, 0.1))


# Updating butyrate dose
def f(variable):
    fb = variable
    first = fb - but_half * bI1 - solution1[1] * bI1
    return (first)


solution3 = fsolve(f, (0.1))

# evaluated parameters
gamma = solution[2]  # activity
Gamma2 = solution[1]  # Migration of regulatory T cells from the intestine T to the blood
Gamma3 = solution[1]  # Migration of regulatory T cells from the blood T to the bone
Fb = solution3[0]  # constant formation rate of butyrate
Ab = solution1[1]  # Absorption of butyrate in blood
Ab1 = solution1[2]  # Absorption of butyrate in bone
with_but1 = solution2[1]  # rate constant parameter for butyrate

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
# x9 = q4*x8/Eta4                       # TGF-beta
x9 = 1
# x9 = q4*math.pow(x8,6)/Eta4                        # TGF-beta
x10 = x9


# sum = x8
def diff(x, T, bI1, null):
    # Updating butyrate dose
    def f(variable):
        fb = variable
        first = fb - but_half * bI1 - solution1[1] * bI1
        return (first)

    solution3 = fsolve(f, (0.1))

    Fb = solution3[0]  # constant formation rate of butyrate

    # formation of naive CD4+ T cell without and with butyrate
    # without butyrate
    formationI = without_but1 * gamma * x3 + Eta1 * x3
    formationB = without_but1 * gamma * x4 + Eta1 * x4
    formationb = without_but1 * gamma * x5 + Eta1 * x5

    # with butyrate
    bB1 = (solution1[1] * bI1) / (but_half + solution1[2])
    bb1 = bB1 * solution1[2] / but_half
    # print(bB1,bb1)

    if solution3[0] != 0:
        formationI1 = with_but1 * x3 * bI1
        formationB1 = with_but1 * x4 * bB1
        formationb1 = with_but1 * x5 * bB1
    else:
        formationI1 = 0
        formationB1 = 0
        formationb1 = 0

    # Differential equations

    y = [0.0 for i in range(len(x))]

    # butyrate in intestine
    y[0] = Fb - but_half * x[0] - Ab * x[0]

    # butyrate in blood
    y[1] = Ab * x[0] - but_half * x[1] - Ab1 * x[1]

    # butyrate in bone
    y[2] = Ab1 * x[1] - but_half * x[2]

    # naive T cells in intestine
    y[3] = formationI - without_but1 * gamma * x[3] - Eta1 * x[3] - with_but1 * x[3] * x[0] + formationI1

    # naive T cells in Blood
    y[4] = formationB - without_but1 * gamma * x[4] - Eta1 * x[4] - with_but1 * x[4] * x[1] + formationB1

    # naive T cells in Bone
    y[5] = formationb - without_but1 * gamma * x[5] - Eta1 * x[5] - with_but1 * x[5] * x[2] + formationb1

    # regulatory T cells in intestine
    y[6] = without_but1 * gamma * x[3] + with_but1 * x[3] * x[0] - Gamma2 * x[6] - Eta2 * x[6]

    # regulatory T cells in blood
    y[7] = Gamma2 * x[6] - Eta2 * x[7] - Gamma3 * x[7] + without_but1 * gamma * x[4] + with_but1 * x[4] * x[1]

    # regulatory T cells in bone
    y[8] = Gamma3 * x[7] - Eta3 * x[8] + with_but1 * x[5] * x[2] + without_but1 * gamma * x[5]

    return y


store = []
N = 280  # Total number of days
# x = (x0,x1,x2,x3,x4,x5)
x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

T = np.arange(0.0, N, 0.001)
result = np.array(odeint(diff, x, T, args=(0.18, 0), atol=1e-22))
result1 = np.array(odeint(diff, x, T, args=(-0.13, 0), atol=1e-22))

count = int(N / 0.001)
count1 = 20
count2 = 0
count3 = 0
tp = 0
tn = 0
tss = 0
for i in range(count - 1):
    for j in range(20):
        ss = result[(i + j + 1), 8] - result[(i + j), 8]
        # print (ss)
        if ss <= 1e-7:
            count2 = count2 + 1
        # ss = result[i, 8] - result[i + 1, 8]
    if count2 == count1:
        tp = i
        print(i * 0.001, 'success')
        break
    count2 = 0

for i in range(count - 1):
    for j in range(20):
        ss = result1[(i + j), 8] - result1[(i + j + 1), 8]
        # print (ss)
        if ss <= 1e-7:
            count2 = count2 + 1
        # ss = result[i, 8] - result[i + 1, 8]
    if count2 == count1:
        tn = i
        print(i * 0.001, 'success')
        break
    count2 = 0
print(tp,tn)
input()
if tp > tn:
    tss = tp
else:
    tss = tn

print('tss', tss)

x_guess = 0
x_fit = 0

rhoTV = .5
rhoTk = 0.33
rms = 1
rms1 = 1
break_check1 = 1
break_check2 = 0
break_check_counter = 0
break_check_counter1 = 20
# while (x_guess - x_fit >= 1e-15):
counter = 0
while (break_check_counter != break_check_counter1):

    x_guess = 1 + (rhoTV * (result[tss, 8] / x8 - 1)) / (Eta4 * (rhoTk - (result[tss, 8] / x8 - 1)))
    x_guess1 = 1 + (rhoTV * (result1[tss, 8] / x8 - 1)) / (Eta4 * (rhoTk - (result1[tss, 8] / x8 - 1)))

    # curve fit
    y_data = np.zeros((2, 3))
    c_sol1 = np.zeros((2, 3))

    tspan = np.array([0, int(28), int(tss * 0.001)])
    C_data = [1, 3.4, x_guess]
    C_data1 = [1, 0.9, x_guess1]

    y_data[0] = np.array(C_data)
    y_data[1] = np.array(C_data1)
    C_data3 = np.ravel((y_data))
    print('ravel', C_data3)


    def fitfunc(t, k1, k2):
        print(t)

        def myode(x, t, bI1, null):
            # Updating butyrate dose
            def f(variable):
                fb = variable
                first = fb - but_half * bI1 - solution1[1] * bI1
                return (first)

            solution3 = fsolve(f, (0.1))

            Fb = solution3[0]  # constant formation rate of butyrate

            # formation of naive CD4+ T cell without and with butyrate
            # without butyrate
            formationI = without_but1 * gamma * x3 + Eta1 * x3
            formationB = without_but1 * gamma * x4 + Eta1 * x4
            formationb = without_but1 * gamma * x5 + Eta1 * x5

            # with butyrate
            bB1 = (solution1[1] * bI1) / (but_half + solution1[2])
            bb1 = bB1 * solution1[2] / but_half
            # print(bB1,bb1)

            if solution3[0] != 0:
                formationI1 = with_but1 * x3 * bI1
                formationB1 = with_but1 * x4 * bB1
                formationb1 = with_but1 * x5 * bB1
            else:
                formationI1 = 0
                formationB1 = 0
                formationb1 = 0

            y = [0.0 for i in range(len(x))]
            # butyrate in intestine
            y[0] = Fb - but_half * x[0] - Ab * x[0]

            # butyrate in blood
            y[1] = Ab * x[0] - but_half * x[1] - Ab1 * x[1]

            # butyrate in bone
            y[2] = Ab1 * x[1] - but_half * x[2]

            # naive T cells in intestine
            y[3] = formationI - without_but1 * gamma * x[3] - Eta1 * x[3] - with_but1 * x[3] * x[0] + formationI1

            # naive T cells in Blood
            y[4] = formationB - without_but1 * gamma * x[4] - Eta1 * x[4] - with_but1 * x[4] * x[1] + formationB1

            # naive T cells in Bone
            y[5] = formationb - without_but1 * gamma * x[5] - Eta1 * x[5] - with_but1 * x[5] * x[2] + formationb1

            # regulatory T cells in intestine
            y[6] = without_but1 * gamma * x[3] + with_but1 * x[3] * x[0] - Gamma2 * x[6] - Eta2 * x[6]

            # regulatory T cells in blood
            y[7] = Gamma2 * x[6] - Eta2 * x[7] - Gamma3 * x[7] + without_but1 * gamma * x[4] + with_but1 * x[4] * x[1]

            # regulatory T cells in bone
            y[8] = Gamma3 * x[7] - Eta3 * x[8] + with_but1 * x[5] * x[2] + without_but1 * gamma * x[5]
            y[9] = Eta4 * (1 - x[9]) + (k1 * (x[8] / x8 - 1)) / (k2 - (x[8] / x8 - 1))
            return y

        x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, C_data[0])  # initial value
        C_sol = np.array(odeint(myode, x, t, args=(0.18, 0)))
        # c_sol1.append(C_sol)
        c_sol1[0] = C_sol[:, 9]
        x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, C_data1[0])  # initial value
        C_sol = np.array(odeint(myode, x, t, args=(-0.11, 0)))
        # c_sol1.append(C_sol)
        c_sol1[1] = C_sol[:, 9]

        print('thissssss', c_sol1)
        return np.ravel(c_sol1)


    k_fit, kcov = curve_fit(fitfunc, tspan, C_data3, p0=([rhoTV,rhoTk]))

    rhoTV = k_fit[0]
    rhoTk = k_fit[1]
    x_fit = 1 + (rhoTV * (result[tss, 8] / x8 - 1)) / (Eta4 * (rhoTk - (result[tss, 8] / x8 - 1)))
    # x_fit = 1 + (rhoTV * (result1[tss, 8] / x8 - 1)) / (Eta4 * (rhoTk - (result1[tss, 8] / x8 - 1)))

    print('thatttttt', k_fit[0], k_fit[1])


    def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())


    break_check2 = break_check1
    fit = fitfunc(tspan, k_fit[0], k_fit[1])
    rms = rmse(fit, C_data3)
    # rms1 = rmse(fit1, C_data1)
    print('rmssss', rms)
    print('constantssss', k_fit[0], k_fit[1])
    break_check1 = rms


    print('break check', break_check1, break_check2)
    if break_check1 == break_check2:
        break_check_counter +=1

    def fitfunc1(t, k1, k2):
        print(t)
        def myode(x, t, bI1, null):
            # Updating butyrate dose
            def f(variable):
                fb = variable
                first = fb - but_half * bI1 - solution1[1] * bI1
                return (first)

            solution3 = fsolve(f, (0.1))

            Fb = solution3[0]  # constant formation rate of butyrate

            # formation of naive CD4+ T cell without and with butyrate
            # without butyrate
            formationI = without_but1 * gamma * x3 + Eta1 * x3
            formationB = without_but1 * gamma * x4 + Eta1 * x4
            formationb = without_but1 * gamma * x5 + Eta1 * x5

            # with butyrate
            bB1 = (solution1[1] * bI1) / (but_half + solution1[2])
            bb1 = bB1 * solution1[2] / but_half
            # print(bB1,bb1)

            if solution3[0] != 0:
                formationI1 = with_but1 * x3 * bI1
                formationB1 = with_but1 * x4 * bB1
                formationb1 = with_but1 * x5 * bB1
            else:
                formationI1 = 0
                formationB1 = 0
                formationb1 = 0

            y = [0.0 for i in range(len(x))]
            # butyrate in intestine
            y[0] = Fb - but_half * x[0] - Ab * x[0]

            # butyrate in blood
            y[1] = Ab * x[0] - but_half * x[1] - Ab1 * x[1]

            # butyrate in bone
            y[2] = Ab1 * x[1] - but_half * x[2]

            # naive T cells in intestine
            y[3] = formationI - without_but1 * gamma * x[3] - Eta1 * x[3] - with_but1 * x[3] * x[0] + formationI1

            # naive T cells in Blood
            y[4] = formationB - without_but1 * gamma * x[4] - Eta1 * x[4] - with_but1 * x[4] * x[1] + formationB1

            # naive T cells in Bone
            y[5] = formationb - without_but1 * gamma * x[5] - Eta1 * x[5] - with_but1 * x[5] * x[2] + formationb1

            # regulatory T cells in intestine
            y[6] = without_but1 * gamma * x[3] + with_but1 * x[3] * x[0] - Gamma2 * x[6] - Eta2 * x[6]

            # regulatory T cells in blood
            y[7] = Gamma2 * x[6] - Eta2 * x[7] - Gamma3 * x[7] + without_but1 * gamma * x[4] + with_but1 * x[4] * x[1]

            # regulatory T cells in bone
            y[8] = Gamma3 * x[7] - Eta3 * x[8] + with_but1 * x[5] * x[2] + without_but1 * gamma * x[5]
            y[9] =  Eta4*(1-x[9]) + (k1*(x[8]/x8-1))/(k2 - (x[8]/x8-1))
            return y

        x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,C_data[0])     #initial value
        C_sol = odeint(myode, x, t,args=(0.18,0))
        #print('thissssss', C_sol)
        return C_sol[:,9]

    def fitfunc2(t, k1, k2):
        print(t)
        def myode(x, t, bI1, null):
            # Updating butyrate dose
            def f(variable):
                fb = variable
                first = fb - but_half * bI1 - solution1[1] * bI1
                return (first)

            solution3 = fsolve(f, (0.1))

            Fb = solution3[0]  # constant formation rate of butyrate

            # formation of naive CD4+ T cell without and with butyrate
            # without butyrate
            formationI = without_but1 * gamma * x3 + Eta1 * x3
            formationB = without_but1 * gamma * x4 + Eta1 * x4
            formationb = without_but1 * gamma * x5 + Eta1 * x5

            # with butyrate
            bB1 = (solution1[1] * bI1) / (but_half + solution1[2])
            bb1 = bB1 * solution1[2] / but_half
            # print(bB1,bb1)

            if solution3[0] != 0:
                formationI1 = with_but1 * x3 * bI1
                formationB1 = with_but1 * x4 * bB1
                formationb1 = with_but1 * x5 * bB1
            else:
                formationI1 = 0
                formationB1 = 0
                formationb1 = 0

            y = [0.0 for i in range(len(x))]
            # butyrate in intestine
            y[0] = Fb - but_half * x[0] - Ab * x[0]

            # butyrate in blood
            y[1] = Ab * x[0] - but_half * x[1] - Ab1 * x[1]

            # butyrate in bone
            y[2] = Ab1 * x[1] - but_half * x[2]

            # naive T cells in intestine
            y[3] = formationI - without_but1 * gamma * x[3] - Eta1 * x[3] - with_but1 * x[3] * x[0] + formationI1

            # naive T cells in Blood
            y[4] = formationB - without_but1 * gamma * x[4] - Eta1 * x[4] - with_but1 * x[4] * x[1] + formationB1

            # naive T cells in Bone
            y[5] = formationb - without_but1 * gamma * x[5] - Eta1 * x[5] - with_but1 * x[5] * x[2] + formationb1

            # regulatory T cells in intestine
            y[6] = without_but1 * gamma * x[3] + with_but1 * x[3] * x[0] - Gamma2 * x[6] - Eta2 * x[6]

            # regulatory T cells in blood
            y[7] = Gamma2 * x[6] - Eta2 * x[7] - Gamma3 * x[7] + without_but1 * gamma * x[4] + with_but1 * x[4] * x[1]

            # regulatory T cells in bone
            y[8] = Gamma3 * x[7] - Eta3 * x[8] + with_but1 * x[5] * x[2] + without_but1 * gamma * x[5]
            y[9] =  Eta4*(1-x[9]) + (k1*(x[8]/x8-1))/(k2 - (x[8]/x8-1))
            return y

        x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,C_data[0])     #initial value
        C_sol = odeint(myode, x, t,args=(-0.11,0))
        #print('thissssss', C_sol)
        return C_sol[:,9]



    # tfit = np.linspace(1, 50)
    tspan1 = np.arange(0.0, 80.736, 0.001)
    fit1 = fitfunc1(tspan1, k_fit[0], k_fit[1])
    fit2 = fitfunc2(tspan1, k_fit[0], k_fit[1])


    plt.rcParams.update({'font.size': 25,'legend.fontsize': 20})
    #plt.rcParams.update({'font.size': 25})
    # plt.xticks(np.arange(0, 30, step=7))
    #plt.yticks(np.arange(0, 8, step=1))
    plt.plot(tspan[0:2], C_data[0:2], 'ro', label='Data', markersize=12)
    plt.plot(tspan[0:2], C_data1[0:2], 'ro', markersize=12)
    plt.plot(tspan[2], C_data[2], 'ro', fillstyle = 'none', label='S.S.', markersize=12)
    plt.plot(tspan[2], C_data1[2], 'ro', fillstyle = 'none', markersize=12)



    plt.plot(tspan1, fit1, 'b-', label='Fitted f(t)', linewidth=3)
    plt.plot(tspan1, fit2, 'b-', linewidth=3)

    #plt.plot(T1, result111[:, 9], color='orange', label='Fitted f(t)', linewidth=3)
    #plt.plot(T1, result222[:, 9], color='orange', linewidth=3)
    #plt.legend(loc='upper left', framealpha=0.05)
    plt.legend(loc='upper left')
    plt.xlabel('Time (days)')
    plt.ylabel('TGF-β fold change')
    plt.title('RMSE= %f' %rms)
    plt.ylim([0,7.3])
    #plt.savefig('images/smb/' + str(counter) + '.png', dpi=300, bbox_inches='tight')
    counter += 1
    # plt.hold(False)
    # plt.show()
    plt.show(block=False)
    plt.pause(1)
    plt.close()

print(rhoTV, rhoTk)
print('check1', result[-1, 8], result[-1, 8] / x8 - 1)
print('check2', result1[-1, 8], result1[-1, 8] / x8 - 1)
# print((result[-1,8]/ x8 - 1)/(rhoT + (result[-1,8]/ x8 - 1)))
# print(q4*x8/Eta4)
# dWnt10/dt = k TGF-beta (fold change and change from baseline)
# plt.plot(T, result[:,0]/volLymphT, T, result[:,1]/volLymphT,T, result[:,2]/volBlood, T, result[:,3]/volBone, T, result[:,4]/volBone)
# plt.plot(T, result[:,0], T, result[:,1], T, result[:,2])
plt.rcParams.update({'font.size': 25})
# plt.xticks(np.arange(0, 30, step=7))
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# plt.yticks(np.arange(0, 0.07, step=0.02))
# plt.plot(T, result[:,7]/(result[:,4] + result[:,7]))
# plt.plot(T, result[:,6]/(result[:,3] + result[:,6]))
# plt.plot(T, result[:,8]/(result[:,5] + result[:,8]), linewidth=3)
# plt.plot(T, TGF)
plt.plot(T, result[:, 8], linewidth=3)
# plt.plot(T, result[:,10], linewidth=3)
# plt.legend(['Bone Treg'])
plt.legend(['with butyrate'])
plt.xlabel('Time (days)')
plt.ylabel('TGF-beta fold change')
# plt.ylabel('Fraction of Naive Tregs in bone')
# plt.plot(T, result[:,3])
# plt.ylim([0.55,0.74])
plt.show()

# pickle.dump(result, open('result1.p', 'wb'))
# pickle.dump(T, open('time.p', 'wb'))
# dfT = pd.DataFrame(result[:,6]/(result[:,3] + result[:,6]))
# dfT = pd.DataFrame(result[:,0])
# dfT.to_excel("Treg_mean6.xlsx")