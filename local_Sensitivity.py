from scipy.optimize import fsolve, least_squares
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random,pickle
import numpy as np


k =8
store = []

for j in range(k):
    A = [1.0 for i in range(k)]
    A[j] = 1.1

    # Input values

    # butyrate dose in the intestine
    bI1 = 0.18         # with antibiotic -0.11, -0.093, -0.13


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
    without_but1 = 0.1718       # without butyrate naive T cell differentiation
    Eta1 = 0.02                 # Half-life of naive T cells
    Eta2 = 0.064                # Half-life of T lymphocytes
    Eta3 = 0.064                # Half-life of T lymphocytes
    #q4 = 5.0e-4                 # TGF-beta production by Tregs in the target organ
    q4 = 2.762                     # TGF-beta production by Tregs in the target organ
    Eta4 = 2.0                  # Half-life of TGF-beta
    but_half = 166.3            # butyrate half life day-1


    # evaluating migration rate, activity and intestine percentage of Tregs
    # x = intestine content, y = migration rate, z = activity
    def f(variable):
        x, y, z = variable
        first = without_but1*z - y*x - Eta2*x
        second = y * x - Eta2 * x7 - y * x7  + without_but1*z
        third = y * x7 - Eta3 * x8 + without_but1*z
        return (first, second, third)

    solution =  fsolve(f, (0.1, 0.1, 0.1))
    #solution = least_squares(f, (0.1, 0.1, 0.1), bounds = ((0, 0, 0), (1, 1, 1)))
    #print(solution)

    # evaluating constant formation and absorption rate for distribution of butyrate
    def f(variable):
        fb , mb, m1b = variable
        first = fb - but_half * bI - mb * bI
        second = mb * bI - but_half * bB - m1b * bB
        third = m1b * bB - but_half * bb
        return (first, second, third)

    solution1 =  fsolve(f, (0.1,0.1,0.1))
    #print(solution1)


    # evaluating rate constant parameter for butyrate
    def f(variable):
        xnew, b = variable
        first = without_but1*solution[2] - xnew*solution[1] - Eta2*xnew + b*bI
        second = xnew*solution[1] - Eta2 * x71 - solution[1]*x71  + without_but1*solution[2] + b*bB
        return (first, second)

    solution2 =  fsolve(f, (0.1,0.1))


    # Updating butyrate dose
    def f(variable):
        fb = variable
        first = fb - but_half * bI1 -  solution1[1] * bI1
        return (first)

    solution3 =  fsolve(f, (0.1))

    # evaluated parameters
    gamma = solution[2]*A[2]           # activity
    Gamma2 = solution[1]*A[6]       # Migration of regulatory T cells from the intestine T to the blood
    Gamma3 = solution[1]*A[7]           # Migration of regulatory T cells from the blood T to the bone
    Fb = solution3[0]*A[3]              # constant formation rate of butyrate
    Ab = solution1[1]*A[4]              # Absorption of butyrate in blood
    Ab1 = solution1[2]*A[5]              # Absorption of butyrate in bone
    with_but1 = solution2[1]*A[1]      # rate constant parameter for butyrate

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
    x9 = q4*x8/Eta4                        # TGF-beta
    #x9 = q4*math.pow(x8,6)/Eta4                        # TGF-beta
    x10 = x9

    without_but1 = 0.1718 * A[0]
    # formation of naive CD4+ T cell without and with butyrate
    # without butyrate
    formationI = without_but1*gamma*x3 + Eta1*x3
    formationB = without_but1*gamma*x4 + Eta1*x4
    formationb = without_but1*gamma*x5 + Eta1*x5

    #with butyrate
    #bB1 = (solution1[1] * bI1)/(but_half + solution1[2])
    #bb1 = bB1*solution1[2]/but_half
    #print(bB1,bb1)
    bB1 = (Ab * bI1)/(but_half + Ab1)

    if solution3[0] != 0:
        formationI1 = with_but1*x3*bI1
        formationB1 = with_but1*x4*bB1
        formationb1 = with_but1*x5*bB1
    else:
        formationI1 = 0
        formationB1 = 0
        formationb1 = 0


    def diff(x,T):

        # Differential equations

        y = [0.0 for i in range(len(x))]

        # butyrate in intestine
        y[0] = Fb - but_half * x[0] - Ab * x[0]

        # butyrate in blood
        y[1] = Ab * x[0] - but_half * x[1] - Ab1 * x[1]

        # butyrate in bone
        y[2] = Ab1 * x[1] - but_half * x[2]

        # naive T cells in intestine
        y[3] =  formationI - without_but1*gamma*x[3] - Eta1*x[3] - with_but1*x[3]*x[0] + formationI1

        # naive T cells in Blood
        y[4] = formationB - without_but1 * gamma * x[4] - Eta1 * x[4] - with_but1 * x[4] * x[1] + formationB1

        # naive T cells in Bone
        y[5] = formationb - without_but1 * gamma * x[5] - Eta1 * x[5] - with_but1 * x[5] * x[2] + formationb1

        # regulatory T cells in intestine
        y[6] = without_but1*gamma*x[3] + with_but1*x[3]*x[0]- Gamma2 * x[6] - Eta2 * x[6]

        # regulatory T cells in blood
        y[7] = Gamma2 * x[6] - Eta2 * x[7] - Gamma3 * x[7] + without_but1*gamma*x[4] + with_but1*x[4]*x[1]

        # regulatory T cells in bone
        y[8] = Gamma3 * x[7] - Eta3 * x[8] + with_but1 * x[5] * x[2] + without_but1*gamma*x[5]

        # TGF-beta production and decay
        #y[9] = q4 * x[8] - Eta4 * x[9]
        #y[9] = q4 * x8 + 10.3*q4*(x[8]-x8) - Eta4 * x[9]
        y[9] = q4 * x8 + 1.41 * q4 * ((x[8]-x8)/ x8) - Eta4 * x[9]

        # Wnt10b formation
        #y[10] = q4 * x8 + 1.7*(x[9]-x9) - Eta4 * x[10]
        #y[10] = 0.1*x[9] - Eta4 * x[10]
        y[10] = q4 * x8 + 1.2 * ((x[9]-x9)/ x9) - Eta4 * x[10]

        return y


    N = 28        # Total number of days
    #x = (x0,x1,x2,x3,x4,x5)
    x = (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)


    T = np.arange(0.0, N, 0.001)
    result = np.array(odeint(diff, x, T))
    xx = int((28-0.001)/0.001)
    store.append(result[xx,8]/(result[xx,5] + result[xx,8]))

    #plt.plot(T, result[:,0]/volLymphT, T, result[:,1]/volLymphT,T, result[:,2]/volBlood, T, result[:,3]/volBone, T, result[:,4]/volBone)
    #plt.plot(T, result[:,0], T, result[:,1], T, result[:,2])
    #plt.rcParams.update({'font.size': 25})
    #plt.xticks(np.arange(0, 30, step=7))
    #plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #plt.yticks(np.arange(0, 0.07, step=0.02))
    #plt.plot(T, result[:,7]/(result[:,4] + result[:,7]))
    #plt.plot(T, result[:,6]/(result[:,3] + result[:,6]))
    #plt.plot(T, result[:,8]/(result[:,5] + result[:,8]))
    #plt.plot(T, result[:,8])
    #plt.plot(T, result[:,9], linewidth=3)
    #plt.plot(T, result[:,10], linewidth=3)
    #plt.legend(['Bone Treg'])
    #plt.legend(['with butyrate'])
    #plt.xlabel('TIME (days)')
    #plt.ylabel('Wnt10b')
    #plt.plot(T, result[:,3])
    #plt.ylim([0,1])
    #plt.show()

    #pickle.dump(result, open('resultb.p', 'wb'))
    #pickle.dump(T, open('time.p', 'wb'))
    #print((result[xx,8]/(result[xx,5] + result[xx,8])))

print(store)

resultb = np.array(pickle.load(open('resultb.p', 'rb')))
xx1 = int((28-0.001)/0.001)
con = resultb[xx1,8]/(resultb[xx1,5] + resultb[xx1,8])

multiplier = 0.1

sensitivity =  (store - con)/(con*multiplier)
print(sensitivity)

r1 = np.arange(8)
print(r1)

plt.figure(figsize=(20.0, 15.0))
plt.rcParams.update({'font.size': 25})
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.yscale('log')
#time = np.arange(0.0, k, 1)
#print(time)
#time = ['b-','b+','γ','FB1', 'AB12', 'AB23', 'δTR12', 'δTR23']
#matplotlib.pyplot.bar(time, sensitivity, width=0.5, bottom=None, align='center', data=None)
matplotlib.pyplot.bar(r1, sensitivity, width=0.5, bottom=None, align='center', data=None)
matplotlib.pyplot.bar(4, abs(sensitivity[4]), width=0.5, bottom=None, align='center', data=None, color = 'r')
plt.xticks([r for r in range(8)], ['b-','b+','γ','FB1', 'AB12', 'AB23', 'δTR12', 'δTR23'])
plt.ylabel('Local sensitivity of bone Tregs')
plt.ylim([0.0001,1])
plt.savefig("1.png", dpi = 300, bbox_inches='tight')
plt.show()
