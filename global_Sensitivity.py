from scipy.optimize import fsolve, least_squares
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random,pickle
import numpy as np
import SALib
from SALib.sample import saltelli
from SALib.sample import latin
from SALib.analyze import sobol

# butyrate dose in the intestine
bI1 = 0.18  # with antibiotic -0.11, -0.093, -0.13

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


# evaluating migration rate, activity and intestine percentage of Tregs
# x = intestine content, y = migration rate, z = activity
def f(variable):
    x, y, z = variable
    first = without_but1 * z - y * x - Eta2 * x
    second = y * x - Eta2 * x7 - y * x7 + without_but1 * z
    third = y * x7 - Eta3 * x8 + without_but1 * z
    return (first, second, third)


solution = np.array(fsolve(f, (0.1, 0.1, 0.1)))


# solution = least_squares(f, (0.1, 0.1, 0.1), bounds = ((0, 0, 0), (1, 1, 1)))
# print(solution)

# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb, mb, m1b = variable
    first = fb - but_half * bI - mb * bI
    second = mb * bI - but_half * bB - m1b * bB
    third = m1b * bB - but_half * bb
    return (first, second, third)


solution1 = np.array(fsolve(f, (0.1, 0.1, 0.1)))


# print(solution1)


# evaluating rate constant parameter for butyrate
def f(variable):
    xnew, b = variable
    first = without_but1 * solution[2] - xnew * solution[1] - Eta2 * xnew + b * bI
    second = xnew * solution[1] - Eta2 * x71 - solution[1] * x71 + without_but1 * solution[2] + b * bB
    return (first, second)


solution2 = np.array(fsolve(f, (0.1, 0.1)))


# Updating butyrate dose
def f(variable):
    fb = variable
    first = fb - but_half * bI1 - solution1[1] * bI1
    return (first)


solution3 = np.array(fsolve(f, (0.1)))

# evaluated parameters
gamma = solution[2]  # activity
Gamma2 = solution[1]  # Migration of regulatory T cells from the intestine T to the blood
Gamma3 = solution[1]  # Migration of regulatory T cells from the blood T to the bone
Fb = solution3[0]  # constant formation rate of butyrate
Ab = solution1[1]  # Absorption of butyrate in blood
Ab1 = solution1[2]  # Absorption of butyrate in bone
with_but1 = solution2[1]  # rate constant parameter for butyrate

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


# formation of naive CD4+ T cell without and with butyrate
# without butyrate
formationI = without_but1 * gamma * x3 + Eta1 * x3
formationB = without_but1 * gamma * x4 + Eta1 * x4
formationb = without_but1 * gamma * x5 + Eta1 * x5

# with butyrate
# bB1 = (solution1[1] * bI1)/(but_half + solution1[2])
# bb1 = bB1*solution1[2]/but_half
bB1 = (Ab * bI1) / (but_half + Ab1)
# print(bB1,bb1)

if solution3[0] != 0:
    formationI1 = with_but1 * x3 * bI1
    formationB1 = with_but1 * x4 * bB1
    formationb1 = with_but1 * x5 * bB1
else:
    formationI1 = 0
    formationB1 = 0
    formationb1 = 0


# print(solution1[0], solution1[1], solution1[2], gamma, Gamma2, with_but1)
# print(formationI, formationB,formationb,formationI1,formationB1,formationb1)

# Input values
def diff(x,T, without_but1, with_but1,gamma,Fb,Ab,Ab1,Gamma2,Gamma3):

    # formation of naive CD4+ T cell without and with butyrate
    # without butyrate
    formationI = without_but1*gamma*x3 + Eta1*x3
    formationB = without_but1*gamma*x4 + Eta1*x4
    formationb = without_but1*gamma*x5 + Eta1*x5


    #with butyrate
    #bB1 = (solution1[1] * bI1)/(but_half + solution1[2])
    #bb1 = bB1*solution1[2]/but_half
    bB1 = (Ab * bI1)/(but_half + Ab1)

    formationI1 = with_but1*x3*bI1
    formationB1 = with_but1*x4*bB1
    formationb1 = with_but1*x5*bB1


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


    return y

store = []
N = 28        # Total number of days
x = (x0,x1,x2,x3,x4,x5,x6,x7,x8)


T = np.arange(0.0, N, 0.1)
result = np.array(odeint(diff, x, T, args= (without_but1, with_but1,gamma,Fb,Ab,Ab1,Gamma2,Gamma3)))


# defining problem
problem = {
  'num_vars': 8, #a's, b's and initial condition
  'names': ['b-','b+','γ','FB1', 'AB12', 'AB23', 'δTR12', 'δTR23'],
  'bounds':  np.column_stack((np.array([without_but1, with_but1,gamma,Fb,Ab,Ab1,Gamma2,Gamma3])*1,np.array([without_but1, with_but1,gamma,Fb,Ab,Ab1,Gamma2,Gamma3])*1.1))
}
# 0.1 and 1.9
# Generate samples
vals = saltelli.sample(problem, 500)



Y = np.zeros([len(vals),9])

for i in range(len(vals)):
  Y[i][:] = odeint(diff,x,T,args=(vals[i][0],vals[i][1],vals[i][2],vals[i][3],vals[i][4],vals[i][5],vals[i][6],vals[i][7]))[len(T)-1]



Si = sobol.analyze(problem, Y[:,8], print_to_console=True)
print(Si['S1'])
print(Si['ST'])


# set width of bar
barWidth = 0.25

# Set position of bar on X axis
r1 = np.arange(8)
print(r1)
r2 = [x + barWidth for x in r1]
print(r2)


#plt.figure(figsize=(20.0, 15.0))
plt.rcParams.update({'font.size': 15})
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.yscale('log')
matplotlib.pyplot.bar(r1, Si['S1'], width=0.25, bottom=None, align='center', data=None)
matplotlib.pyplot.bar(r2, Si['ST'], width=0.25, bottom=None, align='center', data=None)
plt.ylabel('Global sensitivity of bone Tregs')
plt.xticks([r + barWidth/2 for r in range(8)], ['b-','b+','γ','FB1', 'AB12', 'AB23', 'δT12', 'δT23'])
plt.ylim([0.0001,1])
#plt.savefig("images/IECR/" + "11.png", dpi = 300, bbox_inches='tight')
plt.show()
