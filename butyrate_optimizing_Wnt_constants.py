from scipy.optimize import fsolve, least_squares


# Input values

bI1 = 0.18      # butyrate dose in the intestine
N = 28         # Number of days

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
muT_beta = 2.0         # Half-life of TGF-beta
muW = 2.0              # Half-life of Wnt10b
muB = 166.3            # butyrate half life day-1


# estimated parameters from both positive and negative butyrate dose
VT_beta = 6                # rate constant for TGF-beta production from Tregs
kT_beta = 0.489               # rate constant for TGF-beta production from Tregs

rhoW = 1.64                          # rate constant for Wnt10b production induced by TGF-beta

# evaluating migration rate, activity and intestine percentage of Tregs
# x = intestine content, y = migration rate, z = activity
def f(variable):
    x, y, z = variable
    first = b_minus*z - y*x - muT*x
    second = y * x - muT * x7 - y * x7  + b_minus*z
    third = y * x7 - muT * x8 + b_minus*z
    return (first, second, third)

solution =  fsolve(f, (0.1, 0.1, 0.1))
#solution = least_squares(f, (0.1, 0.1, 0.1), bounds = ((0, 0, 0), (1, 1, 1)))
print(solution)

# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb , mb, m1b = variable
    first = fb - muB * bI - mb * bI
    second = mb * bI - muB * bB - m1b * bB
    third = m1b * bB - muB * bb
    return (first, second, third)

solution1 =  fsolve(f, (0.1,0.1,0.1))
print(solution1)


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
gamma = solution[2]           # activity
deltaT12 = solution[1]          # Migration of regulatory T cells from the intestine T to the blood
deltaT23 = solution[1]          # Migration of regulatory T cells from the blood T to the bone
FB1 = solution3[0]             # constant formation rate of butyrate
AB12 = solution1[1]             # Absorption of butyrate in blood
AB23 = solution1[2]             # Absorption of butyrate in bone
b_plus = solution2[1]      # rate constant parameter for butyrate

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
x9 = 1                        # TGF-beta in bone
x10 = 1                       # Wnt10b in bone

# formation of naive CD4+ T cell without and with butyrate
# without butyrate
FN1_minus = b_minus*gamma*x3 + muN*x3
FN2_minus = b_minus*gamma*x4 + muN*x4
FN3_minus = b_minus*gamma*x5 + muN*x5

# check for influx of naive CD4+ T cells
bB1 = (solution1[1] * bI1)/(muB + solution1[2])        # blood butyrate
bb1 = bB1*solution1[2]/muB                             # bone butyrate

if solution3[0] != 0:
    FN1_plus = b_plus*x3*bI1
    FN2_plus = b_plus*x4*bB1
    FN3_plus = b_plus*x5*bb1
else:
    FN1_plus = 0
    FN2_plus = 0
    FN3_plus = 0
