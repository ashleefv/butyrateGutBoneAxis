from scipy.optimize import fsolve, least_squares


# Input values

# butyrate dose in the intestine
bI1 = 0.18      # with antibiotic -0.11

# input parameters for bone metabolism
Nc = 2                       # number of cycle (current code works for 1, 2, 6, and 12 cycle)
cyclelength= 14              # length of each remodeling cycle
case = 2                     # select case (for case 1 set bI1 = 0)

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
VT_beta = 2.86         # rate constant for TGF-beta production from Tregs
kT_beta = 0.373        # rate constant for TGF-beta production from Tregs

rhoW = 1.72            # rate constant for Wnt10b production induced by TGF-beta

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



# parameters for bone metabolism

sc = cyclelength       # scaling factor to convert 100 days cycle

# estimated parameters
kop = 0.051
kod = 0.196
kodw = 1.045
kca = 0.072
kcd = 0.0034
kzo = 0.00897
kzc = 0.4084



# initial value
S0=180
P0=0
B0=0
C0=0
z0=100

# default parameters
kM=25
kg_1= 2.15e-1
kg_2= 2.98e-2
kg_3= 1.28e-3
kg_4= 8.42
Bone=1
alpha_1=0.5
alpha_2=0.1
alpha_3=0.1
beta_1=0.1
delta=0.1
beta_2=0.1
alpha_4=0.1
K_S=200
k1 = .698331
k2=0.015445
g_31=1
g_21=2
g_22=1
g_32=1
g_41=1
g_42=1
g_43=-1
g_44=1
f_12=1
f_14=1
f_23=1
f_34=1
epsilon=1
beta_3=0.1
rho=20

