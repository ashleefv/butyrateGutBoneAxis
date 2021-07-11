from scipy.optimize import fsolve, least_squares


# Input values

# butyrate dose in the intestine
bI1 = 0.18      # with antibiotic -0.11


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
q4 = 2.762                  # TGF-beta production by Tregs in the target organ
Eta4 = 2.0                  # Half-life of TGF-beta
but_half = 166.3            # butyrate half life day-1


# estimated parameters from both positive and negative butyrate dose
rhoTV = 2.659825171419453           # rate constant for TGF-beta production from Tregs
rhoTk = 0.36332299958391195         # rate constant for TGF-beta production from Tregs

rhoW = 1.7                          # rate constant for Wnt10b production induced by TGF-beta


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
print(solution)

# evaluating constant formation and absorption rate for distribution of butyrate
def f(variable):
    fb , mb, m1b = variable
    first = fb - but_half * bI - mb * bI
    second = mb * bI - but_half * bB - m1b * bB
    third = m1b * bB - but_half * bb
    return (first, second, third)

solution1 =  fsolve(f, (0.1,0.1,0.1))
print(solution1)


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
gamma = solution[2]           # activity
Gamma2 = solution[1]          # Migration of regulatory T cells from the intestine T to the blood
Gamma3 = solution[1]          # Migration of regulatory T cells from the blood T to the bone
Fb = solution3[0]             # constant formation rate of butyrate
Ab = solution1[1]             # Absorption of butyrate in blood
Ab1 = solution1[2]             # Absorption of butyrate in bone
with_but1 = solution2[1]      # rate constant parameter for butyrate

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
formationI = without_but1*gamma*x3 + Eta1*x3
formationB = without_but1*gamma*x4 + Eta1*x4
formationb = without_but1*gamma*x5 + Eta1*x5

#with butyrate
bB1 = (solution1[1] * bI1)/(but_half + solution1[2])
bb1 = bB1*solution1[2]/but_half
#print(bB1,bb1)

if solution3[0] != 0:
    formationI1 = with_but1*x3*bI1
    formationB1 = with_but1*x4*bB1
    formationb1 = with_but1*x5*bB1
else:
    formationI1 = 0
    formationB1 = 0
    formationb1 = 0



# parameters for bone metabolism
S0=180
P0=0
B0=0
C0=0
z0=100

sc = 14  # scaling factor

kM=25
Nc=2
cyclelength=sc
tlag=14*(sc/100)
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

kop = 0.051
kod = 0.196
kodw = 1.045
kca = 0.072
kcd = 0.0034
kzo = 0.00897
kzc = 0.4084

case = 2          # select case